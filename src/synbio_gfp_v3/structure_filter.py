from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


def parse_ca_atoms(pdb_path: str | Path) -> dict[int, np.ndarray]:
    coords = {}
    for line in Path(pdb_path).read_text(errors="ignore").splitlines():
        if not line.startswith("ATOM"):
            continue
        atom = line[12:16].strip()
        if atom != "CA":
            continue
        try:
            resi = int(line[22:26])
            x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            coords[resi] = np.array([x, y, z], dtype=float)
        except Exception:
            continue
    return coords


def rmsd(a: np.ndarray, b: np.ndarray) -> float:
    if len(a) == 0 or len(a) != len(b):
        return float("nan")
    # Kabsch alignment for robust C-alpha RMSD.
    a0 = a - a.mean(axis=0)
    b0 = b - b.mean(axis=0)
    u, _, vt = np.linalg.svd(a0.T @ b0)
    r = u @ vt
    ar = a0 @ r
    return float(np.sqrt(((ar - b0) ** 2).sum() / len(a)))


def compute_pdb_metrics(reference_pdb: str | Path, candidate_pdb: str | Path) -> dict[str, float | str]:
    ref = parse_ca_atoms(reference_pdb)
    cand = parse_ca_atoms(candidate_pdb)
    common = sorted(set(ref) & set(cand))
    barrel = [i for i in common if 1 <= i <= 239]
    pocket = [i for i in common if i in {39, 63, 64, 65, 66, 67, 68, 69, 96, 145, 148, 203, 205, 222}]
    ion = [i for i in common if i in {17, 30, 32, 115, 122}]
    def subset_rmsd(indices):
        return rmsd(np.vstack([ref[i] for i in indices]), np.vstack([cand[i] for i in indices])) if len(indices) >= 3 else float("nan")
    return {
        "barrel_rmsd": subset_rmsd(barrel),
        "pocket_rmsd": subset_rmsd(pocket),
        "ion_network_rmsd": subset_rmsd(ion),
        "structure_status": "pdb_checked",
    }


def proxy_structure_metrics(row: pd.Series) -> dict[str, float | str]:
    protected_count = float(row.get("protected_mutation_count", 0))
    pocket_count = float(row.get("pocket_mutation_count", 0))
    ion_count = float(row.get("ion_network_mutation_count", 0))
    score = max(0.0, 1.0 - 0.22 * protected_count - 0.18 * pocket_count - 0.28 * ion_count)
    return {
        "barrel_rmsd": np.nan,
        "pocket_rmsd": np.nan,
        "ion_network_rmsd": np.nan,
        "structure_score": score,
        "structure_status": "proxy_only_no_pdb",
    }


def run_top200_structure_filter(df: pd.DataFrame, config: dict[str, Any], out_dir: str | Path, logger: logging.Logger) -> pd.DataFrame:
    cfg = config.get("structure", {})
    enabled = bool(cfg.get("enabled", True))
    top_n = int(cfg.get("top_n", 200))
    out = df.copy()
    out["structure_score"] = 0.5
    out["structure_status"] = "not_checked"
    if not enabled:
        logger.info("STRUCTURE FILTER SKIPPED | disabled by config")
        return out
    top_idx = out.sort_values("objective", ascending=False).head(top_n).index
    ref_pdb = cfg.get("reference_pdb")
    pred_dir = cfg.get("predicted_pdb_dir")
    logger.info("STRUCTURE FILTER START | top_n=%d | reference_pdb=%s | predicted_pdb_dir=%s", len(top_idx), ref_pdb, pred_dir)
    records = []
    for idx in top_idx:
        row = out.loc[idx]
        metrics = None
        if ref_pdb and pred_dir:
            seq_id = str(row.get("Seq_ID", f"cand_{idx}"))
            candidate_paths = [Path(pred_dir) / f"{seq_id}.pdb", Path(pred_dir) / f"candidate_{idx}.pdb"]
            pdb = next((p for p in candidate_paths if p.exists()), None)
            if pdb:
                metrics = compute_pdb_metrics(ref_pdb, pdb)
                score = 1.0
                if not np.isnan(metrics["barrel_rmsd"]):
                    score -= min(0.45, metrics["barrel_rmsd"] / 5.0)
                if not np.isnan(metrics["pocket_rmsd"]):
                    score -= min(0.35, metrics["pocket_rmsd"] / 4.0)
                if not np.isnan(metrics["ion_network_rmsd"]):
                    score -= min(0.20, metrics["ion_network_rmsd"] / 4.0)
                metrics["structure_score"] = float(max(0.0, score))
        if metrics is None:
            metrics = proxy_structure_metrics(row)
        for k, v in metrics.items():
            out.loc[idx, k] = v
        records.append({"index": int(idx), "sequence": row["sequence"], **metrics})
    pd.DataFrame(records).to_csv(Path(out_dir) / "structure_filter_top200.csv", index=False)
    logger.info("STRUCTURE FILTER DONE | checked=%d | output=%s", len(records), Path(out_dir) / "structure_filter_top200.csv")
    return out
