#!/usr/bin/env python3
"""Augment V7.2 structure metrics with global-fit GFP geometry and PDB provenance.

The primary V7.1 local RMSD columns remain untouched. This script computes an additional
single-transform audit: a Kabsch transform fitted on the beta-barrel is reused for the
chromophore pocket and ion-network subsets. The new values are advisory metrics for the
matrix report and never alter the primary final-selection gate.
"""
from __future__ import annotations

import argparse
import hashlib
import math
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

DEFAULT_BARREL = "1-239"
DEFAULT_POCKET = "65,66,67,68,69,70,95,96,145,146,147,148,203,204"
DEFAULT_ION = "39,40,41,42,43,44,96,97,98,99,164,165,166,167"


def parse_positions(spec: str) -> list[int]:
    """Parse compact residue ranges such as ``1-239,245`` into sorted unique positions."""
    values: list[int] = []
    for token in str(spec).split(","):
        token = token.strip()
        if not token:
            continue
        if "-" in token:
            left, right = token.split("-", 1)
            values.extend(range(int(left), int(right) + 1))
        else:
            values.append(int(token))
    return sorted(set(values))


def parse_ca_atoms(path: Path) -> dict[int, tuple[np.ndarray, float]]:
    """Read CA coordinates and pLDDT-like B-factors from the first compatible chain."""
    atoms: dict[int, tuple[np.ndarray, float]] = {}
    active_chain: str | None = None
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.startswith("ATOM") or line[12:16].strip() != "CA":
                continue
            chain = line[21].strip()
            if active_chain is None:
                active_chain = chain
            if chain != active_chain:
                continue
            try:
                resseq = int(line[22:26])
                xyz = np.array(
                    [float(line[30:38]), float(line[38:46]), float(line[46:54])],
                    dtype=np.float64,
                )
                bfactor = float(line[60:66]) if len(line) >= 66 else math.nan
            except ValueError:
                continue
            atoms[resseq] = (xyz, bfactor)
    return atoms


def fit_kabsch(reference: np.ndarray, mobile: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return row-vector rotation and translation mapping mobile coordinates to reference."""
    if reference.shape != mobile.shape or reference.ndim != 2 or reference.shape[0] < 3:
        raise ValueError("Kabsch fit requires at least three paired coordinates.")
    ref_center = reference.mean(axis=0)
    mob_center = mobile.mean(axis=0)
    ref0 = reference - ref_center
    mob0 = mobile - mob_center
    u, _, vt = np.linalg.svd(mob0.T @ ref0)
    correction = np.eye(3)
    correction[-1, -1] = np.sign(np.linalg.det(u @ vt))
    rotation = u @ correction @ vt
    translation = ref_center - mob_center @ rotation
    return rotation, translation


def rmsd_after_transform(
    ref_atoms: dict[int, tuple[np.ndarray, float]],
    pred_atoms: dict[int, tuple[np.ndarray, float]],
    positions: Iterable[int],
    rotation: np.ndarray,
    translation: np.ndarray,
) -> float:
    """Measure a residue subset using the barrel-fitted transform without local realignment."""
    common = [pos for pos in positions if pos in ref_atoms and pos in pred_atoms]
    if len(common) < 3:
        return math.nan
    reference = np.stack([ref_atoms[pos][0] for pos in common])
    mobile = np.stack([pred_atoms[pos][0] for pos in common])
    aligned = mobile @ rotation + translation
    return float(np.sqrt(np.mean(np.sum((aligned - reference) ** 2, axis=1))))


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def clip_score(value: float, limit: float) -> float:
    if not math.isfinite(value) or limit <= 0:
        return 0.0
    return float(np.clip(1.0 - value / limit, 0.0, 1.0))


def main() -> None:
    parser = argparse.ArgumentParser(description="Add V7.2 global-fit RMSD and PDB hashes to a metrics CSV.")
    parser.add_argument("--metrics-csv", required=True)
    parser.add_argument("--reference-pdb", required=True)
    parser.add_argument("--out-csv", required=True)
    parser.add_argument("--barrel-positions", default=DEFAULT_BARREL)
    parser.add_argument("--pocket-positions", default=DEFAULT_POCKET)
    parser.add_argument("--ion-network-positions", default=DEFAULT_ION)
    parser.add_argument("--max-barrel-rmsd", type=float, default=1.50)
    parser.add_argument("--max-pocket-rmsd", type=float, default=0.85)
    parser.add_argument("--max-ion-network-rmsd", type=float, default=0.75)
    args = parser.parse_args()

    metrics = pd.read_csv(args.metrics_csv)
    required = {"sequence_sha1", "pdb_path"}
    missing = required.difference(metrics.columns)
    if missing:
        raise ValueError(f"Metrics CSV is missing required columns: {sorted(missing)}")

    ref_path = Path(args.reference_pdb)
    if not ref_path.is_file():
        raise FileNotFoundError(ref_path)
    ref_atoms = parse_ca_atoms(ref_path)
    if len(ref_atoms) < 3:
        raise ValueError(f"Reference PDB has too few CA atoms: {ref_path}")

    barrel_positions = parse_positions(args.barrel_positions)
    pocket_positions = parse_positions(args.pocket_positions)
    ion_positions = parse_positions(args.ion_network_positions)

    rows: list[dict[str, object]] = []
    for _, row in metrics.iterrows():
        out = row.to_dict()
        path_value = str(row.get("pdb_path", "") or "").strip()
        pdb_path = Path(path_value) if path_value else None
        out.update(
            {
                "pdb_exists_v72": False,
                "pdb_sha256_v72": "",
                "globalfit_status_v72": "missing_pdb",
                "barrel_rmsd_globalfit_v72": math.nan,
                "pocket_rmsd_globalfit_v72": math.nan,
                "ion_network_rmsd_globalfit_v72": math.nan,
                "globalfit_geometry_score_v72": 0.0,
                "globalfit_local_gap_v72": math.nan,
            }
        )
        if pdb_path is None or not pdb_path.is_file():
            rows.append(out)
            continue
        try:
            pred_atoms = parse_ca_atoms(pdb_path)
            common_barrel = [pos for pos in barrel_positions if pos in ref_atoms and pos in pred_atoms]
            if len(common_barrel) < 3:
                out["globalfit_status_v72"] = "insufficient_barrel_atoms"
                out["pdb_exists_v72"] = True
                out["pdb_sha256_v72"] = sha256_file(pdb_path)
                rows.append(out)
                continue

            reference = np.stack([ref_atoms[pos][0] for pos in common_barrel])
            mobile = np.stack([pred_atoms[pos][0] for pos in common_barrel])
            rotation, translation = fit_kabsch(reference, mobile)
            barrel_rmsd = rmsd_after_transform(ref_atoms, pred_atoms, barrel_positions, rotation, translation)
            pocket_rmsd = rmsd_after_transform(ref_atoms, pred_atoms, pocket_positions, rotation, translation)
            ion_rmsd = rmsd_after_transform(ref_atoms, pred_atoms, ion_positions, rotation, translation)
            scores = [
                clip_score(barrel_rmsd, args.max_barrel_rmsd),
                clip_score(pocket_rmsd, args.max_pocket_rmsd),
                clip_score(ion_rmsd, args.max_ion_network_rmsd),
            ]
            local_values = pd.to_numeric(
                pd.Series(
                    [
                        row.get("barrel_rmsd", math.nan),
                        row.get("pocket_rmsd", math.nan),
                        row.get("ion_network_rmsd", math.nan),
                    ]
                ),
                errors="coerce",
            ).to_numpy(dtype=float)
            global_values = np.array([barrel_rmsd, pocket_rmsd, ion_rmsd], dtype=float)
            finite = np.isfinite(local_values) & np.isfinite(global_values)
            gap = float(np.mean(np.maximum(0.0, global_values[finite] - local_values[finite]))) if finite.any() else math.nan
            out.update(
                {
                    "pdb_exists_v72": True,
                    "pdb_sha256_v72": sha256_file(pdb_path),
                    "globalfit_status_v72": "ok",
                    "barrel_rmsd_globalfit_v72": barrel_rmsd,
                    "pocket_rmsd_globalfit_v72": pocket_rmsd,
                    "ion_network_rmsd_globalfit_v72": ion_rmsd,
                    "globalfit_geometry_score_v72": float(np.mean(scores)),
                    "globalfit_local_gap_v72": gap,
                }
            )
        except Exception as exc:
            out["globalfit_status_v72"] = f"error:{type(exc).__name__}"
        rows.append(out)

    out_df = pd.DataFrame(rows)
    Path(args.out_csv).parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(args.out_csv, index=False)
    print(
        f"Wrote V7.2 global-fit metrics: {args.out_csv} | "
        f"rows={len(out_df)} | ok={(out_df['globalfit_status_v72'] == 'ok').sum()}"
    )


if __name__ == "__main__":
    main()
