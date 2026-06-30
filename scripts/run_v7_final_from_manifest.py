#!/usr/bin/env python3
"""V7 final selector from a frozen Stage-1 candidate universe.

V7 keeps the V6.2 evidence-safe contract and adds optional OpenComplex2-style
orthogonal soft validation. Soft validation can improve confidence and ranking,
but it must never replace the primary real-structure + structure_pass gate.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from run_v62_complete import (  # noqa: E402
    add_epistasis_and_ddg_proxies,
    add_v5_gates,
    add_v5_objective,
    add_v6_attribution,
    add_v62_objective,
    add_v62_physical_thermal_proxy,
    diagnostics_frame,
    ensure_jsonable,
    load_config,
    minmax01,
    mutation_frequency,
    safe_num,
    select_portfolio,
    setup_logger,
    structure_geometry_from_metrics,
    write_v5_reports,
    write_v62_survival_reports,
)

AA20_RE = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]+$")


def sequence_sha1(seq: str) -> str:
    return hashlib.sha1(str(seq).strip().upper().encode("utf-8")).hexdigest()


def write_submission_v7(selected: pd.DataFrame, team_name: str, out_dir: Path) -> Path:
    """Write the sole canonical V7 submission artifact used by the final audit."""
    submission = pd.DataFrame(columns=["Team_Name", "Seq_ID", "Sequence"])
    if not selected.empty:
        submission = pd.DataFrame({
            "Team_Name": [team_name] * len(selected),
            "Seq_ID": selected["Seq_ID"].tolist(),
            "Sequence": selected["sequence"].tolist(),
        })
    output_path = out_dir / "submission_v7_for_upload.csv"
    submission.to_csv(output_path, index=False)
    return output_path


def load_table(path: str | Path) -> pd.DataFrame:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(p)
    return pd.read_excel(p) if p.suffix.lower() in {".xlsx", ".xls"} else pd.read_csv(p)


def ensure_keys(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    if "sequence" not in out.columns:
        raise ValueError("Frozen candidate table must contain sequence.")
    out["sequence"] = out["sequence"].astype(str).str.strip().str.upper()
    if "sequence_sha1" not in out.columns:
        out["sequence_sha1"] = out["sequence"].map(sequence_sha1)
    if "candidate_id" not in out.columns:
        out["candidate_id"] = [f"cand_{i:06d}" for i in range(len(out))]
    if "Seq_ID" not in out.columns:
        out["Seq_ID"] = out["candidate_id"]
    return out


def load_exclusion_sequences(data_dir: str | Path | None, exclusion_csv: str | Path | None) -> set[str]:
    paths: list[Path] = []
    if exclusion_csv:
        paths.append(Path(exclusion_csv))
    if data_dir:
        paths.append(Path(data_dir) / "Exclusion_List.csv")
    for p in paths:
        if p.exists():
            df = load_table(p)
            col = "Sequence" if "Sequence" in df.columns else ("sequence" if "sequence" in df.columns else None)
            if col is None:
                continue
            return set(df[col].astype(str).str.strip().str.upper())
    return set()


def merge_structure_metrics(cand: pd.DataFrame, metrics_csv: str | Path) -> pd.DataFrame:
    out = cand.copy()
    metrics = ensure_keys(load_table(metrics_csv))
    keep = [c for c in [
        "sequence_sha1", "candidate_id", "Seq_ID", "sequence", "barrel_rmsd", "pocket_rmsd",
        "ion_network_rmsd", "mean_plddt", "pdb_path", "has_real_structure", "structure_status", "structure_pass"
    ] if c in metrics.columns]
    metrics = metrics[keep].copy()
    # Prefer exact-sequence rows with explicit pass evidence and complete RMSD values.
    metrics["_pass_rank"] = metrics.get("structure_pass", pd.Series(False, index=metrics.index)).astype(str).str.lower().eq("true").astype(int)
    metrics["_rmsd_rank"] = metrics[[c for c in ["barrel_rmsd", "pocket_rmsd", "ion_network_rmsd"] if c in metrics.columns]].notna().sum(axis=1)
    metrics = metrics.sort_values(["sequence_sha1", "_pass_rank", "_rmsd_rank"], ascending=[True, False, False]).drop_duplicates("sequence_sha1")
    drop_cols = [c for c in ["barrel_rmsd", "pocket_rmsd", "ion_network_rmsd", "mean_plddt", "pdb_path", "has_real_structure", "structure_status", "structure_pass"] if c in out.columns]
    out = out.drop(columns=drop_cols, errors="ignore").merge(metrics.drop(columns=["candidate_id", "Seq_ID", "sequence", "_pass_rank", "_rmsd_rank"], errors="ignore"), on="sequence_sha1", how="left")
    out["structure_status"] = out.get("structure_status", pd.Series("not_checked", index=out.index)).fillna("not_checked")

    rmsd_cols = ["barrel_rmsd", "pocket_rmsd", "ion_network_rmsd"]
    for col in rmsd_cols:
        out[col] = pd.to_numeric(out.get(col, pd.Series(np.nan, index=out.index)), errors="coerce")

    status_real = out["structure_status"].astype(str).eq("rmsd_checked")
    rmsd_complete = out[rmsd_cols].notna().all(axis=1)
    imported_real = out.get("has_real_structure", pd.Series(False, index=out.index)).astype(str).str.lower().eq("true")
    # This is the hard evidence gate: confidence proxies cannot become final real structure evidence.
    out["has_real_structure"] = status_real & rmsd_complete & imported_real
    out["structure_pass"] = out.get("structure_pass", pd.Series(False, index=out.index)).astype(str).str.lower().eq("true") & out["has_real_structure"].astype(bool)
    out["structure_geometry_score"] = structure_geometry_from_metrics(out)
    out.loc[~out["has_real_structure"].astype(bool), "structure_geometry_score"] = 0.05
    return out


def merge_proteinmpnn(cand: pd.DataFrame, proteinmpnn_csv: str | Path | None) -> pd.DataFrame:
    out = cand.copy()
    out["proteinmpnn_status"] = out.get("proteinmpnn_status", "not_provided")
    out["proteinmpnn_score_raw"] = pd.to_numeric(out.get("proteinmpnn_score_raw", pd.Series(np.nan, index=out.index)), errors="coerce")
    out["proteinmpnn_score_norm"] = pd.to_numeric(out.get("proteinmpnn_score_norm", pd.Series(0.5, index=out.index)), errors="coerce").fillna(0.5)
    if not proteinmpnn_csv or not Path(proteinmpnn_csv).exists():
        return out
    pm = load_table(proteinmpnn_csv)
    if pm.empty:
        return out
    if "sequence" in pm.columns and "sequence_sha1" not in pm.columns:
        pm["sequence_sha1"] = pm["sequence"].astype(str).str.strip().str.upper().map(sequence_sha1)
    score_col = next((c for c in ["proteinmpnn_score", "proteinmpnn_score_raw", "score"] if c in pm.columns), None)
    if score_col is None:
        return out
    if "sequence_sha1" in pm.columns:
        key = "sequence_sha1"
    else:
        key = next((c for c in ["candidate_id", "Seq_ID", "seq_id"] if c in pm.columns and c in out.columns), None)
    if key is None:
        return out
    small = pm[[key, score_col]].drop_duplicates(key).rename(columns={score_col: "proteinmpnn_score_raw_imported"})
    out = out.merge(small, on=key, how="left")
    imported = pd.to_numeric(out["proteinmpnn_score_raw_imported"], errors="coerce")
    out["proteinmpnn_score_raw"] = imported.combine_first(out["proteinmpnn_score_raw"])
    out["proteinmpnn_score_norm"] = minmax01(out["proteinmpnn_score_raw"], higher_is_better=False, default=0.5)
    out["proteinmpnn_status"] = np.where(out["proteinmpnn_score_raw"].notna(), "score_imported", out["proteinmpnn_status"])
    out = out.drop(columns=["proteinmpnn_score_raw_imported"], errors="ignore")
    return out


def merge_softcheck(cand: pd.DataFrame, softcheck_csv: str | Path | None) -> pd.DataFrame:
    out = cand.copy()
    out["softcheck_available"] = False
    out["oc2_softcheck_score"] = 0.50
    out["oc2_pass"] = False
    out["oc2_status"] = "not_provided"
    if not softcheck_csv or not Path(softcheck_csv).exists():
        return out
    oc2 = load_table(softcheck_csv)
    if oc2.empty:
        return out
    if "sequence" in oc2.columns and "sequence_sha1" not in oc2.columns:
        oc2["sequence_sha1"] = oc2["sequence"].astype(str).str.strip().str.upper().map(sequence_sha1)
    if "sequence_sha1" not in oc2.columns:
        return out
    rename = {
        "barrel_rmsd": "oc2_barrel_rmsd",
        "pocket_rmsd": "oc2_pocket_rmsd",
        "ion_network_rmsd": "oc2_ion_network_rmsd",
        "mean_plddt": "oc2_mean_plddt",
        "pdb_path": "oc2_pdb_path",
        "structure_status": "oc2_status",
        "structure_pass": "oc2_pass",
        "softcheck_score": "oc2_softcheck_score",
    }
    oc2 = oc2.rename(columns={k: v for k, v in rename.items() if k in oc2.columns and v not in oc2.columns})
    keep = [c for c in [
        "sequence_sha1", "oc2_status", "oc2_pdb_path", "oc2_barrel_rmsd", "oc2_pocket_rmsd",
        "oc2_ion_network_rmsd", "oc2_mean_plddt", "oc2_softcheck_score", "oc2_pass", "oc2_reason"
    ] if c in oc2.columns]
    oc2 = oc2[keep].copy().drop_duplicates("sequence_sha1")
    # Replace placeholder softcheck columns before merging imported OC2/OpenComplex2 evidence.
    out = out.drop(columns=[c for c in oc2.columns if c != "sequence_sha1" and c in out.columns], errors="ignore")
    out = out.merge(oc2, on="sequence_sha1", how="left")
    out["softcheck_available"] = out.get("oc2_status", pd.Series("", index=out.index)).astype(str).eq("rmsd_checked")
    out["oc2_softcheck_score"] = pd.to_numeric(out.get("oc2_softcheck_score", pd.Series(0.5, index=out.index)), errors="coerce").fillna(0.5).clip(0.0, 1.0)
    out["oc2_pass"] = out.get("oc2_pass", pd.Series(False, index=out.index)).astype(str).str.lower().eq("true")
    for col in ["oc2_barrel_rmsd", "oc2_pocket_rmsd", "oc2_ion_network_rmsd", "oc2_mean_plddt"]:
        out[col] = pd.to_numeric(out.get(col, pd.Series(np.nan, index=out.index)), errors="coerce")
    return out



HYDROPHOBIC_AA = set("AVILMFWY")
AROMATIC_AA = set("FWY")
ACIDIC_AA = set("DE")
BASIC_AA = set("KR")
FRAGILE_ALT_AA = set("CPGW")


def _clip01(x: pd.Series | np.ndarray | float) -> pd.Series:
    """Clamp numeric values into [0, 1] while preserving row index when possible."""
    return pd.Series(x).clip(0.0, 1.0)


def _seq_net_charge(seq: str) -> float:
    """Approximate pH-7 net charge for sequence-only solubility screening."""
    seq = str(seq or "").upper()
    return float(seq.count("K") + seq.count("R") + 0.10 * seq.count("H") - seq.count("D") - seq.count("E"))


def _longest_run(seq: str, alphabet: set[str]) -> int:
    best = cur = 0
    for aa in str(seq or "").upper():
        if aa in alphabet:
            cur += 1
            best = max(best, cur)
        else:
            cur = 0
    return best


def _reference_like_sequence(df: pd.DataFrame) -> str:
    """Pick the lowest-mutation sequence as an sfGFP-like baseline for charge-shift estimates."""
    if df.empty or "sequence" not in df.columns:
        return ""
    if "mutation_count" in df.columns:
        idx = pd.to_numeric(df["mutation_count"], errors="coerce").fillna(9999).idxmin()
        return str(df.loc[idx, "sequence"]).strip().upper()
    return str(df.iloc[0]["sequence"]).strip().upper()


def _mutation_fragility_components(muts: pd.Series) -> pd.DataFrame:
    rows: list[dict[str, float]] = []
    for text in muts.astype(str):
        parsed = re.findall(r"([ACDEFGHIKLMNPQRSTVWY])(\d{1,3})([ACDEFGHIKLMNPQRSTVWY])", text.upper())
        denom = max(1, len(parsed))
        acidic_gain = sum(1 for _wt, _pos, alt in parsed if alt in ACIDIC_AA) / denom
        hydrophobic_alt = sum(1 for _wt, _pos, alt in parsed if alt in HYDROPHOBIC_AA) / denom
        aromatic_or_cys_alt = sum(1 for _wt, _pos, alt in parsed if alt in AROMATIC_AA or alt == "C") / denom
        fragile_alt = sum(1 for _wt, _pos, alt in parsed if alt in FRAGILE_ALT_AA) / denom
        charge_flip = sum(1 for wt, _pos, alt in parsed if (wt in ACIDIC_AA and alt in BASIC_AA) or (wt in BASIC_AA and alt in ACIDIC_AA)) / denom
        rows.append({
            "v7_mut_acidic_gain_fraction": acidic_gain,
            "v7_mut_hydrophobic_alt_fraction": hydrophobic_alt,
            "v7_mut_aromatic_or_cys_alt_fraction": aromatic_or_cys_alt,
            "v7_mut_fragile_alt_fraction": fragile_alt,
            "v7_mut_charge_flip_fraction": charge_flip,
        })
    return pd.DataFrame(rows, index=muts.index)


def _add_v7_mutation_neighborhood(df: pd.DataFrame, cfg: dict[str, Any]) -> pd.DataFrame:
    out = df.copy()
    neigh_cfg = cfg.get("v7", {}).get("neighborhood", {})
    if not bool(neigh_cfg.get("enabled", True)) or "mutations" not in out.columns:
        out["v7_neighborhood_support"] = 0.5
        out["v7_neighborhood_fragility"] = 0.5
        return out
    max_positions = int(neigh_cfg.get("max_positions_for_signature", 8))
    signatures: list[tuple[int, ...]] = []
    reduced_counts: dict[tuple[int, ...], int] = {}
    for text in out["mutations"].astype(str):
        positions = sorted({int(pos) for _wt, pos, _alt in re.findall(r"([ACDEFGHIKLMNPQRSTVWY])(\d{1,3})([ACDEFGHIKLMNPQRSTVWY])", text.upper())})
        if len(positions) > max_positions:
            positions = positions[:max_positions]
        sig = tuple(positions)
        signatures.append(sig)
        reduced_counts[sig] = reduced_counts.get(sig, 0) + 1
        # Removing one position creates a cheap local-neighborhood index without O(N^2) Hamming scans.
        for i in range(len(sig)):
            rsig = sig[:i] + sig[i+1:]
            reduced_counts[rsig] = reduced_counts.get(rsig, 0) + 1
    supports = []
    for sig in signatures:
        max_count = reduced_counts.get(sig, 0)
        for i in range(len(sig)):
            max_count = max(max_count, reduced_counts.get(sig[:i] + sig[i+1:], 0))
        supports.append(min(1.0, np.log1p(max_count) / np.log1p(12.0)))
    out["v7_neighborhood_support"] = pd.Series(supports, index=out.index).fillna(0.5).clip(0.0, 1.0)
    mut_count = safe_num(out, "mutation_count", 0.0)
    risk = safe_num(out, "low_brightness_risk", 0.5).clip(0.0, 1.0)
    # Unique high-mutation candidates are not rejected, but they receive an uncertainty penalty.
    out["v7_neighborhood_fragility"] = (0.55 * (1.0 - out["v7_neighborhood_support"]) + 0.30 * risk + 0.15 * np.maximum(0.0, mut_count - 4.0) / 4.0).clip(0.0, 1.0)
    return out


def add_v7_biology_soft_risk(df: pd.DataFrame, cfg: dict[str, Any]) -> pd.DataFrame:
    """Add sequence-level biology priors that complement structure metrics.

    These features are deliberately soft. They estimate expression/aggregation/thermal-retention
    risk from sequence composition and mutation pattern, because no wet-lab feedback is available.
    They must not override the hard ColabFold/PDB structure gate.
    """
    out = df.copy()
    bio_cfg = cfg.get("v7", {}).get("biology_soft_risk", {})
    if not bool(bio_cfg.get("enabled", True)):
        out["v7_biology_validity_score"] = 0.5
        out["v7_soft_risk_penalty"] = 0.0
        return out

    seqs = out["sequence"].astype(str).str.strip().str.upper()
    ref_seq = str(bio_cfg.get("reference_sequence") or _reference_like_sequence(out))
    ref_charge = _seq_net_charge(ref_seq) if ref_seq else 0.0
    lengths = seqs.str.len().replace(0, np.nan)
    aa_counts = {aa: seqs.str.count(aa) for aa in "ACDEFGHIKLMNPQRSTVWY"}
    hydrophobic_frac = sum(aa_counts[aa] for aa in HYDROPHOBIC_AA) / lengths
    aromatic_frac = sum(aa_counts[aa] for aa in AROMATIC_AA) / lengths
    acidic_frac = sum(aa_counts[aa] for aa in ACIDIC_AA) / lengths
    basic_frac = sum(aa_counts[aa] for aa in BASIC_AA) / lengths
    cys_count = aa_counts["C"]
    pro_gly_frac = (aa_counts["P"] + aa_counts["G"]) / lengths
    charges = seqs.map(_seq_net_charge)
    charge_shift = (charges - ref_charge).abs()
    hydrophobic_run = seqs.map(lambda s: _longest_run(s, HYDROPHOBIC_AA))
    repeat_run = seqs.map(lambda s: max((_longest_run(s, {aa}) for aa in "ACDEFGHIKLMNPQRSTVWY"), default=0))

    mut_components = _mutation_fragility_components(out.get("mutations", pd.Series("", index=out.index)))
    out = pd.concat([out, mut_components], axis=1)
    out["v7_net_charge_pH7_approx"] = charges
    out["v7_charge_shift_from_reference"] = charge_shift
    out["v7_hydrophobic_fraction"] = hydrophobic_frac.fillna(0.0)
    out["v7_aromatic_fraction"] = aromatic_frac.fillna(0.0)
    out["v7_acidic_fraction"] = acidic_frac.fillna(0.0)
    out["v7_basic_fraction"] = basic_frac.fillna(0.0)
    out["v7_cys_count"] = cys_count
    out["v7_pro_gly_fraction"] = pro_gly_frac.fillna(0.0)
    out["v7_longest_hydrophobic_run"] = hydrophobic_run
    out["v7_longest_repeat_run"] = repeat_run

    charge_shift_penalty = (charge_shift / float(bio_cfg.get("charge_shift_soft_limit", 8.0))).clip(0.0, 1.0)
    hydrophobic_penalty = ((hydrophobic_frac - float(bio_cfg.get("hydrophobic_fraction_center", 0.39))).abs() / 0.16).clip(0.0, 1.0).fillna(0.5)
    aromatic_penalty = (aromatic_frac / float(bio_cfg.get("aromatic_fraction_soft_limit", 0.12))).clip(0.0, 1.0).fillna(0.0)
    cys_penalty = (cys_count / float(bio_cfg.get("cys_soft_limit", 2.0))).clip(0.0, 1.0)
    run_penalty = ((hydrophobic_run - 4).clip(lower=0) / 4.0).clip(0.0, 1.0)
    repeat_penalty = ((repeat_run - 4).clip(lower=0) / 4.0).clip(0.0, 1.0)
    pro_gly_penalty = ((pro_gly_frac - 0.16).clip(lower=0) / 0.12).clip(0.0, 1.0).fillna(0.0)

    out["v7_sequence_aggregation_risk"] = (
        0.22 * hydrophobic_penalty
        + 0.16 * aromatic_penalty
        + 0.18 * cys_penalty
        + 0.16 * run_penalty
        + 0.10 * repeat_penalty
        + 0.10 * pro_gly_penalty
        + 0.08 * out["v7_mut_hydrophobic_alt_fraction"]
    ).clip(0.0, 1.0)
    # A mildly negative/polar surface is useful, but excessive charge drift from the sfGFP-like baseline is risky.
    out["v7_solubility_score"] = (
        0.40 * (1.0 - out["v7_sequence_aggregation_risk"])
        + 0.25 * (1.0 - charge_shift_penalty)
        + 0.15 * (acidic_frac.fillna(0.0) / 0.12).clip(0.0, 1.0)
        + 0.10 * safe_num(out, "surface_charge_score", 0.5).clip(0.0, 1.0)
        + 0.10 * safe_num(out, "aggregation_safety_score", 0.5).clip(0.0, 1.0)
    ).clip(0.0, 1.0)
    mut_count = safe_num(out, "mutation_count", 0.0)
    out["v7_sequence_fragility_score"] = (
        0.30 * safe_num(out, "low_brightness_risk", 0.5).clip(0.0, 1.0)
        + 0.20 * safe_num(out, "epistasis_penalty", 0.0).clip(0.0, 1.0)
        + 0.15 * np.maximum(0.0, mut_count - 4.0) / 5.0
        + 0.15 * out["v7_mut_fragile_alt_fraction"]
        + 0.10 * out["v7_mut_charge_flip_fraction"]
        + 0.10 * out["v7_sequence_aggregation_risk"]
    ).clip(0.0, 1.0)
    out = _add_v7_mutation_neighborhood(out, cfg)
    risk_inconsistency = (safe_num(out, "brightness_lcb_norm", 0.5).clip(0.0, 1.0) - safe_num(out, "low_risk_score", 0.5).clip(0.0, 1.0)).abs().clip(0.0, 1.0)
    out["v7_surrogate_disagreement_penalty"] = (0.60 * risk_inconsistency + 0.40 * out["v7_neighborhood_fragility"]).clip(0.0, 1.0)
    out["v7_thermal_validity_score"] = (
        0.34 * safe_num(out, "thermal_proxy_v62", 0.5).clip(0.0, 1.0)
        + 0.20 * out["v7_solubility_score"]
        + 0.16 * (1.0 - out["v7_sequence_aggregation_risk"])
        + 0.14 * safe_num(out, "structure_geometry_score", 0.05).clip(0.0, 1.0)
        + 0.10 * safe_num(out, "proteinmpnn_score_norm", 0.5).clip(0.0, 1.0)
        + 0.06 * (1.0 - out["v7_sequence_fragility_score"])
    ).clip(0.0, 1.0)
    out["v7_biology_validity_score"] = (
        0.45 * out["v7_thermal_validity_score"]
        + 0.25 * (1.0 - out["v7_sequence_fragility_score"])
        + 0.20 * out["v7_solubility_score"]
        + 0.10 * out["v7_neighborhood_support"]
    ).clip(0.0, 1.0)
    out["v7_soft_risk_penalty"] = (
        0.36 * out["v7_sequence_fragility_score"]
        + 0.28 * out["v7_sequence_aggregation_risk"]
        + 0.20 * out["v7_surrogate_disagreement_penalty"]
        + 0.16 * charge_shift_penalty
    ).clip(0.0, 1.0)
    return out

def add_v7_confidence_and_objective(df: pd.DataFrame, cfg: dict[str, Any]) -> pd.DataFrame:
    out = df.copy()
    soft_cfg = cfg.get("v7", {}).get("softcheck", {})
    bio_cfg = cfg.get("v7", {}).get("biology_soft_risk", {})
    soft_weight = float(soft_cfg.get("objective_weight", 0.12))
    penalty_weight = float(soft_cfg.get("disagreement_penalty_weight", 0.06))
    biology_weight = float(bio_cfg.get("objective_weight", 0.14))
    biology_penalty_weight = float(bio_cfg.get("risk_penalty_weight", 0.07))
    available = out.get("softcheck_available", pd.Series(False, index=out.index)).astype(bool)
    soft_score = safe_num(out, "oc2_softcheck_score", 0.5).clip(0.0, 1.0)
    main_score = safe_num(out, "v62_structure_score", 0.5).clip(0.0, 1.0)
    pmpnn = safe_num(out, "proteinmpnn_score_norm", 0.5).clip(0.0, 1.0)
    low_risk = safe_num(out, "low_risk_score", 0.5).clip(0.0, 1.0)
    brightness = safe_num(out, "brightness_lcb_norm", 0.5).clip(0.0, 1.0)
    biology_score = safe_num(out, "v7_biology_validity_score", 0.5).clip(0.0, 1.0)
    biology_penalty = safe_num(out, "v7_soft_risk_penalty", 0.0).clip(0.0, 1.0)

    # Softcheck disagreement is intentionally a small penalty: orthogonal models are noisy, but consistent failure is a warning.
    main_pass = out.get("structure_pass", pd.Series(False, index=out.index)).astype(bool)
    soft_pass = out.get("oc2_pass", pd.Series(False, index=out.index)).astype(bool)
    rmsd_delta = (
        (safe_num(out, "oc2_barrel_rmsd", np.nan) - safe_num(out, "barrel_rmsd", np.nan)).abs().fillna(0.0) / 2.5
        + (safe_num(out, "oc2_pocket_rmsd", np.nan) - safe_num(out, "pocket_rmsd", np.nan)).abs().fillna(0.0) / 1.2
        + (safe_num(out, "oc2_ion_network_rmsd", np.nan) - safe_num(out, "ion_network_rmsd", np.nan)).abs().fillna(0.0) / 1.2
    ) / 3.0
    out["v7_softcheck_disagreement_penalty"] = np.where(available & main_pass & ~soft_pass, 1.0, 0.0)
    out["v7_softcheck_disagreement_penalty"] = np.maximum(out["v7_softcheck_disagreement_penalty"], np.where(available, rmsd_delta.clip(0.0, 1.0), 0.0))
    out["v7_orthogonal_structure_score"] = np.where(available, 0.65 * main_score + 0.35 * soft_score, main_score)
    out["v7_surrogate_confidence_score"] = (
        0.30 * out["v7_orthogonal_structure_score"]
        + 0.20 * biology_score
        + 0.16 * pmpnn
        + 0.14 * safe_num(out, "v7_thermal_validity_score", 0.5).clip(0.0, 1.0)
        + 0.12 * low_risk
        + 0.08 * brightness
    ).clip(0.0, 1.0)
    out["v7_confidence_score"] = out["v7_surrogate_confidence_score"]
    base = safe_num(out, "v62_objective", 0.0)
    out["v7_objective"] = (
        (1.0 - soft_weight - biology_weight) * base
        + soft_weight * out["v7_orthogonal_structure_score"]
        + biology_weight * biology_score
        - penalty_weight * out["v7_softcheck_disagreement_penalty"]
        - biology_penalty_weight * biology_penalty
        - 0.10 * out.get("in_exclusion_list", pd.Series(False, index=out.index)).astype(bool).astype(float)
    ).clip(0.0, 1.0)
    if bool(soft_cfg.get("require_for_final", False)):
        out["v5_final_eligible"] = out.get("v5_final_eligible", pd.Series(False, index=out.index)).astype(bool) & available & soft_pass
    # Downstream portfolio selector sorts v6_objective, so alias V7 objective after all V6.2 computations are complete.
    out["v6_objective"] = out["v7_objective"]
    out["v5_objective"] = out["v7_objective"]
    return out

def write_underfilled_recovery(ranked: pd.DataFrame, out_dir: Path, k: int) -> None:
    need = ranked.copy()
    need["needs_structure_recovery"] = ~need.get("structure_pass", pd.Series(False, index=need.index)).astype(bool)
    need = need[need["needs_structure_recovery"] & need.get("v5_basic_gate", pd.Series(False, index=need.index)).astype(bool)]
    need = need.sort_values("v7_objective" if "v7_objective" in need.columns else "v6_objective", ascending=False).head(max(40, k * 20))
    cols = [c for c in ["candidate_id", "Seq_ID", "sequence", "sequence_sha1", "source_branch", "mutations", "v7_objective", "v6_objective", "low_brightness_risk", "mutation_count", "structure_status", "structure_pass", "softcheck_available", "oc2_pass"] if c in need.columns]
    need[cols].to_csv(out_dir / "final_underfilled_needs_recovery.csv", index=False)


def write_v7_summary(selected: pd.DataFrame, ranked: pd.DataFrame, args: argparse.Namespace, out_dir: Path) -> None:
    summary = {
        "frozen_candidates": args.frozen_candidates_csv,
        "structure_metrics": args.structure_metrics_csv,
        "proteinmpnn_scores": args.proteinmpnn_score_csv,
        "softcheck_scores": args.softcheck_csv,
        "ranked_rows": int(len(ranked)),
        "final_rows": int(len(selected)),
        "final_all_real_structure": bool(len(selected) > 0 and selected.get("has_real_structure", pd.Series(False)).astype(bool).all()),
        "final_all_structure_pass": bool(len(selected) > 0 and selected.get("structure_pass", pd.Series(False)).astype(bool).all()),
        "final_softcheck_available": int(selected.get("softcheck_available", pd.Series(False, index=selected.index)).astype(bool).sum()) if not selected.empty else 0,
        "final_softcheck_pass": int(selected.get("oc2_pass", pd.Series(False, index=selected.index)).astype(bool).sum()) if not selected.empty else 0,
        "ranked_softcheck_available": int(ranked.get("softcheck_available", pd.Series(False, index=ranked.index)).astype(bool).sum()),
        "final_mean_biology_validity": float(safe_num(selected, "v7_biology_validity_score", 0.0).mean()) if not selected.empty else 0.0,
        "final_min_thermal_validity": float(safe_num(selected, "v7_thermal_validity_score", 0.0).min()) if not selected.empty else 0.0,
        "final_max_sequence_aggregation_risk": float(safe_num(selected, "v7_sequence_aggregation_risk", 0.0).max()) if not selected.empty else 0.0,
        "final_max_sequence_fragility": float(safe_num(selected, "v7_sequence_fragility_score", 0.0).max()) if not selected.empty else 0.0,
        "final_min_neighborhood_support": float(safe_num(selected, "v7_neighborhood_support", 0.0).min()) if not selected.empty else 0.0,
    }
    (out_dir / "v7_final_summary.json").write_text(json.dumps(ensure_jsonable(summary), indent=2, ensure_ascii=False), encoding="utf-8")
    # Keep a V6.2-compatible summary name for existing downstream notebooks.
    (out_dir / "v62_final_summary.json").write_text(json.dumps(ensure_jsonable(summary), indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(summary, indent=2, ensure_ascii=False))


def main() -> None:
    p = argparse.ArgumentParser(description="Select V7 final Top-6 from frozen Stage-1 candidates with optional soft validation.")
    p.add_argument("--config", default=str(ROOT / "configs" / "v7_adaptive.yaml"))
    p.add_argument("--frozen-candidates-csv", required=True)
    p.add_argument("--structure-metrics-csv", required=True)
    p.add_argument("--proteinmpnn-score-csv", default="")
    p.add_argument("--softcheck-csv", default="")
    p.add_argument("--data-dir", default="")
    p.add_argument("--exclusion-list-csv", default="")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--team-name", default="YourTeamName")
    args = p.parse_args()

    cfg = load_config(args.config)
    cfg.setdefault("v5", {})["require_real_structure_for_final"] = True
    cfg.setdefault("paths", {})["out_dir"] = args.out_dir
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    logger = setup_logger(out_dir, name="synbio_gfp_v7_final")

    cand = ensure_keys(load_table(args.frozen_candidates_csv))
    exclusion = load_exclusion_sequences(args.data_dir, args.exclusion_list_csv)
    cand["in_exclusion_list"] = cand["sequence"].isin(exclusion) if exclusion else False
    cand = merge_structure_metrics(cand, args.structure_metrics_csv)
    cand = merge_proteinmpnn(cand, args.proteinmpnn_score_csv)
    cand = merge_softcheck(cand, args.softcheck_csv)
    cand = add_epistasis_and_ddg_proxies(cand) if "epistasis_penalty" not in cand.columns else cand
    cand = add_v5_gates(cand, cfg)
    cand.loc[cand["in_exclusion_list"].astype(bool), ["v5_basic_gate", "v5_final_eligible"]] = False
    cand = add_v5_objective(cand)
    cand = add_v62_physical_thermal_proxy(cand)
    cand = add_v62_objective(cand)
    cand = add_v7_biology_soft_risk(cand, cfg)
    cand = add_v7_confidence_and_objective(cand, cfg)
    cand = add_v6_attribution(cand)
    ranked = cand.sort_values("v7_objective", ascending=False).reset_index(drop=True)
    ranked.to_csv(out_dir / "all_ranked_candidates_v7.csv", index=False)
    ranked.to_csv(out_dir / "all_ranked_candidates_v62.csv", index=False)
    ranked.to_csv(out_dir / "all_ranked_candidates_v6.csv", index=False)
    risk_cols = [c for c in [
        "candidate_id", "Seq_ID", "sequence_sha1", "source_branch", "mutations", "v7_objective",
        "v7_confidence_score", "v7_biology_validity_score", "v7_thermal_validity_score",
        "v7_solubility_score", "v7_sequence_aggregation_risk", "v7_sequence_fragility_score",
        "v7_neighborhood_support", "v7_soft_risk_penalty", "low_brightness_risk", "mutation_count",
        "structure_pass", "softcheck_available", "oc2_pass",
    ] if c in ranked.columns]
    ranked[risk_cols].head(500).to_csv(out_dir / "v7_soft_risk_top500.csv", index=False)

    selected = select_portfolio(ranked, cfg)
    selected.to_csv(out_dir / "final_top6_v7.csv", index=False)
    selected[[c for c in risk_cols if c in selected.columns]].to_csv(out_dir / "v7_final_soft_risk_report.csv", index=False)
    selected.to_csv(out_dir / "final_top6_v62.csv", index=False)
    selected.to_csv(out_dir / "final_top6_v61.csv", index=False)
    selected.to_csv(out_dir / "final_top6_v6.csv", index=False)
    write_submission_v7(selected, args.team_name, out_dir)

    diagnostics = diagnostics_frame(ranked)
    diagnostics["softcheck_available_pass"] = int(ranked.get("softcheck_available", pd.Series(False, index=ranked.index)).astype(bool).sum())
    diagnostics["softcheck_pass_true"] = int(ranked.get("oc2_pass", pd.Series(False, index=ranked.index)).astype(bool).sum())
    diagnostics["in_exclusion_list"] = int(ranked.get("in_exclusion_list", pd.Series(False, index=ranked.index)).astype(bool).sum())
    diagnostics.to_csv(out_dir / "v7_gate_diagnostics.csv", index=False)
    diagnostics.to_csv(out_dir / "v6_gate_diagnostics.csv", index=False)
    mutation_frequency(ranked).to_csv(out_dir / "mutation_frequency_by_seed_v6.csv", index=False)
    write_v62_survival_reports(ranked, out_dir)
    write_underfilled_recovery(ranked, out_dir, int(cfg.get("v5", {}).get("final_k", 6)))
    priority = ranked.head(int(cfg.get("v5", {}).get("structure_top_n", 200))).copy()
    write_v5_reports(selected, ranked, priority, priority, diagnostics, {"stage4_from_frozen": True, "v7_softcheck": bool(args.softcheck_csv)}, cfg, out_dir)
    write_v7_summary(selected, ranked, args, out_dir)
    logger.info("V7 frozen final done | final=%d", len(selected))


if __name__ == "__main__":
    main()
