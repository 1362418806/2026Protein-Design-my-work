from __future__ import annotations

import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

MUT_RE = re.compile(r"([ACDEFGHIKLMNPQRSTVWY])(\d{1,3})([ACDEFGHIKLMNPQRSTVWY])")


def parse_mutation_text(text: Any) -> list[str]:
    if not isinstance(text, str):
        return []
    return [f"{a}{int(pos)}{b}" for a, pos, b in MUT_RE.findall(text.upper())]


def mutation_position(mut: str) -> int | None:
    m = MUT_RE.fullmatch(str(mut))
    return int(m.group(2)) if m else None


def _safe_num(df: pd.DataFrame, col: str, default: float = 0.0) -> pd.Series:
    if col not in df.columns:
        return pd.Series(default, index=df.index, dtype=float)
    return pd.to_numeric(df[col], errors="coerce").fillna(default)


def _score_column(df: pd.DataFrame) -> str | None:
    for col in ["v6_objective", "v5_objective", "objective", "brightness_lcb_norm"]:
        if col in df.columns:
            return col
    return None


def _normalize_counter(counter: Counter[str] | Counter[int]) -> dict[str, float]:
    if not counter:
        return {}
    max_v = max(float(v) for v in counter.values())
    if max_v <= 0:
        return {}
    return {str(k): float(v) / max_v for k, v in counter.items()}


def _load_table(path: str | None) -> pd.DataFrame:
    if not path:
        return pd.DataFrame()
    p = Path(path)
    if not p.exists():
        return pd.DataFrame()
    return pd.read_excel(p) if p.suffix.lower() in {".xlsx", ".xls"} else pd.read_csv(p)


def build_feedback_policy_from_table(df: pd.DataFrame, cfg: dict[str, Any]) -> dict[str, Any]:
    """Convert ranked candidates with real metrics into causal mutation feedback.

    The policy is deliberately low-dimensional and auditable: branch success rates,
    mutation log-odds, position log-odds, and explicit blacklists.  V6 uses this
    policy both to bias new candidate generation and to score candidates before
    expensive structure prediction.
    """
    if df is None or df.empty:
        return {"enabled": False, "reason": "empty_feedback_table"}

    work = df.copy()
    score_col = _score_column(work)
    if score_col is None:
        work["_score"] = 0.0
        score_col = "_score"
    work["_score_norm"] = _safe_num(work, score_col, 0.0)
    if work["_score_norm"].max() > work["_score_norm"].min():
        work["_score_norm"] = (work["_score_norm"] - work["_score_norm"].min()) / (work["_score_norm"].max() - work["_score_norm"].min())
    else:
        work["_score_norm"] = 0.5

    has_real = work.get("has_real_structure", pd.Series(False, index=work.index)).astype(bool)
    structure_pass = work.get("structure_pass", pd.Series(False, index=work.index)).astype(bool)
    final_eligible = work.get("v5_final_eligible", pd.Series(False, index=work.index)).astype(bool)
    low_risk = _safe_num(work, "low_brightness_risk", 1.0) <= float(cfg.get("v6", {}).get("feedback", {}).get("max_success_low_risk", 0.35))
    pmpnn = _safe_num(work, "proteinmpnn_score_norm", 0.5)
    pmpnn_cut = float(cfg.get("v6", {}).get("feedback", {}).get("min_success_pmpnn_norm", 0.35))
    success = final_eligible | (has_real & structure_pass & low_risk & (pmpnn >= pmpnn_cut))

    # A failure is a measured structure or brightness-risk failure, not merely a proxy miss.
    barrel_fail = _safe_num(work, "barrel_rmsd", 0.0) > float(cfg.get("v5", {}).get("structure_gate", {}).get("max_barrel_rmsd", 1.5))
    pocket_fail = _safe_num(work, "pocket_rmsd", 0.0) > float(cfg.get("v5", {}).get("structure_gate", {}).get("max_pocket_rmsd", 0.85))
    ion_fail = _safe_num(work, "ion_network_rmsd", 0.0) > float(cfg.get("v5", {}).get("structure_gate", {}).get("max_ion_network_rmsd", 0.75))
    risk_fail = _safe_num(work, "low_brightness_risk", 0.0) > float(cfg.get("v5", {}).get("basic_gate", {}).get("max_low_brightness_risk", 0.35))
    failure = (has_real & (~structure_pass | barrel_fail | pocket_fail | ion_fail)) | risk_fail

    good_rows = work[success].sort_values(score_col, ascending=False)
    bad_rows = work[failure].copy()
    if good_rows.empty:
        good_rows = work.sort_values(score_col, ascending=False).head(int(cfg.get("v6", {}).get("feedback", {}).get("fallback_good_rows", 80)))

    pos_good: Counter[int] = Counter()
    pos_bad: Counter[int] = Counter()
    mut_good: Counter[str] = Counter()
    mut_bad: Counter[str] = Counter()
    branch_good: Counter[str] = Counter()
    branch_total: Counter[str] = Counter()
    seed_rows: list[dict[str, Any]] = []

    for _, row in work.iterrows():
        branch = str(row.get("source_branch", row.get("source", "unknown")))
        branch_total[branch] += 1
    for _, row in good_rows.iterrows():
        branch = str(row.get("source_branch", row.get("source", "unknown")))
        branch_good[branch] += 1
        seed_rows.append({
            "candidate_id": row.get("candidate_id", ""),
            "sequence": row.get("sequence", ""),
            "mutations": row.get("mutations", row.get("mutation_list", "")),
            "source_branch": branch,
            "score": float(row.get(score_col, 0.0) or 0.0),
        })
        weight = 1.0 + float(row.get("_score_norm", 0.5) or 0.5)
        for mut in parse_mutation_text(row.get("mutations", row.get("mutation_list", ""))):
            mut_good[mut] += weight
            pos = mutation_position(mut)
            if pos is not None:
                pos_good[pos] += weight

    for _, row in bad_rows.iterrows():
        # Geometry-specific failures get stronger negative weight because they indicate structural fragility.
        weight = 1.0
        if bool(barrel_fail.loc[row.name]) or bool(pocket_fail.loc[row.name]) or bool(ion_fail.loc[row.name]):
            weight += 1.0
        if bool(risk_fail.loc[row.name]):
            weight += 0.5
        for mut in parse_mutation_text(row.get("mutations", row.get("mutation_list", ""))):
            mut_bad[mut] += weight
            pos = mutation_position(mut)
            if pos is not None:
                pos_bad[pos] += weight

    all_pos = sorted(set(pos_good) | set(pos_bad))
    position_log_odds: dict[str, float] = {}
    for pos in all_pos:
        # Additive smoothing prevents single unlucky measurements from permanently banning a site.
        g = float(pos_good.get(pos, 0.0)) + 1.0
        b = float(pos_bad.get(pos, 0.0)) + 1.0
        position_log_odds[str(pos)] = math.log(g / b)

    all_mut = sorted(set(mut_good) | set(mut_bad))
    mutation_log_odds: dict[str, float] = {}
    for mut in all_mut:
        g = float(mut_good.get(mut, 0.0)) + 0.5
        b = float(mut_bad.get(mut, 0.0)) + 0.5
        mutation_log_odds[mut] = math.log(g / b)

    branch_success_rate: dict[str, float] = {}
    for branch, total in branch_total.items():
        # Beta smoothing keeps rare branches alive but prevents them from dominating after one success.
        branch_success_rate[branch] = (float(branch_good.get(branch, 0.0)) + 1.0) / (float(total) + 2.0)

    bad_pos_threshold = float(cfg.get("v6", {}).get("feedback", {}).get("blacklist_position_log_odds", -0.65))
    good_pos_threshold = float(cfg.get("v6", {}).get("feedback", {}).get("whitelist_position_log_odds", 0.35))
    bad_mut_threshold = float(cfg.get("v6", {}).get("feedback", {}).get("blacklist_mutation_log_odds", -0.65))
    good_mut_threshold = float(cfg.get("v6", {}).get("feedback", {}).get("whitelist_mutation_log_odds", 0.25))

    policy = {
        "enabled": True,
        "source_rows": int(len(work)),
        "success_rows": int(success.sum()),
        "failure_rows": int(failure.sum()),
        "score_column": score_col,
        "seed_rows": seed_rows[: int(cfg.get("v6", {}).get("feedback", {}).get("max_seed_rows", 160))],
        "branch_success_rate": branch_success_rate,
        "position_log_odds": position_log_odds,
        "mutation_log_odds": mutation_log_odds,
        "positive_positions": [int(p) for p, v in position_log_odds.items() if float(v) >= good_pos_threshold],
        "negative_positions": [int(p) for p, v in position_log_odds.items() if float(v) <= bad_pos_threshold],
        "positive_mutations": [m for m, v in mutation_log_odds.items() if float(v) >= good_mut_threshold],
        "negative_mutations": [m for m, v in mutation_log_odds.items() if float(v) <= bad_mut_threshold],
        "mut_good_norm": _normalize_counter(mut_good),
        "mut_bad_norm": _normalize_counter(mut_bad),
        "pos_good_norm": _normalize_counter(pos_good),
        "pos_bad_norm": _normalize_counter(pos_bad),
    }
    return policy


def build_feedback_policy_from_sources(current: pd.DataFrame, cfg: dict[str, Any]) -> dict[str, Any]:
    """Merge feedback sources while preventing proxy self-training by default.

    V6.2 only lets previous measured evidence affect generation/scoring unless
    include_current_run_feedback is explicitly enabled for diagnostics. This avoids
    the current candidate table reinforcing its own proxy ranking noise.
    """
    fb_cfg = cfg.get("v6", {}).get("feedback", {})
    tables: list[pd.DataFrame] = []
    for key in ["previous_ranked_csv", "feedback_csv"]:
        p = fb_cfg.get(key)
        table = _load_table(p)
        if not table.empty:
            table = table.copy()
            table["_feedback_source_file"] = str(p)
            tables.append(table)

    include_current = bool(fb_cfg.get("include_current_run_feedback", False))
    if include_current and current is not None and not current.empty:
        current = current.copy()
        current["_feedback_source_file"] = "current_run"
        tables.append(current)

    if not tables:
        return {"enabled": False, "reason": "no_measured_feedback_sources"}

    merged = pd.concat(tables, ignore_index=True, sort=False)
    if bool(fb_cfg.get("measured_only_policy", True)):
        has_real = merged.get("has_real_structure", pd.Series(False, index=merged.index)).astype(bool)
        structure_pass = merged.get("structure_pass", pd.Series(False, index=merged.index)).astype(bool)
        final_eligible = merged.get("v5_final_eligible", pd.Series(False, index=merged.index)).astype(bool)
        measured = has_real | structure_pass | final_eligible
        merged = merged[measured].copy()
        if merged.empty:
            return {"enabled": False, "reason": "no_measured_rows_after_filter"}

    return build_feedback_policy_from_table(merged, cfg)


def _mutation_feedback_score(mutations: list[str], policy: dict[str, Any]) -> tuple[float, int]:
    if not mutations:
        return 0.5, 0
    mut_log = {str(k): float(v) for k, v in policy.get("mutation_log_odds", {}).items()}
    pos_log = {str(k): float(v) for k, v in policy.get("position_log_odds", {}).items()}
    negative_muts = set(policy.get("negative_mutations", []))
    negative_positions = {int(x) for x in policy.get("negative_positions", [])}
    vals = []
    blacklist_hits = 0
    for mut in mutations:
        pos = mutation_position(mut)
        raw = 0.0
        if mut in mut_log:
            raw += 0.65 * mut_log[mut]
        if pos is not None and str(pos) in pos_log:
            raw += 0.35 * pos_log[str(pos)]
        if mut in negative_muts or (pos is not None and pos in negative_positions):
            blacklist_hits += 1
            raw -= 0.9
        vals.append(raw)
    mean_raw = float(np.mean(vals)) if vals else 0.0
    # Logistic transform maps log-odds evidence into a stable 0-1 score.
    score = 1.0 / (1.0 + math.exp(-mean_raw))
    return max(0.0, min(1.0, score)), blacklist_hits


def add_v6_feedback_scores(df: pd.DataFrame, policy: dict[str, Any]) -> pd.DataFrame:
    out = df.copy()
    if not policy.get("enabled", False):
        out["feedback_policy_score"] = 0.5
        out["feedback_mutation_score"] = 0.5
        out["feedback_branch_score"] = 0.5
        out["motif_blacklist_hits"] = 0
        out["feedback_enabled"] = False
        return out
    branch_rates = {str(k): float(v) for k, v in policy.get("branch_success_rate", {}).items()}
    mut_scores = []
    branch_scores = []
    blacklist_hits = []
    for _, row in out.iterrows():
        muts = parse_mutation_text(row.get("mutations", row.get("mutation_list", "")))
        m_score, hits = _mutation_feedback_score(muts, policy)
        branch = str(row.get("source_branch", row.get("source", "unknown")))
        b_score = branch_rates.get(branch, 0.5)
        mut_scores.append(m_score)
        branch_scores.append(max(0.0, min(1.0, b_score)))
        blacklist_hits.append(hits)
    out["feedback_mutation_score"] = mut_scores
    out["feedback_branch_score"] = branch_scores
    out["motif_blacklist_hits"] = blacklist_hits
    # Mutation evidence gets the highest weight because it is the direct causal feedback signal.
    out["feedback_policy_score"] = (
        0.62 * out["feedback_mutation_score"]
        + 0.28 * out["feedback_branch_score"]
        + 0.10 * (1.0 - np.minimum(out["motif_blacklist_hits"], 3) / 3.0)
    ).clip(0.0, 1.0)
    out["feedback_enabled"] = True
    return out


def write_feedback_reports(policy: dict[str, Any], out_dir: str | Path) -> None:
    out = Path(out_dir)
    fb_dir = out / "feedback"
    fb_dir.mkdir(parents=True, exist_ok=True)
    (fb_dir / "v6_feedback_policy.json").write_text(json.dumps(policy, indent=2, ensure_ascii=False), encoding="utf-8")

    pd.DataFrame([
        {"source_branch": k, "success_rate_smoothed": v}
        for k, v in sorted(policy.get("branch_success_rate", {}).items(), key=lambda x: x[1], reverse=True)
    ]).to_csv(fb_dir / "branch_feedback_v6.csv", index=False)
    pd.DataFrame([
        {"position": int(k), "position_log_odds": v}
        for k, v in sorted(policy.get("position_log_odds", {}).items(), key=lambda x: float(x[1]), reverse=True)
    ]).to_csv(fb_dir / "position_feedback_v6.csv", index=False)
    pd.DataFrame([
        {"mutation": k, "mutation_log_odds": v}
        for k, v in sorted(policy.get("mutation_log_odds", {}).items(), key=lambda x: float(x[1]), reverse=True)
    ]).to_csv(fb_dir / "mutation_feedback_v6.csv", index=False)
