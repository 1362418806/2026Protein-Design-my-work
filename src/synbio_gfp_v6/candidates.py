from __future__ import annotations

import random
import re
from pathlib import Path
from typing import Any

import pandas as pd

from .feedback import build_feedback_policy_from_table, mutation_position, parse_mutation_text
from .site_policy import SitePolicy

SAFE_SUBS = {
    "K": "E", "R": "E", "N": "D", "Q": "E",
    "L": "V", "I": "V", "V": "I", "F": "Y", "Y": "F",
    "S": "T", "T": "S", "A": "S", "M": "L",
}
SURFACE_ACID_SUBS = {"K": ["E", "D"], "R": ["E", "D"], "N": ["D"], "Q": ["E"], "A": ["S", "T"], "L": ["S", "T"], "V": ["T", "S"]}
# V6.2 keeps thermal exploration conservative: surface-patch substitutions are
# restricted to polar/acidic replacements that should not invade the barrel core.
SURFACE_PATCH_SUBS = {"K": ["E", "D"], "R": ["E", "D"], "N": ["D", "S"], "Q": ["E", "T"], "A": ["S", "T"], "L": ["S"], "V": ["T"], "S": ["D", "N"], "T": ["E", "Q"]}
# Loop-rigidity branch intentionally avoids global Pro/Gly enrichment. Pro is
# only allowed on a tiny explicit turn whitelist handled inside the generator.
LOOP_RIGIDITY_SUBS = {"G": ["S", "A", "T"], "S": ["N", "T"], "T": ["Q", "S"], "A": ["S", "T"], "N": ["D", "S"], "Q": ["E", "T"], "D": ["N"], "E": ["Q"]}
LOOP_TURN_PRO_WHITELIST = {24, 55, 86, 113, 140, 173, 220}
AA20 = "ACDEFGHIKLMNPQRSTVWY"
MUT_RE = re.compile(r"([ACDEFGHIKLMNPQRSTVWY])(\d{1,3})([ACDEFGHIKLMNPQRSTVWY])")


def mutate(seq: str, changes: list[tuple[int, str]]) -> str:
    arr = list(seq)
    for idx0, aa in changes:
        arr[idx0] = aa
    return "".join(arr)


def _load_optional_table(path: str | None) -> pd.DataFrame:
    if not path:
        return pd.DataFrame()
    p = Path(path)
    if not p.exists():
        return pd.DataFrame()
    return pd.read_excel(p) if p.suffix.lower() in {".xlsx", ".xls"} else pd.read_csv(p)


def _mutation_to_change(scaffold: str, mut: str, policy: SitePolicy) -> tuple[int, str] | None:
    m = MUT_RE.fullmatch(str(mut))
    if not m:
        return None
    wt, pos_text, aa = m.groups()
    pos1 = int(pos_text)
    idx0 = pos1 - 1
    if idx0 < 0 or idx0 >= len(scaffold):
        return None
    if policy.is_protected(pos1):
        return None
    if scaffold[idx0] == aa:
        return None
    # Do not require the previous-run WT letter to match. GFP homologous tables can contain tolerant numbering.
    return idx0, aa


def _weighted_choice(rng: random.Random, items: list[Any], weights: list[float]) -> Any:
    if not items:
        raise ValueError("weighted choice requires non-empty items")
    clean = [max(0.0, float(w)) for w in weights]
    if sum(clean) <= 0:
        clean = [1.0] * len(items)
    return rng.choices(items, weights=clean, k=1)[0]


def _build_prior_from_config(config: dict[str, Any]) -> dict[str, Any]:
    fb_cfg = config.get("v6", {}).get("feedback", {})
    tables = []
    for key in ["previous_ranked_csv", "feedback_csv"]:
        df = _load_optional_table(fb_cfg.get(key))
        if not df.empty:
            tables.append(df)
    if not tables:
        return {"enabled": False, "reason": "no_previous_feedback"}
    return build_feedback_policy_from_table(pd.concat(tables, ignore_index=True, sort=False), config)


def _baseline_candidate(scaffold: str, policy: SitePolicy, rng: random.Random, source: str) -> str:
    allowed = policy.allowed_positions()
    surface = [p - 1 for p in policy.surface_positions]
    loops = [p - 1 for p in policy.loop_positions]

    if source == "conservative_random":
        k = rng.choice([2, 3, 4, 5])
        pool = allowed
    elif source == "surface_patch":
        # V6.2 surface-patch branch targets sparse high-SASA charge/polar edits.
        # It replaces the older coarse acidification behavior with lower mutation burden.
        k = rng.choice([2, 3, 4])
        pool = surface or allowed
    elif source == "loop_rigidity":
        # V6.2 loop-rigidity branch explicitly samples loop priors and avoids beta-strand Pro.
        k = rng.choice([2, 3, 4, 5])
        pool = loops or surface or allowed
    elif source == "tgp_surface":
        k = rng.choice([3, 4, 5])
        pool = surface or allowed
    elif source == "balanced":
        k = rng.choice([4, 5, 6])
        pool = allowed
    else:
        k = rng.choice([6, 7, 8])
        pool = allowed
    if len(pool) < k:
        pool = allowed
    positions = rng.sample(pool, k)
    changes = []
    acid_edits = 0
    for idx0 in positions:
        wt = scaffold[idx0]
        pos1 = idx0 + 1
        if source == "surface_patch" and wt in SURFACE_PATCH_SUBS:
            choices = SURFACE_PATCH_SUBS[wt]
            # Keep acidic edits distributed and bounded to avoid over-acidifying the protein surface.
            if acid_edits >= 2:
                choices = [a for a in choices if a not in {"D", "E"}] or choices
            aa = rng.choice(choices)
            acid_edits += int(aa in {"D", "E"})
        elif source == "loop_rigidity":
            choices = LOOP_RIGIDITY_SUBS.get(wt, [SAFE_SUBS.get(wt, "S")])
            if pos1 in LOOP_TURN_PRO_WHITELIST and rng.random() < 0.08 and wt != "P":
                # Proline is used only as a rare loop-turn option, never as a broad rigidity shortcut.
                choices = choices + ["P"]
            aa = rng.choice([a for a in choices if a != wt] or [SAFE_SUBS.get(wt, rng.choice([a for a in AA20 if a != wt]))])
        elif source == "tgp_surface" and wt in SURFACE_ACID_SUBS:
            aa = rng.choice(SURFACE_ACID_SUBS[wt])
        elif wt in SAFE_SUBS and rng.random() < 0.78:
            aa = SAFE_SUBS[wt]
        else:
            aa = rng.choice([a for a in AA20 if a != wt])
        changes.append((idx0, aa))
    return mutate(scaffold, changes)


def _feedback_candidate(scaffold: str, policy: SitePolicy, rng: random.Random, policy_dict: dict[str, Any], source: str) -> tuple[str, dict[str, Any]]:
    positive_muts = list(policy_dict.get("positive_mutations", []))
    negative_muts = set(policy_dict.get("negative_mutations", []))
    negative_pos = {int(x) for x in policy_dict.get("negative_positions", [])}
    positive_pos = {int(x) for x in policy_dict.get("positive_positions", [])}
    mut_good = {str(k): float(v) for k, v in policy_dict.get("mut_good_norm", {}).items()}
    seed_rows = [r for r in policy_dict.get("seed_rows", []) if isinstance(r, dict) and isinstance(r.get("sequence"), str)]
    changes: dict[int, str] = {}
    parent = ""

    if source in {"feedback_seed_repair", "measured_feedback_repair"} and seed_rows:
        # Start from a measured-good parent, remove failed motifs, then add a small number of successful motifs.
        seed = rng.choice(seed_rows[: max(1, min(len(seed_rows), 80))])
        parent = str(seed.get("candidate_id", ""))
        for mut in parse_mutation_text(seed.get("mutations", "")):
            pos = mutation_position(mut)
            if pos is None or pos in negative_pos or mut in negative_muts:
                continue
            change = _mutation_to_change(scaffold, mut, policy)
            if change is not None:
                changes[change[0]] = change[1]
        add_k = rng.choice([1, 2, 3])
    elif source == "feedback_motif_recombine":
        # Recombine independently successful motifs while respecting the protected-site policy.
        add_k = rng.choice([3, 4, 5, 6])
    elif source == "feedback_surface_amplify":
        # TGP-inspired surface design constrained by observed safe positions from feedback.
        surface = [p for p in policy.surface_positions if p not in negative_pos]
        if positive_pos:
            surface = [p for p in surface if p in positive_pos] or surface
        add_k = rng.choice([4, 5, 6, 7])
        for pos1 in rng.sample(surface, min(add_k, len(surface))):
            idx0 = pos1 - 1
            wt = scaffold[idx0]
            aa = rng.choice(SURFACE_ACID_SUBS.get(wt, [SAFE_SUBS.get(wt, rng.choice([a for a in AA20 if a != wt]))]))
            changes[idx0] = aa
        add_k = rng.choice([0, 1, 2])
    else:
        # Guarded exploration keeps mutation burden moderate and excludes failed positions.
        add_k = rng.choice([4, 5, 6, 7, 8])

    valid_positive = []
    for mut in positive_muts:
        pos = mutation_position(mut)
        if pos is None or pos in negative_pos or mut in negative_muts:
            continue
        ch = _mutation_to_change(scaffold, mut, policy)
        if ch is not None:
            valid_positive.append((mut, ch))
    if valid_positive and add_k > 0:
        sampled = rng.sample(valid_positive, min(add_k, len(valid_positive))) if len(valid_positive) <= add_k else rng.choices(
            valid_positive,
            weights=[1.0 + mut_good.get(m[0], 0.0) for m in valid_positive],
            k=add_k,
        )
        for mut, (idx0, aa) in sampled:
            if (idx0 + 1) not in negative_pos:
                changes[idx0] = aa

    if not changes:
        # Fallback to conservative safe substitutions outside the failed positions.
        allowed = [idx0 for idx0 in policy.allowed_positions() if (idx0 + 1) not in negative_pos]
        for idx0 in rng.sample(allowed, min(4, len(allowed))):
            wt = scaffold[idx0]
            changes[idx0] = SAFE_SUBS.get(wt, rng.choice([a for a in AA20 if a != wt]))
    seq = mutate(scaffold, sorted(changes.items()))
    return seq, {"seed_parent_id": parent, "feedback_mutation_count": len(changes)}


def generate_candidates(scaffold: str, policy: SitePolicy, config: dict[str, Any]) -> pd.DataFrame:
    cg = config.get("candidate_generation", {})
    n = int(cg.get("n_candidates", 20000))
    seed = int(config.get("seed", 42))
    rng = random.Random(seed)
    feedback_policy = _build_prior_from_config(config)
    feedback_enabled = bool(feedback_policy.get("enabled", False)) and bool(config.get("v6", {}).get("feedback", {}).get("enabled", True))

    feedback_fraction = float(config.get("v6", {}).get("feedback", {}).get("feedback_candidate_fraction", 0.55 if feedback_enabled else 0.0))
    n_feedback = int(round(n * feedback_fraction)) if feedback_enabled else 0
    n_baseline = max(0, n - n_feedback)

    rows: list[dict[str, Any]] = []
    seen = {scaffold}

    baseline_sources = ["conservative_random", "surface_patch", "loop_rigidity", "tgp_surface", "balanced", "exploratory"]
    baseline_weights = [0.28, 0.18, 0.14, 0.16, 0.18, 0.06] if feedback_enabled else [0.30, 0.20, 0.16, 0.16, 0.14, 0.04]
    attempts = 0
    while len(rows) < n_baseline and attempts < n_baseline * 40 + 1000:
        attempts += 1
        source = rng.choices(baseline_sources, weights=baseline_weights, k=1)[0]
        seq = _baseline_candidate(scaffold, policy, rng, source)
        if seq in seen:
            continue
        seen.add(seq)
        rows.append({"sequence": seq, "source": source, "generator_mode": "baseline"})

    feedback_sources = ["measured_feedback_repair", "feedback_seed_repair", "feedback_motif_recombine", "feedback_surface_amplify", "feedback_guarded_explore"]
    feedback_weights = [0.30, 0.24, 0.22, 0.16, 0.08]
    attempts = 0
    while feedback_enabled and len(rows) < n and attempts < n_feedback * 80 + 3000:
        attempts += 1
        source = rng.choices(feedback_sources, weights=feedback_weights, k=1)[0]
        seq, meta = _feedback_candidate(scaffold, policy, rng, feedback_policy, source)
        if seq in seen:
            continue
        seen.add(seq)
        rows.append({"sequence": seq, "source": source, "generator_mode": "feedback_conditioned", **meta})

    # If feedback generation underfills because the measured motif set is small, preserve the requested candidate count with safe baseline candidates.
    attempts = 0
    while len(rows) < n and attempts < n * 40 + 1000:
        attempts += 1
        source = rng.choices(baseline_sources, weights=baseline_weights, k=1)[0]
        seq = _baseline_candidate(scaffold, policy, rng, source)
        if seq in seen:
            continue
        seen.add(seq)
        rows.append({"sequence": seq, "source": source, "generator_mode": "baseline_fill"})

    return pd.DataFrame(rows)
