from __future__ import annotations

from typing import Any

import pandas as pd

from .utils import validate_sequence


def contest_filter(df: pd.DataFrame, exclusion: set[str], config: dict[str, Any]) -> tuple[pd.DataFrame, pd.DataFrame]:
    sel_cfg = config.get("selection", {})
    min_len = int(sel_cfg.get("min_len", 220)); max_len = int(sel_cfg.get("max_len", 250))
    risk_max = float(sel_cfg.get("low_brightness_risk_max", 0.55))
    rows = []
    keep = []
    for idx, row in df.iterrows():
        seq = row["sequence"]
        v = validate_sequence(seq, min_len, max_len)
        in_exclusion = seq in exclusion
        risk_ok = float(row.get("low_brightness_risk", 0.0)) <= risk_max
        ok = all(v.values()) and not in_exclusion and risk_ok
        rows.append({"index": idx, **v, "in_exclusion_list": in_exclusion, "risk_ok": risk_ok, "passed": ok})
        if ok:
            keep.append(idx)
    diag = pd.DataFrame(rows)
    return df.loc[keep].copy(), diag


def _take_best(pool: pd.DataFrame, selected: list[int], n: int) -> list[int]:
    ids = []
    for idx in pool.sort_values("objective", ascending=False).index:
        if idx not in selected and idx not in ids:
            ids.append(idx)
        if len(ids) >= n:
            break
    return ids


def layered_top6_selection(df: pd.DataFrame, config: dict[str, Any]) -> pd.DataFrame:
    if df.empty:
        return df.copy()
    selected: list[int] = []
    # Layered strategy: keep Best Top-1 probability by mixing conservative,
    # balanced, thermal-prior, structure-prior, and exploratory mechanisms.
    layers = [
        ("conservative_high_brightness", df[df["mutation_count"].between(2, 4)], 2),
        ("balanced_brightness_stability", df[df["mutation_count"].between(5, 7)], 2),
        ("tgp_thermal_surface", df[df["tgp_surface_prior"] >= df["tgp_surface_prior"].quantile(0.70)], 1),
        ("structure_guided_explorer", df[df["mutation_count"].between(8, 12)], 1),
    ]
    labels = {}
    for label, pool, n in layers:
        take = _take_best(pool, selected, n)
        selected.extend(take)
        for idx in take:
            labels[idx] = label
    if len(selected) < 6:
        for idx in _take_best(df, selected, 6 - len(selected)):
            selected.append(idx)
            labels[idx] = "objective_fallback"
    out = df.loc[selected[:6]].copy()
    out["selection_layer"] = [labels.get(i, "objective_fallback") for i in out.index]
    out = out.sort_values(["selection_layer", "objective"], ascending=[True, False]).reset_index(drop=True)
    out["Seq_ID"] = [f"Seq_{i+1}" for i in range(len(out))]
    return out
