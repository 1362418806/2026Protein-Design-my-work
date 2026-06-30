from __future__ import annotations

import numpy as np
import pandas as pd


def rank01(values: pd.Series, higher_is_better: bool = True) -> pd.Series:
    ranks = values.rank(method="average", ascending=not higher_is_better)
    if len(values) <= 1:
        return pd.Series([1.0] * len(values), index=values.index)
    return 1.0 - (ranks - 1) / (len(values) - 1)


def add_objective_scores(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["brightness_lcb"] = out["brightness_pred"] - 0.75 * out.get("brightness_uncertainty", 0.0)
    out["brightness_rank"] = rank01(out["brightness_lcb"], True)
    out["risk_rank"] = rank01(out["low_brightness_risk"], False)
    out["stability_rank"] = rank01(out["stability_proxy"], True)
    out["prior_rank"] = rank01(out["literature_prior"], True)
    out["structure_rank"] = rank01(out.get("structure_score", pd.Series(0.5, index=out.index)), True)
    out["objective"] = (
        0.38 * out["brightness_rank"]
        + 0.18 * out["risk_rank"]
        + 0.16 * out["stability_rank"]
        + 0.18 * out["prior_rank"]
        + 0.10 * out["structure_rank"]
    )
    return out


def add_proxy_scores(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    # Stability proxy rewards TGP-like surface stabilization and moderate mutation load.
    out["mutation_count"] = out["mutation_count"].astype(int)
    mut_penalty = np.clip((out["mutation_count"] - 6) / 12.0, 0.0, 1.0)
    out["stability_proxy"] = np.clip(
        0.55 + 0.40 * out["tgp_surface_prior"] + 0.20 * out["superfolder_protection_score"] - 0.18 * mut_penalty,
        0.0,
        1.5,
    )
    out["foldability_proxy"] = np.clip(
        0.40 + 0.35 * out["fitness_landscape_prior"] + 0.25 * out["superfolder_protection_score"] - 0.10 * mut_penalty,
        0.0,
        1.0,
    )
    return out
