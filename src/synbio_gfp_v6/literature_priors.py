from __future__ import annotations

import math
from typing import Any

import numpy as np
import pandas as pd

from .site_policy import SitePolicy, CHROMOPHORE_POCKET_8A, ION_NETWORK, SUPERFOLDER_MOTIFS

ACIDIC = set("DE")
BASIC = set("KRH")
HYDROPHOBIC = set("AILMFWVY")
POLAR = set("STNQ")


def _mutations(wt: str, seq: str) -> list[tuple[int, str, str]]:
    return [(i, a, b) for i, (a, b) in enumerate(zip(wt, seq), start=1) if a != b]


def score_one(seq: str, wt: str, policy: SitePolicy) -> dict[str, Any]:
    muts = _mutations(wt, seq)
    n_mut = len(muts)
    protected = [f"{a}{i}{b}" for i, a, b in muts if i in policy.protected_positions]
    pocket = [f"{a}{i}{b}" for i, a, b in muts if i in CHROMOPHORE_POCKET_8A]
    ion = [f"{a}{i}{b}" for i, a, b in muts if i in ION_NETWORK]
    superfolder = [f"{a}{i}{b}" for i, a, b in muts if i in SUPERFOLDER_MOTIFS]

    # Fitness landscape prior: GFP tolerates few mutations; high mutation burden
    # and protected-site changes are penalized to avoid negative epistasis collapse.
    mutation_burden_penalty = max(0.0, (n_mut - 5) / 10.0)
    epistasis_penalty = max(0.0, (n_mut * (n_mut - 1) / 2 - 15) / 60.0)
    protected_penalty = 2.0 * len(protected) + 1.25 * len(pocket) + 1.5 * len(ion)
    fitness_landscape_prior = math.exp(-(0.45 * mutation_burden_penalty + 0.35 * epistasis_penalty + protected_penalty))

    # TGP-inspired surface prior: reward safe surface K/R/N/Q -> D/E or hydrophobic
    # to polar/acidic changes, penalize new positive or hydrophobic surface patches.
    tgp_bonus = 0.0
    tgp_penalty = 0.0
    tgp_events = []
    for i, a, b in muts:
        if i not in policy.surface_positions:
            continue
        if a in set("KRNQ") and b in ACIDIC:
            tgp_bonus += 1.0
            tgp_events.append(f"{a}{i}{b}:surface_negative_charge")
        if a in HYDROPHOBIC and b in (POLAR | ACIDIC):
            tgp_bonus += 0.6
            tgp_events.append(f"{a}{i}{b}:surface_hydrophilic")
        if b in BASIC:
            tgp_penalty += 0.5
        if b in HYDROPHOBIC and a not in HYDROPHOBIC:
            tgp_penalty += 0.4
    tgp_surface_prior = float(1 / (1 + math.exp(-(tgp_bonus - tgp_penalty))))

    # StayGold/mBaoJin-inspired pocket prior: preserve pocket rigidity, avoid charged
    # changes and disruptive Pro/Gly inside chromophore pocket.
    pocket_charge_mut = [m for m in muts if m[0] in CHROMOPHORE_POCKET_8A and (m[2] in ACIDIC or m[2] in BASIC)]
    pocket_pg_mut = [m for m in muts if m[0] in CHROMOPHORE_POCKET_8A and m[2] in set("PG")]
    pocket_score = 1.0 - min(1.0, 0.35 * len(pocket_charge_mut) + 0.30 * len(pocket_pg_mut) + 0.20 * len(pocket))

    superfolder_protection_score = 1.0 - min(1.0, 0.45 * len(ion) + 0.30 * len(superfolder))

    literature_prior = float(np.clip(
        0.34 * fitness_landscape_prior
        + 0.24 * tgp_surface_prior
        + 0.24 * pocket_score
        + 0.18 * superfolder_protection_score,
        0.0,
        1.0,
    ))
    return {
        "fitness_landscape_prior": float(fitness_landscape_prior),
        "tgp_surface_prior": float(tgp_surface_prior),
        "pocket_prior": float(pocket_score),
        "superfolder_protection_score": float(superfolder_protection_score),
        "literature_prior": literature_prior,
        "literature_events": ";".join(tgp_events),
        "protected_mutation_count": len(protected),
        "pocket_mutation_count": len(pocket),
        "ion_network_mutation_count": len(ion),
        "risk_explanation": _risk_text(n_mut, protected, pocket, ion, tgp_events),
    }


def _risk_text(n_mut: int, protected: list[str], pocket: list[str], ion: list[str], tgp_events: list[str]) -> str:
    parts = []
    if n_mut <= 4:
        parts.append("low mutation burden; conservative brightness risk")
    elif n_mut <= 8:
        parts.append("medium mutation burden; balanced exploration")
    else:
        parts.append("high mutation burden; negative epistasis risk should be checked")
    if protected:
        parts.append(f"protected-site mutations detected: {','.join(protected[:5])}")
    if pocket:
        parts.append("chromophore pocket perturbation risk")
    if ion:
        parts.append("S30R ion-network perturbation risk")
    if tgp_events:
        parts.append("TGP-inspired safe surface stabilization events present")
    return " | ".join(parts)


def add_literature_priors(df: pd.DataFrame, wt: str, policy: SitePolicy, sequence_col: str = "sequence") -> pd.DataFrame:
    rows = [score_one(seq, wt, policy) for seq in df[sequence_col].astype(str)]
    prior_df = pd.DataFrame(rows)
    return pd.concat([df.reset_index(drop=True), prior_df], axis=1)
