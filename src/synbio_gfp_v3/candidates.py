from __future__ import annotations

import random
from typing import Any

import pandas as pd

from .site_policy import SitePolicy

SAFE_SUBS = {
    "K": "E", "R": "E", "N": "D", "Q": "E",
    "L": "V", "I": "V", "V": "I", "F": "Y", "Y": "F",
    "S": "T", "T": "S", "A": "S", "M": "L",
}
SURFACE_ACID_SUBS = {"K": ["E", "D"], "R": ["E", "D"], "N": ["D"], "Q": ["E"], "A": ["S", "T"], "L": ["S", "T"], "V": ["T", "S"]}
AA20 = "ACDEFGHIKLMNPQRSTVWY"


def mutate(seq: str, changes: list[tuple[int, str]]) -> str:
    arr = list(seq)
    for idx0, aa in changes:
        arr[idx0] = aa
    return "".join(arr)


def generate_candidates(scaffold: str, policy: SitePolicy, config: dict[str, Any]) -> pd.DataFrame:
    cg = config.get("candidate_generation", {})
    n = int(cg.get("n_candidates", 20000))
    seed = int(config.get("seed", 42))
    rng = random.Random(seed)
    allowed = policy.allowed_positions()
    surface = [p - 1 for p in policy.surface_positions]
    rows = []
    seen = {scaffold}
    while len(rows) < n and len(seen) < n * 5:
        source = rng.choices(
            ["conservative_random", "tgp_surface", "balanced", "exploratory"],
            weights=[0.34, 0.28, 0.28, 0.10],
            k=1,
        )[0]
        if source == "conservative_random":
            k = rng.choice([2, 3, 4, 5])
            pool = allowed
        elif source == "tgp_surface":
            k = rng.choice([3, 4, 5, 6, 7])
            pool = surface or allowed
        elif source == "balanced":
            k = rng.choice([5, 6, 7, 8])
            pool = allowed
        else:
            k = rng.choice([8, 9, 10, 11, 12])
            pool = allowed
        if len(pool) < k:
            pool = allowed
        positions = rng.sample(pool, k)
        changes = []
        for idx0 in positions:
            wt = scaffold[idx0]
            if source == "tgp_surface" and wt in SURFACE_ACID_SUBS:
                aa = rng.choice(SURFACE_ACID_SUBS[wt])
            elif wt in SAFE_SUBS and rng.random() < 0.78:
                aa = SAFE_SUBS[wt]
            else:
                aa = rng.choice([a for a in AA20 if a != wt])
            changes.append((idx0, aa))
        seq = mutate(scaffold, changes)
        if seq in seen:
            continue
        seen.add(seq)
        rows.append({"sequence": seq, "source": source})
    return pd.DataFrame(rows)
