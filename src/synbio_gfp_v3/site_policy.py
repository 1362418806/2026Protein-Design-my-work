from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable

# One-based sfGFP positions. These are intentionally conservative and should be
# treated as hard mutation-protection defaults unless the user explicitly edits them.
CHROMOPHORE = {65, 66, 67}
SUPERFOLDER_MOTIFS = {30, 39, 64, 65, 99, 105, 145, 153, 163, 171, 206}
ION_NETWORK = {17, 30, 32, 115, 122}
CHROMOPHORE_POCKET_8A = {
    39, 63, 64, 65, 66, 67, 68, 69, 94, 96, 145, 148, 165, 167, 203, 205, 222,
}
# Approximate high-risk buried/core positions in GFP beta barrel. This is a safe
# initial prior and is refined by structure_filter when actual structures exist.
BARREL_CORE_PRIOR = {
    8, 10, 12, 14, 16, 18, 20, 22, 43, 45, 47, 49, 51, 53, 70, 71, 73, 75, 77,
    79, 81, 83, 100, 102, 104, 106, 108, 110, 128, 130, 132, 134, 136, 138,
    155, 157, 159, 161, 163, 181, 183, 185, 187, 189, 191, 207, 209, 211, 213,
    215, 217,
}
# Surface-like positions used for TGP-inspired conservative charge engineering.
SURFACE_PRIOR = {
    23, 24, 25, 28, 29, 33, 34, 35, 36, 40, 41, 42, 54, 55, 56, 57, 58, 59,
    60, 85, 86, 87, 88, 89, 90, 91, 112, 113, 114, 116, 117, 118, 119, 120,
    121, 123, 124, 139, 140, 141, 142, 143, 144, 146, 147, 149, 150, 151, 152,
    172, 173, 174, 175, 176, 177, 178, 179, 180, 192, 193, 194, 195, 196, 197,
    198, 199, 200, 201, 202, 218, 219, 220, 221, 223, 224, 225, 226, 227, 228,
    229, 230, 231, 232, 233, 234, 235, 236,
}
LOOP_PRIOR = {23, 24, 25, 54, 55, 56, 85, 86, 87, 112, 113, 114, 139, 140, 141, 172, 173, 174, 218, 219, 220, 221, 222, 223}


@dataclass
class SitePolicy:
    length: int
    protected_positions: set[int] = field(default_factory=set)
    surface_positions: set[int] = field(default_factory=set)
    loop_positions: set[int] = field(default_factory=set)
    notes: dict[str, list[int]] = field(default_factory=dict)

    def is_protected(self, pos1: int) -> bool:
        return pos1 in self.protected_positions

    def allowed_positions(self) -> list[int]:
        # Return zero-based mutable positions for Python sequence indexing.
        return [i - 1 for i in range(1, self.length + 1) if i not in self.protected_positions]

    def explain_mutations(self, wt: str, seq: str) -> dict:
        mutations = []
        protected_hits = []
        surface_hits = []
        loop_hits = []
        for i, (a, b) in enumerate(zip(wt, seq), start=1):
            if a == b:
                continue
            item = f"{a}{i}{b}"
            mutations.append(item)
            if i in self.protected_positions:
                protected_hits.append(item)
            if i in self.surface_positions:
                surface_hits.append(item)
            if i in self.loop_positions:
                loop_hits.append(item)
        return {
            "mutation_count": len(mutations),
            "mutation_list": ";".join(mutations),
            "protected_site_mutations": ";".join(protected_hits),
            "surface_mutations": ";".join(surface_hits),
            "loop_mutations": ";".join(loop_hits),
        }


def build_site_policy(scaffold: str, extra_protected: Iterable[int] | None = None) -> SitePolicy:
    length = len(scaffold)
    protected = set()
    protected |= {p for p in CHROMOPHORE if p <= length}
    protected |= {p for p in SUPERFOLDER_MOTIFS if p <= length}
    protected |= {p for p in ION_NETWORK if p <= length}
    protected |= {p for p in CHROMOPHORE_POCKET_8A if p <= length}
    protected |= {p for p in BARREL_CORE_PRIOR if p <= length}
    if extra_protected:
        protected |= {p for p in extra_protected if 1 <= p <= length}
    surface = {p for p in SURFACE_PRIOR if p <= length and p not in protected}
    loop = {p for p in LOOP_PRIOR if p <= length and p not in protected}
    return SitePolicy(
        length=length,
        protected_positions=protected,
        surface_positions=surface,
        loop_positions=loop,
        notes={
            "chromophore": sorted(CHROMOPHORE),
            "superfolder_motifs": sorted(SUPERFOLDER_MOTIFS),
            "ion_network": sorted(ION_NETWORK),
            "chromophore_pocket_8A": sorted(CHROMOPHORE_POCKET_8A),
            "barrel_core_prior": sorted(BARREL_CORE_PRIOR),
        },
    )
