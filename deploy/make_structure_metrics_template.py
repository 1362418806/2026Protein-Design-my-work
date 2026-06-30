#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser(description="Create a structure metrics CSV template from V4 stage-1 priority candidates")
    parser.add_argument("--priority-csv", required=True, help="Path to structure_priority_top200_v4.csv")
    parser.add_argument("--out", required=True, help="Output CSV path")
    args = parser.parse_args()

    src = pd.read_csv(args.priority_csv)
    keep = [c for c in ["candidate_id", "Seq_ID", "sequence", "mutations", "v4_objective"] if c in src.columns]
    out = src[keep].copy()
    # These columns are intentionally blank: they must be filled by real AF/PDB analysis or an external metrics script.
    out["barrel_rmsd"] = ""
    out["pocket_rmsd"] = ""
    out["ion_network_rmsd"] = ""
    out["structure_status"] = "pending_real_structure"
    out["structure_pass"] = ""
    out.to_csv(args.out, index=False)
    print(f"Wrote template: {Path(args.out).resolve()}")


if __name__ == "__main__":
    main()
