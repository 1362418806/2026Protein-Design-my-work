from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from synbio_gfp_v3.config import load_config, resolve_config
from synbio_gfp_v3.pipeline import run_pipeline


def parse_args():
    p = argparse.ArgumentParser(description="Run SynBio GFP V3 pipeline")
    p.add_argument("--config", required=True)
    p.add_argument("--data-dir", required=True)
    p.add_argument("--team-name", required=True)
    p.add_argument("--feature-mode", choices=["simple", "esm"], default="simple")
    p.add_argument("--esm-model", default="esm2_t30_150M_UR50D")
    p.add_argument("--max-train-samples", type=int, default=20000)
    p.add_argument("--n-candidates", type=int, default=20000)
    p.add_argument("--out-dir", default="outputs/run_v3")
    p.add_argument("--cache-dir", default="cache")
    p.add_argument("--references", default=None)
    p.add_argument("--exclusion-list", default=None)
    p.add_argument("--gfp-data", default=None)
    p.add_argument("--reference-pdb", default=None)
    p.add_argument("--predicted-pdb-dir", default=None)
    p.add_argument("--skip-structure", action="store_true")
    return p.parse_args()


def main():
    args = parse_args()
    cfg = resolve_config(load_config(args.config), args)
    outputs = run_pipeline(cfg)
    print("\nV3 outputs:")
    for k, v in outputs.items():
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
