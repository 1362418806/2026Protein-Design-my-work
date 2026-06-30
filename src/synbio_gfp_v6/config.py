from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml


def deep_update(base: dict[str, Any], updates: dict[str, Any]) -> dict[str, Any]:
    out = dict(base)
    for k, v in updates.items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = deep_update(out[k], v)
        else:
            out[k] = v
    return out


def load_config(path: str | Path) -> dict[str, Any]:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def resolve_config(config: dict[str, Any], args: Any) -> dict[str, Any]:
    cfg = dict(config)
    data_dir = Path(args.data_dir).expanduser().resolve()
    out_dir = Path(args.out_dir).expanduser()
    cfg.setdefault("paths", {})
    cfg["paths"]["data_dir"] = str(data_dir)
    cfg["paths"]["out_dir"] = str(out_dir)
    cfg["paths"]["cache_dir"] = str(Path(args.cache_dir).expanduser())
    cfg["team_name"] = args.team_name
    cfg.setdefault("features", {})["mode"] = args.feature_mode
    cfg["features"]["esm_model"] = args.esm_model
    cfg.setdefault("training", {})["max_train_samples"] = args.max_train_samples
    cfg.setdefault("candidate_generation", {})["n_candidates"] = args.n_candidates
    if args.references:
        cfg["paths"]["reference_sequences_file"] = args.references
    if args.exclusion_list:
        cfg["paths"]["exclusion_list_file"] = args.exclusion_list
    if args.gfp_data:
        cfg["paths"]["gfp_data_file"] = args.gfp_data
    if args.reference_pdb:
        cfg.setdefault("structure", {})["reference_pdb"] = args.reference_pdb
    if args.predicted_pdb_dir:
        cfg.setdefault("structure", {})["predicted_pdb_dir"] = args.predicted_pdb_dir
    if args.skip_structure:
        cfg.setdefault("structure", {})["enabled"] = False
    return cfg
