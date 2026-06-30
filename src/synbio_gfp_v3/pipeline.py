from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .candidates import generate_candidates
from .data import load_exclusion, load_training_table, parse_fasta_like
from .features import featurize
from .literature_priors import add_literature_priors
from .model_cache import ModelCache
from .models import predict_bundle, train_predictive_models
from .report import write_reports, write_submission
from .scoring import add_objective_scores, add_proxy_scores
from .selection import contest_filter, layered_top6_selection
from .site_policy import build_site_policy
from .structure_filter import run_top200_structure_filter
from .utils import ensure_jsonable, setup_logger, set_seed


def _resolve_file(data_dir: Path, filename: str) -> Path:
    p = data_dir / filename
    if p.exists():
        return p
    candidates = list(data_dir.glob(filename))
    if candidates:
        return candidates[0]
    raise FileNotFoundError(f"Cannot find {filename} under {data_dir}")


def run_pipeline(config: dict[str, Any]) -> dict[str, str]:
    out_dir = Path(config["paths"]["out_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)
    logger = setup_logger(out_dir)
    set_seed(int(config.get("seed", 42)))
    t_all = time.time()
    logger.info("V3 PIPELINE START")
    (out_dir / "config_resolved.json").write_text(json.dumps(ensure_jsonable(config), indent=2, ensure_ascii=False), encoding="utf-8")

    data_dir = Path(config["paths"]["data_dir"])
    refs_path = _resolve_file(data_dir, config["paths"].get("reference_sequences_file", "AAseqs of 5 GFP proteins.txt"))
    exclusion_path = _resolve_file(data_dir, config["paths"].get("exclusion_list_file", "Exclusion_List.csv"))
    gfp_data_path = _resolve_file(data_dir, config["paths"].get("gfp_data_file", "GFP_data.xlsx"))

    refs = parse_fasta_like(refs_path)
    scaffold_key = next((k for k in refs if "sfgfp" in k.lower() or "sf" in k.lower()), None) or sorted(refs)[0]
    scaffold = refs[scaffold_key]
    logger.info("REFERENCE LOADED | scaffold=%s | length=%d | refs=%s", scaffold_key, len(scaffold), sorted(refs))
    exclusion = load_exclusion(exclusion_path)
    logger.info("EXCLUSION LOADED | n=%d", len(exclusion))

    policy = build_site_policy(scaffold, config.get("site_policy", {}).get("extra_protected_positions"))
    logger.info("SITE POLICY | protected=%d | surface=%d | loops=%d", len(policy.protected_positions), len(policy.surface_positions), len(policy.loop_positions))

    train_df, meta = load_training_table(
        gfp_data_path,
        scaffold=scaffold,
        max_rows=config.get("training", {}).get("max_train_samples"),
        seed=int(config.get("seed", 42)),
    )
    brightness_col = meta["brightness_col"]
    logger.info("TRAIN TABLE | rows=%d | brightness_col=%s | mutation_col=%s", len(train_df), brightness_col, meta.get("mutation_col"))

    cache = ModelCache(config["paths"].get("cache_dir", "cache"), logger=logger)
    feature_cfg = dict(config.get("features", {}))
    X_train, train_feature_key = cache.get_or_compute_features(
        train_df["full_sequence"].tolist(), feature_cfg, lambda: featurize(train_df["full_sequence"].tolist(), feature_cfg, logger)
    )
    model_key = cache.model_key(train_feature_key, config.get("training", {}), meta)
    bundle = cache.load_model_bundle(model_key)
    if bundle is None:
        bundle = train_predictive_models(X_train, train_df[brightness_col].astype(float).to_numpy(), config, meta, logger)
        cache.save_model_bundle(model_key, bundle)
        cache.save_metrics(model_key, bundle["metrics"])
    metrics = bundle.get("metrics") or cache.load_metrics(model_key) or {}

    logger.info("CANDIDATE GENERATION START")
    cand = generate_candidates(scaffold, policy, config)
    logger.info("CANDIDATE GENERATION DONE | n=%d | by_source=%s", len(cand), cand["source"].value_counts().to_dict())
    explain = [policy.explain_mutations(scaffold, s) for s in cand["sequence"]]
    cand = pd.concat([cand.reset_index(drop=True), pd.DataFrame(explain)], axis=1)
    cand = add_literature_priors(cand, scaffold, policy)
    cand = add_proxy_scores(cand)

    X_cand, cand_feature_key = cache.get_or_compute_features(
        cand["sequence"].tolist(), feature_cfg, lambda: featurize(cand["sequence"].tolist(), feature_cfg, logger)
    )
    pred = predict_bundle(bundle, X_cand)
    for k, v in pred.items():
        cand[k] = v
    cand = add_objective_scores(cand)
    cand["Seq_ID"] = [f"candidate_{i}" for i in range(len(cand))]
    cand = run_top200_structure_filter(cand, config, out_dir, logger)
    cand = add_objective_scores(cand)
    ranked = cand.sort_values("objective", ascending=False).reset_index(drop=True)
    ranked.to_csv(out_dir / "all_ranked_candidates.csv", index=False)

    filtered, diag = contest_filter(ranked, exclusion, config)
    diag.to_csv(out_dir / "filter_diagnostics.csv", index=False)
    logger.info("FILTER DONE | ranked=%d | passed=%d", len(ranked), len(filtered))
    if filtered.empty:
        logger.warning("No candidates passed risk filter; relaxing risk gate to contest-only valid candidates.")
        filtered = ranked.copy()
    selected = layered_top6_selection(filtered, config)
    selected.to_csv(out_dir / "selected_candidates.csv", index=False)
    write_submission(selected, config.get("team_name", "YourTeamName"), out_dir)
    write_reports(selected, ranked, metrics, config, out_dir)
    logger.info("V3 PIPELINE DONE | seconds=%.2f | selected=%d | submission=%s", time.time() - t_all, len(selected), out_dir / "submission.csv")
    return {
        "submission": str(out_dir / "submission.csv"),
        "selected_candidates": str(out_dir / "selected_candidates.csv"),
        "ranked_candidates": str(out_dir / "all_ranked_candidates.csv"),
        "report": str(out_dir / "v3_pipeline_report.md"),
        "metrics": str(out_dir / "training_metrics.json"),
    }
