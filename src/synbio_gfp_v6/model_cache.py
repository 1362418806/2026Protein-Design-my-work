from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any, Callable

import joblib
import numpy as np

from .utils import ensure_jsonable, sequence_hash, stable_hash


class ModelCache:
    """Cache ESM features, trained sklearn models, and training metrics.

    The cache key intentionally includes sequence hashes, feature mode, model name,
    training metadata, and hyperparameter blocks so that cached artifacts are reused
    only when they are logically equivalent to the current run.
    """

    def __init__(self, cache_dir: str | Path, logger: logging.Logger | None = None):
        self.cache_dir = Path(cache_dir)
        self.feature_dir = self.cache_dir / "features"
        self.model_dir = self.cache_dir / "models"
        self.metric_dir = self.cache_dir / "metrics"
        for d in [self.feature_dir, self.model_dir, self.metric_dir]:
            d.mkdir(parents=True, exist_ok=True)
        self.logger = logger or logging.getLogger(__name__)

    def feature_key(self, seqs: list[str], feature_cfg: dict[str, Any]) -> str:
        return stable_hash({"seq_hash": sequence_hash(seqs), "feature_cfg": feature_cfg, "n": len(seqs)})

    def load_features(self, key: str) -> np.ndarray | None:
        path = self.feature_dir / f"features_{key}.npy"
        if path.exists():
            arr = np.load(path)
            self.logger.info("CACHE HIT | ESM/simple features | key=%s | shape=%s", key, arr.shape)
            return arr
        self.logger.info("CACHE MISS | ESM/simple features | key=%s", key)
        return None

    def save_features(self, key: str, features: np.ndarray) -> None:
        path = self.feature_dir / f"features_{key}.npy"
        np.save(path, features)
        self.logger.info("CACHE SAVE | features | key=%s | shape=%s | path=%s", key, features.shape, path)

    def model_key(self, train_feature_key: str, model_cfg: dict[str, Any], train_meta: dict[str, Any]) -> str:
        return stable_hash({"train_feature_key": train_feature_key, "model_cfg": model_cfg, "train_meta": train_meta})

    def load_model_bundle(self, key: str) -> dict[str, Any] | None:
        path = self.model_dir / f"model_bundle_{key}.joblib"
        if path.exists():
            self.logger.info("CACHE HIT | trained model bundle | key=%s", key)
            return joblib.load(path)
        self.logger.info("CACHE MISS | trained model bundle | key=%s", key)
        return None

    def save_model_bundle(self, key: str, bundle: dict[str, Any]) -> None:
        path = self.model_dir / f"model_bundle_{key}.joblib"
        joblib.dump(bundle, path)
        self.logger.info("CACHE SAVE | trained model bundle | key=%s | path=%s", key, path)

    def load_metrics(self, key: str) -> dict[str, Any] | None:
        path = self.metric_dir / f"metrics_{key}.json"
        if path.exists():
            self.logger.info("CACHE HIT | metrics | key=%s", key)
            return json.loads(path.read_text(encoding="utf-8"))
        return None

    def save_metrics(self, key: str, metrics: dict[str, Any]) -> None:
        path = self.metric_dir / f"metrics_{key}.json"
        path.write_text(json.dumps(ensure_jsonable(metrics), indent=2, ensure_ascii=False), encoding="utf-8")
        self.logger.info(
            "CACHE SAVE | metrics | key=%s | brightness_r2=%s | low_brightness_auc=%s | path=%s",
            key,
            metrics.get("brightness_r2"),
            metrics.get("low_brightness_auc"),
            path,
        )

    def get_or_compute_features(self, seqs: list[str], feature_cfg: dict[str, Any], compute_fn: Callable[[], np.ndarray]) -> tuple[np.ndarray, str]:
        key = self.feature_key(seqs, feature_cfg)
        cached = self.load_features(key)
        if cached is not None:
            return cached, key
        features = compute_fn()
        self.save_features(key, features)
        return features, key
