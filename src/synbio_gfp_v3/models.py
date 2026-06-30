from __future__ import annotations

import logging
from typing import Any

import numpy as np
from sklearn.ensemble import ExtraTreesClassifier, ExtraTreesRegressor, RandomForestRegressor
from sklearn.metrics import roc_auc_score, r2_score
from sklearn.model_selection import train_test_split
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge


def train_predictive_models(X: np.ndarray, y: np.ndarray, config: dict[str, Any], meta: dict[str, Any], logger: logging.Logger) -> dict[str, Any]:
    seed = int(config.get("seed", 42))
    threshold_q = float(config.get("training", {}).get("low_brightness_quantile", 0.20))
    low_threshold = float(np.quantile(y, threshold_q))
    y_low = (y <= low_threshold).astype(int)
    X_train, X_val, y_train, y_val, low_train, low_val = train_test_split(
        X, y, y_low, test_size=0.2, random_state=seed, stratify=y_low if len(set(y_low)) > 1 else None
    )
    reg_rf = RandomForestRegressor(n_estimators=260, max_features="sqrt", random_state=seed, n_jobs=-1)
    reg_et = ExtraTreesRegressor(n_estimators=260, max_features="sqrt", random_state=seed + 1, n_jobs=-1)
    reg_ridge = make_pipeline(StandardScaler(with_mean=False), Ridge(alpha=10.0, random_state=seed))
    clf_low = ExtraTreesClassifier(n_estimators=260, max_features="sqrt", random_state=seed + 2, n_jobs=-1, class_weight="balanced")
    logger.info("MODEL TRAIN START | n_train=%d | n_val=%d | dim=%d", len(X_train), len(X_val), X.shape[1])
    reg_rf.fit(X_train, y_train)
    reg_et.fit(X_train, y_train)
    reg_ridge.fit(X_train, y_train)
    clf_low.fit(X_train, low_train)
    pred_rf = reg_rf.predict(X_val)
    pred_et = reg_et.predict(X_val)
    pred_ridge = reg_ridge.predict(X_val)
    pred_mean = (pred_rf + pred_et + pred_ridge) / 3.0
    low_prob = clf_low.predict_proba(X_val)[:, 1]
    metrics = {
        "brightness_r2": float(r2_score(y_val, pred_mean)),
        "brightness_r2_rf": float(r2_score(y_val, pred_rf)),
        "brightness_r2_et": float(r2_score(y_val, pred_et)),
        "brightness_r2_ridge": float(r2_score(y_val, pred_ridge)),
        "low_brightness_auc": float(roc_auc_score(low_val, low_prob)) if len(set(low_val)) > 1 else None,
        "low_brightness_proxy_threshold": low_threshold,
        "n_train": int(len(X_train)),
        "n_val": int(len(X_val)),
        **meta,
    }
    logger.info("MODEL TRAIN DONE | brightness_r2=%.4f | low_auc=%s | threshold=%.4f", metrics["brightness_r2"], metrics["low_brightness_auc"], low_threshold)
    return {"regressors": [reg_rf, reg_et, reg_ridge], "low_classifier": clf_low, "metrics": metrics}


def predict_bundle(bundle: dict[str, Any], X: np.ndarray) -> dict[str, np.ndarray]:
    preds = np.vstack([m.predict(X) for m in bundle["regressors"]])
    brightness = preds.mean(axis=0)
    uncertainty = preds.std(axis=0)
    low_risk = bundle["low_classifier"].predict_proba(X)[:, 1]
    return {"brightness_pred": brightness, "brightness_uncertainty": uncertainty, "low_brightness_risk": low_risk}
