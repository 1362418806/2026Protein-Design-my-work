from __future__ import annotations

import hashlib
import json
import logging
import random
from pathlib import Path
from typing import Any, Iterable

import numpy as np

AA20 = set("ACDEFGHIKLMNPQRSTVWY")


def setup_logger(out_dir: str | Path, name: str = "synbio_gfp_v3") -> logging.Logger:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    fmt = logging.Formatter("%(asctime)s | %(levelname)s | %(message)s")
    stream = logging.StreamHandler()
    stream.setFormatter(fmt)
    file_handler = logging.FileHandler(out / "run.log", encoding="utf-8")
    file_handler.setFormatter(fmt)
    logger.addHandler(stream)
    logger.addHandler(file_handler)
    return logger


def stable_hash(obj: Any) -> str:
    payload = json.dumps(obj, sort_keys=True, ensure_ascii=False, default=str).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()[:16]


def sequence_hash(seqs: Iterable[str]) -> str:
    h = hashlib.sha256()
    for seq in seqs:
        h.update(seq.encode("utf-8"))
        h.update(b"\n")
    return h.hexdigest()[:16]


def validate_sequence(seq: str, min_len: int = 220, max_len: int = 250) -> dict[str, bool]:
    return {
        "valid_length": min_len <= len(seq) <= max_len,
        "starts_with_M": seq.startswith("M"),
        "valid_aa20": all(c in AA20 for c in seq),
    }


def set_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)


def ensure_jsonable(value: Any) -> Any:
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, dict):
        return {str(k): ensure_jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [ensure_jsonable(v) for v in value]
    return value
