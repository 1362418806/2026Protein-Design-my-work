from __future__ import annotations

import logging
import time
from typing import Any

import numpy as np

AA_ORDER = "ACDEFGHIKLMNPQRSTVWY"
AA_TO_IDX = {a: i for i, a in enumerate(AA_ORDER)}

HYDRO = {
    "A": 1.8, "C": 2.5, "D": -3.5, "E": -3.5, "F": 2.8, "G": -0.4,
    "H": -3.2, "I": 4.5, "K": -3.9, "L": 3.8, "M": 1.9, "N": -3.5,
    "P": -1.6, "Q": -3.5, "R": -4.5, "S": -0.8, "T": -0.7, "V": 4.2,
    "W": -0.9, "Y": -1.3,
}
CHARGE = {"D": -1, "E": -1, "K": 1, "R": 1, "H": 0.1}


def simple_features(seqs: list[str]) -> np.ndarray:
    rows = []
    for seq in seqs:
        n = max(1, len(seq))
        comp = np.zeros(len(AA_ORDER), dtype=np.float32)
        for aa in seq:
            if aa in AA_TO_IDX:
                comp[AA_TO_IDX[aa]] += 1.0
        comp /= n
        hyd = np.mean([HYDRO.get(a, 0.0) for a in seq])
        charge = np.sum([CHARGE.get(a, 0.0) for a in seq]) / n
        aromatic = sum(a in "FWY" for a in seq) / n
        polar = sum(a in "STNQ" for a in seq) / n
        pro_gly = sum(a in "PG" for a in seq) / n
        rows.append(np.concatenate([comp, np.array([n, hyd, charge, aromatic, polar, pro_gly], dtype=np.float32)]))
    return np.vstack(rows).astype(np.float32)


def esm_features(seqs: list[str], model_name: str, batch_size: int = 16, device: str | None = None, logger: logging.Logger | None = None) -> np.ndarray:
    logger = logger or logging.getLogger(__name__)
    t0 = time.time()
    import torch
    try:
        from esm import pretrained as esm_pretrained
    except Exception as exc:
        raise ImportError("FAIR ESM is required for --feature-mode esm. Install with `pip install fair-esm`.") from exc
    logger.info("ESM LOAD START | model=%s", model_name)
    model, alphabet = esm_pretrained.load_model_and_alphabet(model_name)
    batch_converter = alphabet.get_batch_converter()
    if device is None:
        device = "cuda" if torch.cuda.is_available() else "cpu"
    model = model.eval().to(device)
    logger.info("ESM LOAD DONE | model=%s | device=%s | seconds=%.2f", model_name, device, time.time() - t0)
    reps = []
    with torch.no_grad():
        for start in range(0, len(seqs), batch_size):
            batch = seqs[start:start + batch_size]
            labels = [(f"seq_{start+i}", s) for i, s in enumerate(batch)]
            _, _, tokens = batch_converter(labels)
            tokens = tokens.to(device)
            out = model(tokens, repr_layers=[model.num_layers], return_contacts=False)
            token_reps = out["representations"][model.num_layers]
            for i, seq in enumerate(batch):
                reps.append(token_reps[i, 1:len(seq)+1].mean(0).detach().cpu().numpy())
            if (start // batch_size) % 25 == 0:
                logger.info("ESM PROGRESS | %d/%d sequences", min(start + batch_size, len(seqs)), len(seqs))
    arr = np.vstack(reps).astype(np.float32)
    logger.info("ESM FEATURE DONE | shape=%s | seconds=%.2f", arr.shape, time.time() - t0)
    return arr


def featurize(seqs: list[str], feature_cfg: dict[str, Any], logger: logging.Logger | None = None) -> np.ndarray:
    mode = feature_cfg.get("mode", "simple")
    if mode == "simple":
        return simple_features(seqs)
    if mode == "esm":
        return esm_features(
            seqs,
            model_name=feature_cfg.get("esm_model", "esm2_t30_150M_UR50D"),
            batch_size=int(feature_cfg.get("batch_size", 16)),
            device=feature_cfg.get("device"),
            logger=logger,
        )
    raise ValueError(f"Unsupported feature mode: {mode}")
