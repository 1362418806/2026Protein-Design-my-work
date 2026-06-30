from __future__ import annotations

import re
from pathlib import Path
from typing import Iterable

import pandas as pd

AA3_TO_AA1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def parse_fasta_like(path: str | Path) -> dict[str, str]:
    text = Path(path).read_text(encoding="utf-8", errors="ignore")
    refs: dict[str, str] = {}
    name = None
    chunks: list[str] = []
    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name and chunks:
                refs[name] = "".join(chunks).replace(" ", "").upper()
            name = line[1:].strip().split()[0]
            chunks = []
        elif re.fullmatch(r"[A-Za-z*\-\s]+", line):
            chunks.append(re.sub(r"[^A-Za-z]", "", line).upper())
    if name and chunks:
        refs[name] = "".join(chunks).replace(" ", "").upper()
    if not refs:
        # Support simple name sequence tables.
        for raw in text.splitlines():
            parts = re.split(r"[,\t ]+", raw.strip())
            if len(parts) >= 2 and re.fullmatch(r"[A-Za-z]+", parts[-1]):
                refs[parts[0]] = parts[-1].upper()
    return refs


def load_exclusion(path: str | Path) -> set[str]:
    p = Path(path)
    if not p.exists():
        return set()
    if p.suffix.lower() in {".csv", ".tsv"}:
        sep = "\t" if p.suffix.lower() == ".tsv" else ","
        df = pd.read_csv(p, sep=sep)
        seq_col = next((c for c in df.columns if "seq" in c.lower()), df.columns[-1])
        return set(df[seq_col].astype(str).str.strip().str.upper())
    return {x.strip().upper() for x in p.read_text(errors="ignore").splitlines() if x.strip()}


def apply_mutations(scaffold: str, mutation_text: str) -> str | None:
    seq = list(scaffold)
    if not isinstance(mutation_text, str) or not mutation_text.strip():
        return scaffold
    # Supports S30R, S30 R, Ser30Arg-like text fragments.
    pattern = re.compile(r"([ACDEFGHIKLMNPQRSTVWY])\s*(\d{1,3})\s*([ACDEFGHIKLMNPQRSTVWY])")
    matches = pattern.findall(mutation_text.upper())
    if not matches:
        return scaffold
    for wt, pos, mut in matches:
        idx = int(pos) - 1
        if idx < 0 or idx >= len(seq):
            return None
        # The official table can mix homologous numbering; keep reconstruction tolerant.
        seq[idx] = mut
    return "".join(seq)


def load_training_table(gfp_data_path: str | Path, scaffold: str, max_rows: int | None = None, seed: int = 42) -> tuple[pd.DataFrame, dict]:
    p = Path(gfp_data_path)
    df = pd.read_excel(p) if p.suffix.lower() in {".xlsx", ".xls"} else pd.read_csv(p)
    brightness_col = next((c for c in df.columns if "brightness" in c.lower() or "fluores" in c.lower()), None)
    sequence_col = next((c for c in df.columns if "sequence" in c.lower() or c.lower() in {"seq", "aa_seq"}), None)
    mutation_col = next((c for c in df.columns if "mutation" in c.lower() or "aamut" in c.lower()), None)
    gfp_type_col = next((c for c in df.columns if "gfp type" in c.lower() or "type" == c.lower()), None)
    if brightness_col is None:
        raise ValueError(f"Cannot find brightness/fluorescence column in {p}")
    if sequence_col is not None:
        df["full_sequence"] = df[sequence_col].astype(str).str.strip().str.upper()
    elif mutation_col is not None:
        df["full_sequence"] = [apply_mutations(scaffold, str(x)) for x in df[mutation_col]]
    else:
        raise ValueError("Need either sequence column or mutation column to reconstruct training sequences.")
    df = df.dropna(subset=["full_sequence", brightness_col]).copy()
    df = df[df["full_sequence"].str.len().between(200, 260)]
    if max_rows and len(df) > max_rows:
        df = df.sample(max_rows, random_state=seed).reset_index(drop=True)
    meta = {
        "brightness_col": brightness_col,
        "sequence_col": sequence_col,
        "mutation_col": mutation_col,
        "gfp_type_col": gfp_type_col,
        "n_training_rows": int(len(df)),
    }
    return df, meta
