"""Readers for RLM entropy output.

Full format (RLM ``-s entropy``/``all``, file given by ``-e``): 23 tab-separated
columns with a ``#chr`` header::

    #chr start end entropy epipolymorphism gggg gggG ... GGGG mean_methylation coverage

Legacy format (what v1/v2 consumed): 5 headerless columns::

    chrom start end entropy methylation

The legacy format has no read depth or epiallele counts, so only the pooled
null at an assumed depth can be computed for it.
"""

from __future__ import annotations

import gzip
from typing import Iterator

import numpy as np
import pandas as pd

from .core import EPIALLELE_NAMES

RLM_COLUMNS = ["chrom", "start", "end", "rlm_entropy", "epipolymorphism",
               *EPIALLELE_NAMES, "rlm_mean_methylation", "coverage"]
LEGACY_COLUMNS = ["chrom", "start", "end", "rlm_entropy", "rlm_mean_methylation"]


def _open_text(path):
    path = str(path)
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def detect_format(path) -> str:
    """Return 'rlm' or 'legacy' from the first non-empty line."""
    with _open_text(path) as handle:
        for line in handle:
            if not line.strip():
                continue
            if line.startswith("#chr") or line.startswith("chr\tstart"):
                return "rlm"
            n = len(line.rstrip("\n").split("\t"))
            if n == len(RLM_COLUMNS):
                return "rlm"
            if n == len(LEGACY_COLUMNS):
                return "legacy"
            raise ValueError(
                f"{path}: first line has {n} tab-separated columns; expected 23 "
                "(RLM entropy output) or 5 (legacy chrom/start/end/entropy/methylation)")
    raise ValueError(f"{path} is empty")


def read_chunks(path, fmt: str, chunksize: int) -> Iterator[pd.DataFrame]:
    names = RLM_COLUMNS if fmt == "rlm" else LEGACY_COLUMNS
    dtypes = {"chrom": str, "start": np.int64, "end": np.int64}
    reader = pd.read_csv(path, sep="\t", header=None, names=names, comment="#",
                         dtype=dtypes, chunksize=chunksize)
    for chunk in reader:
        yield chunk.reset_index(drop=True)


def epiallele_matrix(chunk: pd.DataFrame) -> np.ndarray:
    return chunk[EPIALLELE_NAMES].to_numpy(dtype=np.int64)
