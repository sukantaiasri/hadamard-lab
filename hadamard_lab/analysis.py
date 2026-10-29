"""Verification, normalisation, export and descriptive properties of Hadamard matrices."""
from __future__ import annotations

import io
import json
import math
from dataclasses import dataclass, asdict
from typing import Dict, List

import numpy as np
import pandas as pd

FORMATS = {
    "pm1": "+1 / −1",
    "sign": "+ / −",
    "binary": "1 / 0",
}
NORMALISATIONS = {
    "none": "As constructed",
    "semi": "Semi-normalised (first row all +1)",
    "full": "Normalised (first row and column all +1)",
}


# --------------------------------------------------------------------------
def normalise(H: np.ndarray, mode: str) -> np.ndarray:
    """Equivalence-preserving sign changes: negating rows/columns keeps H Hadamard."""
    H = H.astype(np.int8).copy()
    if mode == "none":
        return H
    if mode not in ("semi", "full"):
        raise ValueError(mode)
    H *= H[0:1, :]                      # make first row +1 (negate offending columns)
    if mode == "full":
        H *= H[:, 0:1]                  # make first column +1 (negate offending rows)
    return H


def is_semi_normalised(H) -> bool:
    return bool((H[0] == 1).all() or (H[:, 0] == 1).all())


def is_normalised(H) -> bool:
    return bool((H[0] == 1).all() and (H[:, 0] == 1).all())


# --------------------------------------------------------------------------
@dataclass
class Verification:
    order: int
    entries_pm1: bool
    row_orthogonal: bool
    column_orthogonal: bool
    max_row_deviation: int
    max_col_deviation: int
    passed: bool
    log10_determinant: float          # = (n/2) log10 n, attained by Hadamard matrices

    def as_dict(self) -> dict:
        return asdict(self)


def verify(H: np.ndarray) -> Verification:
    n = H.shape[0]
    pm1 = bool(H.shape == (n, n) and np.isin(H, (-1, 1)).all())
    Hf = H.astype(np.float64)
    I = n * np.eye(n)
    dr = int(np.abs(Hf @ Hf.T - I).max())
    dc = int(np.abs(Hf.T @ Hf - I).max())
    return Verification(n, pm1, dr == 0, dc == 0, dr, dc, pm1 and dr == 0 and dc == 0,
                        (n / 2) * math.log10(n) if n > 0 else 0.0)


# --------------------------------------------------------------------------
@dataclass
class Properties:
    order: int
    plus_entries: int
    minus_entries: int
    excess: int
    max_excess_bound: float
    symmetric: bool
    skew_type: bool
    constant_diagonal: bool
    normalised: bool
    design_parameters: str

    def as_rows(self) -> List[tuple]:
        return [
            ("Order n", self.order),
            ("Entries +1 / −1", f"{self.plus_entries} / {self.minus_entries}"),
            ("Excess (sum of entries)", self.excess),
            ("Excess bound n^{3/2}", f"{self.max_excess_bound:.1f}"),
            ("Symmetric", "yes" if self.symmetric else "no"),
            ("Skew type (H + Hᵀ = 2I)", "yes" if self.skew_type else "no"),
            ("Constant diagonal", "yes" if self.constant_diagonal else "no"),
            ("Normalised (row & column)", "yes" if self.normalised else "no"),
            ("Associated block design", self.design_parameters),
        ]


def properties(H: np.ndarray) -> Properties:
    n = H.shape[0]
    plus = int((H == 1).sum())
    I = np.eye(n, dtype=np.int8)
    skew = n > 1 and bool(np.array_equal(H.astype(int) + H.T.astype(int), 2 * np.eye(n, dtype=int)))
    if n >= 8 and n % 4 == 0:
        t = n // 4
        design = f"2-({4*t-1}, {2*t-1}, {t-1}) design (from the normalised form)"
    elif n == 4:
        design = "2-(3, 1, 0) (trivial)"
    else:
        design = "—"
    return Properties(n, plus, n * n - plus, int(H.astype(int).sum()), float(n ** 1.5),
                      bool(np.array_equal(H, H.T)), skew, bool(len(set(np.diag(H).tolist())) == 1),
                      is_normalised(H), design)


def gram(H: np.ndarray) -> np.ndarray:
    Hf = H.astype(np.float64)
    return Hf @ Hf.T


def row_sums(H: np.ndarray) -> np.ndarray:
    return H.astype(int).sum(axis=1)


# --------------------------------------------------------------------------
def to_text(H: np.ndarray, fmt: str = "pm1") -> str:
    if fmt == "pm1":
        w = 2
        return "\n".join(" ".join(f"{v:>{w}d}" for v in row) for row in H.astype(int))
    if fmt == "sign":
        return "\n".join(" ".join("+" if v == 1 else "-" for v in row) for row in H)
    if fmt == "binary":
        return "\n".join(" ".join("1" if v == 1 else "0" for v in row) for row in H)
    raise ValueError(fmt)


def to_csv(H: np.ndarray, fmt: str = "pm1") -> str:
    buf = io.StringIO()
    for row in H:
        if fmt == "pm1":
            buf.write(",".join(str(int(v)) for v in row))
        elif fmt == "sign":
            buf.write(",".join("+" if v == 1 else "-" for v in row))
        elif fmt == "binary":
            buf.write(",".join("1" if v == 1 else "0" for v in row))
        else:
            raise ValueError(fmt)
        buf.write("\n")
    return buf.getvalue()


def metadata_json(meta: Dict) -> str:
    return json.dumps(meta, indent=2, ensure_ascii=False, default=str)


# --------------------------------------------------------------------------
def plackett_burman(H: np.ndarray) -> pd.DataFrame:
    """Two-level orthogonal main-effect design with N = n runs and n-1 factors.

    Take the form of H whose first column is +1 and delete that column; every
    remaining column is balanced and pairwise orthogonal (Plackett & Burman, 1946).
    """
    Hn = normalise(H, "full")
    D = Hn[:, 1:]
    n = D.shape[1]
    return pd.DataFrame(D.astype(int), index=pd.RangeIndex(1, D.shape[0] + 1, name="Run"),
                        columns=[f"X{i}" for i in range(1, n + 1)])


def design_is_orthogonal(D: pd.DataFrame) -> bool:
    X = D.values.astype(float)
    return bool(np.array_equal(X.T @ X, len(X) * np.eye(X.shape[1])))
