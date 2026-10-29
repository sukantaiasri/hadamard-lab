"""Mutually orthogonal Latin squares, orthogonal arrays OA(s^2, s+1, s, 2) and lattice designs."""
from __future__ import annotations

from typing import List, Optional, Tuple

import numpy as np
import pandas as pd

from .gf import GF, field, is_prime_power, prime_power

S_LIMIT = 997                 # largest prime power below 1000
MAX_MATERIALISE = 20_000_000  # cells; larger arrays are streamed / sampled instead


def validate_s(s) -> Tuple[bool, str]:
    try:
        s = int(s)
    except (TypeError, ValueError):
        return False, "Enter an integer."
    if s < 2:
        return False, "The order must be at least 2."
    if s > S_LIMIT:
        return False, f"Orders above {S_LIMIT} are not supported."
    if not is_prime_power(s):
        nxt = next(q for q in range(s + 1, S_LIMIT + 2) if is_prime_power(q))
        prv = max(q for q in range(2, s) if is_prime_power(q))
        return False, (f"A complete set of MOLS from finite fields needs a prime power; "
                       f"{s} is not one (nearest: {prv} and {nxt}).")
    return True, ""


def latin_square(s: int, k: int) -> np.ndarray:
    """``L_k[i, j] = k*i + j`` in GF(s), for k = 1 .. s-1 (symbols 0 .. s-1)."""
    if not 1 <= k <= s - 1:
        raise ValueError("k must be in 1..s-1")
    f: GF = field(s)
    i = np.arange(s)
    ki = f.mul(np.full(s, k), i)                    # k*i for each row i
    return f.add(ki[:, None], i[None, :]).astype(np.int16)


def oa_columns(s: int) -> int:
    return s + 1


def oa_rows(s: int, start: int = 0, stop: Optional[int] = None) -> np.ndarray:
    """Rows ``start:stop`` of OA(s^2, s+1, s, 2); row r = (i, j) with r = i*s + j.

    Columns: i, j, then k*i + j for k = 1..s-1.
    """
    n = s * s
    stop = n if stop is None else min(stop, n)
    f = field(s)
    r = np.arange(start, stop)
    i, j = r // s, r % s
    cols = [i, j]
    if s > 2:
        K = np.arange(1, s)
        ki = f.mul(K[:, None], i[None, :])          # (s-1, rows)
        cols.extend(f.add(ki, j[None, :]))
    else:
        cols.append(f.add(f.mul(np.ones_like(i), i), j))
    return np.stack(cols, axis=1).astype(np.int16)


def oa_column(s: int, c: int) -> np.ndarray:
    """Column ``c`` (0 .. s) of OA(s^2, s+1, s, 2) over all s^2 rows, without building the array."""
    r = np.arange(s * s)
    i, j = r // s, r % s
    if c == 0:
        return i.astype(np.int16)
    if c == 1:
        return j.astype(np.int16)
    f = field(s)
    return f.add(f.mul(np.full_like(i, c - 1), i), j).astype(np.int16)


def oa_matrix(s: int) -> np.ndarray:
    if s * s * (s + 1) > MAX_MATERIALISE:
        raise MemoryError(f"OA({s}^2, {s}+1, {s}, 2) has {s*s*(s+1):,} cells; use oa_rows() in blocks.")
    return oa_rows(s)


# ---- verification --------------------------------------------------------
def is_latin(L: np.ndarray) -> bool:
    s = L.shape[0]
    full = np.arange(s)
    return bool(all((np.sort(L[i]) == full).all() for i in range(s)) and
                all((np.sort(L[:, j]) == full).all() for j in range(s)))


def are_orthogonal(L1: np.ndarray, L2: np.ndarray) -> bool:
    s = L1.shape[0]
    code = L1.astype(np.int64) * s + L2.astype(np.int64)
    return bool(np.bincount(code.ravel(), minlength=s * s).min() == 1 and code.max() < s * s)


def _pair_ok(a: np.ndarray, b: np.ndarray, s: int) -> bool:
    return bool((np.bincount(a.astype(np.int64) * s + b, minlength=s * s) == 1).all())


def verify_oa(s: int, max_pairs: Optional[int] = None, seed: int = 0) -> Tuple[bool, str]:
    """Check strength 2: every column pair shows each of the s^2 symbol pairs exactly once.

    Exhaustive for small s, random sample of column pairs for large s (stated in the message).
    """
    ncols = s + 1
    all_pairs = [(a, b) for a in range(ncols) for b in range(a + 1, ncols)]
    full_cost = len(all_pairs) * s * s
    exhaustive = max_pairs is None and full_cost <= 4e8
    if exhaustive:
        pairs = all_pairs
    else:
        rng = np.random.default_rng(seed)
        m = max_pairs or 300
        idx = rng.choice(len(all_pairs), size=min(m, len(all_pairs)), replace=False)
        pairs = [all_pairs[t] for t in idx]
    cache = {}
    col = lambda c: cache.setdefault(c, oa_column(s, c)) if len(cache) < 64 else oa_column(s, c)
    ok = all(_pair_ok(col(a), col(b), s) for a, b in pairs)
    scope = (f"all {len(pairs):,} column pairs" if exhaustive
             else f"{len(pairs):,} randomly sampled of {len(all_pairs):,} column pairs")
    return ok, f"strength-2 balance {'verified' if ok else 'FAILED'} for {scope}"


def verify_mols_set(s: int, sample: int = 200, seed: int = 0) -> Tuple[bool, str]:
    """Latin property of each square and pairwise orthogonality (exhaustive for small s)."""
    ks = list(range(1, s))
    pairs = [(a, b) for a in ks for b in ks if a < b]
    rng = np.random.default_rng(seed)
    exhaustive = len(ks) * s * s * 2 + len(pairs) * s * s <= 3e8
    if not exhaustive:
        ks = sorted(rng.choice(ks, size=min(sample, len(ks)), replace=False).tolist())
        pairs = [(ks[a], ks[b]) for a in range(len(ks)) for b in range(a + 1, min(len(ks), a + 6))]
    ok = all(is_latin(latin_square(s, k)) for k in ks)
    ok = ok and all(are_orthogonal(latin_square(s, a), latin_square(s, b)) for a, b in pairs)
    scope = (f"all {len(ks)} squares and {len(pairs):,} pairs" if exhaustive
             else f"{len(ks)} sampled squares and {len(pairs):,} sampled pairs")
    return ok, f"Latin and orthogonality properties {'verified' if ok else 'FAILED'} for {scope}"


def superimposed(s: int, ks: List[int]) -> List[List[str]]:
    Ls = [latin_square(s, k) for k in ks]
    return [["".join(f"{L[i, j]}" if s <= 10 else f"{L[i, j]}," for L in Ls).rstrip(",")
             for j in range(s)] for i in range(s)]


# ==========================================================================
# lattice designs
# ==========================================================================
def validate_lattice(v) -> Tuple[bool, str]:
    try:
        v = int(v)
    except (TypeError, ValueError):
        return False, "Enter an integer."
    k = int(round(v ** 0.5))
    if k * k != v:
        lo, hi = int(v ** 0.5) ** 2, (int(v ** 0.5) + 1) ** 2
        return False, f"{v} is not a perfect square (nearest: {lo} and {hi})."
    if k < 3:
        return False, "Use at least 9 treatments (k ≥ 3)."
    if k > 31:
        return False, "At most 961 treatments (k ≤ 31) are supported."
    return True, ""


def max_replicates(k: int) -> int:
    return k + 1 if is_prime_power(k) else 3


def lattice_blocks(k: int, r: int) -> List[List[List[int]]]:
    """``blocks[rep][block] = [treatments]`` (treatments numbered 1 .. k^2)."""
    if not 2 <= r <= max_replicates(k):
        raise ValueError(f"replicates must be between 2 and {max_replicates(k)} for k = {k}")
    T = np.arange(k * k).reshape(k, k) + 1
    reps = [[T[i, :].tolist() for i in range(k)],            # rows
            [T[:, j].tolist() for j in range(k)]]            # columns
    for c in range(1, r - 1):                                 # each MOLS-class
        if is_prime_power(k):
            L = latin_square(k, c) if k > 2 else (np.arange(k)[:, None] + np.arange(k)[None, :]) % k
        else:
            L = (np.arange(k)[:, None] + np.arange(k)[None, :]) % k
        reps.append([T[L == sym].tolist() for sym in range(k)])
    return reps[:r]


def lattice_incidence(k: int, r: int) -> np.ndarray:
    reps = lattice_blocks(k, r)
    v = k * k
    rows = []
    for rep in reps:
        for blk in rep:
            row = np.zeros(v, dtype=np.int8)
            row[np.array(blk) - 1] = 1
            rows.append(row)
    return np.array(rows).T            # v x (r*k)


def lattice_field_book(k: int, r: int) -> pd.DataFrame:
    recs = []
    for ri, rep in enumerate(lattice_blocks(k, r), 1):
        for bi, blk in enumerate(rep, 1):
            for pi, t in enumerate(blk, 1):
                recs.append((ri, bi, pi, t))
    return pd.DataFrame(recs, columns=["Replicate", "Block", "Plot", "Treatment"])


def lattice_summary(k: int, r: int) -> dict:
    N = lattice_incidence(k, r).astype(np.int64)
    v = k * k
    conc = N @ N.T
    off = conc[~np.eye(v, dtype=bool)]
    lam_counts = {int(a): int(b) for a, b in zip(*np.unique(off, return_counts=True))}
    summary = {
        "treatments": v, "block_size": k, "replicates": r, "blocks": r * k,
        "plots": v * r, "concurrence_counts": lam_counts,
        "each_treatment_once_per_replicate": bool((N.sum(axis=1) == r).all()),
        "balanced": len(lam_counts) == 1,
    }
    Cmat = r * np.eye(v) - conc / k                    # information matrix
    w = np.linalg.eigvalsh(Cmat)
    nz = w[w > 1e-8]
    eff = nz / r
    summary["average_efficiency"] = float(len(eff) / np.sum(1.0 / eff)) if len(eff) == v - 1 else float("nan")
    summary["canonical_efficiencies"] = {round(float(a), 6): int(b) for a, b in
                                         zip(*np.unique(np.round(eff, 6), return_counts=True))}
    return summary
