"""Constructions of Hadamard matrices.

Every construction is a small, self-contained function that returns an
``int8`` matrix with entries +-1.  The *planner* (see :mod:`planner`) decides
which construction applies to a requested order; this module only knows how to
carry one out and how to describe it.

Notation follows the classical literature (Hedayat & Wallis 1978,
Seberry & Yamada 1992, Craigen 1996, Seberry & Yamada 2020).

Version 2.0 adds the sequence-based and difference-set constructions implemented in
:mod:`sequences` (Golay pairs, Legendre pairs and quadruples, cyclic Hadamard difference
sets, Scarpis–Đoković) and the data behind the orders 668, 716, 876 and 892.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np

from .gf import field as gf_field, is_prime, is_prime_power, modulus_string, prime_power
from . import sequences as SQ

DATA_FILE = Path(__file__).with_name("data") / "sequences.json"

# Baumert-Hall array of order 12 (Hall 1986): an orthogonal design OD(12; 3,3,3,3).
BAUMERT_HALL_12 = (
    'A A A B -B C -C -D B C -D -D',
    'A -A B -A -B -D D -C -B -D -C -C',
    'A -B -A A -D D -B B -C -D C -C',
    'B A -A -A D D D C C -B -B -C',
    'B -D D D A A A C -C B -C B',
    'B C -D D A -A C -A -D C B -B',
    'D -C B -B A -C -A A B C D -D',
    '-C -D -C -D C A -A -A -D B -B -B',
    'D -C -B -B -B C C -D A A A D',
    '-D -B C C C B B -D A -A D -A',
    'C -B -C C D -B -D -B A -D -A A',
    '-C -D -D C -C -B B B D A -A -A',
)


# ==========================================================================
# data access
# ==========================================================================
@lru_cache(maxsize=1)
def sequence_data() -> dict:
    with open(DATA_FILE, "r", encoding="utf-8") as fh:
        raw = json.load(fh)
    return {
        "williamson": {int(k): v for k, v in raw["williamson"].items()},
        "base": {int(k): v for k, v in raw["base"].items()},
        "turyn": {int(k): v for k, v in raw["turyn"].items()},
        "sds": {int(k): v for k, v in raw["sds"].items()},
        "tseq": {int(k): v for k, v in raw.get("tseq", {}).items()},
        "lpairs": {int(k): v for k, v in raw.get("legendre_pairs", {}).items()},
        "lquads": {int(k): v for k, v in raw.get("legendre_quadruples", {}).items()},
        "meta": raw["_meta"],
    }


T_LENGTH_LIMIT = 250        # T-sequences longer than this are never needed for n <= 1000


def t_sequence_sources() -> Dict[int, List[Tuple[str, int]]]:
    """``length -> [(source, key), ...]`` for every way of obtaining T-sequences.

    Sources: ``base`` (stored base sequences), ``turyn`` (stored Turyn-type sequences),
    ``tseq`` (T-sequences stored directly), ``golay`` (``((a+b)/2, (a-b)/2, 0, 0)`` of a Golay
    pair of length g) and ``golay2`` (length ``2g+1`` from the base sequences ``(a,1),(a,-1),b,b``).
    """
    d = sequence_data()
    out: Dict[int, List[Tuple[str, int]]] = {}
    for n, seqs in d["base"].items():
        p = len(seqs[0]) - len(seqs[2])
        out.setdefault(2 * n + p, []).append(("base", n))
    for t in d["turyn"]:
        out.setdefault(3 * t - 1, []).append(("turyn", t))
    for t in d["tseq"]:
        out.setdefault(t, []).append(("tseq", t))
    for g in SQ.golay_lengths(T_LENGTH_LIMIT):
        out.setdefault(g, []).append(("golay", g))
        if 2 * g + 1 <= T_LENGTH_LIMIT:
            out.setdefault(2 * g + 1, []).append(("golay2", g))
    return dict(sorted(out.items()))


def describe_t_source(source: str, key: int) -> str:
    return {"base": f"base sequences, m = {key}", "turyn": f"Turyn-type sequences, t = {key}",
            "tseq": f"stored T-sequences, t = {key}", "golay": f"Golay pair of length {key}",
            "golay2": f"Golay pair of length {key} (Turyn: length 2g + 1)"}[source]


def t_sequence_lengths() -> Dict[int, str]:
    """Lengths of the T-sequences obtainable from any source, with the preferred source."""
    return {L: describe_t_source(*srcs[0]) for L, srcs in t_sequence_sources().items()}


# ==========================================================================
# elementary matrix tools
# ==========================================================================
def circulant(first_row: Sequence[int]) -> np.ndarray:
    """Circulant matrix with the given first row: ``C[i, j] = row[(j - i) mod n]``."""
    r = np.asarray(first_row, dtype=np.int32)
    n = len(r)
    idx = (np.arange(n)[None, :] - np.arange(n)[:, None]) % n
    return r[idx]


def back_identity(n: int) -> np.ndarray:
    return np.eye(n, dtype=np.int32)[::-1]


def _i8(m: np.ndarray) -> np.ndarray:
    return np.asarray(m, dtype=np.int8)


def is_hadamard(h: np.ndarray) -> bool:
    n = h.shape[0]
    if h.shape != (n, n) or not np.isin(h, (-1, 1)).all():
        return False
    g = h.astype(np.float64) @ h.T.astype(np.float64)
    return bool(np.array_equal(g, n * np.eye(n)))


def aperiodic_autocorrelation(a: Sequence[int]) -> np.ndarray:
    a = np.asarray(a, dtype=np.int64)
    n = len(a)
    return np.array([int(np.dot(a[: n - s], a[s:])) for s in range(n)])


def periodic_autocorrelation(a: Sequence[int]) -> np.ndarray:
    a = np.asarray(a, dtype=np.int64)
    return np.array([int(np.dot(a, np.roll(a, -s))) for s in range(len(a))])


# ==========================================================================
# arrays
# ==========================================================================
def williamson_array(A, B, C, D) -> np.ndarray:
    return np.block([[A, B, C, D],
                     [-B, A, -D, C],
                     [-C, D, A, -B],
                     [-D, -C, B, A]])


def goethals_seidel_array(A, B, C, D) -> np.ndarray:
    """Goethals-Seidel array with ``R`` the back-diagonal identity."""
    R = back_identity(A.shape[0])
    return np.block([[A, B @ R, C @ R, D @ R],
                     [-B @ R, A, D.T @ R, -C.T @ R],
                     [-C @ R, -D.T @ R, A, B.T @ R],
                     [-D @ R, C.T @ R, -B.T @ R, A]])


def baumert_hall_array(A, B, C, D) -> np.ndarray:
    """Baumert-Hall array of order 12 filled with four commuting +-1 matrices."""
    blocks = {"A": A, "B": B, "C": C, "D": D}
    rows = []
    for line in BAUMERT_HALL_12:
        row = []
        for tok in line.split():
            m = blocks[tok[-1]]
            row.append(-m if tok.startswith("-") else m)
        rows.append(row)
    return np.block(rows)


# ==========================================================================
# T-sequences
# ==========================================================================
def t_sequences_from_base(A, B, C, D) -> List[np.ndarray]:
    """Four T-sequences of length ``2n+p`` from base sequences of lengths n+p, n+p, n, n."""
    A, B, C, D = (np.asarray(x, dtype=np.int64) for x in (A, B, C, D))
    n, p = len(C), len(A) - len(C)
    z_n = np.zeros(n, dtype=np.int64)
    z_np = np.zeros(n + p, dtype=np.int64)
    t1 = np.concatenate([(A + B) // 2, z_n])
    t2 = np.concatenate([(A - B) // 2, z_n])
    t3 = np.concatenate([z_np, (C + D) // 2])
    t4 = np.concatenate([z_np, (C - D) // 2])
    return [t1, t2, t3, t4]


def base_from_turyn(X, Y, Z, W) -> List[np.ndarray]:
    """Base sequences (2n-1, 2n-1, n, n) from Turyn-type sequences (n, n, n, n-1)."""
    X, Y, Z, W = (np.asarray(x, dtype=np.int64) for x in (X, Y, Z, W))
    return [np.concatenate([Z, W]), np.concatenate([Z, -W]), X, Y]


def is_t_sequence_set(ts: Sequence[np.ndarray]) -> bool:
    ts = [np.asarray(t) for t in ts]
    n = len(ts[0])
    if any(len(t) != n for t in ts):
        return False
    support = sum(np.abs(t) for t in ts)
    if not (support == 1).all():
        return False
    tot = sum(aperiodic_autocorrelation(t) for t in ts)
    return bool((tot[1:] == 0).all())


def _gs_matrices_from_t(ts) -> List[np.ndarray]:
    T1, T2, T3, T4 = (circulant(t) for t in ts)
    return [T1 + T2 + T3 + T4, -T1 + T2 + T3 - T4, -T1 - T2 + T3 + T4, -T1 + T2 - T3 + T4]


def _t_sequences(source: str, key: int) -> List[np.ndarray]:
    d = sequence_data()
    if source == "base":
        return t_sequences_from_base(*d["base"][key])
    if source == "turyn":
        return t_sequences_from_base(*base_from_turyn(*d["turyn"][key]))
    if source == "tseq":
        return [np.asarray(t, dtype=np.int64) for t in d["tseq"][key]]
    if source == "golay":
        return SQ.t_sequences_from_golay(*SQ.golay_pair(key))
    if source == "golay2":
        return t_sequences_from_base(*SQ.base_sequences_from_golay(*SQ.golay_pair(key)))
    raise KeyError(source)


# ==========================================================================
# the constructions proper
# ==========================================================================
def sylvester(order: int) -> np.ndarray:
    k = order.bit_length() - 1
    if order != 2 ** k:
        raise ValueError("Sylvester construction needs a power of two")
    h = np.ones((1, 1), dtype=np.int32)
    h2 = np.array([[1, 1], [1, -1]], dtype=np.int32)
    for _ in range(k):
        h = np.kron(h, h2)
    return _i8(h)


def kronecker(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return _i8(np.kron(a.astype(np.int32), b.astype(np.int32)))


def conference_matrix(q: int) -> np.ndarray:
    """Conference matrix ``S = [[0, e'], [e, Q]]`` of order q+1 (q = 1 mod 4: symmetric,
    q = 3 mod 4 the same layout with ``-e`` gives the skew matrix, see :func:`paley1`)."""
    Q = gf_field(q).jacobsthal().astype(np.int32)
    S = np.zeros((q + 1, q + 1), dtype=np.int32)
    S[0, 1:] = 1
    S[1:, 0] = 1
    S[1:, 1:] = Q
    return S


def paley1(q: int) -> np.ndarray:
    """Paley I: order q+1 for a prime power ``q = 3 (mod 4)``; ``H = S + I`` is skew-type."""
    if q % 4 != 3 or not is_prime_power(q):
        raise ValueError("Paley I needs a prime power q = 3 (mod 4)")
    Q = gf_field(q).jacobsthal().astype(np.int32)
    S = np.zeros((q + 1, q + 1), dtype=np.int32)
    S[0, 1:] = 1
    S[1:, 0] = -1
    S[1:, 1:] = Q
    return _i8(S + np.eye(q + 1, dtype=np.int32))


def paley2_general(q: int, h: np.ndarray) -> np.ndarray:
    """Williamson/Paley II: from H_h (h even) and q = 1 (mod 4) get order h(q+1).
    ``H = A (x) I + B (x) S`` with ``B = U A``, ``U = I_{h/2} (x) [[0,1],[-1,0]]``."""
    if q % 4 != 1 or not is_prime_power(q):
        raise ValueError("Paley II needs a prime power q = 1 (mod 4)")
    hh = h.shape[0]
    if hh % 2:
        raise ValueError("h must be even")
    S = conference_matrix(q)
    A = h.astype(np.int32)
    U = np.kron(np.eye(hh // 2, dtype=np.int32), np.array([[0, 1], [-1, 0]], dtype=np.int32))
    B = U @ A
    return _i8(np.kron(A, np.eye(q + 1, dtype=np.int32)) + np.kron(B, S))


def paley2(q: int) -> np.ndarray:
    return paley2_general(q, np.array([[1, 1], [1, -1]], dtype=np.int8))


def williamson_matrices(w: int) -> List[np.ndarray]:
    rows = sequence_data()["williamson"][w]
    return [circulant(r) for r in rows]


def williamson(w: int) -> np.ndarray:
    A, B, C, D = williamson_matrices(w)
    return _i8(williamson_array(A, B, C, D))


def goethals_seidel_from_t(ts) -> np.ndarray:
    A, B, C, D = _gs_matrices_from_t(ts)
    return _i8(goethals_seidel_array(A, B, C, D))


def goethals_seidel_turyn(t: int) -> np.ndarray:
    return goethals_seidel_from_t(_t_sequences("turyn", t))


def goethals_seidel_base(n: int) -> np.ndarray:
    return goethals_seidel_from_t(_t_sequences("base", n))


def baumert_hall(w: int) -> np.ndarray:
    if w == 1:
        one = np.ones((1, 1), dtype=np.int32)
        A = B = C = D = one
    else:
        A, B, C, D = williamson_matrices(w)
    return _i8(baumert_hall_array(A, B, C, D))


def cooper_wallis(w: int, source: str, key: int) -> np.ndarray:
    """Cooper-Wallis: Williamson matrices of order w and T-sequences of length L.

    The four blocks are the rows of a coefficient array whose *columns* must be
    orthogonal.  A transcription with ``W = ... - T4 (x) A`` breaks that orthogonality; the
    last term must be ``+ T4 (x) A`` (checked in the tests).
    """
    A, B, C, D = williamson_matrices(w)
    T1, T2, T3, T4 = (circulant(t) for t in _t_sequences(source, key))
    k = np.kron
    X = k(T1, A) + k(T2, B) + k(T3, C) + k(T4, D)
    Y = -k(T1, B) + k(T2, A) + k(T3, D) - k(T4, C)
    Z = -k(T1, C) - k(T2, D) + k(T3, A) + k(T4, B)
    W = -k(T1, D) + k(T2, C) - k(T3, B) + k(T4, A)
    return _i8(goethals_seidel_array(X, Y, Z, W))


def ehlich(k: int) -> np.ndarray:
    """Ehlich: order (k-1)^2 from the skew-type Paley matrix of order k+1 (k = 3 mod 4)
    and a Paley II core of order m = k-2 = 1 (mod 4).

    ``H = [[1, e'], [e, K]]`` with ``K = A (x) G - I_m (x) (J_k - I_k) + J_m (x) I_k``.
    NOTE: the formula is sometimes transcribed with ``+ I_m (x) (J_k - I_k)``; with that sign
    ``H`` is *not* orthogonal (the row sums of K must equal -1 and the J(x)I
    coefficient of K K' must vanish, which forces the opposite sign).  The corrected
    form is verified for every admissible k <= 31 in the test-suite.
    """
    m = k - 2
    if k % 4 != 3 or not is_prime_power(k) or m % 4 != 1 or not is_prime_power(m):
        raise ValueError("Ehlich needs prime powers k = 3 (mod 4) and k-2 = 1 (mod 4)")
    A = gf_field(m).jacobsthal().astype(np.int32)          # AA' = mI - J
    G = gf_field(k).jacobsthal().astype(np.int32)          # skew, zero diagonal
    Ik, Im = np.eye(k, dtype=np.int32), np.eye(m, dtype=np.int32)
    Jk, Jm = np.ones((k, k), dtype=np.int32), np.ones((m, m), dtype=np.int32)
    K = np.kron(A, G) - np.kron(Im, Jk - Ik) + np.kron(Jm, Ik)
    N = m * k
    H = np.ones((N + 1, N + 1), dtype=np.int32)
    H[1:, 1:] = K
    return _i8(H)


# Sign table for Miyamoto's construction.  Entry (i, j) = (a, b, alpha, beta, gamma):
#   X_ij = [[alpha, beta e'], [gamma e, T_ij]],   T_ij = a U_ij (x) J2 + b V_ij (x) M2.
# The block skeleton U, V follows Miyamoto (1991); the border signs as commonly transcribed do
# not give an orthogonal matrix, so this table was determined by exhaustive
# row-by-row search at q = 5 and then verified for every admissible q <= 249 (see tests).
MIYAMOTO_SIGNS = (
    ((1, 1, 1, 1, 1), (1, 1, 1, 1, -1), (1, 1, 1, 1, -1), (1, 1, 1, 1, -1)),
    ((1, 1, 1, -1, -1), (1, 1, 1, 1, -1), (1, 1, -1, 1, 1), (1, 1, -1, -1, -1)),
    ((1, 1, 1, -1, 1), (1, 1, -1, 1, 1), (-1, 1, -1, -1, 1), (1, 1, 1, 1, -1)),
    ((1, 1, 1, -1, 1), (1, 1, -1, -1, -1), (1, 1, 1, 1, 1), (-1, 1, -1, 1, 1)),
)


def miyamoto(q: int, K: np.ndarray) -> np.ndarray:
    """Miyamoto: order 4q from a prime power q = 1 (mod 4) and an Hadamard matrix K of order q-1."""
    if q % 4 != 1 or not is_prime_power(q):
        raise ValueError("Miyamoto needs a prime power q = 1 (mod 4)")
    if K.shape[0] != q - 1:
        raise ValueError("K must have order q-1")
    m = (q - 1) // 2
    Q = gf_field(q).jacobsthal().astype(np.int32)
    P = np.array([y for y in range(1, q) if Q[0, y] == 1])      # nonzero squares
    Nn = np.array([y for y in range(1, q) if Q[0, y] == -1])    # non-squares
    C1 = -Q[np.ix_(P, P)]
    C2 = Q[np.ix_(P, Nn)]
    C4 = Q[np.ix_(Nn, Nn)]
    K = K.astype(np.int32)
    K1, K2 = K[:m, :m], K[:m, m:]
    K3, K4 = -K[m:, :m], K[m:, m:]
    Z = np.zeros((m, m), dtype=np.int32)
    I = np.eye(m, dtype=np.int32)
    U = np.block([[C1, C2, Z, Z],
                  [-C2.T, C4, Z, Z],
                  [Z, Z, C1, C2],
                  [Z, Z, -C2.T, C4]])
    V = np.block([[I, Z, K1, K2],
                  [Z, I, K3, -K4],
                  [-K1.T, -K3.T, I, Z],
                  [-K2.T, K4.T, Z, I]])
    J2 = np.ones((2, 2), dtype=np.int32)
    M2 = np.array([[1, -1], [-1, 1]], dtype=np.int32)
    ones_row = np.ones((1, 2 * m), dtype=np.int32)

    def X(i, j):
        a, b, al, be, ga = MIYAMOTO_SIGNS[i][j]
        sl = (slice(i * m, (i + 1) * m), slice(j * m, (j + 1) * m))
        T = a * np.kron(U[sl], J2) + b * np.kron(V[sl], M2)
        top = np.concatenate([np.array([[al]], dtype=np.int32), be * ones_row], axis=1)
        bot = np.concatenate([ga * ones_row.T, T], axis=1)
        return np.vstack([top, bot])

    return _i8(np.block([[X(i, j) for j in range(4)] for i in range(4)]))


def sds_circulant(n: int, subset: Sequence[int]) -> np.ndarray:
    """``A_S[x, y] = -1`` if ``y - x`` lies in S, ``+1`` otherwise (group Z_n)."""
    mask = np.zeros(n, dtype=bool)
    mask[list(subset)] = True
    idx = (np.arange(n)[None, :] - np.arange(n)[:, None]) % n
    return np.where(mask[idx], -1, 1).astype(np.int32)


def supplementary_difference_sets(n: int) -> np.ndarray:
    sets = sequence_data()["sds"][n]
    A, B, C, D = (sds_circulant(n, s) for s in sets)
    return _i8(goethals_seidel_array(A, B, C, D))


# ---- version 2.0: sequence and difference-set constructions (see sequences.py) ----
def goethals_seidel_tseq(source: str, key: int) -> np.ndarray:
    """Goethals–Seidel array from T-sequences of any source (Golay-derived or stored)."""
    return goethals_seidel_from_t(_t_sequences(source, key))


def golay_two_circulant(g: int) -> np.ndarray:
    """Order 2g from a Golay pair of length g: ``[[A, B], [-B', A']]``."""
    return SQ.golay_two_circulant(*SQ.golay_pair(g))


def legendre_pair_sources(limit: int = 1000) -> Dict[int, List[dict]]:
    return SQ.legendre_pair_sources(limit, sequence_data()["lpairs"])


def legendre_pair(spec: dict):
    return SQ.legendre_pair(spec, sequence_data()["lpairs"])


def legendre_two_circulant(spec: dict) -> np.ndarray:
    """Order 2l + 2 from a Legendre pair of length l (two-circulant core)."""
    return SQ.two_circulant_core(*legendre_pair(spec))


def legendre_quadruple(m: int) -> np.ndarray:
    """Order 4m + 4 from a stored Legendre quadruple of even length m (bordered Goethals–Seidel)."""
    return SQ.legendre_quadruple_matrix(sequence_data()["lquads"][m])


def difference_set_core(spec: dict) -> np.ndarray:
    """Order l + 1 from a cyclic/abelian Hadamard difference set of a non-Paley type."""
    return SQ.difference_set_core_matrix(spec)


def scarpis(q: int, A: np.ndarray) -> np.ndarray:
    """Order q(q+1) from a Hadamard matrix A of order q+1, q = 3 (mod 4) a prime power."""
    return SQ.scarpis(q, A)


# ==========================================================================
# registry: every method knows how to enumerate its parameter sets for an order
# ==========================================================================
@dataclass(frozen=True)
class Method:
    key: str
    name: str
    short: str                      # one-line description
    reference: str
    priority: int                   # lower = preferred by the automatic planner


METHODS: Dict[str, Method] = {m.key: m for m in [
    Method("trivial", "Trivial", "Orders 1 and 2.", "", 0),
    Method("sylvester", "Sylvester (Kronecker power of H2)",
           "n = 2^k from the k-fold Kronecker power of the order-2 matrix.",
           "Sylvester (1867)", 1),
    Method("paley1", "Paley I",
           "n = q + 1 with q a prime power, q ≡ 3 (mod 4). Skew-Hadamard type.",
           "Paley (1933)", 2),
    Method("paley2", "Paley II",
           "n = 2(q + 1) with q a prime power, q ≡ 1 (mod 4).",
           "Paley (1933)", 3),
    Method("golay2c", "Golay pair – two circulants",
           "n = 2g from a Golay complementary pair of length g = 2^a·10^b·26^c in [[A, B], [−Bᵀ, Aᵀ]].",
           "Golay (1961); Turyn (1974); Borwein & Ferguson (2004)", 4),
    Method("williamson", "Williamson",
           "n = 4w from four symmetric circulant matrices of order w in the Williamson array.",
           "Williamson (1944)", 5),
    Method("gs_turyn", "Goethals–Seidel (Turyn-type sequences)",
           "n = 4(3t − 1) from Turyn-type sequences of length t via base and T-sequences.",
           "Goethals & Seidel (1967); Turyn (1974); London & Kotsireas (2025)", 6),
    Method("gs_base", "Goethals–Seidel (base sequences)",
           "n = 4(2m + p) from base sequences of lengths m+p, m+p, m, m via T-sequences.",
           "Goethals & Seidel (1967); Seberry & Yamada (1992); Đoković (2010)", 7),
    Method("gs_golay", "Goethals–Seidel (Golay-derived T-sequences)",
           "n = 4g or n = 4(2g + 1) from the T-sequences of a Golay pair of length g.",
           "Golay (1961); Turyn (1974)", 8),
    Method("gs_tseq", "Goethals–Seidel (stored T-sequences)",
           "n = 4t from T-sequences of length t stored directly (t = 73).",
           "Cooper & Wallis (1972); Egan, Kotsireas, London & Ó Catháin (2026)", 9),
    Method("legendre", "Legendre pair – two-circulant core",
           "n = 2ℓ + 2 from a Legendre pair of odd length ℓ (periodic autocorrelations summing to −2).",
           "Fletcher, Gysin & Seberry (2001); Kotsireas & Koutschan (2021)", 10),
    Method("dscore", "Cyclic Hadamard difference set – circulant core",
           "n = ℓ + 1 from a Singer (m-sequence), twin-prime-power or Hall sextic-residue difference set.",
           "Singer (1938); Stanton & Sprott (1958); Hall (1956); Whiteman (1962)", 11),
    Method("lquad", "Legendre quadruple – bordered Goethals–Seidel",
           "n = 4m + 4 from four sequences of even length m with periodic autocorrelations summing to −4.",
           "Egan, Kotsireas, London & Ó Catháin (2026); Alpöge, Voinov & Reynolds-Haertle (2026)", 12),
    Method("scarpis", "Scarpis – Đoković",
           "n = q(q + 1) from an Hadamard matrix of order q + 1, q ≡ 3 (mod 4) a prime power.",
           "Scarpis (1898); Đoković (2017)", 13),
    Method("paley2_general", "Paley II – generalised (Williamson)",
           "n = h(q + 1) from an Hadamard matrix of order h > 2 and q ≡ 1 (mod 4).",
           "Williamson (1944)", 14),
    Method("baumert_hall", "Baumert–Hall",
           "n = 12w from the order-12 Baumert–Hall array filled with Williamson matrices.",
           "Baumert & Hall (1965)", 15),
    Method("ehlich", "Ehlich",
           "n = (k − 1)² from a skew Paley matrix of order k + 1 and a Paley II core of order k − 2.",
           "Ehlich (1965)", 16),
    Method("miyamoto", "Miyamoto",
           "n = 4q for a prime power q ≡ 1 (mod 4) when an Hadamard matrix of order q − 1 exists.",
           "Miyamoto (1991)", 17),
    Method("cooper_wallis", "Cooper–Wallis",
           "n = 4wL from Williamson matrices of order w and T-sequences of length L.",
           "Cooper & Wallis (1972)", 18),
    Method("sds", "Supplementary difference sets",
           "n = 4N from four supplementary difference sets in Z_N in the Goethals–Seidel array.",
           "Wallis & Whiteman (1972); Đoković (1992–2008); Alpöge, Voinov & Reynolds-Haertle (2026)", 19),
    Method("kronecker", "Kronecker product",
           "n = a·b from Hadamard matrices of orders a and b.",
           "Sylvester (1867)", 20),
]}

# methods new in version 2.0 (shown in the "Construction lab" tab)
NEW_IN_V2 = ("golay2c", "gs_golay", "gs_tseq", "legendre", "dscore", "lquad", "scarpis")


def is_hadamard_order(n: int) -> bool:
    return n in (1, 2) or (n > 0 and n % 4 == 0)


def candidate_params(method: str, n: int) -> List[dict]:
    """All parameter sets for which ``method`` *could* produce order ``n``.

    ``deps`` lists the orders of Hadamard matrices that must themselves be
    constructible; the planner resolves those recursively.
    """
    d = sequence_data()
    out: List[dict] = []
    if method == "trivial":
        if n in (1, 2):
            out.append({"deps": []})
    elif method == "sylvester":
        if n >= 2 and n & (n - 1) == 0:
            out.append({"k": n.bit_length() - 1, "deps": []})
    elif method == "paley1":
        q = n - 1
        if n >= 4 and q % 4 == 3 and is_prime_power(q):
            out.append({"q": q, "deps": []})
    elif method == "paley2":
        if n % 2 == 0:
            q = n // 2 - 1
            if q >= 5 and q % 4 == 1 and is_prime_power(q):
                out.append({"q": q, "deps": []})
    elif method == "paley2_general":
        for h in range(4, n, 4):
            if n % h == 0:
                q = n // h - 1
                if q >= 5 and q % 4 == 1 and is_prime_power(q):
                    out.append({"q": q, "h": h, "deps": [h]})
    elif method == "williamson":
        if n % 4 == 0 and n // 4 in d["williamson"] and n // 4 >= 3:
            out.append({"w": n // 4, "deps": []})
    elif method == "gs_turyn":
        if n % 4 == 0:
            for t in d["turyn"]:
                if n == 4 * (3 * t - 1):
                    out.append({"t": t, "deps": []})
    elif method == "gs_base":
        if n % 4 == 0:
            for m, seqs in d["base"].items():
                p = len(seqs[0]) - len(seqs[2])
                if n == 4 * (2 * m + p):
                    out.append({"m": m, "p": p, "deps": []})
    elif method == "baumert_hall":
        if n % 12 == 0:
            w = n // 12
            if w == 1 or w in d["williamson"]:
                out.append({"w": w, "deps": []})
    elif method == "ehlich":
        for k in range(7, 40):
            if (k - 1) ** 2 == n and k % 4 == 3 and is_prime_power(k) \
                    and (k - 2) % 4 == 1 and is_prime_power(k - 2):
                out.append({"k": k, "m": k - 2, "deps": []})
    elif method == "miyamoto":
        if n % 4 == 0:
            q = n // 4
            if q >= 5 and q % 4 == 1 and is_prime_power(q):
                out.append({"q": q, "deps": [q - 1]})
    elif method == "cooper_wallis":
        if n % 4 == 0:
            quarter = n // 4
            for w in d["williamson"]:
                if quarter % w == 0 and w >= 3:
                    L = quarter // w
                    if L in t_sequence_lengths():
                        src, key = _t_source(L)
                        out.append({"w": w, "L": L, "source": src, "key": key, "deps": []})
    elif method == "sds":
        if n % 4 == 0 and n // 4 in d["sds"]:
            out.append({"N": n // 4, "deps": []})
    elif method == "golay2c":
        if n % 2 == 0 and n // 2 in SQ.golay_lengths(n // 2):
            out.append({"g": n // 2, "deps": []})
    elif method == "gs_golay":
        if n % 4 == 0:
            for src, key in t_sequence_sources().get(n // 4, []):
                if src in ("golay", "golay2"):
                    out.append({"t": n // 4, "g": key, "source": src, "key": key, "deps": []})
    elif method == "gs_tseq":
        if n % 4 == 0 and n // 4 in d["tseq"]:
            out.append({"t": n // 4, "source": "tseq", "key": n // 4, "deps": []})
    elif method == "legendre":
        if n % 4 == 0:
            L = (n - 2) // 2
            for spec in legendre_pair_sources(n).get(L, []):
                out.append({"l": L, "spec": _freeze_spec(spec), "deps": []})
    elif method == "dscore":
        L = n - 1
        for spec in SQ.two_level_sources(n).get(L, []):
            out.append({"l": L, "spec": _freeze_spec(spec), "deps": []})
    elif method == "lquad":
        if n % 4 == 0 and (n - 4) // 4 in d["lquads"]:
            out.append({"m": (n - 4) // 4, "deps": []})
    elif method == "scarpis":
        for q in range(3, 40, 4):
            if q * (q + 1) == n and is_prime_power(q):
                out.append({"q": q, "deps": [q + 1]})
    elif method == "kronecker":
        for a in range(2, int(n ** 0.5) + 1):
            if n % a == 0:
                b = n // a
                if is_hadamard_order(a) and is_hadamard_order(b):
                    out.append({"a": a, "b": b, "deps": [a, b]})
    else:
        raise KeyError(method)
    return out


def _t_source(L: int) -> Tuple[str, int]:
    return t_sequence_sources()[L][0]


def _freeze_spec(spec: dict) -> str:
    """Specs are stored in the (hashable) planner as a compact string such as ``twin:q=11``."""
    return spec["kind"] + ":" + ",".join(f"{k}={v}" for k, v in sorted(spec.items()) if k != "kind")


def thaw_spec(s: str) -> dict:
    kind, _, rest = s.partition(":")
    out: dict = {"kind": kind}
    for item in rest.split(","):
        if not item:
            continue
        k, v = item.split("=")
        out[k] = (v == "True") if v in ("True", "False") else int(v)
    return out


# ==========================================================================
# building (given already-built dependencies)
# ==========================================================================
def build_direct(method: str, n: int, params: dict, deps: Dict[int, np.ndarray]) -> np.ndarray:
    if method == "trivial":
        return _i8(np.ones((1, 1))) if n == 1 else _i8([[1, 1], [1, -1]])
    if method == "sylvester":
        return sylvester(n)
    if method == "paley1":
        return paley1(params["q"])
    if method == "paley2":
        return paley2(params["q"])
    if method == "paley2_general":
        return paley2_general(params["q"], deps[params["h"]])
    if method == "williamson":
        return williamson(params["w"])
    if method == "gs_turyn":
        return goethals_seidel_turyn(params["t"])
    if method == "gs_base":
        return goethals_seidel_base(params["m"])
    if method == "baumert_hall":
        return baumert_hall(params["w"])
    if method == "ehlich":
        return ehlich(params["k"])
    if method == "miyamoto":
        return miyamoto(params["q"], deps[params["q"] - 1])
    if method == "cooper_wallis":
        return cooper_wallis(params["w"], params["source"], params["key"])
    if method == "sds":
        return supplementary_difference_sets(params["N"])
    if method == "golay2c":
        return golay_two_circulant(params["g"])
    if method in ("gs_golay", "gs_tseq"):
        return goethals_seidel_tseq(params["source"], params["key"])
    if method == "legendre":
        return legendre_two_circulant(thaw_spec(params["spec"]))
    if method == "dscore":
        return difference_set_core(thaw_spec(params["spec"]))
    if method == "lquad":
        return legendre_quadruple(params["m"])
    if method == "scarpis":
        return scarpis(params["q"], deps[params["q"] + 1])
    if method == "kronecker":
        return kronecker(deps[params["a"]], deps[params["b"]])
    raise KeyError(method)
