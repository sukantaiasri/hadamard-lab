"""Complementary sequences, difference-set sequences and the arrays that turn them
into Hadamard matrices (added in version 2.0).

Four families of ingredients are handled here:

* **Golay complementary pairs** of length ``g = 2^a 10^b 26^c`` (Golay 1961; Turyn 1974),
  generated from the seeds of lengths 2, 10 and 26 by Turyn's product.
* **Sequences with two-level periodic autocorrelation** (cyclic Hadamard difference sets):
  m-sequences / Singer sets (length ``2^k - 1``), twin-prime-power sets (Stanton & Sprott 1958;
  Whiteman 1962) and Hall's sextic-residue sets (Hall 1956).  Over a non-cyclic group the
  matrix is *group-developed* rather than circulant.
* **Legendre pairs** (Fletcher, Gysin & Seberry 2001; Kotsireas & Koutschan 2021; Egan,
  Kotsireas, London & Ó Catháin 2026): two ``+-1`` sequences of odd length ``l`` whose periodic
  autocorrelations sum to ``-2`` at every non-zero shift; they give order ``2l + 2`` through the
  two-circulant core.
* **Legendre quadruples** (Egan et al. 2026, Theorem 5): four sequences of even length ``m``
  whose periodic autocorrelations sum to ``-4``; a Goethals–Seidel array with a border of
  width four gives order ``4m + 4``.  The matrices of orders 668 and 716 announced by
  Alpöge, Voinov and Reynolds-Haertle (2026) are of this type.

Also here: the Scarpis construction in Đoković's generalisation (2017) and the recovery of
sequences stored in ``data/sequences.json`` (searched quadruples / pairs, T-sequences).

All functions return NumPy integer arrays; nothing here depends on Shiny.
"""
from __future__ import annotations

from functools import lru_cache
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .gf import field as gf_field, is_prime, is_prime_power, prime_power

# ==========================================================================
# autocorrelations and small helpers
# ==========================================================================
def paf(a: Sequence[int]) -> np.ndarray:
    """Periodic autocorrelation ``P_a(s) = sum_i a_i a_{i+s}`` (indices mod length)."""
    a = np.asarray(a, dtype=np.int64)
    return np.array([int(np.dot(a, np.roll(a, -s))) for s in range(len(a))])


def napf(a: Sequence[int]) -> np.ndarray:
    """Non-periodic (aperiodic) autocorrelation ``N_a(s) = sum_{i<n-s} a_i a_{i+s}``."""
    a = np.asarray(a, dtype=np.int64)
    n = len(a)
    return np.array([int(np.dot(a[: n - s], a[s:])) for s in range(n)])


def circulant(first_row: Sequence[int]) -> np.ndarray:
    r = np.asarray(first_row, dtype=np.int32)
    n = len(r)
    idx = (np.arange(n)[None, :] - np.arange(n)[:, None]) % n
    return r[idx]


def back_identity(n: int) -> np.ndarray:
    return np.eye(n, dtype=np.int32)[::-1]


def _i8(m) -> np.ndarray:
    return np.asarray(m, dtype=np.int8)


def is_hadamard(h: np.ndarray) -> bool:
    n = h.shape[0]
    if h.shape != (n, n) or not np.isin(h, (-1, 1)).all():
        return False
    g = h.astype(np.int64) @ h.T.astype(np.int64)
    return bool(np.array_equal(g, n * np.eye(n, dtype=np.int64)))


def is_golay_pair(a, b) -> bool:
    a, b = np.asarray(a), np.asarray(b)
    if len(a) != len(b) or not (np.isin(a, (-1, 1)).all() and np.isin(b, (-1, 1)).all()):
        return False
    t = napf(a) + napf(b)
    return bool((t[1:] == 0).all())


def is_legendre_pair(a, b) -> bool:
    a, b = np.asarray(a), np.asarray(b)
    if len(a) != len(b) or len(a) % 2 == 0:
        return False
    t = paf(a) + paf(b)
    return bool((t[1:] == -2).all())


def is_legendre_quadruple(seqs) -> bool:
    seqs = [np.asarray(s) for s in seqs]
    m = len(seqs[0])
    if len(seqs) != 4 or any(len(s) != m for s in seqs) or m % 2:
        return False
    t = sum(paf(s) for s in seqs)
    return bool((t[1:] == -4).all())


def has_two_level_paf(a) -> bool:
    """Periodic autocorrelation equal to ``-1`` at every non-zero shift (Hadamard difference set)."""
    p = paf(a)
    return bool((p[1:] == -1).all())


# ==========================================================================
# Golay complementary pairs
# ==========================================================================
# Seeds (verified by exhaustive search; lengths 2, 10 and 26 are Golay's original lengths).
GOLAY_SEEDS: Dict[int, Tuple[Tuple[int, ...], Tuple[int, ...]]] = {
    1: ((1,), (1,)),
    2: ((1, 1), (1, -1)),
    10: ((1, 1, 1, 1, 1, -1, 1, -1, -1, 1), (1, -1, 1, -1, -1, -1, 1, 1, -1, -1)),
    26: ((1, 1, 1, 1, -1, 1, 1, -1, -1, 1, -1, 1, -1, 1, -1, -1, 1, -1, 1, 1, 1, -1, -1, 1, 1, 1),
         (1, 1, 1, -1, -1, 1, 1, 1, -1, 1, -1, -1, -1, -1, -1, 1, -1, 1, 1, -1, -1, 1, -1, -1, -1, -1)),
}


def golay_product(a, b, c, d) -> Tuple[np.ndarray, np.ndarray]:
    """Turyn's product: Golay pairs of lengths ``m`` and ``n`` give one of length ``mn``.

    ``x = a (x) (c+d)/2 + b~ (x) (c-d)/2``,  ``y = b (x) (c+d)/2 - a~ (x) (c-d)/2``
    where ``~`` reverses a sequence (checked for every combination of the seeds in the tests).
    """
    a, b, c, d = (np.asarray(s, dtype=np.int64) for s in (a, b, c, d))
    cp, cm = (c + d) // 2, (c - d) // 2
    x = np.kron(a, cp) + np.kron(b[::-1], cm)
    y = np.kron(b, cp) - np.kron(a[::-1], cm)
    return x, y


def golay_lengths(limit: int) -> List[int]:
    """All ``g = 2^a 10^b 26^c`` with ``2 <= g <= limit``."""
    out = set()
    g26 = 1
    while g26 <= limit:
        g10 = g26
        while g10 <= limit:
            g2 = g10
            while g2 <= limit:
                if g2 >= 2:
                    out.add(g2)
                g2 *= 2
            g10 *= 10
        g26 *= 26
    return sorted(out)


@lru_cache(maxsize=None)
def golay_pair(g: int) -> Tuple[np.ndarray, np.ndarray]:
    """A Golay pair of length ``g`` (``g = 2^a 10^b 26^c``), built from the seeds by products."""
    if g in GOLAY_SEEDS:
        a, b = GOLAY_SEEDS[g]
        return np.array(a, dtype=np.int64), np.array(b, dtype=np.int64)
    for f in (26, 10, 2):
        if g % f == 0 and g // f >= 1:
            a, b = golay_pair(f)
            c, d = golay_pair(g // f)
            return golay_product(a, b, c, d)
    raise ValueError(f"{g} is not of the form 2^a 10^b 26^c")


def golay_two_circulant(a, b) -> np.ndarray:
    """``H = [[A, B], [-B^T, A^T]]`` with circulants ``A, B`` of a Golay pair: order ``2g``."""
    A, B = circulant(a), circulant(b)
    return _i8(np.block([[A, B], [-B.T, A.T]]))


def t_sequences_from_golay(a, b) -> List[np.ndarray]:
    """T-sequences of length ``g``: ``((a+b)/2, (a-b)/2, 0, 0)``."""
    a, b = np.asarray(a, dtype=np.int64), np.asarray(b, dtype=np.int64)
    z = np.zeros(len(a), dtype=np.int64)
    return [(a + b) // 2, (a - b) // 2, z, z]


def base_sequences_from_golay(a, b) -> List[np.ndarray]:
    """Base sequences ``BS(g+1, g) = ((a,1), (a,-1), b, b)`` from a Golay pair (Turyn 1974)."""
    a, b = np.asarray(a, dtype=np.int64), np.asarray(b, dtype=np.int64)
    return [np.append(a, 1), np.append(a, -1), b.copy(), b.copy()]


# ==========================================================================
# sequences with two-level periodic autocorrelation (cyclic Hadamard difference sets)
# ==========================================================================
def m_sequence(k: int) -> np.ndarray:
    """Maximal-length shift-register sequence of period ``2^k - 1`` (Singer difference set).

    The recurrence uses the primitive polynomial over GF(2) found by :mod:`gf`; the output is
    ``+-1`` with ``PAF(s) = -1`` for every ``s != 0``.
    """
    if k < 2:
        raise ValueError("k must be at least 2")
    f = gf_field(2 ** k)
    coeffs = f.modulus[:k]            # x^k = sum_j coeffs[j] x^j  (mod 2)
    n = 2 ** k - 1
    state = [1] + [0] * (k - 1)
    bits = []
    for _ in range(n):
        bits.append(state[0])
        new = sum(c * s for c, s in zip(coeffs, state)) % 2
        state = state[1:] + [new]
    seq = np.where(np.array(bits) == 1, -1, 1).astype(np.int64)
    return seq


def legendre_symbol_sequence(p: int, a0: int = 1) -> np.ndarray:
    """``(a0, chi(1), ..., chi(p-1))`` for a prime ``p`` (the Legendre sequence)."""
    if not is_prime(p) or p == 2:
        raise ValueError("p must be an odd prime")
    chi = gf_field(p).chi().astype(np.int64)
    chi[0] = a0
    return chi


def hall_sextic_sequence(p: int) -> np.ndarray:
    """Hall's sextic-residue difference set for a prime ``p = 4x^2 + 27`` (Hall 1956).

    The set is the union of three of the six cyclotomic classes of index 6; the right union
    (``C_0 + C_1 + C_3`` for a suitable primitive root) is found by testing the twenty candidates.
    """
    if not is_prime(p) or p % 6 != 1 or not _is_hall_prime(p):
        raise ValueError("p must be a prime of the form 4x^2 + 27")
    f = gf_field(p)
    log = f.log
    cls = np.full(p, -1, dtype=np.int64)
    cls[1:] = log[1:] % 6
    from itertools import combinations
    for trip in combinations(range(6), 3):
        seq = np.ones(p, dtype=np.int64)
        seq[np.isin(cls, trip)] = -1
        seq[0] = 1
        for s in (seq, -seq):
            if has_two_level_paf(s) and s.sum() == -1:
                return s
    raise RuntimeError("no sextic residue difference set found")  # pragma: no cover


def _is_hall_prime(p: int) -> bool:
    x = 0
    while 4 * x * x + 27 <= p:
        if 4 * x * x + 27 == p:
            return True
        x += 1
    return False


def twin_prime_power_function(q: int):
    """Stanton–Sprott set in ``GF(q)^+ x GF(q+2)^+``; returns ``(f, add_table)`` of the group.

    ``f`` is a ``+-1`` array over the group ``G`` (elements ``0 .. q(q+2)-1``, element ``x*(q+2)+y``)
    with ``f(g) = -1`` exactly on the complement of the difference set ``D``, and
    ``add_table[g, h] = g + h``.  When ``q`` and ``q+2`` are both primes ``G`` is cyclic and the
    function is also returned in the cyclic (CRT) labelling by :func:`twin_prime_sequence`.
    """
    if not (is_prime_power(q) and is_prime_power(q + 2)) or q % 2 == 0:
        raise ValueError("q and q+2 must be odd prime powers")
    f1, f2 = gf_field(q), gf_field(q + 2)
    c1, c2 = f1.chi().astype(np.int64), f2.chi().astype(np.int64)
    n2 = q + 2
    vals = np.ones(q * n2, dtype=np.int64)
    for x in range(q):
        for y in range(n2):
            g = x * n2 + y
            if y == 0:
                vals[g] = 1                     # (x, 0) in D
            elif x == 0:
                vals[g] = -1                    # (0, y), y != 0: not in D
            else:
                vals[g] = 1 if c1[x] * c2[y] == 1 else -1
    # normalise so that the row sum of the core is -1
    if vals.sum() != -1:
        vals = -vals
    add = (f1.add_table[:, None, :, None] * n2 + f2.add_table[None, :, None, :]).reshape(q * n2, q * n2)
    return vals, add.astype(np.int32)


def twin_prime_sequence(q: int) -> np.ndarray:
    """The twin-prime sequence of length ``q(q+2)`` in cyclic (CRT) labelling; needs both primes."""
    if not (is_prime(q) and is_prime(q + 2)):
        raise ValueError("cyclic labelling needs q and q+2 prime")
    vals, _ = twin_prime_power_function(q)
    n = q * (q + 2)
    seq = np.zeros(n, dtype=np.int64)
    for z in range(n):
        seq[z] = vals[(z % q) * (q + 2) + (z % (q + 2))]
    return seq


def group_developed(vals: np.ndarray, add: np.ndarray) -> np.ndarray:
    """``M[g, h] = f(h - g)`` over an abelian group given by its addition table."""
    n = len(vals)
    neg = np.zeros(n, dtype=np.int64)
    for g in range(n):
        neg[g] = int(np.where(add[g] == 0)[0][0])
    diff = add[neg[:, None], np.arange(n)[None, :]]   # (-g) + h
    return vals[diff]


def bordered_core(M: np.ndarray) -> np.ndarray:
    """``H = [[1, e'], [e, M]]`` for a core with ``M M' = (l+1) I - J`` and ``M J = -J``."""
    n = M.shape[0] + 1
    H = np.ones((n, n), dtype=np.int32)
    H[1:, 1:] = M
    return _i8(H)


# sources of two-level sequences / difference sets, keyed by length l (order l+1 for the core)
def two_level_sources(limit: int) -> Dict[int, List[dict]]:
    """``length -> [{"kind", ...}]`` for every l < limit with a known cyclic/abelian Hadamard
    difference set other than Paley's quadratic residues (those are the Paley I construction)."""
    out: Dict[int, List[dict]] = {}
    k = 2
    while 2 ** k - 1 < limit:
        out.setdefault(2 ** k - 1, []).append({"kind": "singer", "k": k, "cyclic": True})
        k += 1
    q = 3
    while q * (q + 2) < limit:
        if q % 2 and is_prime_power(q) and is_prime_power(q + 2):
            out.setdefault(q * (q + 2), []).append({"kind": "twin", "q": q,
                                                    "cyclic": is_prime(q) and is_prime(q + 2)})
        q += 2
    x = 1
    while 4 * x * x + 27 < limit:
        p = 4 * x * x + 27
        if is_prime(p):
            out.setdefault(p, []).append({"kind": "hall", "p": p, "cyclic": True})
        x += 1
    return dict(sorted(out.items()))


def two_level_sequence(spec: dict) -> np.ndarray:
    """The cyclic sequence of a source (``cyclic`` must be true)."""
    if spec["kind"] == "singer":
        s = m_sequence(spec["k"])
    elif spec["kind"] == "twin":
        s = twin_prime_sequence(spec["q"])
    elif spec["kind"] == "hall":
        s = hall_sextic_sequence(spec["p"])
    elif spec["kind"] == "paley":
        s = legendre_symbol_sequence(spec["p"], 1)
    else:
        raise KeyError(spec["kind"])
    if s.sum() != -1:
        s = -s
    return s


def difference_set_core_matrix(spec: dict) -> np.ndarray:
    """Hadamard matrix of order ``l + 1`` from a two-level source (circulant or group-developed core)."""
    if spec["kind"] == "twin" and not spec.get("cyclic", True):
        vals, add = twin_prime_power_function(spec["q"])
        return bordered_core(group_developed(vals, add))
    return bordered_core(circulant(two_level_sequence(spec)))


def describe_two_level(spec: dict) -> str:
    if spec["kind"] == "singer":
        return f"m-sequence (Singer set) over GF(2^{spec['k']}), length {2 ** spec['k'] - 1}"
    if spec["kind"] == "twin":
        q = spec["q"]
        grp = f"Z_{q * (q + 2)}" if spec["cyclic"] else f"GF({q})⁺ × GF({q + 2})⁺ (non-cyclic)"
        return f"twin prime powers {q}, {q + 2}; group {grp}"
    if spec["kind"] == "hall":
        return f"Hall sextic residues, p = {spec['p']} = 4x² + 27"
    if spec["kind"] == "paley":
        return f"quadratic residues modulo {spec['p']}"
    return spec["kind"]


# ==========================================================================
# Legendre pairs and the two-circulant core
# ==========================================================================
def two_circulant_core(a, b) -> np.ndarray:
    """Fletcher–Gysin–Seberry array of order ``2l + 2`` from a Legendre pair ``(a, b)``.

    ``A`` and ``B`` are the circulants, normalised so that every row sums to ``+1``::

        [[-, -,  j,   j ],
         [-, +,  j,  -j ],
         [j', j', A,   B ],
         [j',-j', B', -A']]
    """
    a, b = np.asarray(a, dtype=np.int64), np.asarray(b, dtype=np.int64)
    if a.sum() < 0:
        a = -a
    if b.sum() < 0:
        b = -b
    L = len(a)
    A, B = circulant(a), circulant(b)
    j = np.ones((1, L), dtype=np.int32)
    top = np.block([[np.array([[-1, -1]], dtype=np.int32), j, j],
                    [np.array([[-1, 1]], dtype=np.int32), j, -j]])
    bot = np.block([[j.T, j.T, A, B],
                    [j.T, -j.T, B.T, -A.T]])
    return _i8(np.vstack([top, bot]))


def legendre_pair_from_two_level(seq) -> Tuple[np.ndarray, np.ndarray]:
    """A sequence with ``PAF = -1`` off-peak gives the Legendre pair ``(s, s)``."""
    s = np.asarray(seq, dtype=np.int64)
    return s.copy(), s.copy()


def legendre_pair_paley(p: int) -> Tuple[np.ndarray, np.ndarray]:
    """Legendre pair of prime length ``p``: for ``p = 3 (mod 4)`` two copies of the Legendre
    sequence; for ``p = 1 (mod 4)`` the two choices ``a_0 = +1`` and ``a_0 = -1`` (their
    autocorrelations are ``-1 + 2chi(s)`` and ``-1 - 2chi(s)``)."""
    if p % 4 == 3:
        s = legendre_symbol_sequence(p, 1)
        return s, s.copy()
    return legendre_symbol_sequence(p, 1), legendre_symbol_sequence(p, -1)


def legendre_pair_sources(limit: int, stored: Dict[int, list]) -> Dict[int, List[dict]]:
    """``length -> [{"kind", ...}]`` for odd lengths ``l`` with ``2l + 2 <= limit``."""
    out: Dict[int, List[dict]] = {}
    for L in range(3, (limit - 2) // 2 + 1, 2):
        srcs = []
        if is_prime(L):
            srcs.append({"kind": "paley", "p": L})
        if L in stored:
            srcs.append({"kind": "stored", "l": L})
        for spec in two_level_sources(L + 1).get(L, []):
            if spec["cyclic"]:
                srcs.append(dict(spec))
        if srcs:
            out[L] = srcs
    return out


def legendre_pair(spec: dict, stored: Dict[int, list]) -> Tuple[np.ndarray, np.ndarray]:
    if spec["kind"] == "paley":
        return legendre_pair_paley(spec["p"])
    if spec["kind"] == "stored":
        a, b = stored[spec["l"]]
        return np.array(a, dtype=np.int64), np.array(b, dtype=np.int64)
    return legendre_pair_from_two_level(two_level_sequence(spec))


def describe_legendre_source(spec: dict) -> str:
    if spec["kind"] == "paley":
        return f"Legendre symbol sequences modulo {spec['p']}"
    if spec["kind"] == "stored":
        return f"stored pair of length {spec['l']} (computer search)"
    return "two copies of a two-level sequence: " + describe_two_level(spec)


# ==========================================================================
# Legendre quadruples and the bordered Goethals–Seidel array (Egan et al. 2026, Thm 5)
# ==========================================================================
def goethals_seidel_array(A, B, C, D) -> np.ndarray:
    R = back_identity(A.shape[0])
    return np.block([[A, B @ R, C @ R, D @ R],
                     [-B @ R, A, D.T @ R, -C.T @ R],
                     [-C @ R, -D.T @ R, A, B.T @ R],
                     [-D @ R, C.T @ R, -B.T @ R, A]])


def normalise_quadruple(seqs) -> List[np.ndarray]:
    """Order and sign the four sequences so that the row sums are ``(2, 0, 0, 0)``."""
    seqs = [np.asarray(s, dtype=np.int64) for s in seqs]
    sums = [int(s.sum()) for s in seqs]
    if sorted(abs(x) for x in sums) != [0, 0, 0, 2]:
        raise ValueError("a Legendre quadruple has row sums (±2, 0, 0, 0)")
    k = [abs(x) for x in sums].index(2)
    first = seqs[k] if sums[k] > 0 else -seqs[k]
    return [first] + [s for i, s in enumerate(seqs) if i != k]


def bordered_goethals_seidel(A, B, C, D) -> np.ndarray:
    """Theorem 5 of Egan, Kotsireas, London and Ó Catháin (2026): order ``4m + 4`` from circulants
    of even order ``m`` with ``AA' + BB' + CC' + DD' = (4m+4) I - 4J`` and row sums ``2, 0, 0, 0``."""
    m = A.shape[0]
    j = np.ones((1, m), dtype=np.int32)
    jt = j.T
    K = np.array([[1, -1, -1, -1], [-1, 1, -1, -1], [-1, -1, 1, -1], [-1, -1, -1, 1]], dtype=np.int32)
    E = np.block([[j, j, j, j], [j, -j, j, -j], [j, -j, -j, j], [j, j, -j, -j]])
    F = np.block([[jt, jt, jt, jt], [-jt, jt, jt, -jt], [-jt, -jt, jt, jt], [-jt, jt, -jt, jt]])
    M = goethals_seidel_array(A, B, C, D)
    return _i8(np.block([[K, E], [F, M]]))


def legendre_quadruple_matrix(seqs) -> np.ndarray:
    a, b, c, d = normalise_quadruple(seqs)
    return bordered_goethals_seidel(circulant(a), circulant(b), circulant(c), circulant(d))


# ==========================================================================
# Scarpis (1898) in the generalisation of Đoković (2017)
# ==========================================================================
def normalise_hadamard(A: np.ndarray) -> np.ndarray:
    A = A.astype(np.int32).copy()
    if A[0, 0] == -1:
        A = -A
    A = A * A[:, :1]          # rows: make first column +1
    A = A * A[:1, :]          # columns: make first row +1
    return A


def scarpis(q: int, A: np.ndarray) -> np.ndarray:
    """Order ``q(q+1)`` from a Hadamard matrix ``A`` of order ``q + 1``, ``q = 3 (mod 4)`` a prime power.

    Following Đoković (2017): normalise ``A``; let ``C`` be its core with rows ``c(x)``, ``x`` in
    ``GF(q)``.  ``B_0 = A' (x) j`` (``A'`` = ``A`` without its first row) and for ``r`` in
    ``GF(q)`` the block row ``B_r = [j' (x) c(r) | B_{r,x}, x in GF(q)]`` where row ``k`` of
    ``B_{r,x}`` is ``c(x r + k)``.  Only the bijection ``{1..q} -> GF(q)`` is a choice.
    """
    if q % 4 != 3 or not is_prime_power(q):
        raise ValueError("Scarpis needs a prime power q = 3 (mod 4)")
    n = q + 1
    if A.shape != (n, n):
        raise ValueError("A must have order q + 1")
    A = normalise_hadamard(A)
    Cc = A[1:, 1:]                            # core, rows indexed by GF(q) elements 0..q-1
    f = gf_field(q)
    mul, add = f.mul_table, f.add_table
    j = np.ones((1, q), dtype=np.int32)
    B0 = np.kron(A[1:, :], j)                 # q x nq
    rows = [B0]
    for r in range(q):
        blocks = [np.kron(j.T, Cc[r:r + 1, :])]         # q x q, every row c(r)
        for x in range(q):
            idx = add[mul[x, r], np.arange(q)]           # x r + k for k = 0..q-1
            blocks.append(Cc[idx, :])
        rows.append(np.hstack(blocks))
    return _i8(np.vstack(rows))


# ==========================================================================
# T-sequences (shared helpers used by constructions.py)
# ==========================================================================
def is_t_sequence_set(ts) -> bool:
    ts = [np.asarray(t) for t in ts]
    n = len(ts[0])
    if any(len(t) != n for t in ts):
        return False
    support = sum(np.abs(t) for t in ts)
    if not (support == 1).all():
        return False
    tot = sum(napf(t) for t in ts)
    return bool((tot[1:] == 0).all())


def t_sequences_from_base(A, B, C, D) -> List[np.ndarray]:
    A, B, C, D = (np.asarray(x, dtype=np.int64) for x in (A, B, C, D))
    n, p = len(C), len(A) - len(C)
    z_n = np.zeros(n, dtype=np.int64)
    z_np = np.zeros(n + p, dtype=np.int64)
    return [np.concatenate([(A + B) // 2, z_n]), np.concatenate([(A - B) // 2, z_n]),
            np.concatenate([z_np, (C + D) // 2]), np.concatenate([z_np, (C - D) // 2])]
