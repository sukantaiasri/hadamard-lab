"""Finite-field arithmetic GF(p^k) for the combinatorial constructions.

Elements are the integers ``0 .. q-1``.  For a prime field the integer *is* the
residue.  For an extension field the base-``p`` digits of the integer are the
coefficients (constant term first) of a polynomial over GF(p), reduced modulo a
primitive polynomial that is found by direct search (no lookup table needed).

Everything needed by the Paley-type constructions, MOLS and orthogonal arrays
is exposed as dense NumPy tables, which is fast for ``q <= 1000``.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from itertools import product
from typing import Optional, Tuple

import numpy as np


# --------------------------------------------------------------------------
# number-theory helpers
# --------------------------------------------------------------------------
def is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    i = 3
    while i * i <= n:
        if n % i == 0:
            return False
        i += 2
    return True


def prime_power(q: int) -> Optional[Tuple[int, int]]:
    """Return ``(p, k)`` with ``q == p**k`` for a prime ``p``, else ``None``."""
    if q < 2:
        return None
    for p in range(2, int(q ** 0.5) + 2):
        if q % p == 0:
            if not is_prime(p):
                return None
            k = 0
            m = q
            while m % p == 0:
                m //= p
                k += 1
            return (p, k) if m == 1 else None
    return (q, 1)  # q itself is prime


def is_prime_power(q: int) -> bool:
    return prime_power(q) is not None


def prime_powers_up_to(limit: int):
    return [q for q in range(2, limit + 1) if is_prime_power(q)]


# --------------------------------------------------------------------------
# the field
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class GF:
    """A finite field GF(q), q = p**k, with precomputed arithmetic tables."""

    q: int
    p: int
    k: int
    modulus: Tuple[int, ...]      # monic primitive polynomial, low -> high (length k+1)
    exp: np.ndarray               # exp[i] = g**i, length q-1 (g = x, a primitive element)
    log: np.ndarray               # log[a] = i with g**i = a  (log[0] = -1)
    digits: np.ndarray            # (q, k) base-p digits of each element
    _add: np.ndarray              # (q, q) addition table
    _neg: np.ndarray              # (q,)   additive inverses

    # ---- arithmetic -------------------------------------------------------
    def add(self, a, b):
        return self._add[a, b]

    def neg(self, a):
        return self._neg[a]

    def sub(self, a, b):
        return self._add[a, self._neg[b]]

    def mul(self, a, b):
        a = np.asarray(a)
        b = np.asarray(b)
        out = self.exp[(self.log[a] + self.log[b]) % (self.q - 1)]
        return np.where((a == 0) | (b == 0), 0, out)

    @property
    def add_table(self) -> np.ndarray:
        return self._add

    @property
    def sub_table(self) -> np.ndarray:
        """``sub_table[a, b] = a - b``."""
        return self._add[:, self._neg]

    @property
    def mul_table(self) -> np.ndarray:
        r = np.arange(self.q)
        return self.mul(r[:, None], r[None, :]).astype(np.int32)

    # ---- quadratic character & Jacobsthal matrix -------------------------
    def chi(self) -> np.ndarray:
        """Quadratic character as a length-q vector with values 0, +1, -1."""
        v = np.zeros(self.q, dtype=np.int8)
        if self.p == 2:
            raise ValueError("The quadratic character needs an odd characteristic.")
        nz = np.arange(1, self.q)
        v[nz] = np.where(self.log[nz] % 2 == 0, 1, -1)
        return v

    def jacobsthal(self) -> np.ndarray:
        """Jacobsthal matrix ``Q[x, y] = chi(x - y)`` (order q)."""
        chi = self.chi()
        if self.k == 1:  # circulant shortcut, identical result
            idx = (np.arange(self.q)[:, None] - np.arange(self.q)[None, :]) % self.q
            return chi[idx]
        return chi[self.sub_table]


def _primitive_modulus(p: int, k: int) -> Tuple[int, ...]:
    """Search for a monic primitive polynomial of degree ``k`` over GF(p)."""
    q = p ** k
    if k == 1:
        # x - g for a primitive root g of the prime field; represent modulus as (-g, 1)
        for g in range(1, p) if p > 2 else [1]:
            x, order = g, 1
            while x != 1:
                x = (x * g) % p
                order += 1
            if order == p - 1:
                return ((-g) % p, 1)
        return (1, 1)  # p == 2
    for tail in product(range(p), repeat=k):
        if tail[0] == 0:
            continue                      # constant term must be non-zero
        coeffs = list(tail) + [1]         # monic
        # multiply by x repeatedly, starting from 1
        cur = [1] + [0] * (k - 1)
        order = 0
        while True:
            top = cur[-1]
            cur = [0] + cur[:-1]
            if top:
                cur = [(c - top * m) % p for c, m in zip(cur, coeffs[:k])]
            order += 1
            if cur == [1] + [0] * (k - 1) or order > q:
                break
        if order == q - 1:
            return tuple(coeffs)
    raise RuntimeError(f"no primitive polynomial found for GF({p}^{k})")


@lru_cache(maxsize=64)
def field(q: int) -> GF:
    """Return (and cache) the field GF(q)."""
    pk = prime_power(q)
    if pk is None:
        raise ValueError(f"{q} is not a prime power")
    p, k = pk
    modulus = _primitive_modulus(p, k)

    digits = np.zeros((q, k), dtype=np.int64)
    for a in range(q):
        m = a
        for j in range(k):
            digits[a, j] = m % p
            m //= p
    weights = p ** np.arange(k)

    # exp / log tables through repeated multiplication by x
    exp = np.zeros(q - 1, dtype=np.int64)
    if k == 1:
        g = (-modulus[0]) % p if p > 2 else 1
        x = 1
        for i in range(q - 1):
            exp[i] = x
            x = (x * g) % p
    else:
        cur = [1] + [0] * (k - 1)
        for i in range(q - 1):
            exp[i] = int(np.dot(cur, weights))
            top = cur[-1]
            cur = [0] + cur[:-1]
            if top:
                cur = [(c - top * m) % p for c, m in zip(cur, modulus[:k])]
    log = np.full(q, -1, dtype=np.int64)
    log[exp] = np.arange(q - 1)
    if len(set(exp.tolist())) != q - 1:
        raise RuntimeError("generator does not have full order")  # defensive

    add = ((digits[:, None, :] + digits[None, :, :]) % p) @ weights
    neg = ((-digits) % p) @ weights
    return GF(q, p, k, modulus, exp, log, digits, add.astype(np.int32), neg.astype(np.int32))


def modulus_string(f: GF) -> str:
    """Human-readable defining polynomial, e.g. ``x^2 + x + 2 over GF(3)``."""
    if f.k == 1:
        return f"prime field Z/{f.p}"
    terms = []
    for power in range(f.k, -1, -1):
        c = f.modulus[power]
        if c == 0:
            continue
        mon = "" if power == 0 else ("x" if power == 1 else f"x^{power}")
        coef = "" if (c == 1 and power > 0) else str(c)
        terms.append(f"{coef}{mon}")
    return " + ".join(terms) + f" over GF({f.p})"
