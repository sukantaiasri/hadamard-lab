"""Planning and executing constructions.

The planner answers three questions

* *which methods can produce order n?*      (:func:`feasible_methods`)
* *which order is constructible at all?*    (:func:`existence_table`)
* *how exactly is it built, step by step?*  (:func:`construct`, a recursion tree)

Feasibility is computed by dynamic programming over increasing orders, because
several methods (Kronecker, generalised Paley II, Miyamoto) need an Hadamard
matrix of a smaller order.
"""
from __future__ import annotations

import time
from dataclasses import dataclass, field
from functools import lru_cache
from typing import Dict, List, Optional, Tuple

import numpy as np

from . import constructions as C
from .gf import field as gf_field, is_prime_power, modulus_string, prime_power

ORDER_LIMIT = 1000
KNOWN_OPEN = (668, 716, 892)          # no Hadamard matrix known (n <= 1000)
KNOWN_BUT_NOT_IMPLEMENTED = (876,)  # exists in the literature (Djokovic); sequence data not available here


# ==========================================================================
# validity of a request
# ==========================================================================
def validate_order(n) -> Tuple[bool, str]:
    try:
        n_int = int(n)
    except (TypeError, ValueError):
        return False, "Enter an integer order."
    if n_int != n:
        return False, "The order must be an integer."
    if n_int < 1:
        return False, "The order must be at least 1."
    if n_int > ORDER_LIMIT:
        return False, f"Orders above {ORDER_LIMIT} are not supported."
    if n_int not in (1, 2) and n_int % 4 != 0:
        return False, (f"No Hadamard matrix of order {n_int} exists: the order must be 1, 2 "
                       "or a multiple of 4 (a necessary condition proved by Hadamard).")
    return True, ""


# ==========================================================================
# feasibility DP
# ==========================================================================
@dataclass(frozen=True)
class Option:
    method: str
    params: Tuple[Tuple[str, object], ...]     # hashable form of the parameter dict
    deps: Tuple[int, ...]
    depth: int                                   # recursion depth (0 = direct construction)

    @property
    def p(self) -> dict:
        return dict(self.params)


def _freeze(d: dict) -> Tuple[Tuple[str, object], ...]:
    return tuple(sorted((k, v) for k, v in d.items() if k != "deps"))


@lru_cache(maxsize=1)
def _dp(limit: int = ORDER_LIMIT) -> Dict[int, List[Option]]:
    """For every admissible order up to ``limit`` the feasible options, best first."""
    table: Dict[int, List[Option]] = {}
    ordered = sorted(C.METHODS.values(), key=lambda m: m.priority)
    for n in [1, 2] + list(range(4, limit + 1, 4)):
        opts: List[Option] = []
        for meth in ordered:
            cands = []
            for prm in C.candidate_params(meth.key, n):
                deps = tuple(prm["deps"])
                if all(d in table and table[d] for d in deps):
                    depth = 0 if not deps else 1 + max(min(o.depth for o in table[d]) for d in deps)
                    cands.append(Option(meth.key, _freeze(prm), deps, depth))
            cands.sort(key=lambda o: (o.depth, max(o.deps) if o.deps else 0))
            opts.extend(cands)
        table[n] = opts
    return table


def feasible_options(n: int) -> List[Option]:
    return list(_dp().get(n, []))


def feasible_methods(n: int) -> List[str]:
    seen, out = set(), []
    for o in feasible_options(n):
        if o.method not in seen:
            seen.add(o.method)
            out.append(o.method)
    return out


def is_constructible(n: int) -> bool:
    return bool(_dp().get(n))


def existence_table(limit: int = ORDER_LIMIT):
    """Rows (order, status, first_method, methods) for every multiple of 4 up to ``limit``."""
    rows = []
    for n in range(4, limit + 1, 4):
        ms = feasible_methods(n)
        if ms:
            status = "constructible"
        elif n in KNOWN_OPEN:
            status = "open (no matrix known)"
        elif n in KNOWN_BUT_NOT_IMPLEMENTED:
            status = "exists in literature, not implemented"
        else:  # pragma: no cover  (cannot happen with the shipped data)
            status = "not implemented"
        rows.append((n, status, ms[0] if ms else "", ms))
    return rows


# ==========================================================================
# human-readable reasons (why a method does / does not apply)
# ==========================================================================
def _pp_text(q: int) -> str:
    pk = prime_power(q)
    if pk is None:
        return f"{q} is not a prime power"
    p, k = pk
    return f"{q} = {p}" + (f"^{k}" if k > 1 else "") + " is a prime power"


def explain_method(method: str, n: int) -> Tuple[bool, str]:
    """Return (applicable, one-line explanation) for a method and order."""
    d = C.sequence_data()
    feas = {o.method for o in feasible_options(n)}
    ok = method in feas
    if method == "trivial":
        return (n in (1, 2)), "Only orders 1 and 2."
    if method == "sylvester":
        return ok, ("n = 2^%d" % (n.bit_length() - 1)) if ok else f"{n} is not a power of 2."
    if method == "paley1":
        q = n - 1
        if ok:
            return True, f"q = n − 1 = {q}: {_pp_text(q)} with q ≡ 3 (mod 4)."
        return False, (f"q = n − 1 = {q} is not a prime power ≡ 3 (mod 4)." if q > 0 else "n too small.")
    if method == "paley2":
        q = n // 2 - 1
        if ok:
            return True, f"q = n/2 − 1 = {q}: {_pp_text(q)} with q ≡ 1 (mod 4)."
        return False, f"q = n/2 − 1 = {q} is not a prime power ≡ 1 (mod 4)."
    if method == "paley2_general":
        opts = [o for o in feasible_options(n) if o.method == method]
        if opts:
            o = opts[0].p
            return True, (f"n = h(q+1) with h = {o['h']} and q = {o['q']} ≡ 1 (mod 4); "
                          f"an Hadamard matrix of order {o['h']} is available.")
        return False, "No factorisation n = h(q+1) with h ≥ 4 Hadamard and q ≡ 1 (mod 4) a prime power."
    if method == "williamson":
        w = n // 4
        if ok:
            return True, f"Williamson matrices of order {w} are available."
        known = sorted(d["williamson"])
        return False, f"Needs Williamson matrices of order n/4 = {w}; stored orders: {', '.join(map(str, known))}."
    if method == "gs_turyn":
        if ok:
            t = [o.p["t"] for o in feasible_options(n) if o.method == method][0]
            return True, f"Turyn-type sequences of length {t}: n = 4(3·{t} − 1)."
        return False, ("Needs n = 4(3t − 1) with Turyn-type sequences of length t ∈ "
                       + ", ".join(map(str, sorted(d["turyn"]))) + ".")
    if method == "gs_base":
        if ok:
            o = [o.p for o in feasible_options(n) if o.method == method][0]
            return True, f"Base sequences with m = {o['m']}, p = {o['p']}: n = 4(2m + p)."
        return False, "No stored base sequences give this order."
    if method == "baumert_hall":
        w = n // 12 if n % 12 == 0 else None
        if ok:
            return True, f"n = 12·{w} with Williamson matrices of order {w}."
        return False, (f"n is not 12w for a stored Williamson order w." if w is None
                       else f"No Williamson matrices of order {w} stored.")
    if method == "ehlich":
        if ok:
            o = [o.p for o in feasible_options(n) if o.method == method][0]
            return True, f"k = {o['k']} (≡ 3 mod 4), k − 2 = {o['m']} (≡ 1 mod 4); n = (k−1)²."
        return False, "Needs n = (k−1)² with k, k−2 prime powers ≡ 3, 1 (mod 4) respectively."
    if method == "miyamoto":
        q = n // 4
        if ok:
            return True, f"q = n/4 = {q}: {_pp_text(q)}, q ≡ 1 (mod 4); an Hadamard matrix of order {q-1} is available."
        if n % 4 == 0 and q % 4 == 1 and is_prime_power(q):
            return False, f"q = {q} qualifies, but no Hadamard matrix of order {q-1} is available."
        return False, f"q = n/4 = {q} is not a prime power ≡ 1 (mod 4)."
    if method == "cooper_wallis":
        if ok:
            o = [o.p for o in feasible_options(n) if o.method == method][0]
            return True, f"n = 4·{o['w']}·{o['L']}: Williamson order {o['w']}, T-sequences of length {o['L']}."
        return False, "Needs n = 4wL with stored Williamson order w and T-sequence length L."
    if method == "sds":
        if ok:
            return True, f"4 supplementary difference sets in Z_{n // 4}."
        return False, "Needs stored supplementary difference sets with N = n/4 ∈ " + ", ".join(map(str, sorted(d["sds"]))) + "."
    if method == "kronecker":
        if ok:
            o = [o.p for o in feasible_options(n) if o.method == method][0]
            return True, f"n = {o['a']} × {o['b']} (both orders constructible)."
        return False, "No factorisation into two constructible Hadamard orders."
    return False, ""


def alternatives(n: int):
    """One row per method: (key, name, applicable, explanation)."""
    rows = []
    for m in sorted(C.METHODS.values(), key=lambda m: m.priority):
        if m.key == "trivial" and n not in (1, 2):
            continue
        ok, why = explain_method(m.key, n)
        rows.append((m.key, m.name, ok, why))
    return rows


# ==========================================================================
# construction (recursive, with trace)
# ==========================================================================
@dataclass
class Node:
    order: int
    method: str
    params: dict
    children: List["Node"] = field(default_factory=list)
    notes: List[str] = field(default_factory=list)
    seconds: float = 0.0

    @property
    def label(self) -> str:
        return C.METHODS[self.method].name

    def walk(self, depth=0):
        yield depth, self
        for c in self.children:
            yield from c.walk(depth + 1)


@dataclass
class Construction:
    matrix: np.ndarray
    tree: Node
    seconds: float

    @property
    def order(self) -> int:
        return self.matrix.shape[0]

    @property
    def method(self) -> str:
        return self.tree.method


def _notes(method: str, n: int, p: dict) -> List[str]:
    out: List[str] = []
    if method in ("paley1", "paley2", "paley2_general", "miyamoto", "ehlich"):
        qs = {"paley1": [p.get("q")], "paley2": [p.get("q")], "paley2_general": [p.get("q")],
              "miyamoto": [p.get("q")], "ehlich": [p.get("k"), p.get("m")]}[method]
        for q in qs:
            f = gf_field(q)
            out.append(f"GF({q}): {modulus_string(f)}; generator x (primitive).")
    if method == "williamson":
        out.append(f"Four symmetric circulant ±1 matrices of order {p['w']}; "
                   f"A²+B²+C²+D² = {4 * p['w']}·I (verified on import).")
    if method in ("gs_turyn", "gs_base", "cooper_wallis"):
        src = "Turyn-type" if (method == "gs_turyn" or p.get("source") == "turyn") else "base"
        out.append(f"T-sequences derived from {src} sequences; circulant blocks fed into the Goethals–Seidel array.")
    if method == "sds":
        sets = C.sequence_data()["sds"][p["N"]]
        out.append(f"Subset sizes {[len(s) for s in sets]}, λ = {sum(len(s) for s in sets) - p['N']} in Z_{p['N']}.")
    return out


def _build_node(opt: Option, n: int, cache: Dict[int, Tuple[np.ndarray, Node]]) -> Tuple[np.ndarray, Node]:
    t0 = time.perf_counter()
    p = opt.p
    deps: Dict[int, np.ndarray] = {}
    children: List[Node] = []
    for d in opt.deps:
        if d not in cache:
            cache[d] = _build_best(d, cache)
        deps[d] = cache[d][0]
        if cache[d][1] not in children:
            children.append(cache[d][1])
    mat = C.build_direct(opt.method, n, p, deps)
    node = Node(n, opt.method, dict(p), children, _notes(opt.method, n, p), time.perf_counter() - t0)
    return mat, node


def _build_best(n: int, cache) -> Tuple[np.ndarray, Node]:
    opts = feasible_options(n)
    if not opts:
        raise ValueError(f"No construction available for order {n}.")
    return _build_node(opts[0], n, cache)


def construct(n: int, method: Optional[str] = None, option_index: int = 0) -> Construction:
    """Build an Hadamard matrix of order ``n``.

    ``method=None`` picks the preferred applicable method; otherwise the named
    method is used (``option_index`` selects among several parameter sets).
    Raises ``ValueError`` with an explanatory message if impossible.
    """
    ok, msg = validate_order(n)
    if not ok:
        raise ValueError(msg)
    n = int(n)
    t0 = time.perf_counter()
    cache: Dict[int, Tuple[np.ndarray, Node]] = {}
    opts = feasible_options(n)
    if method:
        opts = [o for o in opts if o.method == method]
        if not opts:
            _, why = explain_method(method, n)
            raise ValueError(f"{C.METHODS[method].name} cannot produce order {n}. {why}")
    if not opts:
        if n in KNOWN_OPEN:
            raise ValueError(f"Order {n} is an open case of the Hadamard conjecture: "
                             "no Hadamard matrix of this order is currently known.")
        raise ValueError(f"No construction for order {n} is implemented in this version.")
    opt = opts[min(option_index, len(opts) - 1)]
    mat, node = _build_node(opt, n, cache)
    return Construction(mat, node, time.perf_counter() - t0)


def describe_node(node: Node) -> str:
    """Compact one-line description of a node, e.g. ``Paley I (q = 11)``."""
    p = node.params
    show = {k: v for k, v in p.items() if k not in ("source", "key")}
    inner = ", ".join(f"{k} = {v}" for k, v in show.items())
    return f"{node.label}" + (f" ({inner})" if inner else "")
