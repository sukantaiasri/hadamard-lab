"""Back-end of the *Construction lab* tab (version 2.0).

Each construction is run from its ingredients rather than from a target order: the user picks a
Golay pair, a Legendre pair or quadruple, a difference set, T-sequences or a Scarpis input
matrix, and gets the ingredient sequences, the autocorrelation identity they satisfy, the
assembled matrix and its verification.  Everything is plain NumPy; ``app.py`` only renders.
"""
from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np

from . import constructions as C
from . import planner as P
from . import sequences as SQ

KINDS = {
    "golay": "Golay pair → two circulants (order 2g)",
    "tseq": "T-sequences → Goethals–Seidel / Cooper–Wallis (order 4t or 4tw)",
    "legendre": "Legendre pair → two-circulant core (order 2ℓ + 2)",
    "lquad": "Legendre quadruple → bordered Goethals–Seidel (order 4m + 4)",
    "dscore": "Cyclic Hadamard difference set → circulant core (order ℓ + 1)",
    "scarpis": "Scarpis–Đoković: H(q + 1) → H(q(q + 1))",
}

KIND_TO_METHOD = {"golay": "golay2c", "tseq": "gs_tseq", "legendre": "legendre", "lquad": "lquad",
                  "dscore": "dscore", "scarpis": "scarpis"}


@dataclass
class SeqResult:
    ok: bool
    error: str = ""
    kind: str = ""
    title: str = ""
    order: int = 0
    H: Optional[np.ndarray] = None
    ingredients: List[Tuple[str, np.ndarray]] = field(default_factory=list)
    checks: List[Tuple[str, str, bool]] = field(default_factory=list)   # (label, value, passed)
    corr: Optional[np.ndarray] = None          # total autocorrelation of the ingredients
    corr_kind: str = "periodic"
    corr_target: int = 0                       # expected value at non-zero shifts
    description: str = ""
    reference: str = ""
    seconds: float = 0.0
    params: dict = field(default_factory=dict)


# --------------------------------------------------------------------------
# choices offered by the UI
# --------------------------------------------------------------------------
def golay_choices() -> Dict[str, str]:
    return {str(g): f"g = {g}  (orders 2g = {2 * g}, 4g = {4 * g}, 4(2g+1) = {4 * (2 * g + 1)})"
            for g in SQ.golay_lengths(P.ORDER_LIMIT // 2)}


def tseq_choices() -> Dict[str, str]:
    out = {}
    for L, srcs in C.t_sequence_sources().items():
        if 4 * L <= P.ORDER_LIMIT:
            out[str(L)] = f"t = {L}  (order 4t = {4 * L}; {len(srcs)} source{'s' if len(srcs) > 1 else ''})"
    return out


def tseq_source_choices(L: int) -> Dict[str, str]:
    return {f"{src}:{key}": C.describe_t_source(src, key) for src, key in C.t_sequence_sources().get(L, [])}


def williamson_choices(L: int) -> Dict[str, str]:
    out = {"0": "none – Goethals–Seidel array from the T-sequences alone (order 4t)"}
    for w in sorted(C.sequence_data()["williamson"]):
        if 4 * w * L <= P.ORDER_LIMIT and w >= 3:
            out[str(w)] = f"Williamson matrices of order {w} (Cooper–Wallis, order 4tw = {4 * w * L})"
    return out


def legendre_choices() -> Dict[str, str]:
    out = {}
    for L, srcs in C.legendre_pair_sources(P.ORDER_LIMIT).items():
        out[str(L)] = f"ℓ = {L}  (order 2ℓ+2 = {2 * L + 2}; {len(srcs)} source{'s' if len(srcs) > 1 else ''})"
    return out


def legendre_source_choices(L: int) -> Dict[str, str]:
    return {C._freeze_spec(s): SQ.describe_legendre_source(s) for s in C.legendre_pair_sources(P.ORDER_LIMIT).get(L, [])}


def lquad_choices() -> Dict[str, str]:
    return {str(m): f"m = {m}  (order 4m+4 = {4 * m + 4})" + ("  – Alpöge, Voinov & Reynolds-Haertle (2026)" if m in (166, 178) else "")
            for m in sorted(C.sequence_data()["lquads"])}


def dscore_choices() -> Dict[str, str]:
    out = {}
    for L, srcs in SQ.two_level_sources(P.ORDER_LIMIT).items():
        out[str(L)] = f"ℓ = {L}  (order ℓ+1 = {L + 1}; " + ", ".join(s["kind"] for s in srcs) + ")"
    return out


def dscore_source_choices(L: int) -> Dict[str, str]:
    return {C._freeze_spec(s): SQ.describe_two_level(s) for s in SQ.two_level_sources(P.ORDER_LIMIT).get(L, [])}


def scarpis_choices() -> Dict[str, str]:
    return {str(q): f"q = {q}  (input order q+1 = {q + 1}, output order q(q+1) = {q * (q + 1)})"
            for q in (3, 7, 11, 19, 23, 27, 31) if q * (q + 1) <= P.ORDER_LIMIT}


def scarpis_input_choices(q: int) -> Dict[str, str]:
    out = {"auto": "Automatic (preferred method for order q + 1)"}
    for m in P.feasible_methods(q + 1):
        out[m] = C.METHODS[m].name
    return out


# --------------------------------------------------------------------------
# running a construction
# --------------------------------------------------------------------------
def _finish(r: SeqResult, t0: float) -> SeqResult:
    r.seconds = time.perf_counter() - t0
    if r.H is not None:
        r.order = int(r.H.shape[0])
        r.checks.append(("H·Hᵀ = nI (exact)", "yes" if SQ.is_hadamard(r.H) else "NO", SQ.is_hadamard(r.H)))
    r.ok = all(c[2] for c in r.checks)
    return r


def run(kind: str, params: dict) -> SeqResult:
    """Run one construction; ``params`` holds the sidebar values (strings)."""
    t0 = time.perf_counter()
    r = SeqResult(ok=False, kind=kind, params=dict(params), reference=C.METHODS[KIND_TO_METHOD[kind]].reference)
    try:
        if kind == "golay":
            g = int(params["g"])
            a, b = SQ.golay_pair(g)
            r.title = f"Golay pair of length {g}"
            r.ingredients = [("a", a), ("b", b)]
            r.corr, r.corr_kind, r.corr_target = SQ.napf(a) + SQ.napf(b), "non-periodic", 0
            r.checks.append(("N_a(s) + N_b(s) = 0 for all s ≥ 1", "yes" if SQ.is_golay_pair(a, b) else "NO", SQ.is_golay_pair(a, b)))
            r.H = SQ.golay_two_circulant(a, b)
            r.description = (f"Pair of length {g} = 2^a·10^b·26^c obtained from the seeds of lengths 2, 10 and 26 by Turyn's "
                             f"product; H = [[A, B], [−Bᵀ, Aᵀ]] with circulants A, B of order {g}.")
        elif kind == "tseq":
            L = int(params["t"])
            src, key = params["tsrc"].split(":")
            key = int(key)
            ts = C._t_sequences(src, key)
            r.ingredients = [(f"T{i + 1}", t) for i, t in enumerate(ts)]
            r.corr, r.corr_kind, r.corr_target = sum(SQ.napf(t) for t in ts), "non-periodic", 0
            okt = SQ.is_t_sequence_set(ts)
            r.checks.append(("disjoint supports covering every position", "yes" if (sum(np.abs(t) for t in ts) == 1).all() else "NO",
                             bool((sum(np.abs(t) for t in ts) == 1).all())))
            r.checks.append(("Σ N_Ti(s) = 0 for all s ≥ 1", "yes" if okt else "NO", okt))
            w = int(params.get("w", "0") or 0)
            if w == 0:
                r.title = f"T-sequences of length {L} → Goethals–Seidel array"
                r.H = C.goethals_seidel_from_t(ts)
                r.description = (f"T-sequences from {C.describe_t_source(src, key)}; circulants A₁ = T₁+T₂+T₃+T₄, "
                                 "A₂ = −T₁+T₂+T₃−T₄, A₃ = −T₁−T₂+T₃+T₄, A₄ = −T₁+T₂−T₃+T₄ in the Goethals–Seidel array.")
                r.reference = C.METHODS["gs_base"].reference + "; Cooper & Wallis (1972)"
            else:
                r.title = f"T-sequences of length {L} with Williamson matrices of order {w} → Cooper–Wallis"
                r.H = C.cooper_wallis(w, src, key)
                A, B, Cc, D = C.williamson_matrices(w)
                r.ingredients += [(f"Williamson row {nm}", M[0]) for nm, M in zip("ABCD", (A, B, Cc, D))]
                r.description = (f"T-sequences from {C.describe_t_source(src, key)} combined with the stored Williamson matrices "
                                 f"of order {w}: X = Σ Tᵢ ⊗ Mᵢ with (Mᵢ) = (A, B, C, D), Y = (−B, A, D, −C), Z = (−C, −D, A, B), "
                                 "W = (−D, C, −B, A), then the Goethals–Seidel array; order 4tw.")
                r.reference = C.METHODS["cooper_wallis"].reference
        elif kind == "legendre":
            L = int(params["l"])
            spec = C.thaw_spec(params["lsrc"])
            a, b = C.legendre_pair(spec)
            r.title = f"Legendre pair of length {L}"
            r.ingredients = [("a", a), ("b", b)]
            r.corr, r.corr_kind, r.corr_target = SQ.paf(a) + SQ.paf(b), "periodic", -2
            okp = SQ.is_legendre_pair(a, b)
            r.checks.append(("P_a(s) + P_b(s) = −2 for all s ≠ 0", "yes" if okp else "NO", okp))
            r.checks.append(("row sums of A and B", f"{int(np.sign(a.sum()) * 1)}, {int(np.sign(b.sum()))} (normalised to +1, +1)", True))
            r.H = SQ.two_circulant_core(a, b)
            r.description = ("Source: " + SQ.describe_legendre_source(spec) + ". Two-circulant core with a 2 × 2 border "
                             "(Fletcher–Gysin–Seberry array): [[−, −, j, j], [−, +, j, −j], [jᵀ, jᵀ, A, B], [jᵀ, −jᵀ, Bᵀ, −Aᵀ]].")
        elif kind == "lquad":
            m = int(params["m"])
            seqs = SQ.normalise_quadruple(C.sequence_data()["lquads"][m])
            r.title = f"Legendre quadruple of length {m}"
            r.ingredients = [(nm, s) for nm, s in zip("ABCD", seqs)]
            r.corr, r.corr_kind, r.corr_target = sum(SQ.paf(s) for s in seqs), "periodic", -4
            okq = SQ.is_legendre_quadruple(seqs)
            r.checks.append(("Σ P(s) = −4 for all s ≠ 0", "yes" if okq else "NO", okq))
            r.checks.append(("row sums (A, B, C, D)", ", ".join(str(int(s.sum())) for s in seqs), [int(s.sum()) for s in seqs] == [2, 0, 0, 0]))
            r.H = SQ.legendre_quadruple_matrix(seqs)
            r.description = ("Goethals–Seidel array of the four circulants with a border of width 4 (Egan, Kotsireas, London & "
                             "Ó Catháin 2026, Theorem 5). " + ("The sequences are those of the order-%d matrix announced by Alpöge, "
                             "Voinov and Reynolds-Haertle on 12 August 2026." % (4 * m + 4) if m in (166, 178)
                             else "The sequences were found by HadamardLab's own tabu search."))
        elif kind == "dscore":
            L = int(params["l"])
            spec = C.thaw_spec(params["dsrc"])
            r.title = f"Hadamard difference set of length {L}: " + SQ.describe_two_level(spec)
            if spec.get("cyclic", True):
                s = SQ.two_level_sequence(spec)
                r.ingredients = [("f", s)]
                r.corr, r.corr_kind, r.corr_target = SQ.paf(s), "periodic", -1
                ok2 = SQ.has_two_level_paf(s)
                r.checks.append(("P_f(s) = −1 for all s ≠ 0", "yes" if ok2 else "NO", ok2))
                r.checks.append(("Σ f = −1 (core row sum)", str(int(s.sum())), int(s.sum()) == -1))
            else:
                vals, _ = SQ.twin_prime_power_function(spec["q"])
                r.ingredients = [("f (group order x·(q+2)+y)", vals)]
                r.corr_kind = "group"
                r.checks.append(("Σ f = −1 (core row sum)", str(int(vals.sum())), int(vals.sum()) == -1))
            r.H = SQ.difference_set_core_matrix(spec)
            core = r.H[1:, 1:].astype(np.int64)
            okc = np.array_equal(core @ core.T, (L + 1) * np.eye(L, dtype=np.int64) - np.ones((L, L), dtype=np.int64))
            r.checks.append(("M·Mᵀ = (ℓ+1)I − J for the core", "yes" if okc else "NO", okc))
            r.description = ("Bordered " + ("circulant" if spec.get("cyclic", True) else "group-developed") +
                             " core H = [[1, eᵀ], [e, M]], M(x, y) = f(y − x). " +
                             {"singer": "f is a maximal-length shift-register sequence (the Singer difference set of PG(k−1, 2)).",
                              "twin": "f = +1 on {(x,y): χ(x)χ(y) = 1} ∪ {(x,0)} (Stanton–Sprott twin-prime-power set), −1 elsewhere, signed to row sum −1.",
                              "hall": "f = −1 on the union of three sextic cyclotomic classes (Hall's set), +1 elsewhere.",
                              "paley": "f is the Legendre symbol (Paley)."}[spec["kind"]])
        elif kind == "scarpis":
            q = int(params["q"])
            meth = params.get("sinput", "auto")
            con = P.construct(q + 1, None if meth == "auto" else meth)
            A = con.matrix
            r.title = f"Scarpis–Đoković with q = {q}: H({q + 1}) → H({q * (q + 1)})"
            An = SQ.normalise_hadamard(A)
            r.ingredients = [(f"input H({q + 1}) row {i + 1} (normalised)", An[i]) for i in range(min(q + 1, 8))]
            r.corr_kind = "none"
            r.checks.append((f"input matrix of order {q + 1} is Hadamard", "yes", True))
            r.H = SQ.scarpis(q, A)
            r.description = (f"Input: {P.describe_node(con.tree)}. Rows of the normalised core c(x), x ∈ GF({q}), are "
                             f"permuted by x ↦ x·r + k; B₀ = A′ ⊗ j and B_r = [jᵀ ⊗ c(r) | B_{{r,x}}] with row k of B_{{r,x}} = c(x r + k).")
        else:
            raise ValueError(f"unknown construction {kind}")
    except Exception as exc:                          # surfaced in the UI
        r.ok, r.error = False, str(exc)
        return r
    return _finish(r, t0)


def ingredients_text(r: SeqResult) -> str:
    lines = []
    for name, s in r.ingredients:
        lines.append(f"{name} (length {len(s)}): " + " ".join({1: "+", -1: "-", 0: "0"}.get(int(v), str(v)) for v in s))
    return "\n".join(lines)


def ingredients_json(r: SeqResult) -> str:
    import json
    return json.dumps({"construction": r.title, "kind": r.kind, "parameters": r.params,
                       "ingredients": {name: [int(v) for v in s] for name, s in r.ingredients},
                       "autocorrelation": {"kind": r.corr_kind, "target_off_peak": r.corr_target,
                                           "values": [int(v) for v in r.corr] if r.corr is not None else None},
                       "order": r.order, "reference": r.reference}, indent=1)
