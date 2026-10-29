"""The bundled sequence data must satisfy the defining equations (guards against import/parse errors)."""
import numpy as np
from hadamard_lab import constructions as C


def test_williamson_rows():
    for n, rows in C.sequence_data()["williamson"].items():
        rows = [np.array(r) for r in rows]
        assert all(len(r) == n and set(r) <= {-1, 1} for r in rows)
        assert all((r == np.roll(r[::-1], 1)).all() for r in rows), f"order {n} not symmetric"
        tot = sum(C.periodic_autocorrelation(r) for r in rows)
        assert tot[0] == 4 * n and not tot[1:].any(), n


def test_base_sequences():
    for n, seqs in C.sequence_data()["base"].items():
        L = [len(s) for s in seqs]
        assert L[0] == L[1] and L[2] == L[3] == n
        m = max(L)
        tot = np.zeros(m, dtype=int)
        for s in seqs:
            a = C.aperiodic_autocorrelation(s)
            tot[: len(a)] += a
        assert not tot[1:].any(), n


def test_turyn_sequences():
    for n, (X, Y, Z, W) in C.sequence_data()["turyn"].items():
        assert (len(X), len(Y), len(Z), len(W)) == (n, n, n, n - 1)
        tot = np.zeros(n, dtype=int)
        for s, wgt in ((X, 1), (Y, 1), (Z, 2), (W, 2)):
            a = C.aperiodic_autocorrelation(s)
            tot[: len(a)] += wgt * a
        assert not tot[1:].any(), n


def test_supplementary_difference_sets():
    for N, sets in C.sequence_data()["sds"].items():
        tot = np.zeros(N, dtype=int)
        for s in sets:
            ind = np.zeros(N, dtype=int)
            ind[s] = 1
            assert len(set(s)) == len(s)
            tot += np.array([int(np.dot(ind, np.roll(ind, -d))) for d in range(N)])
        lam = sum(len(s) for s in sets) - N
        assert (tot[1:] == lam).all(), N


def test_t_sequences():
    d = C.sequence_data()
    for n in d["base"]:
        assert C.is_t_sequence_set(C.t_sequences_from_base(*d["base"][n]))
    for t in d["turyn"]:
        assert C.is_t_sequence_set(C.t_sequences_from_base(*C.base_from_turyn(*d["turyn"][t])))
