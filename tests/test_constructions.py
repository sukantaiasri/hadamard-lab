import numpy as np
import pytest
from hadamard_lab import constructions as C
from hadamard_lab.gf import is_prime_power


def ok(H, n):
    return H.shape == (n, n) and C.is_hadamard(H)


def test_sylvester_paley():
    assert ok(C.sylvester(64), 64)
    for q in (3, 7, 11, 19, 23, 27, 343):
        assert ok(C.paley1(q), q + 1)
    for q in (5, 9, 13, 25, 81):
        assert ok(C.paley2(q), 2 * (q + 1))
    assert ok(C.paley2_general(5, C.sylvester(4)), 24)
    assert ok(C.paley2_general(9, C.sylvester(8)), 80)
    assert ok(C.paley2_general(5, C.paley1(11)), 72)


def test_paley1_is_skew():
    H = C.paley1(19).astype(int)
    assert (H + H.T == 2 * np.eye(20)).all()


@pytest.mark.parametrize("w", sorted(C.sequence_data()["williamson"]))
def test_williamson(w):
    assert ok(C.williamson(w), 4 * w)


def test_baumert_hall():
    assert ok(C.baumert_hall(1), 12)
    for w in sorted(C.sequence_data()["williamson"]):
        if 12 * w <= 1000:
            assert ok(C.baumert_hall(w), 12 * w)


def test_baumert_hall_array_is_orthogonal_design():
    rng = np.random.default_rng(1)
    for _ in range(5):
        v = rng.integers(-9, 10, 4).astype(float)
        M = C.baumert_hall_array(*[np.array([[x]]) for x in v]).astype(float)
        assert np.allclose(M @ M.T, 3 * (v ** 2).sum() * np.eye(12))


@pytest.mark.parametrize("k", [7, 11, 19, 27, 31])
def test_ehlich(k):
    assert ok(C.ehlich(k), (k - 1) ** 2)


def test_ehlich_printed_sign_fails():
    """Documents why the sign in the Ehlich formula must be negative."""
    k, m = 7, 5
    from hadamard_lab.gf import field
    A, G = field(m).jacobsthal().astype(int), field(k).jacobsthal().astype(int)
    K = np.kron(A, G) + np.kron(np.eye(m, dtype=int), np.ones((k, k), int) - np.eye(k, dtype=int)) + np.kron(np.ones((m, m), int), np.eye(k, dtype=int))
    H = np.ones((m * k + 1,) * 2, dtype=int)
    H[1:, 1:] = K
    assert not C.is_hadamard(H)


def test_miyamoto():
    for q in range(5, 250, 4):
        if not is_prime_power(q):
            continue
        n = q - 1
        K = None
        if n & (n - 1) == 0:
            K = C.sylvester(n)
        elif is_prime_power(n - 1) and (n - 1) % 4 == 3:
            K = C.paley1(n - 1)
        elif n % 4 == 0 and n // 4 in C.sequence_data()["williamson"]:
            K = C.williamson(n // 4)
        elif n % 2 == 0 and is_prime_power(n // 2 - 1) and (n // 2 - 1) % 4 == 1:
            K = C.paley2(n // 2 - 1)
        if K is not None:
            assert ok(C.miyamoto(q, K), 4 * q), q


def test_goethals_seidel_and_sds():
    d = C.sequence_data()
    for n, s in d["base"].items():
        assert ok(C.goethals_seidel_base(n), 4 * (2 * n + len(s[0]) - len(s[2])))
    for t in d["turyn"]:
        assert ok(C.goethals_seidel_turyn(t), 4 * (3 * t - 1))
    for N in d["sds"]:
        assert ok(C.supplementary_difference_sets(N), 4 * N)


def test_cooper_wallis():
    for w in C.sequence_data()["williamson"]:
        for L in C.t_sequence_lengths():
            if w >= 3 and 4 * w * L <= 1000:
                src, key = C._t_source(L)
                assert ok(C.cooper_wallis(w, src, key), 4 * w * L), (w, L)
