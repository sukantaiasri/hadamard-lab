import numpy as np
import pytest
from hadamard_lab.gf import field, prime_power, prime_powers_up_to


def test_prime_power():
    assert prime_power(27) == (3, 3) and prime_power(1) is None and prime_power(12) is None
    assert prime_power(997) == (997, 1) and len(prime_powers_up_to(1000)) == 193


@pytest.mark.parametrize("q", prime_powers_up_to(1000)[::7] + [4, 8, 9, 16, 25, 27, 32, 64, 81, 125, 128, 243, 625, 729, 961])
def test_field_axioms(q):
    f = field(q)
    r = np.random.default_rng(q)
    a, b, c = (r.integers(0, q, 100) for _ in range(3))
    assert (f.add(a, b) == f.add(b, a)).all()
    assert (f.mul(a, f.add(b, c)) == f.add(f.mul(a, b), f.mul(a, c))).all()
    assert (f.add(a, f.neg(a)) == 0).all()
    assert sorted(f.exp.tolist()) == list(range(1, q))      # g generates the multiplicative group


@pytest.mark.parametrize("q", [5, 9, 13, 25, 29, 49, 81, 3, 7, 11, 27, 43, 243])
def test_jacobsthal(q):
    Q = field(q).jacobsthal().astype(int)
    assert (Q @ Q.T == q * np.eye(q) - np.ones((q, q))).all()
    assert (Q == (Q.T if q % 4 == 1 else -Q.T)).all()
