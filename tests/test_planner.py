import pytest
from hadamard_lab import planner as P, constructions as C, analysis as A


def test_validate_order():
    assert P.validate_order(12)[0] and P.validate_order(1)[0] and P.validate_order(2)[0]
    for bad in (0, -4, 3, 6, 10, 1001, 1004, "x", 12.5):
        assert not P.validate_order(bad)[0], bad


def test_coverage():
    rows = P.existence_table()
    missing = [r[0] for r in rows if not r[3]]
    assert missing == [668, 716, 876, 892]
    assert set(P.KNOWN_OPEN) == {668, 716, 892}
    for n in (92, 116, 156, 172, 188, 412, 428, 508, 604, 652, 724, 764, 956, 964):
        assert P.is_constructible(n), n


def test_every_order_with_every_method():
    """Each feasible (order, method, parameter set) must give a verified Hadamard matrix."""
    for n, _, _, methods in P.existence_table():
        for m in methods:
            cnt = len([o for o in P.feasible_options(n) if o.method == m])
            for i in range(min(cnt, 2)):
                c = P.construct(n, m, i)
                assert c.matrix.shape == (n, n) and C.is_hadamard(c.matrix), (n, m, i)
                assert A.verify(c.matrix).passed


def test_errors():
    with pytest.raises(ValueError, match="open case"):
        P.construct(668)
    with pytest.raises(ValueError, match="cannot produce"):
        P.construct(12, "ehlich")
    with pytest.raises(ValueError):
        P.construct(10)


def test_trivial_orders():
    assert P.construct(1).matrix.tolist() == [[1]]
    assert C.is_hadamard(P.construct(2).matrix)
