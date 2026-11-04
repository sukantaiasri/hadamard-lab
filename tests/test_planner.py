import pytest
from hadamard_lab import planner as P, constructions as C, analysis as A


def test_validate_order():
    assert P.validate_order(12)[0] and P.validate_order(1)[0] and P.validate_order(2)[0]
    for bad in (0, -4, 3, 6, 10, 1001, 1004, "x", 12.5):
        assert not P.validate_order(bad)[0], bad


def test_coverage():
    rows = P.existence_table()
    missing = [r[0] for r in rows if not r[3]]
    assert missing == []                      # every admissible order <= 1000 (since version 2.0)
    assert P.KNOWN_OPEN == () and P.KNOWN_BUT_NOT_IMPLEMENTED == ()
    for n in (92, 116, 156, 172, 188, 412, 428, 508, 604, 652, 668, 716, 724, 764, 876, 892, 956, 964):
        assert P.is_constructible(n), n


def test_former_open_orders_use_2026_data():
    assert P.construct(668).tree.method == "lquad" and P.construct(668).tree.params["m"] == 166
    assert P.construct(716).tree.method == "lquad" and P.construct(716).tree.params["m"] == 178
    assert "sds" in P.feasible_methods(892) and P.construct(892, "sds").tree.params["N"] == 223
    assert "cooper_wallis" in P.feasible_methods(876)
    c = P.construct(876, "cooper_wallis")
    assert c.tree.params["L"] == 73 and c.tree.params["w"] == 3 and C.is_hadamard(c.matrix)


def test_new_methods_reach_expected_orders():
    rows = {n: ms for n, _, _, ms in P.existence_table()}
    assert [n for n in rows if "scarpis" in rows[n]] == [12, 56, 132, 380, 552, 756, 992]
    assert [n for n in rows if "golay2c" in rows[n]][:6] == [4, 8, 16, 20, 32, 40]
    assert all("legendre" in rows[2 * p + 2] for p in (3, 5, 7, 11, 13) if 2 * p + 2 in rows)
    assert "dscore" in rows[100] and "dscore" in rows[128] and "dscore" in rows[900]
    assert "lquad" in rows[44] and "gs_tseq" in rows[292] and "gs_golay" in rows[84]


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
    with pytest.raises(ValueError, match="cannot produce"):
        P.construct(12, "ehlich")
    with pytest.raises(ValueError):
        P.construct(10)


def test_trivial_orders():
    assert P.construct(1).matrix.tolist() == [[1]]
    assert C.is_hadamard(P.construct(2).matrix)
