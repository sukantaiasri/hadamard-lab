import numpy as np
import pytest
from hadamard_lab import mols, analysis as A, planner as P
from hadamard_lab.gf import prime_powers_up_to


@pytest.mark.parametrize("s", [2, 3, 4, 5, 7, 8, 9, 11, 16, 25, 27, 32, 49])
def test_mols_and_oa(s):
    assert mols.verify_mols_set(s)[0]
    assert mols.verify_oa(s)[0]
    assert (mols.oa_rows(s) == np.stack([mols.oa_column(s, c) for c in range(s + 1)], 1)).all()


def test_large_s_sampled():
    assert mols.verify_oa(997, max_pairs=10)[0]


def test_validate_s():
    assert mols.validate_s(9)[0]
    ok, msg = mols.validate_s(12)
    assert not ok and "11" in msg and "13" in msg
    assert not mols.validate_s(1000)[0] and not mols.validate_s(1)[0]


def test_lattice():
    assert mols.validate_lattice(25)[0] and not mols.validate_lattice(24)[0]
    s = mols.lattice_summary(5, 6)                  # balanced lattice = BIBD, E = k/(k+1)
    assert s["balanced"] and abs(s["average_efficiency"] - 5 / 6) < 1e-9
    s3 = mols.lattice_summary(6, 3)
    assert not s3["balanced"] and s3["each_treatment_once_per_replicate"]
    assert mols.max_replicates(6) == 3 and mols.max_replicates(7) == 8


def test_normalisation_and_properties():
    H = P.construct(52).matrix
    for mode, chk in (("semi", A.is_semi_normalised), ("full", A.is_normalised)):
        N = A.normalise(H, mode)
        assert chk(N) and A.verify(N).passed
    assert (A.normalise(H, "none") == H).all()
    D = A.plackett_burman(H)
    assert D.shape == (52, 51) and A.design_is_orthogonal(D)
    assert A.properties(P.construct(12).matrix).skew_type


def test_verify_detects_errors():
    H = P.construct(12).matrix.copy()
    H[3, 4] *= -1
    v = A.verify(H)
    assert not v.passed and v.max_row_deviation == 2


def test_export_formats():
    H = P.construct(4).matrix
    assert A.to_text(H, "sign").splitlines()[0] == "+ + + +"
    assert A.to_text(H, "binary").splitlines()[1] == "1 0 1 0"
    assert A.to_csv(H, "pm1").splitlines()[1] == "1,-1,1,-1"
