"""Golay pairs, two-level sequences, Legendre pairs/quadruples, Scarpis and the sequence-lab back-end."""
import itertools

import numpy as np
import pytest

from hadamard_lab import constructions as C, planner as P, seqlab as SL, sequences as S


def test_golay_seeds_and_product():
    for g, (a, b) in S.GOLAY_SEEDS.items():
        assert S.is_golay_pair(a, b), g
    for (g1, p1), (g2, p2) in itertools.product(S.GOLAY_SEEDS.items(), repeat=2):
        x, y = S.golay_product(*p1, *p2)
        assert len(x) == g1 * g2 and S.is_golay_pair(x, y), (g1, g2)


@pytest.mark.parametrize("g", S.golay_lengths(500))
def test_golay_lengths(g):
    a, b = S.golay_pair(g)
    assert S.is_golay_pair(a, b)
    assert C.is_hadamard(S.golay_two_circulant(a, b))
    assert S.is_t_sequence_set(S.t_sequences_from_golay(a, b))
    ts = S.t_sequences_from_base(*S.base_sequences_from_golay(a, b))
    assert len(ts[0]) == 2 * g + 1 and S.is_t_sequence_set(ts)


def test_golay_lengths_set():
    assert S.golay_lengths(60) == [2, 4, 8, 10, 16, 20, 26, 32, 40, 52]


def test_two_level_sources():
    src = S.two_level_sources(1000)
    assert {3, 7, 15, 31, 35, 43, 63, 99, 127, 143, 223, 255, 283, 323, 511, 575, 675, 783, 811, 899} == set(src)
    for L, specs in src.items():
        for spec in specs:
            H = S.difference_set_core_matrix(spec)
            assert H.shape == (L + 1, L + 1) and C.is_hadamard(H), spec
            if spec["cyclic"]:
                s = S.two_level_sequence(spec)
                assert S.has_two_level_paf(s) and s.sum() == -1


def test_m_sequence_and_hall():
    assert S.has_two_level_paf(S.m_sequence(6)) and len(S.m_sequence(6)) == 63
    assert S.has_two_level_paf(S.hall_sextic_sequence(43))
    assert S.has_two_level_paf(S.twin_prime_sequence(11))
    with pytest.raises(ValueError):
        S.twin_prime_sequence(9)          # GF(9) x GF(11) is not cyclic


@pytest.mark.parametrize("p", [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97])
def test_legendre_pairs_from_primes(p):
    a, b = S.legendre_pair_paley(p)
    assert S.is_legendre_pair(a, b)
    H = S.two_circulant_core(a, b)
    assert H.shape == (2 * p + 2, 2 * p + 2) and C.is_hadamard(H)


def test_legendre_pair_sources_cover_every_prime_length():
    src = C.legendre_pair_sources(1000)
    for L in range(3, 499, 2):
        if S.is_prime(L):
            assert L in src and any(s["kind"] == "paley" for s in src[L])
    for L, specs in src.items():
        for spec in specs:
            a, b = C.legendre_pair(spec)
            assert S.is_legendre_pair(a, b) and C.is_hadamard(S.two_circulant_core(a, b)), spec


def test_legendre_quadruples():
    for m, seqs in C.sequence_data()["lquads"].items():
        H = S.legendre_quadruple_matrix(seqs)
        assert H.shape == (4 * m + 4, 4 * m + 4) and C.is_hadamard(H), m
    with pytest.raises(ValueError):
        S.normalise_quadruple([[1, 1], [1, 1], [1, 1], [1, 1]])


def test_scarpis_every_q_and_any_input():
    for q in (3, 7, 11, 19, 23, 27, 31):
        H = S.scarpis(q, C.paley1(q))
        assert H.shape == (q * (q + 1), q * (q + 1)) and C.is_hadamard(H), q
    assert C.is_hadamard(S.scarpis(7, C.sylvester(8)))
    assert C.is_hadamard(S.scarpis(3, -C.sylvester(4)[::-1]))      # un-normalised input is normalised internally
    assert C.is_hadamard(S.scarpis(11, C.williamson(3)))
    with pytest.raises(ValueError):
        S.scarpis(5, C.paley2(5)[:6, :6])


def test_t_sequence_sources_include_golay_and_stored():
    src = C.t_sequence_sources()
    assert ("golay", 10) in src[10] and ("golay2", 10) in src[21] and ("tseq", 73) in src[73]
    for L in (10, 21, 73):
        for s, k in src[L]:
            assert S.is_t_sequence_set(C._t_sequences(s, k)), (s, k)


def test_seqlab_runs_every_kind():
    cases = [("golay", {"g": "20"}), ("tseq", {"t": "21", "tsrc": "golay2:10", "w": "0"}),
             ("tseq", {"t": "7", "tsrc": "base:3", "w": "3"}), ("legendre", {"l": "21", "lsrc": "stored:l=21"}),
             ("legendre", {"l": "15", "lsrc": "singer:cyclic=True,k=4"}), ("lquad", {"m": "12"}),
             ("dscore", {"l": "63", "dsrc": "twin:cyclic=False,q=7"}), ("dscore", {"l": "31", "dsrc": "hall:cyclic=True,p=31"}),
             ("scarpis", {"q": "11", "sinput": "paley1"})]
    for kind, prm in cases:
        r = SL.run(kind, prm)
        assert r.ok and C.is_hadamard(r.H), (kind, prm, r.error)
        assert SL.ingredients_text(r) and SL.ingredients_json(r)
    bad = SL.run("golay", {"g": "12"})
    assert not bad.ok and bad.error


def test_seqlab_choices_consistent_with_planner():
    assert set(map(int, SL.golay_choices())) == set(S.golay_lengths(500))
    assert set(map(int, SL.scarpis_choices())) == {3, 7, 11, 19, 23, 27, 31}
    for L in map(int, SL.legendre_choices()):
        assert SL.legendre_source_choices(L)
    for L in map(int, SL.dscore_choices()):
        assert SL.dscore_source_choices(L)
    assert "7" in SL.scarpis_input_choices(7) or "sylvester" in SL.scarpis_input_choices(7)


def test_thaw_freeze_spec_roundtrip():
    for spec in ({"kind": "twin", "q": 11, "cyclic": True}, {"kind": "paley", "p": 13}, {"kind": "stored", "l": 21}):
        assert C.thaw_spec(C._freeze_spec(spec)) == spec
