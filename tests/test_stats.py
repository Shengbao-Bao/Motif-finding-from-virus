"""Self-checks for :mod:`stats` (run: python tests/test_stats.py)."""

import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))

import stats


def brute_poisson_sf(x, lam):
    """P(X >= x) by direct summation (only usable for modest lam)."""
    total = 0.0
    for k in range(int(math.floor(x))):
        total += math.exp(-lam + k * math.log(lam) - math.lgamma(k + 1))
    return 1.0 - total


def check_poisson():
    cases = [(0, 10.0), (1, 10.0), (5, 10.0), (15, 10.0), (30, 10.0), (3, 0.5), (10, 3.0)]
    for x, lam in cases:
        analytic = stats.poisson_sf(x, lam)
        numeric = brute_poisson_sf(x, lam)
        assert abs(analytic - numeric) < 1e-9, (x, lam, analytic, numeric)
    # depletion must NOT be significant: P(X >= observed) ~ 1 when observed << lam
    assert stats.poisson_sf(100, 200) > 0.999, stats.poisson_sf(100, 200)
    assert stats.poisson_sf(200, 100) < 1e-15, stats.poisson_sf(200, 100)
    # large lambda must not overflow
    assert 0.0 < stats.poisson_sf(5000, 4000) < 1.0
    print("poisson_sf OK")


def check_bh():
    pvalues = [1e-8, 1e-4, 0.01, 0.02, 0.5, 0.9]
    qvalues = stats.bh_fdr(pvalues)
    assert all(q >= p for p, q in zip(pvalues, qvalues))
    assert all(q <= 1.0 for q in qvalues)
    # monotone: q-values must be non-decreasing in p (step-up property)
    ordered = sorted(zip(pvalues, qvalues), key=lambda item: item[0])
    for (_, q1), (_, q2) in zip(ordered, ordered[1:]):
        assert q2 >= q1 - 1e-12, ordered
    # all p=1 -> all q=1 (the direction bug used to break this)
    assert all(q == 1.0 for q in stats.bh_fdr([1.0] * 10))
    assert stats.bh_fdr([]) == []
    print("bh_fdr OK")


def check_markov():
    seqs = ["ACGTACGTAC", "ACGTACGTAC"]
    model = stats.MarkovModel.from_sequences(seqs)
    assert abs(sum(model.mono.values()) - 1.0) < 1e-9
    # a homopolymer never seen in the data must still get a non-zero probability
    assert model.prob("AAAAAA") > 0.0
    assert model.prob("ACGTAC") > model.prob("AAAAAA")
    # expected window counts scale linearly with the number of windows
    single = model.expected_windows("ACGT", 7)
    assert abs(model.expected_windows("ACGT", 14) - 2 * single) < 1e-12
    print("MarkovModel OK")


if __name__ == "__main__":
    check_poisson()
    check_bh()
    check_markov()
    print("all stats self-checks passed")
