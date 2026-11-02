"""Core statistics for k-mer enrichment analysis (standard library only).

1. :class:`MarkovModel` -- order-1 (dinucleotide) background probabilities.
2. :func:`poisson_sf`   -- right tail P(X >= x) for X ~ Poisson(lam), via the
   regularised incomplete gamma function (accepts real, i.e. weighted, counts).
3. :func:`bh_fdr`       -- Benjamini-Hochberg FDR with the correct monotone
   step-up direction.
"""

from __future__ import annotations

import math

BASES = "ACGT"
TINY = 1e-300


class MarkovModel:
    """Order-1 Markov background, optionally weighted.

    ``P(motif) = f(motif[0]) * prod_i f(motif[i] | motif[i-1])``

    Dinucleotide frequencies are smoothed with a Dirichlet prior
    (``alpha`` pseudo-counts split evenly over the four bases) so that unseen
    dinucleotides never produce a zero probability.
    """

    def __init__(self, mono_counts, di_counts, total_bases, alpha=1.0):
        self.alpha = float(alpha)
        total = float(total_bases)
        if total <= 0:
            self.mono = {b: 0.25 for b in BASES}
        else:
            self.mono = {b: mono_counts.get(b, 0.0) / total for b in BASES}
        self.cond = {}
        for a in BASES:
            denom = mono_counts.get(a, 0.0) + self.alpha
            for b in BASES:
                num = di_counts.get(a + b, 0.0) + self.alpha * 0.25
                self.cond[(a, b)] = num / denom if denom > 0 else 0.25

    @classmethod
    def from_sequences(cls, sequences, weights=None, alpha=1.0):
        """Estimate the model from an iterable of sequences.

        ``weights`` lets every sequence carry a cluster weight (1/cluster_size)
        so that over-sampled lineages do not dominate the background either.
        """
        weights = list(weights) if weights is not None else [1.0] * len(sequences)
        mono_counts = {b: 0.0 for b in BASES}
        di_counts = {}
        total = 0.0
        for seq, weight in zip(sequences, weights):
            if weight <= 0:
                continue
            for i, base in enumerate(seq):
                mono_counts[base] = mono_counts.get(base, 0.0) + weight
                total += weight
                if i:
                    key = seq[i - 1 : i + 1]
                    di_counts[key] = di_counts.get(key, 0.0) + weight
        return cls(mono_counts, di_counts, total, alpha=alpha)

    def prob(self, motif: str) -> float:
        """Probability of observing ``motif`` under the background model."""
        if not motif:
            return 1.0
        p = self.mono.get(motif[0], 0.25)
        for i in range(1, len(motif)):
            p *= self.cond.get((motif[i - 1], motif[i]), 0.25)
        return p

    def expected_windows(self, motif: str, weighted_slots: float) -> float:
        """Expected number of windows equal to ``motif``."""
        return weighted_slots * self.prob(motif)


# --------------------------------------------------------------------------- #
# Poisson tail via the regularised incomplete gamma function
# --------------------------------------------------------------------------- #

def poisson_sf(x: float, lam: float) -> float:
    """Return ``P(X >= x)`` for ``X ~ Poisson(lam)``.

    Two standard identities connect the Poisson tail to the regularised
    incomplete gamma functions::

        P(X <= k) = Q(k + 1, lam)      (upper regularised gamma)
        P(X >= k) = P(k, lam)          (lower regularised gamma)

    The *lower* branch is the survival function, so ``poisson_sf`` uses
    ``_gammp``. Getting this backwards silently turns "depleted" motifs into
    "highly significant enrichment" -- exactly the failure the self-test in
    ``tests/test_stats.py`` guards against.

    This form is numerically stable for large ``lam`` and also accepts real
    (cluster-weighted) counts, which a Poisson sum cannot.
    """
    if lam <= 0.0:
        return 1.0 if x <= 0 else 0.0
    if x <= 0.0:
        return 1.0
    value = _gammp(x, lam)
    return value if value > 0.0 else 1e-300


def _gser(a: float, x: float) -> float:
    """Series expansion of the regularised lower incomplete gamma P(a, x)."""
    ap = a
    total = 1.0 / a
    delta = total
    for _ in range(1000):
        ap += 1.0
        delta *= x / ap
        total += delta
        if abs(delta) < abs(total) * 1e-15:
            break
    return total * math.exp(-x + a * math.log(x) - math.lgamma(a))


def _gcf(a: float, x: float) -> float:
    """Continued fraction (Lentz) for the regularised upper incomplete gamma."""
    b = x + 1.0 - a
    c = 1.0 / TINY
    d = 1.0 / b if b != 0.0 else 1.0 / TINY
    h = d
    for i in range(1, 1000):
        an = -i * (i - a)
        b += 2.0
        d = an * d + b
        if abs(d) < TINY:
            d = TINY
        c = b + an / c
        if abs(c) < TINY:
            c = TINY
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < 1e-15:
            break
    return math.exp(-x + a * math.log(x) - math.lgamma(a)) * h


def _gammq(a: float, x: float) -> float:
    if x <= 0.0:
        return 1.0
    if x < a + 1.0:
        return 1.0 - _gser(a, x)
    return _gcf(a, x)


def _gammp(a: float, x: float) -> float:
    """Regularised lower incomplete gamma P(a, x) = 1 - Q(a, x)."""
    if x <= 0.0:
        return 0.0
    if x < a + 1.0:
        return _gser(a, x)
    return 1.0 - _gcf(a, x)


# --------------------------------------------------------------------------- #
# Multiple testing
# --------------------------------------------------------------------------- #

def bh_fdr(pvalues):
    """Benjamini-Hochberg FDR, returning a q-value per input p-value.

    ``q_(i) = min_{j >= i} ( p_(j) * m / j )``: the running minimum has to be
    taken from the largest p-value downwards (monotone step-up). Walking the
    other way pins every q-value to 1.0 -- a bug worth remembering.
    """
    m = len(pvalues)
    if m == 0:
        return []
    order = sorted(range(m), key=lambda i: pvalues[i])
    qvalues = [1.0] * m
    prev = 1.0
    for rank in range(m, 0, -1):
        idx = order[rank - 1]
        value = pvalues[idx] * m / rank
        if value > prev:
            value = prev
        prev = value
        qvalues[idx] = min(1.0, value)
    return qvalues
