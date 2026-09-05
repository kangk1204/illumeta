#!/usr/bin/env python3
"""Han-Eskin RE2: a random-effects test whose null also forbids heterogeneity.

The problem it addresses
------------------------
This workflow pools at most five cohorts, and the two references already available are both
unsatisfactory at that size, in opposite directions. The standard normal reference is
anti-conservative: with inputs correct by construction it returns a genomic inflation factor
of 1.21 at k = 5 and 1.41 at k = 3, and calls thousands of null probes per 200,000. The
Knapp-Hartung t(k-1) reference is exactly calibrated everywhere tested -- within 0.02 of unity,
Kolmogorov-Smirnov deviation below 0.005, no false call anywhere -- and has no genome-wide
power at all: on real data it returns 0, 0 and 4 sites against 1,654, 5,130 and 4,650. That is
a property of a t distribution with four degrees of freedom, which cannot produce P-values
small enough for a Benjamini-Hochberg threshold over 300,000 tests.

Han and Eskin (2011) identified why classical random effects is conservative in this setting:
its null hypothesis implicitly allows heterogeneity to exist, so a heterogeneous alternative
is not far from it. RE2 tests the sharper null that the mean effect is zero *and* the
between-study variance is zero, against an alternative where either may differ. Rejecting that
null is a weaker claim than rejecting the classical one -- a significant RE2 result can mean a
real mean effect, real heterogeneity, or both -- and the module's output says so, because a
test that answers a compound question must not be read as if it answered the first part alone.

What is deliberately not claimed
--------------------------------
The asymptotic null of the RE2 statistic is a 50:50 mixture of chi-squared with one and two
degrees of freedom. Han and Eskin note that this approximation is poor for small numbers of
studies, which is precisely the regime here, and supply tabulated corrections instead. This
module computes the mixture P-value and labels it as asymptotic. It does not claim
calibration: the only calibration statement this project makes about any reference is the one
its sign-flipping null measures, and that measurement is run separately.
"""
from __future__ import annotations

import math
from typing import Sequence

# Maximising the profile likelihood over tau^2 is a one-dimensional problem on [0, inf).
# A golden-section search on a bounded interval is used rather than a derivative method
# because the likelihood is flat near zero for homogeneous inputs, which is the common case,
# and a Newton step there is numerically unhelpful.
_GOLDEN = (math.sqrt(5.0) - 1.0) / 2.0
_TOL = 1e-12
_MAX_ITER = 200


def _weighted_mean(effects, ses, valid, tau2):
    num = den = 0.0
    for e, s, ok in zip(effects, ses, valid):
        if not ok:
            continue
        w = 1.0 / (s * s + tau2)
        num += w * e
        den += w
    return (num / den if den > 0 else float("nan")), den


def _loglik(effects, ses, valid, tau2):
    """Profile log-likelihood at a given tau^2, with mu at its conditional maximum."""
    mu, den = _weighted_mean(effects, ses, valid, tau2)
    if not math.isfinite(mu) or den <= 0:
        return float("-inf")
    total = 0.0
    for e, s, ok in zip(effects, ses, valid):
        if not ok:
            continue
        var = s * s + tau2
        if var <= 0:
            return float("-inf")
        total += math.log(2.0 * math.pi * var) + ((e - mu) ** 2) / var
    return -0.5 * total


def _null_loglik(effects, ses, valid):
    """Log-likelihood under mu = 0 and tau^2 = 0."""
    total = 0.0
    for e, s, ok in zip(effects, ses, valid):
        if not ok:
            continue
        var = s * s
        if var <= 0:
            return float("-inf")
        total += math.log(2.0 * math.pi * var) + (e * e) / var
    return -0.5 * total


def _maximise_tau2(effects, ses, valid):
    """Maximum-likelihood tau^2 on [0, upper], where upper brackets any plausible value."""
    spread = [e for e, ok in zip(effects, valid) if ok]
    if len(spread) < 2:
        return 0.0
    mean = sum(spread) / len(spread)
    var = sum((x - mean) ** 2 for x in spread) / max(len(spread) - 1, 1)
    upper = max(var * 10.0, 1e-8)

    at_zero = _loglik(effects, ses, valid, 0.0)
    lo, hi = 0.0, upper
    a = hi - _GOLDEN * (hi - lo)
    b = lo + _GOLDEN * (hi - lo)
    fa, fb = _loglik(effects, ses, valid, a), _loglik(effects, ses, valid, b)
    for _ in range(_MAX_ITER):
        if hi - lo < _TOL:
            break
        if fa > fb:
            hi, b, fb = b, a, fa
            a = hi - _GOLDEN * (hi - lo)
            fa = _loglik(effects, ses, valid, a)
        else:
            lo, a, fa = a, b, fb
            b = lo + _GOLDEN * (hi - lo)
            fb = _loglik(effects, ses, valid, b)
    best = (lo + hi) / 2.0
    # tau^2 is bounded below by zero and the maximum often sits exactly on that boundary;
    # the search cannot represent a boundary solution, so it is checked explicitly.
    return best if _loglik(effects, ses, valid, best) > at_zero else 0.0


def _chi2_sf(stat: float, df: int) -> float:
    """Upper tail of chi-squared for the only two degrees of freedom this test needs."""
    if not math.isfinite(stat) or stat <= 0:
        return 1.0
    if df == 1:
        return math.erfc(math.sqrt(stat / 2.0))
    if df == 2:
        return math.exp(-stat / 2.0)
    raise ValueError(f"unsupported df: {df}")


def re2_one(
    effects: Sequence[float],
    ses: Sequence[float],
    valid: Sequence[bool],
) -> dict[str, float]:
    """RE2 statistic and its asymptotic mixture P-value for one feature.

    Returns ``re2_stat``, ``re2_p`` (asymptotic, see the module docstring), the maximum
    likelihood ``re2_tau2``, the mean under the alternative ``re2_mu``, and ``re2_k``.
    """
    usable = [
        bool(ok) and math.isfinite(e) and math.isfinite(s) and s > 0
        for e, s, ok in zip(effects, ses, valid)
    ]
    k = sum(usable)
    nan = float("nan")
    if k < 2:
        return {"re2_stat": nan, "re2_p": nan, "re2_tau2": nan, "re2_mu": nan, "re2_k": k}

    tau2 = _maximise_tau2(effects, ses, usable)
    alt = _loglik(effects, ses, usable, tau2)
    null = _null_loglik(effects, ses, usable)
    if not (math.isfinite(alt) and math.isfinite(null)):
        return {"re2_stat": nan, "re2_p": nan, "re2_tau2": nan, "re2_mu": nan, "re2_k": k}

    stat = 2.0 * (alt - null)
    # The likelihood under the alternative can only exceed the null's up to numerical noise;
    # a small negative statistic is that noise and is clamped rather than propagated.
    if stat < 0.0:
        stat = 0.0
    p = 0.5 * _chi2_sf(stat, 1) + 0.5 * _chi2_sf(stat, 2)
    mu, _ = _weighted_mean(effects, ses, usable, tau2)
    return {
        "re2_stat": stat,
        "re2_p": min(1.0, max(0.0, p)),
        "re2_tau2": tau2,
        "re2_mu": mu,
        "re2_k": k,
    }
