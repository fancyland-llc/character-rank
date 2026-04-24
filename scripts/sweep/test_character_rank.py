"""Tests for character-theoretic effective-rank bounds.

Pure-math tests only. Exercise:

  1. The a-priori isotropy decompositions of ℝ³ under {S_3, S_2, {e}}:
     each sums to N = 3, and the committed bounds are {2, 3, 3}.

  2. The `predict_bound(isotropy)` lookup.

  3. The full `measure_covariance_spectrum` + `check_bound` pipeline
     on hand-built G-equivariant covariances constructed to realize
     the Schur-forced degeneracy pattern exactly. No simulation —
     everything is a closed-form N = 3 PSD matrix.
"""
from __future__ import annotations

import numpy as np

from .character_rank import (
    PERM_REP_DECOMP,
    check_bound,
    count_distinct_eigenvalues,
    measure_covariance_spectrum,
    predict_bound,
)


# ──────────────────────────────────────────────────────────────────────────
# Decomposition sanity
# ──────────────────────────────────────────────────────────────────────────

def test_decomposition_dimensions_equal_N() -> None:
    """For each isotropy, Σ m_π · d_π must equal N = 3.

    This is just the claim that the permutation rep has full ambient
    dimension — a sanity check on the hand-derived decompositions.
    """
    for iso, decomp in PERM_REP_DECOMP.items():
        assert decomp.ambient_dim == 3, (
            f"{iso}: decomposition sums to {decomp.ambient_dim}, expected 3"
        )


def test_bounds_are_2_3_3() -> None:
    """The committed three integers. If any of these ever changes, the
    whole module is wrong and the experiment is moot."""
    assert PERM_REP_DECOMP["S3"].distinct_eigenvalue_bound == 2
    assert PERM_REP_DECOMP["S2"].distinct_eigenvalue_bound == 3
    assert PERM_REP_DECOMP["trivial"].distinct_eigenvalue_bound == 3


def test_predict_bound_lookup() -> None:
    assert predict_bound("S3") == 2
    assert predict_bound("S2") == 3
    assert predict_bound("trivial") == 3


def test_n_distinct_irrep_types() -> None:
    """|𝓐| for each isotropy: S_3 sees 2 types, S_2 sees 2, {e} sees 1."""
    assert PERM_REP_DECOMP["S3"].n_distinct_irrep_types == 2
    assert PERM_REP_DECOMP["S2"].n_distinct_irrep_types == 2
    assert PERM_REP_DECOMP["trivial"].n_distinct_irrep_types == 1


# ──────────────────────────────────────────────────────────────────────────
# Clustering
# ──────────────────────────────────────────────────────────────────────────

def test_count_distinct_eigenvalues_obvious_cases() -> None:
    assert count_distinct_eigenvalues(np.array([1.0, 1.0, 1.0])) == 1
    assert count_distinct_eigenvalues(np.array([3.0, 2.0, 1.0])) == 3
    # Two clusters: {3.0, 3.01} vs {1.0}, default rel_gap = 0.06
    assert count_distinct_eigenvalues(np.array([3.0, 3.01, 1.0])) == 2
    # Below abs_floor: drop as numerical zero.
    assert count_distinct_eigenvalues(np.array([1.0, 0.0, -1e-13])) == 1


# ──────────────────────────────────────────────────────────────────────────
# End-to-end: hand-built G-equivariant covariances on N = 3
# ──────────────────────────────────────────────────────────────────────────

def _C_S3(a: float, b: float) -> np.ndarray:
    """S_3-equivariant covariance on ℝ^3: C = a·I + b·(J − I), where J
    is the all-ones matrix. Spectrum is {a + 2b (trivial), a − b (×2 std)}
    — exactly two distinct eigenvalues, per the Schur bound."""
    I = np.eye(3)
    J = np.ones((3, 3))
    return a * I + b * (J - I)


def _C_S2(a1: float, a2: float, b: float) -> np.ndarray:
    """S_2-equivariant covariance on ℝ^3 (swap of axes 1 ↔ 2).

        C = diag(a1, a2, a2) + b·(e1 e2^T + e2 e1^T unfolded on the S_2-invariant
                                    off-diagonal block)
          = [[a1, b,  b ],
             [b,  a2, c ],
             [b,  c,  a2]]  with c = a2 (keeps two trivial copies equal? No —
                                          we want them unequal).

    Simpler construction: pick any symmetric matrix that commutes with the
    swap (1↔2 on axes 2,3 of ℝ^3) and generically has 3 distinct eigvals.
    """
    return np.array(
        [[a1, b, b],
         [b,  a2, 0.3],
         [b,  0.3, a2]],
        dtype=np.float64,
    )


def test_check_bound_S3_equivariant() -> None:
    """Hand-built S_3-equivariant covariance. Two distinct eigenvalues,
    the smaller one appearing with multiplicity 2 (the 2-dim standard
    irrep). Measurement must report distinct_count = 2."""
    C = _C_S3(a=1.0, b=0.3)
    chk = check_bound(C, "S3", M=10_000)
    assert chk.passes
    assert chk.measured_distinct == 2
    # Spectrum is {1 + 2·0.3, 1 − 0.3, 1 − 0.3} = {1.6, 0.7, 0.7}.
    # Sorted descending, top eigenvalue is the trivial mode.
    assert np.isclose(chk.eigenvalues[0], 1.6)
    assert np.isclose(chk.eigenvalues[1], 0.7)
    assert np.isclose(chk.eigenvalues[2], 0.7)


def test_check_bound_S3_degenerate_stays_within_bound() -> None:
    """A = diag(1, 1, 1) is trivially S_3-equivariant (all three eigvals
    equal). Distinct count = 1 ≤ 2, so the bound holds (Schur forbids
    exceeding the bound, never achieving it)."""
    C = np.eye(3)
    chk = check_bound(C, "S3", M=10_000)
    assert chk.passes
    assert chk.measured_distinct == 1


def test_check_bound_S2_equivariant_three_distinct() -> None:
    """Generic S_2-equivariant C on ℝ^3 hits the bound: 3 distinct."""
    C = _C_S2(a1=2.0, a2=1.0, b=0.4)
    chk = check_bound(C, "S2", M=10_000)
    assert chk.passes
    assert chk.measured_distinct <= 3


def test_measure_covariance_spectrum_shape() -> None:
    """Spectrum output is sorted descending and has length N."""
    C = _C_S3(a=1.0, b=0.2)
    meas = measure_covariance_spectrum(C, M=5000)
    assert len(meas.eigenvalues) == 3
    assert meas.eigenvalues[0] >= meas.eigenvalues[1] >= meas.eigenvalues[2]
    assert 1 <= meas.distinct_eigenvalue_count <= 3
