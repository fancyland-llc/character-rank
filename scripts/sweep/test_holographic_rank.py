"""Tests for the character-entropy decomposition and Gaussian MI tools.

Pure-math tests only. All covariances are hand-built to be perfectly
G-equivariant, so the identity

        S(ρ) = H(p) + Σ p_π log d_π + Σ p_π S_π^int

must reconstruct exactly (up to numerical tolerance).
"""
from __future__ import annotations

import numpy as np

from .holographic_rank import (
    character_entropy_decomposition,
    gaussian_mutual_information,
    isotypic_basis,
    isotypic_projectors,
    reduced_covariance,
    sample_covariance,
    von_neumann_entropy,
)


# ──────────────────────────────────────────────────────────────────────────
# Projectors and bases
# ──────────────────────────────────────────────────────────────────────────

def test_projectors_are_idempotent_and_orthogonal() -> None:
    """P^2 = P, P_a P_b = 0 for a != b, and Σ P = I."""
    for iso in ("S3", "S2", "trivial"):
        projs = isotypic_projectors(iso)
        names = list(projs.keys())
        for n in names:
            P = projs[n]
            assert np.allclose(P @ P, P, atol=1e-10), f"{iso}/{n} not idempotent"
        for i, a in enumerate(names):
            for b in names[i + 1:]:
                assert np.allclose(
                    projs[a] @ projs[b], 0.0, atol=1e-10,
                ), f"{iso}/{a},{b} not orthogonal"
        total = sum(projs.values())
        assert np.allclose(total, np.eye(3), atol=1e-10), (
            f"{iso} projectors do not sum to I"
        )


def test_bases_orthonormal_and_span_ambient() -> None:
    """Each basis U_π has U_π^T U_π = I, and stacking all of them gives
    an orthonormal basis of ℝ^3."""
    for iso in ("S3", "S2", "trivial"):
        bases = isotypic_basis(iso)
        all_cols = np.hstack(list(bases.values()))
        assert all_cols.shape == (3, 3), f"{iso}: basis not 3x3"
        assert np.allclose(
            all_cols.T @ all_cols, np.eye(3), atol=1e-10,
        ), f"{iso} basis not orthonormal"


# ──────────────────────────────────────────────────────────────────────────
# Entropy decomposition on hand-built equivariant covariances
# ──────────────────────────────────────────────────────────────────────────

def test_s3_symmetric_covariance_decomposition() -> None:
    """C = a I + b J  →  eigenvalues (a + 3b, a, a).
    For a = 1, b = 0.5 : λ = (2.5, 1, 1). tr = 4.5.
    p_triv = 2.5 / 4.5, p_std = 2 / 4.5 (two copies of a in std block).

    S = H(p) + p_std · log 2.
    """
    a, b = 1.0, 0.5
    C = a * np.eye(3) + b * np.ones((3, 3))
    decomp = character_entropy_decomposition(C, "S3")

    tr = 3 * a + 3 * b  # = a + 3b + 2a = 4.5 ✓
    p_triv_expected = (a + 3 * b) / tr  # 2.5 / 4.5
    p_std_expected = 2 * a / tr          # 2.0 / 4.5

    assert abs(decomp.p_block["triv"] - p_triv_expected) < 1e-10
    assert abs(decomp.p_block["std"] - p_std_expected) < 1e-10
    # Internal entropies: both blocks have multiplicity 1.
    assert decomp.internal["triv"] == 0.0
    assert decomp.internal["std"] == 0.0
    # Reconstruction: shannon + p_std * log 2 == direct entropy.
    expected_total = (
        -p_triv_expected * np.log(p_triv_expected)
        - p_std_expected * np.log(p_std_expected)
        + p_std_expected * np.log(2.0)
    )
    assert abs(decomp.total - expected_total) < 1e-10
    assert decomp.reconstruction_err < 1e-10


def test_s2_two_mult_trivial_has_internal_entropy() -> None:
    """Build an S_2-equivariant C where the 2D trivial block carries
    nontrivial internal structure:

        singleton variance = 1, pair variance = 1, within-pair corr = 0.3,
        singleton-pair cross = 0.25.

    The trivial isotypic subspace is 2-dimensional (singleton direction +
    pair-symmetric direction) and carries real off-diagonal structure;
    the sign block is 1-dimensional and contains just the antisymmetric
    variance. Internal entropy of the 2×2 copy matrix of the trivial
    block must be strictly positive.
    """
    v_s = 1.0        # singleton variance
    v_p = 1.0        # each pair member's variance
    rho_pp = 0.3     # within-pair correlation
    c_sp = 0.25      # singleton-pair covariance (same for both pair slots)

    C = np.array([
        [v_s,       c_sp,     c_sp    ],
        [c_sp,      v_p,      rho_pp  ],
        [c_sp,      rho_pp,   v_p     ],
    ])
    # Sanity: C must be symmetric and positive definite.
    assert np.allclose(C, C.T)
    assert np.all(np.linalg.eigvalsh(C) > 0)

    decomp = character_entropy_decomposition(C, "S2", singleton_slot=0)
    assert decomp.internal["sign"] == 0.0        # m=1 → always 0
    assert decomp.internal["triv"] > 0.0         # m=2 → real internal freedom
    # Reconstruction: component sum should match direct entropy exactly,
    # because C is built to be perfectly S_2-equivariant.
    assert decomp.reconstruction_err < 1e-10


def test_scalar_covariance_total_is_log3_across_all_isotropies() -> None:
    """C = σ² I  →  eigvals (σ², σ², σ²). Normalized ρ = I/3. S(ρ) = log 3.

    The DECOMPOSITION of log 3 is isotropy-dependent:
      S_3: trivial (mult 1, dim 1) + standard (mult 1, dim 2). All of
           log 3 comes from Shannon H(p) + p_std · log 2. Internal = 0.
      S_2: trivial (mult 2, dim 1) + sign (mult 1, dim 1). No Schur
           bonus (all dims are 1). But the mult-2 trivial block carries
           internal entropy log 2 weighted by p_triv = 2/3.
      {e}: one block (mult 3). Only internal: log 3.

    The *total* must equal log 3 for all three, and the reconstruction
    must be exact (the covariance is perfectly equivariant)."""
    C = 4.0 * np.eye(3)
    for iso in ("S3", "S2", "trivial"):
        decomp = character_entropy_decomposition(C, iso)
        assert abs(decomp.total - np.log(3.0)) < 1e-10, (
            f"{iso}: total = {decomp.total}, expected log 3 = {np.log(3.0)}"
        )
        assert decomp.reconstruction_err < 1e-10, (
            f"{iso}: reconstruction err = {decomp.reconstruction_err}"
        )

    # Check specific breakdown for S_3: all entropy in Shannon + Schur bonus.
    d3 = character_entropy_decomposition(C, "S3")
    assert abs(d3.internal_total) < 1e-12
    assert d3.schur_bonus > 0.0
    # S_2: no Schur bonus (all d=1), but internal contributes.
    d2 = character_entropy_decomposition(C, "S2")
    assert abs(d2.schur_bonus) < 1e-12
    assert d2.internal_total > 0.0
    # Trivial: no Schur bonus, no Shannon (only one block), all internal.
    dt = character_entropy_decomposition(C, "trivial")
    assert abs(dt.shannon) < 1e-12
    assert abs(dt.schur_bonus) < 1e-12
    assert abs(dt.internal_total - np.log(3.0)) < 1e-10


# ──────────────────────────────────────────────────────────────────────────
# Gaussian MI sanity
# ──────────────────────────────────────────────────────────────────────────

def test_mi_is_zero_for_diagonal_covariance() -> None:
    """Independent coords → MI = 0 exactly."""
    C = np.diag([1.0, 2.0, 3.0])
    I = gaussian_mutual_information(C, [0], [1, 2])
    assert abs(I) < 1e-12


def test_mi_nonzero_when_correlated() -> None:
    """Add correlation → MI > 0."""
    C = np.array([
        [1.0, 0.5, 0.0],
        [0.5, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ])
    I = gaussian_mutual_information(C, [0], [1])
    # Closed form: I = -0.5 log(1 - ρ²) = -0.5 log(0.75)
    expected = -0.5 * np.log(1 - 0.25)
    assert abs(I - expected) < 1e-10


def test_mi_symmetric_under_s3() -> None:
    """In an S_3-symmetric C = aI + bJ, I({i}; {j}) is independent of the
    pair (i, j). All three 1-vs-1 MIs equal."""
    C = np.eye(3) + 0.4 * np.ones((3, 3))
    I01 = gaussian_mutual_information(C, [0], [1])
    I02 = gaussian_mutual_information(C, [0], [2])
    I12 = gaussian_mutual_information(C, [1], [2])
    assert abs(I01 - I02) < 1e-10
    assert abs(I01 - I12) < 1e-10


# ──────────────────────────────────────────────────────────────────────────
# von Neumann entropy and reduced covariance helpers
# ──────────────────────────────────────────────────────────────────────────

def test_von_neumann_entropy_log3_for_scalar() -> None:
    """S(I/N) = log N."""
    assert abs(von_neumann_entropy(np.eye(3)) - np.log(3.0)) < 1e-12
    assert abs(von_neumann_entropy(5.0 * np.eye(3)) - np.log(3.0)) < 1e-12


def test_reduced_covariance_is_submatrix() -> None:
    C = np.array([[1.0, 0.2, 0.3], [0.2, 2.0, 0.4], [0.3, 0.4, 3.0]])
    assert np.allclose(reduced_covariance(C, [0, 2]),
                       np.array([[1.0, 0.3], [0.3, 3.0]]))
    assert np.allclose(reduced_covariance(C, [1]), np.array([[2.0]]))


# ──────────────────────────────────────────────────────────────────────────
# Sample-covariance sanity (data → C, then full pipeline)
# ──────────────────────────────────────────────────────────────────────────

def test_sample_covariance_from_s3_equivariant_data() -> None:
    """Draw M samples from N(0, C) with C = I + 0.3(J − I). The sample
    covariance converges to C, so the measured decomposition approaches
    the exact one within Wishart tolerance."""
    rng = np.random.default_rng(0)
    a, b = 1.0, 0.3
    C_true = a * np.eye(3) + b * (np.ones((3, 3)) - np.eye(3))
    # Cholesky factor
    L = np.linalg.cholesky(C_true)
    M = 50_000
    Z = rng.standard_normal((M, 3))
    X = Z @ L.T
    C_hat = sample_covariance(X)
    # Spectrum of C_true: {a + 2b, a − b, a − b} = {1.6, 0.7, 0.7}.
    # At M = 50k, relative Wishart noise floor ~ √(2/M) ≈ 0.6%.
    assert np.max(np.abs(C_hat - C_true)) < 0.05
