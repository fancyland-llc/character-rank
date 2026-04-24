"""Character-theoretic effective-rank bound for G-equivariant covariances.

═══════════════════════════════════════════════════════════════════
THE CLAIM
═══════════════════════════════════════════════════════════════════

Let a finite group G_α ≤ S_N act on ℝ^N by coordinate permutation,
and let C ∈ ℝ^{N×N} be a symmetric covariance matrix commuting with
the G_α-action (g C g^{-1} = C for every g ∈ G_α). Decompose the
permutation representation of G_α on ℝ^N into irreducibles:

                ℝ^N = ⊕_{π ∈ Ĝ_α} m_π · V_π

with multiplicities m_π and irrep dimensions d_π. Let
𝓐 = {π : m_π > 0} be the set of appearing types.

THEOREM (Character-theoretic effective-rank bound).

    #{distinct eigenvalues of C}  ≤  Σ_{π ∈ 𝓐} m_π,

with each distinct eigenvalue appearing in C's spectrum with
multiplicity d_π · k_π for some k_π ≤ m_π (forced by Schur's
lemma). When every m_π = 1, equality holds generically:
the number of distinct eigenvalues equals |𝓐|.

═══════════════════════════════════════════════════════════════════
THE THREE ISOTROPY TYPES AT N = 3
═══════════════════════════════════════════════════════════════════

    isotropy  ℝ^3 decomposition                              bound
    ───────── ─────────────────────────────────────────────── ─────
    S_3       1·trivial  ⊕  1·standard(d=2)                     2
    S_2       2·trivial  ⊕  1·sign                               3
    {e}       3·trivial                                          3

═══════════════════════════════════════════════════════════════════
SAMPLING REQUIREMENT (Wishart × Nyquist)
═══════════════════════════════════════════════════════════════════

The theorem is an M → ∞ statement. On an M-sample empirical
covariance, the Wishart distribution splits genuinely-equal
eigenvalues with standard deviation Δλ/λ ~ σ² · √(N/M). To resolve
a Schur-forced pair of equal eigenvalues at relative precision δ
against an n_σ-sigma Wishart noise floor requires

                 M  ≥  n_σ²  ·  d_max  /  δ²

where d_max is the largest appearing irrep dimension. For N = 3,
three-sigma (n_σ = 3), d_max = 2 (the standard irrep of S_3), and
1% precision (δ = 0.01): M ≥ 1.8 × 10⁵.
"""
from __future__ import annotations

import dataclasses as dc
from typing import Literal

import numpy as np


# ──────────────────────────────────────────────────────────────────────────
# Isotropy classification
# ──────────────────────────────────────────────────────────────────────────

IsotropyName = Literal["S3", "S2", "trivial"]


@dc.dataclass(frozen=True)
class IrrepDecomp:
    """Decomposition of ℝ^N as a G_α-module by irreducible type.

    irrep_multiplicities : m_π for each appearing irrep π. Σ m_π · d_π = N.
    irrep_dims           : d_π for the same irreps, same order.
    """
    isotropy: IsotropyName
    irrep_multiplicities: tuple[int, ...]
    irrep_dims: tuple[int, ...]

    @property
    def n_distinct_irrep_types(self) -> int:
        """|𝓐|: how many distinct irrep types appear."""
        return len(self.irrep_multiplicities)

    @property
    def ambient_dim(self) -> int:
        return sum(m * d for m, d in zip(
            self.irrep_multiplicities, self.irrep_dims,
        ))

    @property
    def distinct_eigenvalue_bound(self) -> int:
        """Upper bound on the number of distinct eigenvalues of a
        G_α-equivariant covariance: Σ_π m_π. Each distinct eigenvalue
        carries a Schur-forced multiplicity of d_π."""
        return sum(self.irrep_multiplicities)


# Decompositions of ℝ^3 under the three isotropy subgroups of S_3.
#   S_3: trivial (m=1, d=1) ⊕ standard (m=1, d=2). Bound 2.
#   S_2: trivial (m=2, d=1) ⊕ sign (m=1, d=1). Bound 3 (the two trivial
#        copies are not forced equal by Schur).
#   {e}: trivial (m=3, d=1). Bound 3.
PERM_REP_DECOMP: dict[IsotropyName, IrrepDecomp] = {
    "S3":      IrrepDecomp("S3",      (1, 1), (1, 2)),
    "S2":      IrrepDecomp("S2",      (2, 1), (1, 1)),
    "trivial": IrrepDecomp("trivial", (3,),   (1,)),
}


def predict_bound(isotropy: IsotropyName) -> int:
    """A-priori distinct-eigenvalue bound for a given isotropy type.
    Looks up the precomputed decomposition in PERM_REP_DECOMP."""
    return PERM_REP_DECOMP[isotropy].distinct_eigenvalue_bound


# ──────────────────────────────────────────────────────────────────────────
# Wishart-scaled tolerance + distinct-eigenvalue clustering
# ──────────────────────────────────────────────────────────────────────────

def wishart_rel_gap(M: int, n_sigma: float = 3.0) -> float:
    """Heuristic relative-gap tolerance for clustering eigenvalues of an
    M-sample 3×3 covariance under Wishart(M, σ²I_N) noise. The splitting
    of genuinely-equal eigenvalues scales as σ²·√(2/M); we multiply by
    n_sigma for a three-sigma tolerance by default.
    """
    return float(n_sigma * np.sqrt(2.0 / max(1, M)))


def count_distinct_eigenvalues(
    eigvals: np.ndarray,
    *,
    rel_gap: float = 0.06,
    abs_floor: float = 1e-10,
) -> int:
    """Count distinct eigenvalues by relative-gap clustering. Two
    eigenvalues are considered the same cluster if their relative gap
    to the running cluster mean is below `rel_gap`. Values below
    `abs_floor` are treated as numerical zeros and dropped.
    """
    vals = np.sort(np.asarray(eigvals, dtype=np.float64))[::-1]
    vals = vals[vals > abs_floor]
    if vals.size == 0:
        return 0
    clusters: list[list[float]] = [[float(vals[0])]]
    for v in vals[1:]:
        mean_prev = float(np.mean(clusters[-1]))
        if abs(v - mean_prev) / max(mean_prev, abs_floor) > rel_gap:
            clusters.append([float(v)])
        else:
            clusters[-1].append(float(v))
    return len(clusters)


# ──────────────────────────────────────────────────────────────────────────
# Measurement on a raw covariance matrix
# ──────────────────────────────────────────────────────────────────────────

@dc.dataclass(frozen=True)
class RankMeasurement:
    """Empirical spectrum of a G_α-equivariant covariance matrix.

    eigenvalues              : sorted descending, len = N.
    distinct_eigenvalue_count: count from Wishart-scaled clustering.
    effective_rank_09        : |{i : Σ_{j≤i} λ_j² / Σ_j λ_j² ≥ 0.9}|.
    """
    eigenvalues: tuple[float, ...]
    distinct_eigenvalue_count: int
    effective_rank_09: float


def measure_covariance_spectrum(
    C: np.ndarray,
    *,
    M: int | None = None,
    rel_gap: float | None = None,
    abs_floor: float = 1e-10,
) -> RankMeasurement:
    """Return the sorted spectrum of an already-computed covariance C,
    together with a Wishart-scaled distinct-eigenvalue count and the
    effective rank at 90% variance.

    Parameters
    ----------
    C : (N, N) symmetric PSD matrix.
    M : sample size used to estimate C. If provided (and `rel_gap` is not),
        the clustering tolerance defaults to `wishart_rel_gap(M)`.
    rel_gap : override the relative-gap tolerance manually.
    abs_floor : eigenvalues below this magnitude are treated as zero.
    """
    eigvals = np.sort(np.linalg.eigvalsh(C))[::-1]
    if rel_gap is None:
        rel_gap = wishart_rel_gap(M) if M is not None else 0.06
    d = count_distinct_eigenvalues(
        eigvals, rel_gap=rel_gap, abs_floor=abs_floor,
    )
    s2 = eigvals ** 2
    total = float(s2.sum())
    if total <= 0.0:
        eff_09 = float("nan")
    else:
        cum = np.cumsum(s2) / total
        eff_09 = float(int(np.searchsorted(cum, 0.9) + 1))
    return RankMeasurement(
        eigenvalues=tuple(float(v) for v in eigvals),
        distinct_eigenvalue_count=d,
        effective_rank_09=eff_09,
    )


# ──────────────────────────────────────────────────────────────────────────
# Pass/fail check against the Schur bound
# ──────────────────────────────────────────────────────────────────────────

@dc.dataclass(frozen=True)
class BoundCheck:
    isotropy: IsotropyName
    bound: int
    measured_distinct: int
    measured_eff_09: float
    eigenvalues: tuple[float, ...]
    passes: bool


def check_bound(
    C: np.ndarray,
    isotropy: IsotropyName,
    *,
    M: int | None = None,
    rel_gap: float | None = None,
) -> BoundCheck:
    """Measure the spectrum of C and test it against the Schur-forced
    upper bound for the given isotropy. Passes iff
    measured_distinct ≤ PERM_REP_DECOMP[isotropy].distinct_eigenvalue_bound.

    Schur's lemma forces degeneracies; it never forbids them. A
    measurement strictly below the bound is allowed and is a signature
    of either accidental degeneracy or a noise-floor collapse.
    """
    meas = measure_covariance_spectrum(C, M=M, rel_gap=rel_gap)
    bound = PERM_REP_DECOMP[isotropy].distinct_eigenvalue_bound
    return BoundCheck(
        isotropy=isotropy,
        bound=bound,
        measured_distinct=meas.distinct_eigenvalue_count,
        measured_eff_09=meas.effective_rank_09,
        eigenvalues=meas.eigenvalues,
        passes=meas.distinct_eigenvalue_count <= bound,
    )


__all__ = [
    "IsotropyName",
    "IrrepDecomp",
    "PERM_REP_DECOMP",
    "predict_bound",
    "wishart_rel_gap",
    "count_distinct_eigenvalues",
    "RankMeasurement",
    "measure_covariance_spectrum",
    "BoundCheck",
    "check_bound",
]
