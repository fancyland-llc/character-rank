"""Character-entropy decomposition for G-equivariant covariances (N = 3).

═══════════════════════════════════════════════════════════════════
THE THEOREM
═══════════════════════════════════════════════════════════════════

Let G_α ≤ S_N act on ℝ^N by coordinate permutation, and let C be a
G_α-equivariant symmetric PSD matrix. Write the isotypic decomposition

            ℝ^N = ⊕_π V_π^{⊕ m_π},       C = ⊕_π C_π,

and let p_π = tr(C_π) / tr(C) be the Born-rule weight on the π-block,
d_π = dim V_π the irrep dimension, and S_π^int the von-Neumann entropy
of the m_π × m_π "copy matrix" representing C restricted to the m_π
copies of V_π (equivalently, C_π up to the Schur-forced d_π-fold
multiplicity). Define ρ = C / tr(C). Then

    S(ρ) = -Σ_i ρ_i log ρ_i
         = H(p)                      ← Shannon over irrep types
         + Σ_π p_π · log d_π          ← Schur bonus (forced by symmetry)
         + Σ_π p_π · S_π^int          ← dynamics-controlled internal term,

an exact identity on G_α-equivariant C, up to numerical error from the
isotypic projection on an empirical (not perfectly equivariant) matrix.

═══════════════════════════════════════════════════════════════════
HOLOGRAPHIC READING
═══════════════════════════════════════════════════════════════════

    Standard Ryu-Takayanagi          Equivariant-Markov analog
    ─────────────────────────        ──────────────────────────
    boundary region A                coordinate subset A ⊆ {1..N}
    bulk AdS geometry                full trajectory space
    radial coordinate                Markov depth T
    area |γ_A|/4G                    Schur bonus Σ p_π log d_π
    quantum corrections              internal entropy Σ p_π S_π^int

The Schur bonus is depth-independent (it's set by the symmetry, not
the dynamics); internal entropy is the piece the generator can grow
or erase. The central operational claim of the framework is that if
the Markov generator preserves G_α-equivariance and admits a steady
state, then S(ρ) saturates in Markov depth T and each of the three
pieces stabilizes — a discrete area law.

═══════════════════════════════════════════════════════════════════
SCOPE (N = 3)
═══════════════════════════════════════════════════════════════════

The isotypic projectors/bases implemented here are specialized to
N = 3 under the three isotropy subgroups {S_3, S_2, {e}}. Extending
to larger N is mechanical — decompose the permutation representation
once per isotropy type, drop in orthonormal basis columns for each
isotypic subspace, and the rest of the machinery (character_entropy_
decomposition, gaussian_mutual_information) works unchanged.
"""
from __future__ import annotations

import dataclasses as dc
from typing import Sequence

import numpy as np

from .character_rank import IsotropyName


# ──────────────────────────────────────────────────────────────────────────
# Isotypic projectors for N = 3
# ──────────────────────────────────────────────────────────────────────────

def isotypic_projectors(
    isotropy: IsotropyName,
    singleton_slot: int = 0,
) -> dict[str, np.ndarray]:
    """Orthogonal projectors onto the isotypic subspaces of ℝ^3 under the
    given isotropy. Keyed by irrep name.

    S_3:
        P_triv  = (1/3) J_3               rank 1
        P_std   = I_3 - (1/3) J_3         rank 2

    S_2 (with the singleton coordinate at `singleton_slot`, the other two
    forming the exchangeable pair):
        P_triv  = e_s e_s^T + (1/2)(e_p + e_q)(e_p + e_q)^T     rank 2
        P_sign  = (1/2) (e_p - e_q)(e_p - e_q)^T                rank 1

    {e}:
        P_triv  = I_3                     rank 3
    """
    if isotropy == "S3":
        J = np.ones((3, 3)) / 3.0
        return {"triv": J, "std": np.eye(3) - J}

    if isotropy == "S2":
        if singleton_slot not in (0, 1, 2):
            raise ValueError(f"singleton_slot must be 0/1/2, got {singleton_slot}")
        pair = tuple(i for i in range(3) if i != singleton_slot)
        p, q = pair
        P_triv = np.zeros((3, 3))
        P_triv[singleton_slot, singleton_slot] = 1.0
        v_sym = np.zeros(3)
        v_sym[p] = v_sym[q] = 1.0 / np.sqrt(2.0)
        P_triv = P_triv + np.outer(v_sym, v_sym)
        v_anti = np.zeros(3)
        v_anti[p] = 1.0 / np.sqrt(2.0)
        v_anti[q] = -1.0 / np.sqrt(2.0)
        P_sign = np.outer(v_anti, v_anti)
        return {"triv": P_triv, "sign": P_sign}

    if isotropy == "trivial":
        return {"triv": np.eye(3)}

    raise ValueError(f"unknown isotropy: {isotropy}")


def isotypic_basis(
    isotropy: IsotropyName,
    singleton_slot: int = 0,
) -> dict[str, np.ndarray]:
    """Orthonormal basis matrices for each isotypic subspace of ℝ^3.

    Each value is an (N, m_π · d_π) column matrix whose columns span the
    π-isotypic subspace.  Used for computing S_π^int via restriction
    C_π = U_π^T C U_π followed by (conceptual) Schur-trace over the d_π
    copies; in the small cases below every d_π = 1 holds in the mult > 1
    branches, so C_π is directly the m_π × m_π copy matrix.
    """
    if isotropy == "S3":
        U_triv = np.ones((3, 1)) / np.sqrt(3.0)
        U_std = np.array([
            [1.0,  1.0],
            [-1.0, 1.0],
            [0.0, -2.0],
        ])
        U_std[:, 0] /= np.sqrt(2.0)
        U_std[:, 1] /= np.sqrt(6.0)
        return {"triv": U_triv, "std": U_std}

    if isotropy == "S2":
        pair = tuple(i for i in range(3) if i != singleton_slot)
        p, q = pair
        e_sing = np.zeros((3, 1))
        e_sing[singleton_slot, 0] = 1.0
        e_pair_sym = np.zeros((3, 1))
        e_pair_sym[p, 0] = e_pair_sym[q, 0] = 1.0 / np.sqrt(2.0)
        U_triv = np.hstack([e_sing, e_pair_sym])          # (3, 2)
        U_sign = np.zeros((3, 1))
        U_sign[p, 0] = 1.0 / np.sqrt(2.0)
        U_sign[q, 0] = -1.0 / np.sqrt(2.0)
        return {"triv": U_triv, "sign": U_sign}

    if isotropy == "trivial":
        return {"triv": np.eye(3)}

    raise ValueError(f"unknown isotropy: {isotropy}")


# ──────────────────────────────────────────────────────────────────────────
# Gaussian entropy + mutual information on an arbitrary covariance
# ──────────────────────────────────────────────────────────────────────────

def _normalize(C: np.ndarray) -> np.ndarray:
    tr = float(np.trace(C))
    if tr <= 0.0:
        raise ValueError("covariance has non-positive trace")
    return C / tr


def von_neumann_entropy(C: np.ndarray, *, normalize: bool = True) -> float:
    """Density-matrix entropy of a PSD matrix C.

    If normalize=True, treat ρ = C / tr(C) and compute S(ρ) = -Σ λ log λ.
    Zero eigenvalues are dropped (0 log 0 = 0 by convention).
    """
    if normalize:
        C = _normalize(C)
    eigs = np.linalg.eigvalsh(C)
    eigs = eigs[eigs > 1e-12]
    if eigs.size == 0:
        return 0.0
    return float(-np.sum(eigs * np.log(eigs)))


def reduced_covariance(C: np.ndarray, subset: Sequence[int]) -> np.ndarray:
    """Principal submatrix C[subset, subset] — the marginal covariance."""
    idx = np.asarray(subset, dtype=int)
    return C[np.ix_(idx, idx)]


def gaussian_mutual_information(
    C: np.ndarray,
    A: Sequence[int],
    B: Sequence[int],
) -> float:
    """Shannon mutual information of a Gaussian N(0, C) between disjoint
    marginals A and B:

        I(A:B) = (1/2) · log( det(C_A) · det(C_B) / det(C_{A∪B}) ).

    Returns 0 if A or B is empty. A and B must be disjoint.
    """
    A = list(A)
    B = list(B)
    if not A or not B:
        return 0.0
    if set(A) & set(B):
        raise ValueError("A and B must be disjoint")
    AB = sorted(set(A) | set(B))
    _, ld_A = np.linalg.slogdet(reduced_covariance(C, A))
    _, ld_B = np.linalg.slogdet(reduced_covariance(C, B))
    _, ld_AB = np.linalg.slogdet(reduced_covariance(C, AB))
    return float(0.5 * (ld_A + ld_B - ld_AB))


# ──────────────────────────────────────────────────────────────────────────
# Character-decomposed entropy (the main theorem)
# ──────────────────────────────────────────────────────────────────────────

@dc.dataclass(frozen=True)
class CharEntropyDecomp:
    """Per-isotypic-block entropy decomposition of a G_α-equivariant C.

    p_block          : Born weight p_π = tr(C_π) / tr(C), keyed by irrep.
    d_block          : Schur-forced multiplicity d_π of each irrep.
    m_block          : appearance multiplicity m_π in ℝ^N.
    shannon          : H(p).
    schur_bonus      : Σ_π p_π log d_π.
    internal         : dict of S_π^int — zero whenever m_π = 1.
    internal_total   : Σ_π p_π S_π^int.
    total            : shannon + schur_bonus + internal_total.
    direct_entropy   : S(ρ) computed directly from C's eigenvalues.
    reconstruction_err: |total - direct_entropy| — sanity check.
    """
    isotropy: IsotropyName
    blocks: tuple[str, ...]
    p_block: dict[str, float]
    d_block: dict[str, int]
    m_block: dict[str, int]
    shannon: float
    schur_bonus: float
    internal: dict[str, float]
    internal_total: float
    total: float
    direct_entropy: float
    reconstruction_err: float


def character_entropy_decomposition(
    C: np.ndarray,
    isotropy: IsotropyName,
    *,
    singleton_slot: int = 0,
) -> CharEntropyDecomp:
    """Decompose S(C/tr C) into Shannon + Schur-bonus + internal pieces.

    Exact when C is perfectly G-equivariant; on an empirical estimate the
    reconstruction_err field reports the residual between the component
    sum and the direct density-matrix entropy.
    """
    bases = isotypic_basis(isotropy, singleton_slot=singleton_slot)
    tr_C = float(np.trace(C))
    if tr_C <= 0.0:
        raise ValueError("covariance has non-positive trace")

    if isotropy == "S3":
        d = {"triv": 1, "std": 2}
        m = {"triv": 1, "std": 1}
    elif isotropy == "S2":
        d = {"triv": 1, "sign": 1}
        m = {"triv": 2, "sign": 1}
    else:  # trivial
        d = {"triv": 1}
        m = {"triv": 3}

    p_block: dict[str, float] = {}
    internal: dict[str, float] = {}

    for name, U in bases.items():
        C_block = U.T @ C @ U
        tr_block = float(np.trace(C_block))
        p_block[name] = tr_block / tr_C
        if m[name] == 1:
            internal[name] = 0.0
            continue
        # In every case with m > 1 in the N = 3 tables, d = 1, so C_block
        # already is the m × m copy matrix; no Schur-trace needed.
        tr_copy = float(np.trace(C_block))
        if tr_copy <= 0.0:
            internal[name] = 0.0
        else:
            internal[name] = von_neumann_entropy(C_block)

    p_vals = np.array(list(p_block.values()))
    p_vals = p_vals[p_vals > 1e-12]
    shannon = float(-np.sum(p_vals * np.log(p_vals)))
    schur = float(sum(p_block[n] * np.log(d[n]) for n in bases.keys()))
    internal_total = float(sum(
        p_block[n] * internal[n] for n in bases.keys()
    ))
    total = shannon + schur + internal_total
    direct = von_neumann_entropy(C)

    return CharEntropyDecomp(
        isotropy=isotropy,
        blocks=tuple(bases.keys()),
        p_block=p_block,
        d_block=d,
        m_block=m,
        shannon=shannon,
        schur_bonus=schur,
        internal=internal,
        internal_total=internal_total,
        total=total,
        direct_entropy=direct,
        reconstruction_err=abs(total - direct),
    )


# ──────────────────────────────────────────────────────────────────────────
# Convenience: sample covariance from a (M, N) matrix of observations
# ──────────────────────────────────────────────────────────────────────────

def sample_covariance(X: np.ndarray) -> np.ndarray:
    """Unbiased sample covariance of an (M, N) matrix X of observations.

    Centers across the first axis (across the M independent samples),
    then returns X_c^T X_c / (M - 1). No centering across coordinates —
    that would kill the trivial irrep.
    """
    X = np.asarray(X, dtype=np.float64)
    if X.ndim != 2:
        raise ValueError(f"X must be (M, N); got shape {X.shape}")
    M = X.shape[0]
    if M < 2:
        raise ValueError(f"need M ≥ 2 samples; got {M}")
    Xc = X - X.mean(axis=0, keepdims=True)
    return (Xc.T @ Xc) / (M - 1)


__all__ = [
    "isotypic_projectors",
    "isotypic_basis",
    "von_neumann_entropy",
    "reduced_covariance",
    "gaussian_mutual_information",
    "CharEntropyDecomp",
    "character_entropy_decomposition",
    "sample_covariance",
]
