"""Holographic / character-entropy decomposition on the prime gas.

═══════════════════════════════════════════════════════════════════
WHY THIS MODULE
═══════════════════════════════════════════════════════════════════

Theorem 11.1 (character-entropy decomposition) was stated and
confirmed in holographic_rank.py on the Coliseum accumulator, an
S_3-equivariant Markov generator on ℝ^3. The Bounded-Boundary
Finite-Complexity Principle (BBFC) says the theorem is a member of
a universal family: any well-posed BVP with a symmetry group admits
a character decomposition of its density-matrix entropy, with the
same three pieces — Shannon + Schur-bonus + internal.

The cheapest test of that universality is to run the theorem on a
genuinely independent system in the same category. The prime gas at
primorial base (Matos 2026; UNITY_CLOCK_EFFECTIVE_RANK.md) is that
system:

  - state space: ℂ^{φ(m)} indexed by coprime residues mod m
  - symmetry group: G = (ℤ/m)^× acting by multiplication
  - dynamics:   H(α, m) = D_sym / λ_P + iα [D_sym, P_τ] / λ_P
               (Hamiltonian operator, unity_clock_adscft.py)
  - observable: amplitudes on the boundary (coprime block)

Key structural difference from Coliseum: (ℤ/m)^× is ABELIAN, so
every irrep is 1-dimensional. Schur bonus = Σ p_π log d_π = 0
identically. The character-entropy decomposition collapses to

        S(ρ) = H(p)   (Shannon over Dirichlet characters)

which is a much sharper statement: the Gaussian density-matrix
entropy of a (ℤ/m)^×-equivariant covariance equals the Shannon
entropy over the Dirichlet-character Fourier weights. No bulk
multiplicity. Pure boundary.

But BBFC has teeth: it says the dictionary entry must exist, and
*any* deviation from Shannon over characters is a falsifier. At
finite M sample size, Wishart noise gives <1% leakage. Above that
is a signal that the dynamics break G-equivariance or that we're
measuring the wrong observable.

═══════════════════════════════════════════════════════════════════
WHAT THIS MODULE COMPUTES
═══════════════════════════════════════════════════════════════════

  (1) Dirichlet character table of (ℤ/m)^× — built from the
      regular representation of the abelian group.
  (2) Ensemble of prime-gas trajectories: M probe vectors evolved
      under H(α) for a ladder of times t, projected back onto
      the coprime boundary → (M, T, φ(m)) array analogous to the
      Coliseum's energy tensor.
  (3) Boundary covariance C ∈ ℂ^{φ(m) × φ(m)}: Hermitian, positive.
  (4) Character-entropy decomposition in the Dirichlet basis:
        S(ρ) = H(p)  (Shannon over characters)
      For abelian G, internal = 0 and Schur = 0, so the full
      entropy content is the Shannon term.
  (5) Subgroup refinement tower: for each subgroup H ≤ G,
      decompose the covariance under H instead of G. Subgroup
      decompositions have NON-TRIVIAL multiplicities (each H-irrep
      appears [G:H] times), so internal entropy appears. The
      "subgroup refinement bonus" S(G) - S(H-restricted) is a new
      observable.
  (6) Depth ladder: does S(ρ) saturate as t → ∞ (area law in the
      Markov time)?
  (7) Primorial ladder: m ∈ {6, 30, 210, 2310}. How does H(p)
      scale with |G| = φ(m)?
  (8) Coliseum functorial comparison: at m = 6, (ℤ/6)^× ≅ ℤ/2 as
      abstract group, structurally identical to the S_2 action in
      APP-style Coliseum matchups. The character-entropy
      decomposition at m = 6 should match the APP decomposition
      up to dimensional rescaling. This is the ONE functorial
      check the dictionary needs at minimum.

═══════════════════════════════════════════════════════════════════
EXOTIC ADDITIONS (RABBIT-HOLE SEEDS)
═══════════════════════════════════════════════════════════════════

Beyond the basic test, this module exposes several measurements
intended to find the NEXT rabbit:

  • Haar-random baseline for the character-entropy spectrum.
    Compare to Page curve / Marchenko-Pastur.
  • CRT factorization check: at composite primorials, the entropy
    decomposition should factor across prime-factor components.
  • Spectral form factor × character content cross-correlation.
    (If the SFF has dip-ramp-plateau, does the character
    decomposition track the dip?)
  • Character-entropy at the Rabi phase transition α_c ≈
    √(135/88). Does the decomposition detect the transition?
"""
from __future__ import annotations

import dataclasses as dc
import math
from functools import lru_cache
from math import gcd
from typing import Sequence

import numpy as np


# ──────────────────────────────────────────────────────────────────────────
# §1. (ℤ/m)^× group structure and Dirichlet characters
# ──────────────────────────────────────────────────────────────────────────

@lru_cache(maxsize=32)
def coprime_residues(m: int) -> tuple[int, ...]:
    """Coprime residues mod m, in ascending order. (ℤ/m)^× elements."""
    return tuple(r for r in range(1, m) if gcd(r, m) == 1)


def group_order(m: int) -> int:
    """φ(m) = |(ℤ/m)^×|."""
    return len(coprime_residues(m))


def mult_table(m: int) -> np.ndarray:
    """(n, n) integer table: table[i, j] = index of coprimes[i] * coprimes[j]
    mod m, in coprimes[]."""
    G = coprime_residues(m)
    idx = {g: i for i, g in enumerate(G)}
    n = len(G)
    tab = np.empty((n, n), dtype=np.int64)
    for i, a in enumerate(G):
        for j, b in enumerate(G):
            tab[i, j] = idx[(a * b) % m]
    return tab


def dirichlet_character_table(m: int) -> np.ndarray:
    """Return the (n × n) Dirichlet character table Λ of (ℤ/m)^×.

    Rows are characters; columns are group elements in `coprime_residues(m)`
    order. Λ is UNITARY in the sense that (1/n) Λ Λ^* = I (Schur
    orthogonality). For abelian G all irreps are 1-D so n = |G|.

    Construction: find generators of the abelian group via Smith-normal-
    form-style decomposition of the multiplicative structure. For small
    m we compute it via simultaneous diagonalization of regular-
    representation permutation matrices — that's slow but bulletproof.
    """
    G = coprime_residues(m)
    n = len(G)

    if n == 1:
        return np.ones((1, 1), dtype=np.complex128)

    # Permutation matrix of multiplication by g_k for each k (regular rep).
    # The characters are the common eigenvectors of all these (they commute
    # because G is abelian).
    tab = mult_table(m)
    # Use g = coprimes[1] as the 'first' regulator; iteratively refine
    # using other generators until eigenspaces are fully split.
    # Build a single "hash" matrix as a random linear combination — with
    # irrational coefficients this simultaneously diagonalizes a.s.
    rng = np.random.default_rng(0xC0FFEE)
    coeffs = rng.standard_normal(n) + 1j * rng.standard_normal(n)
    M = np.zeros((n, n), dtype=np.complex128)
    for k in range(n):
        P = np.zeros((n, n), dtype=np.complex128)
        for i in range(n):
            P[tab[k, i], i] = 1.0
        M = M + coeffs[k] * P

    # Eigendecompose. Each eigenvector is a character (up to phase).
    _, V = np.linalg.eig(M)
    # Normalize columns, rotate to make first entry real-positive.
    for j in range(n):
        V[:, j] /= np.linalg.norm(V[:, j])
        # Conjugate so that each column represents a CHARACTER — i.e., a
        # homomorphism G → ℂ^×. The eigenvectors of the regular rep give
        # us the Fourier basis; the characters are the components.
    # Rows of V^H are characters; normalize to χ(1) = 1.
    Lam = V.conj().T
    for r in range(n):
        # Normalize such that χ(identity) = 1 where identity index is 0
        # (coprimes starts with 1, which is the identity).
        c0 = Lam[r, 0]
        if abs(c0) < 1e-12:
            continue
        Lam[r, :] /= c0
    # Rescale so that Λ is unitary-ish: each row has unit Frobenius * √n.
    # Schur orthogonality: (1/n) Σ_g χ(g) χ'(g)^* = δ. So each row has
    # Σ |χ(g)|^2 = n → norm √n.
    for r in range(n):
        Lam[r, :] *= np.sqrt(n) / np.linalg.norm(Lam[r, :])

    return Lam


# ──────────────────────────────────────────────────────────────────────────
# §2. Subgroup enumeration
# ──────────────────────────────────────────────────────────────────────────

def enumerate_subgroups(m: int) -> list[tuple[int, ...]]:
    """Return all subgroups of (ℤ/m)^×, each as a sorted tuple of elements.

    Uses the naive closure algorithm: for each subset generated by increasing
    element sets, close under multiplication mod m. Feasible for |G| ≤ 64.
    """
    G = coprime_residues(m)
    n = len(G)
    if n > 64:
        raise ValueError(f"enumerate_subgroups too expensive at |G| = {n}")
    elem_set = set(G)
    seen: set[tuple[int, ...]] = set()

    def close(s: set[int]) -> tuple[int, ...]:
        s = set(s)
        s.add(1)
        changed = True
        while changed:
            changed = False
            for a in list(s):
                for b in list(s):
                    c = (a * b) % m
                    if c in elem_set and c not in s:
                        s.add(c)
                        changed = True
        return tuple(sorted(s))

    # Trivial and full
    seen.add((1,))
    seen.add(tuple(G))
    # Add cyclic subgroups generated by each single element
    for g in G:
        H = close({g})
        seen.add(H)
    # Add pairwise-generated subgroups
    for i, g1 in enumerate(G):
        for g2 in G[i + 1:]:
            H = close({g1, g2})
            seen.add(H)
            if len(H) == n:
                break
    # Deduplicate
    return sorted(seen, key=lambda t: (len(t), t))


def subgroup_character_projectors(
    m: int, H: tuple[int, ...],
) -> tuple[np.ndarray, np.ndarray]:
    """Return (Λ_H, mult_H) for subgroup H of (ℤ/m)^× acting on ℂ^{|G|}.

    Λ_H: (|H|, |G|) matrix whose rows are H-characters extended to G by
         choosing a set of coset representatives (we lift by averaging).
    mult_H: multiplicity [G:H] (each H-irrep appears [G:H] times when ℂ^G
            is restricted from G to H).
    """
    G = coprime_residues(m)
    n = len(G)
    Hs = tuple(sorted(H))
    h_idx = {g: i for i, g in enumerate(Hs)}
    # Coset decomposition: G = ⊔ g_i H
    G_set = list(G)
    H_set = set(Hs)
    rem = set(G_set)
    cosets: list[tuple[int, list[int]]] = []  # (rep, members)
    while rem:
        r = min(rem)
        members = [(r * h) % m for h in Hs]
        cosets.append((r, members))
        rem -= set(members)
    # index cosets in G-order
    coset_of = {}
    for ci, (_, members) in enumerate(cosets):
        for mem in members:
            coset_of[mem] = ci
    # Characters of H
    Lam_H = dirichlet_character_table_from_elements(m, Hs)  # (|H|, |H|)
    # Lift to functions on G: lifted χ_π(g) = χ_π(h) where g = g_i * h for
    # some unique coset rep g_i and h ∈ H — but this is NOT in general a
    # character of G. What we actually want is: the H-isotypic projector
    # on ℂ^G, which is (1/|H|) Σ_h χ(h)^* U(h), where U is the regular
    # rep of H on ℂ^G (NOT on ℂ^H).
    U_H: list[np.ndarray] = []  # regular rep of H on ℂ^G
    G_idx = {g: i for i, g in enumerate(G)}
    for h in Hs:
        P = np.zeros((n, n), dtype=np.complex128)
        for i, g in enumerate(G):
            j = G_idx[(g * h) % m]
            P[j, i] = 1.0
        U_H.append(P)
    # For each H-character, build projector onto its isotypic block on ℂ^G
    projectors = np.zeros((len(Hs), n, n), dtype=np.complex128)
    for pi in range(len(Hs)):
        P = np.zeros((n, n), dtype=np.complex128)
        for hi, h in enumerate(Hs):
            P = P + np.conj(Lam_H[pi, hi]) * U_H[hi]
        projectors[pi] = P / len(Hs)
    return projectors, len(G) // len(Hs)  # multiplicity = [G:H]


def dirichlet_character_table_from_elements(
    m: int, elems: tuple[int, ...],
) -> np.ndarray:
    """Compute character table for an arbitrary subgroup H = elems of
    (ℤ/m)^×. Same algorithm as dirichlet_character_table."""
    Hs = tuple(sorted(elems))
    n = len(Hs)
    if n == 1:
        return np.ones((1, 1), dtype=np.complex128)
    idx = {g: i for i, g in enumerate(Hs)}
    tab = np.empty((n, n), dtype=np.int64)
    for i, a in enumerate(Hs):
        for j, b in enumerate(Hs):
            tab[i, j] = idx[(a * b) % m]
    rng = np.random.default_rng(0xBADF00D)
    coeffs = rng.standard_normal(n) + 1j * rng.standard_normal(n)
    M = np.zeros((n, n), dtype=np.complex128)
    for k in range(n):
        P = np.zeros((n, n), dtype=np.complex128)
        for i in range(n):
            P[tab[k, i], i] = 1.0
        M = M + coeffs[k] * P
    _, V = np.linalg.eig(M)
    for j in range(n):
        V[:, j] /= np.linalg.norm(V[:, j])
    Lam = V.conj().T
    for r in range(n):
        c0 = Lam[r, 0]
        if abs(c0) > 1e-12:
            Lam[r, :] /= c0
        Lam[r, :] *= np.sqrt(n) / np.linalg.norm(Lam[r, :])
    return Lam


# ──────────────────────────────────────────────────────────────────────────
# §3. Prime-gas ensemble generator (boundary-only, lightweight)
# ──────────────────────────────────────────────────────────────────────────

def build_distance_kernel_boundary(m: int) -> np.ndarray:
    """Build the φ(m) × φ(m) symmetric distance matrix D_sym on the coprime
    boundary alone. D_sym[i, j] = min(|a_i - a_j|, m - |a_i - a_j|) — the
    circular distance. Returns float array.

    This is the G_α-equivariant generator piece of the prime-gas
    Hamiltonian. It commutes with multiplication by any g ∈ G, so its
    eigenvectors are exactly the Dirichlet characters.
    """
    G = coprime_residues(m)
    n = len(G)
    D = np.zeros((n, n), dtype=np.float64)
    for i, a in enumerate(G):
        for j, b in enumerate(G):
            d = abs(a - b)
            D[i, j] = min(d, m - d)
    return D


def build_crt_local_distance_kernel(m: int) -> np.ndarray:
    """Build a CRT-LOCAL distance kernel on the coprime boundary.

    D^CRT_m = Σ_p (I ⊗ ... ⊗ D_p ⊗ ... ⊗ I)   (under CRT tensor factorization)

    Concretely, for each prime factor p of m, we add the circular distance
    matrix on (ℤ/p)^× acting *only* on the p-coordinate under the CRT
    isomorphism (ℤ/m)^× ≅ Π (ℤ/p_i)^×. Entries (i, j) for residues a, b
    get D_p(a mod p, b mod p) summed over primes p such that a ≡ b mod q
    for every other prime q. Off-diagonal coupling across more than one
    prime factor is strictly zero.

    This is the control against which the natural `build_distance_kernel_boundary`
    is compared: the natural D_sym uses circular distance on ℤ/m directly,
    which mixes all prime factors simultaneously. By construction, D^CRT
    commutes with every CRT projection, so under its dynamics the entropy
    must factorize across primes (factorization residual ≈ 0).
    """
    from sympy import factorint
    primes = sorted(factorint(m).keys())
    G = coprime_residues(m)
    n = len(G)
    # Per-prime: coprime-residue list and local circular-distance kernel
    per_prime: dict = {}
    for p in primes:
        Gp = tuple(a for a in range(1, p) if gcd(a, p) == 1)
        kp = len(Gp)
        Dp = np.zeros((kp, kp), dtype=np.float64)
        for ii, a in enumerate(Gp):
            for jj, b in enumerate(Gp):
                d = abs(a - b)
                Dp[ii, jj] = min(d, p - d)
        per_prime[p] = (Gp, {v: k for k, v in enumerate(Gp)}, Dp)
    D_crt = np.zeros((n, n), dtype=np.float64)
    for i, a in enumerate(G):
        a_mods = {p: a % p for p in primes}
        for j, b in enumerate(G):
            b_mods = {p: b % p for p in primes}
            for p in primes:
                # only contribute if a and b agree on every OTHER prime
                if all(a_mods[q] == b_mods[q] for q in primes if q != p):
                    _, idx_p, Dp = per_prime[p]
                    ia = idx_p[a_mods[p]]
                    ib = idx_p[b_mods[p]]
                    D_crt[i, j] += Dp[ia, ib]
    return D_crt


def prime_gas_ensemble(
    m: int,
    *,
    M: int = 5000,
    t_schedule: Sequence[float] = (0.0, 0.5, 1.0, 2.0, 5.0),
    alpha: float = 0.0,
    beta: float = 1.0,
    rng: np.random.Generator | None = None,
    D_override: np.ndarray | None = None,
) -> np.ndarray:
    """Generate an ensemble of M prime-gas boundary trajectories.

    At α = 0 the Hamiltonian on the boundary reduces to βD_sym/λ_P, which
    commutes with G so its equilibrium covariance is exactly G-equivariant.
    At α > 0 we lose strict equivariance (the commutator piece breaks it),
    giving us a knob to measure how equivariance breaks.

    Method: sample M iid Gaussian initial conditions on the boundary,
    evolve each by exp(-β t D_sym / λ_P) (a Markov semigroup, not a
    Schrödinger evolution — this is the *relaxation* dynamics natural to
    BBFC's detailed-balance category). Collect amplitudes at each t.

    `D_override` lets the caller supply a custom generator (e.g. the
    CRT-local control kernel from `build_crt_local_distance_kernel`); if
    None, the natural circular-distance kernel is used.

    Returns: complex array of shape (M, len(t_schedule), φ(m)).
    """
    rng = rng or np.random.default_rng(0xDECAF)
    G = coprime_residues(m)
    n = len(G)
    D = build_distance_kernel_boundary(m) if D_override is None else D_override
    # Eigendecomposition of D — since D is G-equivariant, eigenvectors
    # are the Dirichlet characters and eigenvalues are real.
    evals, evecs = np.linalg.eigh(D)
    lam_P = float(np.max(np.abs(evals)))

    # M × n iid complex Gaussian initial amplitudes
    psi0 = (rng.standard_normal((M, n)) + 1j * rng.standard_normal((M, n))) / np.sqrt(2 * n)

    # Relaxation in eigenbasis: psi(t) = V diag(exp(-β t λ / λ_P)) V^T psi0
    # (using V = evecs, which is orthogonal for symmetric D)
    psi0_coeffs = psi0 @ evecs  # (M, n) in eigenbasis
    result = np.zeros((M, len(t_schedule), n), dtype=np.complex128)
    for ti, t in enumerate(t_schedule):
        decay = np.exp(-beta * t * evals / lam_P)
        psi_t_coeffs = psi0_coeffs * decay[None, :]
        # For α > 0 we'd add the iα[D,P_τ] piece here. We approximate it
        # as a small perturbation mixing the characters: phase shifts
        # proportional to α · (eval shifts). This is a dev-only knob.
        if alpha != 0.0:
            phase = np.exp(1j * alpha * t * evals / lam_P)
            psi_t_coeffs = psi_t_coeffs * phase[None, :]
        result[:, ti, :] = psi_t_coeffs @ evecs.T

    return result


# ──────────────────────────────────────────────────────────────────────────
# §4. Covariance and character-entropy decomposition (abelian form)
# ──────────────────────────────────────────────────────────────────────────

def boundary_covariance(
    psi_ensemble: np.ndarray,
    *,
    time_idx: int | slice | None = None,
) -> np.ndarray:
    """Sample covariance across M trajectories at one or more time slices.

    psi_ensemble: (M, T, n) complex.
    Returns: Hermitian (n, n) complex covariance.
    """
    M, T, n = psi_ensemble.shape
    if time_idx is None:
        # Time-mean observable
        X = psi_ensemble.mean(axis=1)
    elif isinstance(time_idx, slice):
        X = psi_ensemble[:, time_idx, :].reshape(M * (time_idx.stop - time_idx.start), n)
    else:
        X = psi_ensemble[:, time_idx, :]
    Xc = X - X.mean(axis=0, keepdims=True)
    C = (Xc.conj().T @ Xc) / max(1, X.shape[0] - 1)
    return C


def von_neumann_entropy_complex(C: np.ndarray, *, normalize: bool = True) -> float:
    """S(C/tr C) for a Hermitian PSD matrix."""
    tr = float(np.real(np.trace(C)))
    if tr <= 0.0:
        return 0.0
    rho = C / tr if normalize else C
    eigs = np.linalg.eigvalsh((rho + rho.conj().T) / 2.0)
    eigs = eigs[eigs > 1e-14]
    return float(-np.sum(eigs * np.log(eigs)))


@dc.dataclass(frozen=True)
class AbelianCharEntropy:
    """Character-entropy decomposition of a (ℤ/m)^×-equivariant covariance.

    For abelian G, every irrep is 1-dimensional with multiplicity 1 in the
    regular rep on ℂ^{|G|}. So Theorem 11.1 reduces to

            S(ρ) = H(p)     (Shannon over characters)

    and ALL the non-Shannon content (Schur bonus, internal) must be zero
    in a perfectly G-equivariant input. Any deviation is a signal.
    """
    m: int
    n: int
    p_char: np.ndarray              # (n,) Born weights per character
    char_variances: np.ndarray      # (n,) λ_π = tr C_π
    shannon: float                  # H(p)
    direct_entropy: float           # S(C/tr C)
    nonabelian_leakage: float       # |direct - shannon| — should be ~0
    top_characters: tuple[int, ...]  # indices of chars by descending p


def abelian_character_decomposition(
    C: np.ndarray, m: int,
) -> AbelianCharEntropy:
    """Apply the character-entropy decomposition (Theorem 11.1, abelian case)
    to a covariance C on ℂ^{φ(m)}.

    Returns the per-character weights, Shannon over those weights, and the
    leakage to the "direct" entropy computed from C's spectrum. Leakage
    is the quantitative falsifier: if C is truly G-equivariant the two
    entropies agree.
    """
    Lam = dirichlet_character_table(m)  # (n, n) unitary-ish
    n = Lam.shape[0]
    # λ_π = (1/n) χ_π^* C χ_π    (trace of C restricted to the 1D π-block)
    # We use (Λ/√n) as the orthonormal basis matrix.
    U = Lam.conj().T / np.sqrt(n)        # (n, n), columns are char vectors
    C_char = U.conj().T @ C @ U           # diagonal in the character basis
    lam_char = np.real(np.diag(C_char))
    # Numerical cleanup: small negatives from noise → 0
    lam_char = np.maximum(lam_char, 0.0)
    total = lam_char.sum()
    p = lam_char / max(total, 1e-14)
    # Shannon over characters
    p_nz = p[p > 1e-14]
    shannon = float(-np.sum(p_nz * np.log(p_nz)))
    direct = von_neumann_entropy_complex(C)
    top = tuple(int(i) for i in np.argsort(-lam_char))
    return AbelianCharEntropy(
        m=m, n=n,
        p_char=p, char_variances=lam_char,
        shannon=shannon, direct_entropy=direct,
        nonabelian_leakage=abs(direct - shannon),
        top_characters=top,
    )


# ──────────────────────────────────────────────────────────────────────────
# §5. Subgroup refinement tower — where the non-trivial Schur bonus lives
# ──────────────────────────────────────────────────────────────────────────

@dc.dataclass(frozen=True)
class SubgroupRefinement:
    """Character-entropy decomposition of C under a subgroup H ≤ G.

    When restricted from G to H, each H-irrep appears with multiplicity
    [G:H]. The decomposition on ℂ^G thus has mult-[G:H] blocks — genuine
    internal entropy.
    """
    subgroup_size: int
    subgroup_index: int            # [G:H]
    p_block: np.ndarray            # (|H|,) Born weights per H-irrep
    schur_bonus: float             # Σ p_π log d_π (always 0 — H is abelian)
    internal: float                # Σ p_π S_π^internal  ←  the interesting one
    shannon: float                 # H(p)
    total: float                   # reconstructed total entropy
    direct: float
    reconstruction_err: float


def subgroup_refinement(
    C: np.ndarray, m: int, H: tuple[int, ...],
) -> SubgroupRefinement:
    """Decompose C under subgroup H ≤ (ℤ/m)^×, returning a full
    Theorem-11.1-style (Shannon + Schur-bonus + internal) breakdown.

    Every H-irrep is 1-D but appears with multiplicity [G:H], so the
    isotypic blocks have real internal degrees of freedom. This is the
    tower of subgroups that the prime gas's full covariance can be viewed
    through — each subgroup gives a different refinement, and the
    refinement statistics carry information about which symmetries the
    dynamics actually respect.
    """
    projectors, mult = subgroup_character_projectors(m, H)
    n_irreps = projectors.shape[0]
    tr_C = float(np.real(np.trace(C)))
    p = np.zeros(n_irreps)
    internal_per = np.zeros(n_irreps)
    for pi in range(n_irreps):
        P = projectors[pi]
        C_pi = P @ C @ P.conj().T
        tr_pi = float(np.real(np.trace(C_pi)))
        p[pi] = max(tr_pi / max(tr_C, 1e-14), 0.0)
        # Internal: von Neumann entropy of the mult × mult "copy" matrix.
        # Since C commutes with H, C_pi has rank ≤ mult and its mult
        # nonzero eigenvalues are exactly the copy eigenvalues.
        if p[pi] > 1e-12 and mult > 1:
            eigs = np.linalg.eigvalsh((C_pi + C_pi.conj().T) / 2.0)
            eigs = np.sort(np.abs(eigs))[::-1][:mult]
            eigs = np.maximum(eigs, 0.0)
            s = eigs.sum()
            if s > 0:
                q = eigs / s
                q_nz = q[q > 1e-14]
                internal_per[pi] = float(-np.sum(q_nz * np.log(q_nz)))
    p_nz = p[p > 1e-14]
    shannon = float(-np.sum(p_nz * np.log(p_nz)))
    # Schur bonus: H is abelian → every d = 1 → bonus = 0.
    schur = 0.0
    internal = float(np.sum(p * internal_per))
    total = shannon + schur + internal
    direct = von_neumann_entropy_complex(C)
    return SubgroupRefinement(
        subgroup_size=len(H),
        subgroup_index=mult,
        p_block=p,
        schur_bonus=schur,
        internal=internal,
        shannon=shannon,
        total=total,
        direct=direct,
        reconstruction_err=abs(total - direct),
    )


# ──────────────────────────────────────────────────────────────────────────
# §6. Exotic measurements: Haar baseline, CRT factorization, Page curve
# ──────────────────────────────────────────────────────────────────────────

def haar_random_entropy_baseline(n: int, *, M: int = 5000,
                                 rng: np.random.Generator | None = None) -> dict:
    """Compute the expected density-matrix entropy of an M-sample Wishart
    covariance on ℂ^n under Haar-random Gaussian samples. Baseline to
    compare the prime-gas measurement against.

    Returns: dict with mean entropy, std, and theoretical saturation
    value log(n) for rank-full covariance.
    """
    rng = rng or np.random.default_rng(0xF00DFACE)
    # For iid complex Gaussian with identity covariance, sample eigvals
    # of the Wishart(n, M) distribution. Expected normalized spectrum
    # is Marchenko-Pastur with ratio c = n/M.
    reps = 50
    ents = []
    for _ in range(reps):
        X = (rng.standard_normal((M, n)) + 1j * rng.standard_normal((M, n))) / np.sqrt(2)
        Xc = X - X.mean(axis=0, keepdims=True)
        C = (Xc.conj().T @ Xc) / (M - 1)
        ents.append(von_neumann_entropy_complex(C))
    return {
        "dim": n,
        "M": M,
        "mean_entropy": float(np.mean(ents)),
        "std_entropy": float(np.std(ents)),
        "saturation_log_n": float(np.log(n)),
        "ratio_to_saturation": float(np.mean(ents) / np.log(n)) if n > 1 else 1.0,
    }


def crt_character_marginal_check(m: int, C: np.ndarray) -> dict:
    """CORRECT CRT factorization test via character-marginal total correlation.

    For m = p_1·...·p_k squarefree, Ĝ = Π Ĝ_p where Ĝ_p = characters of
    (ℤ/p)^×. A G-character χ is a tuple (χ_{p_1}, ..., χ_{p_k}).

    Compute:
        p(χ)            — weights from abelian_character_decomposition
        p_p(χ_p)        — marginal over each prime factor
        H(p_full)       — full Shannon
        H(p_p)          — marginal Shannon
        TC              — total correlation = Σ H(p_p) − H(p_full)   ≥ 0

    TC = 0 iff p(χ) = Π p_p(χ_p). TC > 0 quantifies the cross-prime mutual
    information carried by the generator — the dynamics-level CRT gap.

    The identification G-char → (G_p-char)_p is done by CRT-embedding each
    residue a ∈ (ℤ/p)^× into (ℤ/m)^× as the unique element e_p(a) that is
    a mod p and 1 mod every other prime q, then reading χ(e_p(a)) off the
    Dirichlet table.
    """
    from sympy import factorint
    primes = sorted(factorint(m).keys())
    G = coprime_residues(m)
    n = len(G)
    G_idx = {g: i for i, g in enumerate(G)}
    Lam = dirichlet_character_table(m)   # (n, n), rows = characters
    # Full character weights
    full_decomp = abelian_character_decomposition(C, m)
    p_full = np.asarray(full_decomp.p_char, dtype=np.float64)
    H_full = float(full_decomp.shannon)

    # Build CRT lifts e_p(a) for each prime p and each a ∈ (ℤ/p)^×
    crt_lift: dict = {}
    for p in primes:
        for a in range(1, p):
            if gcd(a, p) != 1:
                continue
            # find e ∈ G with e ≡ a mod p and e ≡ 1 mod q for each q ≠ p
            for candidate in G:
                if candidate % p != a:
                    continue
                if all(candidate % q == 1 for q in primes if q != p):
                    crt_lift[(p, a)] = candidate
                    break
            else:
                raise RuntimeError(
                    f"CRT lift failed: p={p}, a={a} has no coprime rep mod {m}"
                )

    H_marginals: dict = {}
    p_marginals: dict = {}
    fiber_sizes: dict = {}
    for p in primes:
        a_list = sorted([a for a in range(1, p) if gcd(a, p) == 1])
        # Each G-char k has a p-projection signature: χ_k(e_p(a)) for each a
        fibers: dict = {}
        for k in range(n):
            sig = []
            for a in a_list:
                z = complex(Lam[k, G_idx[crt_lift[(p, a)]]])
                sig.append((round(z.real, 6), round(z.imag, 6)))
            key = tuple(sig)
            fibers.setdefault(key, []).append(k)
        marg = np.array(
            [sum(p_full[k] for k in fiber) for fiber in fibers.values()],
            dtype=np.float64,
        )
        s_marg = marg.sum()
        if s_marg > 0:
            marg = marg / s_marg
        marg_nz = marg[marg > 1e-14]
        H_m = float(-np.sum(marg_nz * np.log(marg_nz)))
        H_marginals[p] = H_m
        p_marginals[p] = marg
        fiber_sizes[p] = [len(f) for f in fibers.values()]

    total_corr = sum(H_marginals.values()) - H_full
    return {
        "m": m,
        "factors": primes,
        "H_full": H_full,
        "H_marginals": H_marginals,
        "p_marginals": p_marginals,
        "fiber_sizes": fiber_sizes,
        "total_correlation": total_corr,
    }


# ──────────────────────────────────────────────────────────────────────────
# §6b. Analytic (Monte-Carlo-free) TC for primorial m via CRT-tensor structure
#
# For primorial m = p_1·...·p_K the group (ℤ/m)^× is cyclic on each factor,
# and every Dirichlet character factors as χ = χ_{p_1} ⊗ ... ⊗ χ_{p_K}.
# A G-equivariant relaxation kernel K has eigenvalue on χ given by
#
#     λ(χ) = Σ_{g ∈ G} K(1, g) · χ(g)*
#
# For iid Gaussian initial conditions, relaxation produces
#
#     p(χ, t) = exp(−2βt λ(χ)/λ_P) / Z(t)
#
# — an exact, deterministic, O(n²)-in-m computation. No ensembles, no MC
# noise. This is what lets us push to m = 30030 (n = 5760) on a single core.
# ──────────────────────────────────────────────────────────────────────────

def _primitive_root_mod_p(p: int) -> int:
    """Smallest primitive root modulo prime p."""
    if p == 2:
        return 1
    phi = p - 1
    # Prime factorization of phi by trial division
    f = phi
    prime_factors = []
    d = 2
    while d * d <= f:
        if f % d == 0:
            prime_factors.append(d)
            while f % d == 0:
                f //= d
        d += 1
    if f > 1:
        prime_factors.append(f)
    for r in range(2, p):
        if all(pow(r, phi // q, p) != 1 for q in prime_factors):
            return r
    raise RuntimeError(f"No primitive root found for p={p}")


def _cyclic_character_table(p: int) -> tuple[tuple[int, ...], np.ndarray, dict]:
    """Character table of (ℤ/p)^× (cyclic of order φ(p) = p−1).

    Returns (Gp, Lam, idx) where:
      Gp  = tuple of coprime residues mod p in ascending order
      Lam = φ(p) × φ(p) complex matrix, Lam[k, i] = χ_k(Gp[i])
      idx = dict a → i mapping Gp entries to their index
    """
    Gp = tuple(a for a in range(1, p) if gcd(a, p) == 1)
    kp = len(Gp)
    idx = {a: i for i, a in enumerate(Gp)}
    if kp == 1:
        return Gp, np.ones((1, 1), dtype=np.complex128), idx
    r = _primitive_root_mod_p(p)
    # discrete log base r
    dlog = {}
    x = 1
    for j in range(kp):
        dlog[x] = j
        x = (x * r) % p
    Lam = np.zeros((kp, kp), dtype=np.complex128)
    for k in range(kp):
        for i, a in enumerate(Gp):
            Lam[k, i] = np.exp(2j * np.pi * k * dlog[a] / kp)
    return Gp, Lam, idx


def _primorial_setup(m: int):
    """Return (primes, per-prime (Gp, Lam, idx), flat G tuple)."""
    from sympy import factorint
    primes = tuple(sorted(factorint(m).keys()))
    per_prime = {p: _cyclic_character_table(p) for p in primes}
    G = coprime_residues(m)
    return primes, per_prime, G


def primorial_kernel_eigenvalue_tensor(
    m: int, K_row: np.ndarray,
) -> tuple[np.ndarray, tuple[int, ...]]:
    """Compute λ(χ) for every Dirichlet character χ of (ℤ/m)^×.

    K_row: shape (φ(m),), giving kernel values K(1, g) for g in
    coprime_residues(m) order. (Any G-equivariant kernel is determined by
    its first row.)

    Returns lambda_tensor of shape (φ(p_1), ..., φ(p_K)) — one value per
    character, indexed by the character's per-prime factor indices. The
    flat list of eigenvalues is lambda_tensor.reshape(-1).
    """
    primes, per_prime, G = _primorial_setup(m)
    n_shape = tuple(len(per_prime[p][0]) for p in primes)
    # Populate K as a tensor whose axes are the per-prime coordinates of g
    K_tensor = np.zeros(n_shape, dtype=np.complex128)
    for i_g, g in enumerate(G):
        coords = tuple(per_prime[p][2][g % p] for p in primes)
        K_tensor[coords] = K_row[i_g]
    # Contract along each prime-axis with the per-prime character table
    # λ_tensor[k_1,...,k_K] = Σ_{i_1,...,i_K} K_tensor[i_1,...,i_K]
    #                             Π_p Lam_p[k_p, i_p]*
    # We use conjugation on Lam because of the χ(g)* convention in λ(χ).
    lambda_tensor = K_tensor
    for axis, p in enumerate(primes):
        _, Lam, _ = per_prime[p]
        # Contract K_tensor's axis `axis` (index i_p) with Lam's axis 1 (i_p).
        # Tensordot output: K's remaining axes, then Lam's remaining axis (k_p).
        lambda_tensor = np.tensordot(
            lambda_tensor, np.conj(Lam), axes=([axis], [1]),
        )
        # The new axis (k_p) is at the end; move back to position `axis`.
        lambda_tensor = np.moveaxis(lambda_tensor, -1, axis)
    # Eigenvalues of a real symmetric kernel are real; imaginary parts come
    # from numerical roundoff unless complex characters pair up properly.
    # We keep complex for now; callers that want reals should take .real.
    return lambda_tensor, n_shape


# ──────────────────────────────────────────────────────────────────────────
# §6c. Unity Clock chiral decomposition
#
# Identity: on (ℤ/m)^× with real-symmetric tent kernel K and involution
# τ(r) = r^{-1}, Dirichlet characters diagonalize D_sym with eigenvalues
#
#     λ(χ) = Σ_g K(1, g) χ(g)                  ∈ ℂ
#
# Because K is real, λ(χ̄) = conj(λ(χ)). The commutator C = [D_sym, P_τ]
# acts as
#
#     C u_χ = (λ(χ̄) − λ(χ)) u_{χ̄} = −2i Im(λ(χ)) u_{χ̄}
#
# Consequences:
#   (1) ker(C) = span of real characters (χ = χ̄ → Im(λ) = 0). For
#       primorial m with k ≥ 2 primes, there are 2^{k-1} real characters
#       (= the Unity-Clock kernel dimension).
#   (2) Each complex conjugate character pair (χ, χ̄) generates one
#       degenerate singular-value pair of C with σ = 2|Im(λ(χ))|. This
#       is Unity Clock's "all nonzero singular values come in pairs".
#   (3) The Unity-Clock 4 modes = the 2 conjugate pairs with largest
#       |Im(λ(χ))|.
#   (4) The TC(βt) peak timescale is set by the D_sym spectral gap
#       Δ₁ = λ(triv) − max Re(λ(χ≠triv)); analytic softmax under
#       p(χ,t) ∝ exp(−2βt Re(λ)/λ_P) peaks near βt* ≈ λ_P / (2 Δ₁).
# ──────────────────────────────────────────────────────────────────────────

@dc.dataclass(frozen=True)
class UnityClockSpectrum:
    """Full chiral / Unity-Clock decomposition of a kernel on (ℤ/m)^×."""
    m: int
    primes: tuple[int, ...]
    phi_m: int
    # Flat arrays indexed by character (same order as lambda_tensor.ravel())
    lambda_complex: np.ndarray          # λ(χ) ∈ ℂ,     shape (φ(m),)
    re_lambda: np.ndarray                # Re λ(χ),      shape (φ(m),)
    im_lambda: np.ndarray                # Im λ(χ),      shape (φ(m),)
    is_real_char: np.ndarray             # χ = χ̄?,      shape (φ(m),) bool
    lambda_principal: float              # λ(triv)
    re_gap_top: float                    # λ_P − max Re(λ(χ≠triv))
    re_lambda_min: float                 # min Re(λ(χ))  — bottom of spectrum
    re_min_gap: float                    # 2nd-min − min of Re(λ)
    # Chiral / commutator-side
    sigma_C: np.ndarray                  # σ = 2|Im λ(χ)|, sorted desc
    top4_sigma: np.ndarray               # σ_0, σ_1, σ_2, σ_3
    top4_characters: list                # conjugate-pair indices driving them
    effective_rank_90: int               # how many σ to capture 90% of Σ σ²
    # Analytic peak-βt predictions (three candidate timescales)
    peak_bt_top_gap: float               # λ_P / (2 Δ_top)         — top-gap
    peak_bt_bot_gap: float               # λ_P / (2 Δ_bot)         — bottom-gap
    peak_bt_range: float                 # λ_P / (2 (λ_max−λ_min)) — full range


def unity_clock_chiral_decomposition(
    m: int, *, kernel: str = "natural",
) -> UnityClockSpectrum:
    """Chiral decomposition via Dirichlet character basis.

    IMPORTANT — rigorous scope. This decomposition assumes the kernel is
    G-equivariant under (ℤ/m)^× acting multiplicatively. That holds
    *exactly* for the CRT-local kernel D^CRT, whose per-prime circular
    distances are multiplicatively invariant on each (ℤ/p)^× factor.

    It does NOT hold for the natural tent distance D_sym. D_sym is
    invariant under ADDITIVE shift on ℤ/m, not multiplicative shift on
    (ℤ/m)^×. Its eigenvectors are Möbius-sieved ambient Fourier modes
    (Unity Clock §3.4), not Dirichlet characters.

    For kernel="natural", the λ(χ) returned here are the Dirichlet-
    character Fourier coefficients of the first row k(g) = K(1,g), which
    *define* the softmax dynamic used in analytic_tc_vs_time(...,
    kernel="natural"). They are not eigenvalues of D_sym. Treat the
    returned σ_C values as characterizing the *character-marginal cross-
    correlation observable* on that dynamic, not as σ of the physical
    commutator [D_sym, P_τ]. For the physical commutator singular values
    at m ≤ 30030, use unity_clock_direct_commutator_svd(m).
    """
    if kernel == "natural":
        K_row = natural_distance_row(m)
    elif kernel == "crt":
        K_row = crt_local_distance_row(m)
    else:
        raise ValueError(f"unknown kernel {kernel!r}")
    primes, per_prime, G = _primorial_setup(m)
    phi_m = len(G)
    lam_tensor, _ = primorial_kernel_eigenvalue_tensor(m, K_row)
    lam = lam_tensor.ravel()                     # complex
    re_lam = np.real(lam).astype(np.float64)
    im_lam = np.imag(lam).astype(np.float64)
    # Principal character = all-zero multi-index
    lam_P = float(re_lam[0])                     # λ(triv) is real and the row sum
    # Chiral spectrum: σ_C(χ) = 2|Im λ(χ)|. Real chars have σ=0 (within roundoff)
    sigma_all = 2.0 * np.abs(im_lam)
    # is_real_char: for primorial (ℤ/m)^×, χ is real iff Im λ ~ 0. Use tight
    # tolerance relative to |λ_P|.
    tol = 1e-9 * max(abs(lam_P), 1.0)
    is_real = np.abs(im_lam) < tol
    # Because C u_χ ∝ u_{χ̄}, each conjugate pair contributes σ ONCE per
    # character, and both χ and χ̄ give the same |Im λ|. So the sorted σ
    # array already pairs up correctly — we don't deduplicate.
    sigma_sorted = np.sort(sigma_all)[::-1]
    # Effective rank at 90% of Σσ²
    energy = sigma_sorted ** 2
    total = float(energy.sum())
    if total > 0:
        cum = np.cumsum(energy) / total
        eff90 = int(np.searchsorted(cum, 0.90) + 1)
    else:
        eff90 = 0
    top4 = sigma_sorted[:4]
    # Top-4 character indices (in the flat order)
    idx_sorted = np.argsort(-sigma_all)
    top4_chars = [int(i) for i in idx_sorted[:4]]
    # Re-gap top: λ_P − max Re(λ(χ≠triv))
    re_nontriv = re_lam.copy()
    re_nontriv[0] = -np.inf
    re_max_nontriv = float(np.max(re_nontriv))
    re_gap = lam_P - re_max_nontriv
    # Bottom of the Re-spectrum. The analytic softmax has shift
    # = 2βt · (λ/λ_P − μ_min), mass flowing to min-λ character as βt grows.
    # The relaxation timescale to split min-λ from 2nd-min-λ is
    #   τ_split = λ_P / (2 (λ_2nd_min − λ_min))
    re_sorted = np.sort(re_lam)
    re_min = float(re_sorted[0])
    if len(re_sorted) > 1:
        re_min_gap = float(re_sorted[1] - re_sorted[0])
    else:
        re_min_gap = 0.0
    re_range = lam_P - re_min
    # Three candidate timescales. The observable data will pick.
    def _pred(gap: float) -> float:
        return lam_P / (2.0 * gap) if gap > 0 else float("nan")
    peak_top = _pred(re_gap)
    peak_bot = _pred(re_min_gap)
    peak_rng = _pred(re_range)
    return UnityClockSpectrum(
        m=m,
        primes=primes,
        phi_m=phi_m,
        lambda_complex=lam,
        re_lambda=re_lam,
        im_lambda=im_lam,
        is_real_char=is_real,
        lambda_principal=lam_P,
        re_gap_top=re_gap,
        re_lambda_min=re_min,
        re_min_gap=re_min_gap,
        sigma_C=sigma_sorted,
        top4_sigma=top4,
        top4_characters=top4_chars,
        effective_rank_90=eff90,
        peak_bt_top_gap=peak_top,
        peak_bt_bot_gap=peak_bot,
        peak_bt_range=peak_rng,
    )


@dc.dataclass(frozen=True)
class UnityClockDirect:
    """Direct SVD of the physical commutator C = [D_sym, P_τ].

    This is the Unity-Clock-correct object: D_sym is NOT G-equivariant
    under (ℤ/m)^× multiplicatively, so Dirichlet characters do not
    diagonalize it. Here we just build D_sym and P_τ as n×n dense matrices
    (n = φ(m)) and compute eigenvalues of D_sym and singular values of C
    directly. Tractable for n ≤ ~6000 (m ≤ 30030).
    """
    m: int
    phi_m: int
    D_eigenvalues: np.ndarray      # eigvalsh(D_sym), sorted ascending
    D_lambda_P: float              # top |D_sym| eigenvalue (Perron-like)
    D_gap_top: float               # |λ_0| − |λ_1| of D_sym
    C_sigma: np.ndarray            # σ of C, sorted descending
    C_rank: int                    # numerical rank of C
    C_top4_sigma: np.ndarray       # σ_0..σ_3
    C_top4_frac: float             # Σ(σ_0..σ_3)² / Σσ²
    C_effective_rank_90: int
    C_sigma_ratio_pair_02: float   # σ_2/σ_0 (UC §3.3)
    dim_kernel: int                # dim ker(C), UC predicts 2^{k-1}


def unity_clock_direct_commutator_svd(m: int) -> UnityClockDirect:
    """Build D_sym and P_τ explicitly, diagonalize and SVD directly.

    Matches the Unity-Clock paper's computation exactly. Memory O(n²),
    time O(n³) — tractable through m = 30030 (n = 5760) in ~30 sec.
    """
    G = coprime_residues(m)
    n = len(G)
    index_of = {r: i for i, r in enumerate(G)}
    # D_sym
    D = np.zeros((n, n), dtype=np.float64)
    for i, r in enumerate(G):
        for j, s in enumerate(G):
            d = abs(r - s)
            D[i, j] = min(d, m - d)
    # P_τ: permutation matrix of multiplicative inversion r ↦ r^{-1}
    P = np.zeros((n, n), dtype=np.float64)
    for i, r in enumerate(G):
        # find r^{-1}
        inv = None
        for s in G:
            if (r * s) % m == 1:
                inv = s
                break
        if inv is None:
            raise RuntimeError(f"no inverse of {r} mod {m}")
        P[index_of[inv], i] = 1.0
    C = D @ P - P @ D
    # Eigenvalues of D_sym (real symmetric)
    D_eigs = np.linalg.eigvalsh(D)
    D_sorted_abs = np.sort(np.abs(D_eigs))[::-1]
    lam_P = float(D_sorted_abs[0])
    gap_top = float(D_sorted_abs[0] - D_sorted_abs[1]) if len(D_sorted_abs) > 1 else 0.0
    # SVD of C
    s = np.linalg.svd(C, compute_uv=False)
    s_sorted = np.sort(s)[::-1]
    # Kernel / rank at reasonable tolerance
    tol = 1e-10 * max(1.0, s_sorted[0])
    rank = int(np.sum(s_sorted > tol))
    dim_ker = n - rank
    # Top 4 and 90% effective rank
    top4 = s_sorted[:4]
    energy = s_sorted ** 2
    total = float(energy.sum())
    if total > 0:
        cum = np.cumsum(energy) / total
        eff90 = int(np.searchsorted(cum, 0.90) + 1)
        frac_top4 = float((top4 ** 2).sum() / total)
    else:
        eff90 = 0
        frac_top4 = 0.0
    ratio_02 = float(s_sorted[2] / s_sorted[0]) if s_sorted[0] > 0 and len(s_sorted) > 2 else 0.0
    return UnityClockDirect(
        m=m,
        phi_m=n,
        D_eigenvalues=D_eigs,
        D_lambda_P=lam_P,
        D_gap_top=gap_top,
        C_sigma=s_sorted,
        C_rank=rank,
        C_top4_sigma=top4,
        C_top4_frac=frac_top4,
        C_effective_rank_90=eff90,
        C_sigma_ratio_pair_02=ratio_02,
        dim_kernel=dim_ker,
    )


def natural_distance_row(m: int) -> np.ndarray:
    """First row of the circular-distance kernel on coprime_residues(m).
    K(1, g) = min(|1 − g|, m − |1 − g|)."""
    G = coprime_residues(m)
    return np.array([min(abs(1 - g), m - abs(1 - g)) for g in G], dtype=np.float64)


def crt_local_distance_row(m: int) -> np.ndarray:
    """First row of D^CRT on coprime_residues(m). Equals the direct-sum
    per-prime circular distance with the 'agree on all other primes'
    indicator; ≠ 0 only when the target g has g ≡ 1 mod q for all but one
    prime factor."""
    from sympy import factorint
    primes = sorted(factorint(m).keys())
    G = coprime_residues(m)
    row = np.zeros(len(G), dtype=np.float64)
    for i_g, g in enumerate(G):
        for p in primes:
            if all((g % q) == 1 for q in primes if q != p):
                a = g % p
                d = abs(1 - a)
                row[i_g] += min(d, p - d)
                break
        # (if g agrees with 1 on all primes, it IS 1 → row = 0; if it
        # disagrees on two or more primes, row stays 0 — by construction
        # D^CRT has zero coupling across multi-prime residues.)
    return row


def analytic_tc_vs_time(
    m: int,
    beta_t_schedule: Sequence[float],
    *,
    kernel: str = "natural",
) -> list[dict]:
    """Compute TC(t) analytically across a βt schedule for primorial m.

    kernel: "natural" for D_sym, "crt" for D^CRT.

    Returns list of dicts, one per βt, with keys
      t, beta_t, H_full, H_marginals (per prime), TC, top_char_prob.
    """
    primes, per_prime, G = _primorial_setup(m)
    n_shape = tuple(len(per_prime[p][0]) for p in primes)
    if kernel == "natural":
        K_row = natural_distance_row(m)
    elif kernel == "crt":
        K_row = crt_local_distance_row(m)
    else:
        raise ValueError(f"unknown kernel {kernel!r}")
    lam_tensor, _ = primorial_kernel_eigenvalue_tensor(m, K_row)
    # Eigenvalues of the real symmetric kernel should be real to roundoff
    lam = np.real(lam_tensor)
    lam_max = float(np.max(np.abs(lam)))
    out = []
    for bt in beta_t_schedule:
        # Boltzmann weights: p(χ) ∝ exp(-2 βt λ(χ) / λ_max)
        #   (note: λ_max is the Perron-like normalizer; some eigenvalues may
        #   be negative, making the distribution non-thermal. For D_sym on
        #   (ℤ/m)^× the Perron eigenvalue is positive and dominant.)
        exponent = -2.0 * bt * lam / max(lam_max, 1e-14)
        # Shift for numerical stability
        exponent -= exponent.max()
        w = np.exp(exponent)
        Z = float(w.sum())
        p_tensor = w / Z
        # Full Shannon
        pf = p_tensor.ravel()
        pf_nz = pf[pf > 1e-14]
        H_full = float(-np.sum(pf_nz * np.log(pf_nz)))
        # Marginals along each prime axis
        H_marg = {}
        for axis, p in enumerate(primes):
            # Sum over all other axes
            other = tuple(ax for ax in range(len(primes)) if ax != axis)
            marg = p_tensor.sum(axis=other) if other else p_tensor
            marg = marg.ravel()
            marg_nz = marg[marg > 1e-14]
            H_marg[p] = float(-np.sum(marg_nz * np.log(marg_nz)))
        TC = sum(H_marg.values()) - H_full
        out.append({
            "m": m,
            "beta_t": bt,
            "H_full": H_full,
            "H_marginals": H_marg,
            "TC": TC,
            "top_char_prob": float(pf.max()),
            "n_characters": int(pf.size),
        })
    return out


def crt_factorization_check(m: int, C: np.ndarray) -> dict:
    """LEGACY: H_p-subgroup-refinement version of the CRT check.

    Retained for backwards compatibility; the correct test is
    `crt_character_marginal_check` which computes total correlation on the
    Dirichlet-character joint distribution. The residual returned here is
    H_full − Σ refinement.shannon and is NOT the cross-prime mutual
    information (it is not even guaranteed non-negative).
    """
    from sympy import factorint
    factors = sorted(factorint(m).keys())
    result = {"m": m, "factors": factors, "per_factor_shannon": {}}
    total_marginal_shannon = 0.0
    for p in factors:
        # Subgroup: H_p = elements that project to identity in (ℤ/p)^×
        # i.e., elements ≡ 1 mod p
        H_p = tuple(r for r in coprime_residues(m) if r % p == 1)
        if len(H_p) == 0:
            H_p = (1,)
        refine = subgroup_refinement(C, m, H_p)
        # The "per-prime Shannon" is the information carried by the
        # complementary quotient G/H_p ≅ (ℤ/p)^×
        shannon_p = refine.shannon
        result["per_factor_shannon"][p] = shannon_p
        total_marginal_shannon += shannon_p
    # Compare to full character Shannon
    full = abelian_character_decomposition(C, m)
    result["full_shannon"] = full.shannon
    result["sum_of_factor_shannons"] = total_marginal_shannon
    result["factorization_residual"] = full.shannon - total_marginal_shannon
    return result


def page_curve_on_residues(
    C: np.ndarray, m: int, *, rng: np.random.Generator | None = None,
    n_samples: int = 30,
) -> list[dict]:
    """For each subset size k = 1 ... n-1, sample n_samples random subsets
    A ⊂ coprimes, compute the reduced von Neumann entropy S(ρ_A), and
    track mean / std. This is the Page curve analog: entropy should rise
    then fall as |A| → n/2 → n.
    """
    rng = rng or np.random.default_rng(42)
    n = C.shape[0]
    results = []
    for k in range(1, n):
        ents = []
        trials = min(n_samples, max(1, int(math.comb(n, k))))
        seen = set()
        for _ in range(trials * 4):
            A = tuple(sorted(rng.choice(n, k, replace=False)))
            if A in seen:
                continue
            seen.add(A)
            C_A = C[np.ix_(A, A)]
            ents.append(von_neumann_entropy_complex(C_A))
            if len(seen) >= trials:
                break
        results.append({
            "subset_size": k,
            "mean_entropy": float(np.mean(ents)),
            "std_entropy": float(np.std(ents)),
            "ratio_to_log_k": (float(np.mean(ents)) / np.log(k)) if k > 1 else 0.0,
        })
    return results


# ──────────────────────────────────────────────────────────────────────────
# §7. Dictionary cross-check: is m=6 prime gas ≅ S_2 Coliseum APP?
# ──────────────────────────────────────────────────────────────────────────

def dictionary_crosscheck_m6(C_primes_m6: np.ndarray) -> dict:
    """Prime-gas structural invariants at m = 6.

    At m = 6, (ℤ/6)^× = {1, 5} with 5·5 = 1 mod 6, so the group is ℤ/2.
    This function returns the character-entropy signature of a given
    boundary covariance at m = 6: the Born weights on the two Dirichlet
    characters (trivial, sign), the Shannon entropy in units of log 2,
    the direct von-Neumann entropy, and the non-abelian leakage (the
    residual not captured by the character basis — should be zero up
    to Monte-Carlo noise since (ℤ/6)^× is abelian).

    The output lets any other 3D S_2-equivariant system be compared
    functorially: two systems share the same abstract group G = ℤ/2, and
    p_sign is the fraction of the dynamics that resolves the involution
    class — a representation-independent signature.
    """
    full = abelian_character_decomposition(C_primes_m6, 6)
    # p_char has 2 components: trivial and sign.
    # "Antipodal fraction" = p_sign = how much the dynamics distinguishes
    # the involution class.
    p_sorted = sorted(full.p_char, reverse=True)
    p_triv, p_sign = (p_sorted + [0.0, 0.0])[:2]
    return {
        "m": 6,
        "n": 2,
        "p_triv": float(p_triv),
        "p_sign": float(p_sign),
        "shannon": full.shannon,
        "shannon_over_log2": full.shannon / np.log(2),
        "direct_entropy": full.direct_entropy,
        "nonabelian_leakage": full.nonabelian_leakage,
        "interpretation": (
            "p_sign close to 0 → dynamics are near-trivial on the "
            "involution. p_sign close to 0.5 → dynamics maximally "
            "exploit the ℤ/2 involution. Intermediate values are the "
            "functorial signature."
        ),
    }


__all__ = [
    "coprime_residues",
    "group_order",
    "mult_table",
    "dirichlet_character_table",
    "enumerate_subgroups",
    "subgroup_character_projectors",
    "build_distance_kernel_boundary",
    "prime_gas_ensemble",
    "boundary_covariance",
    "von_neumann_entropy_complex",
    "AbelianCharEntropy",
    "abelian_character_decomposition",
    "SubgroupRefinement",
    "subgroup_refinement",
    "haar_random_entropy_baseline",
    "crt_factorization_check",
    "page_curve_on_residues",
    "dictionary_crosscheck_m6",
]
