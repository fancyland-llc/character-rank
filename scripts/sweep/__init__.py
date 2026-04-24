"""Character-theoretic effective-rank bounds and chiral decomposition for
G-equivariant Markov generators — math-only reproducibility package for
`CHARACTER_RANK.md`.

Modules:

  - character_rank            : Theorem 3.1 — distinct-eigenvalue bound
                                from Schur's lemma on a G_α-equivariant
                                covariance, together with the Wishart-
                                scaled distinct-eigenvalue clustering
                                and the pass/fail check.
  - holographic_rank          : Theorem 11.1 — character-entropy
                                decomposition S(ρ) = H(p) + Σ p_π log d_π
                                + Σ p_π S_π^int on N = 3, with isotypic
                                projectors/bases and Gaussian mutual
                                information.
  - holographic_rank_primes   : Prime-gas Dirichlet-character machinery
                                on (ℤ/m)^×: character tables, CRT tensor
                                kernels, the natural tent distance D_sym
                                and its CRT-local G-equivariant twin
                                D^CRT, the analytic TC(βt) dynamic, and
                                the Unity Clock chiral decompositions.

Tests `test_character_rank.py`, `test_holographic_rank.py` pin the
main-body claims on hand-built equivariant covariances.
"""
