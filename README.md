# Character-Theoretic Effective Rank

**From Schur's Lemma to a Stereoscopic Measurement Apparatus on Continuous Hamiltonian Substrates**

[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.19744573-blue)](https://doi.org/10.5281/zenodo.19744573)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

This repository contains the preprint and the math-only reproducibility bundle for *Character-Theoretic Effective Rank* (Matos, 2026). The paper isolates a representation-theoretic mechanism — the Schur-forced effective-rank ceiling — that constrains the covariance spectrum of any detailed-balance Markov generator on a state space carrying a finite-group action, sharpens it to a closed-form sampling theorem, derives its number-theoretic strengthening on coprime-residue lattices, and demonstrates that the resulting operator class extends to continuous Hamiltonian substrates.

Five named theorems organize the work:

1. **Character-Theoretic Effective Rank Bound** (§3) — distinct-eigenvalue ceiling forced by Schur's lemma on any G-equivariant covariance.
2. **Wishart-Nyquist Sampling Theorem** (§4) — closed-form sample-size requirement for the Schur-forced degeneracies to emerge above the empirical noise floor.
3. **Triangle-Wave Overtone Theorem** (§6) — exact parity-chiral degeneracy, fundamental amplitude *m·φ(m)/π²*, and rank-4 block ratios *(2i−1)²/(2i+1)²* on the commutator *[D_sym, P_τ]* for coprime-residue lattices.
4. **Universal Coupling Invariant** (§7) — two independent character-theoretic derivations agreeing on *κ² = 2/3*.
5. **A-Optimality of Texture and Heisenberg Obstruction** (§§13–14) — harmonic-mean texture ceiling *√2* on the rank-4 subspace of a continuous Hamiltonian substrate, with a Fourier-dual cost converting the ceiling into a manifold-curvature measurement program.

Empirical verification spans a finite *S₃* interacting-particle system (Fancyland Coliseum energy-gauge accumulator at *N = 3* uids), the prime gas at *m₀ = 6*, the coprime-residue lattices for *m* up to 30030, and a 768-dimensional semantic embedding manifold of *N = 353* entities.

## Repository Structure

```
paper/
  CHARACTER_RANK.pdf     ← Preprint PDF (the artifact archived on Zenodo)
  CHARACTER_RANK.md      ← Markdown source
  CHARACTER_RANK.tex     ← LaTeX source (pandoc-generated)
  header.tex             ← pandoc build helper
scripts/
  sweep/                 ← Math-only reproducibility package (pytest-tested)
    character_rank.py            ← Theorem 3.1 support: isotypic projectors,
                                   Wishart-scaled clustering, Schur-bound check
    holographic_rank.py          ← Theorem 5.1 support: character-entropy
                                   decomposition, Gaussian mutual information
    holographic_rank_primes.py   ← Theorem 6.1 / §7 support: CRT-local kernel,
                                   tent distance, multiplicative involution,
                                   direct commutator SVD
    test_character_rank.py       ← 9 math-only tests for Theorem 3.1
    test_holographic_rank.py     ← 11 math-only tests for Theorem 5.1
  run_holographic_primes.py      ← §6 CRT-tower spectrum, §7 κ²=2/3 scan
  run_unity_clock_chirality.py   ← §6 Lemma 6.4.1 reproduction, §6.3 fundamentals
  run_horizon_trajectory.py      ← §6.3 flatline fit on 47 anchors
  run_overtone_check.py          ← §6.3 Table 6.1 overtone ladder
  verify_triangle_wave_theorem.py ← §6 three-witness cross-verification
  probe_cosine_sobolev.py        ← §12.4 Lemma 12.2 — falsifier for the Sobolev
                                   regularity hypothesis (cosine kernel rank-(1+d))
LICENSE                  ← MIT
```

## Reproducing the Theorems

Requires **Python ≥ 3.10**, **NumPy ≥ 2.0**, **SciPy ≥ 1.10**, and **pytest** (for the `sweep/` test suite). No GPU. No external services.

### Schur bound and entropy decomposition (Theorems 3.1, 5.1)

```bash
cd scripts
python -m pytest sweep/ -q
```

Twenty hand-built equivariant covariances with analytically-known spectra are checked against the Schur bound (`test_character_rank.py`) and against the character-entropy decomposition identity (`test_holographic_rank.py`). Wall time: < 1 s.

### Triangle-Wave Overtone Theorem (Theorem 6.1)

```bash
cd scripts
python run_holographic_primes.py --M 5000      # CRT tower scan + κ²=2/3 fit
python run_unity_clock_chirality.py            # Direct commutator SVD reproduction
python run_horizon_trajectory.py               # σ₄/σ₃ flatline fit, 47 anchors
python run_overtone_check.py                   # Table 6.1 overtone ladder
python verify_triangle_wave_theorem.py         # Three-witness cross-verification
```

`verify_triangle_wave_theorem.py` runs three independent code paths against each anchor modulus *m*:

- **Witness (A)** — closed-form rational predictions via `fractions.Fraction`
- **Witness (B)** — from-scratch NumPy re-implementation
- **Witness (C)** — production pipeline from `sweep/holographic_rank_primes.py`

Witnesses (B) and (C) agree to bit-identical IEEE 754 floats on every tested anchor; (B) agrees with the rational prediction (A) within the *O(1/φ(m))* envelope of Theorem 6.1(iv). The final line of the verifier prints `THREE-WAY CROSS-VERIFICATION PASSES.`

End-to-end wall time for the full §6 reproduction: approximately 20 minutes on a Ryzen 9 7950X single core.

### Sobolev-hypothesis falsifier (Lemma 12.2)

```bash
cd scripts
python probe_cosine_sobolev.py
```

Verifies the cosine kernel on *S^(d−1)* has spherical-harmonic content only at degrees 0 and 1 (analytical rank exactly *1 + d*) — establishing the finite-rank truncation mechanism of Lemma 12.2 and falsifying the Sobolev-regularity hypothesis discussed in §17. Wall time: < 5 s.

## Companion Patents

The mathematics in this preprint underlies a trilogy of US provisional patent applications protecting the corresponding engineering apparatus. The mathematics is published openly under MIT license; the apparatus and method claims are reserved.

| Layer | US Provisional | Title | Filed |
|---|---|---|---|
| State storage | **64/031,440** | Fault-Injection-Immune Computational Unit Using Primorial Coprime Residue Topology (the Arithmetic Qubit) | 2026-04-07 |
| Spectral compute | **64/033,689** | Holographic Eigen-Solver Using QM Boundary Projection on Coprime Residue Lattices | 2026-04-08 |
| Stereoscopic measurement | **64/048,617** | Tensegrity Interferometer: Stereoscopic Query Resolution and Manifold-Curvature Measurement on Continuous Hamiltonian Substrates | 2026-04-24 |

Implementations practicing the disclosed mathematics in the manner claimed in the patents above require a license from the assignee (Fancyland LLC).

## Citation

```bibtex
@misc{matos2026character,
  author       = {Matos, Antonio P.},
  title        = {Character-Theoretic Effective Rank: From Schur's Lemma to a Stereoscopic Measurement Apparatus on Continuous Hamiltonian Substrates},
  year         = {2026},
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.19744573},
  url          = {https://doi.org/10.5281/zenodo.19744573}
}
```

## Companion Preprints

This work cites and extends five companion preprints on the same Zenodo channel:

- *The Unity Clock* — [10.5281/zenodo.19478727](https://doi.org/10.5281/zenodo.19478727)
- *The Arithmetic Black Hole* — [10.5281/zenodo.19442006](https://doi.org/10.5281/zenodo.19442006)
- *Active Transport / Arithmetic Qubit* — [10.5281/zenodo.19243258](https://doi.org/10.5281/zenodo.19243258)
- *Universal Two-Prime Formula* — [10.5281/zenodo.19210625](https://doi.org/10.5281/zenodo.19210625)
- *The Patient Compass* — [10.5281/zenodo.19561701](https://doi.org/10.5281/zenodo.19561701)

## License

Code in this repository is released under the [MIT License](LICENSE). The preprint PDF is licensed under [CC-BY 4.0](https://creativecommons.org/licenses/by/4.0/) via the Zenodo deposit.

## Author

**Antonio P. Matos** — Independent Researcher; Fancyland LLC / Lattice OS
ORCID: [0009-0002-0722-3752](https://orcid.org/0009-0002-0722-3752)

---

*Fancyland LLC — Lattice OS research infrastructure.*
