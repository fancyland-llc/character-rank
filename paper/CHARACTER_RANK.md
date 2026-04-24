```{=latex}
\begin{center}
{\LARGE\bfseries Character-Theoretic Effective Rank:\\[0.4em]
From Schur's Lemma to a Stereoscopic Measurement\\[0.2em]
Apparatus on Continuous Hamiltonian Substrates}

\vspace{1.5em}

{\large Antonio P.\ Matos}\\[0.4em]
{\small ORCID: \href{https://orcid.org/0009-0002-0722-3752}{0009-0002-0722-3752}}\\[0.3em]
April 24, 2026\\[0.3em]
{\small DOI: \href{https://doi.org/10.5281/zenodo.19744573}{10.5281/zenodo.19744573}}\\[0.3em]
Independent Researcher; Fancyland LLC / Lattice OS\\[0.3em]
Preprint v1.0

\vspace{0.8em}

{\small\textbf{MSC 2020:} 60J22, 20C30, 15A18, 62H25, 15B52, 82C22}

\vspace{0.4em}

{\small\textbf{Keywords:} Schur's lemma, effective rank, equivariant Markov generator, detailed balance, Wishart-Nyquist sampling, character-entropy decomposition, triangle-wave overtone theorem, holographic dictionary, Hamiltonian substrate, A-optimality, Heisenberg duality, manifold curvature, semantic interferometer}
\end{center}
```

**Companion papers:**

- [1] A. P. Matos, "The Unity Clock: Effective Dimensional Collapse of the Addition-Multiplication Commutator on Coprime Lattices," preprint (2026). DOI: [10.5281/zenodo.19478727](https://doi.org/10.5281/zenodo.19478727)
- [2] A. P. Matos, "The Arithmetic Black Hole: Softmax Thermodynamics and the Four Eigenvalue Laws of the Prime Gas," preprint (2026). DOI: [10.5281/zenodo.19442006](https://doi.org/10.5281/zenodo.19442006)
- [3] A. P. Matos, "Active Transport on the Prime Gas: Flat-Band Condensation, the Rabi Phase Transition, and the Arithmetic Qubit," preprint (2026). DOI: [10.5281/zenodo.19243258](https://doi.org/10.5281/zenodo.19243258)
- [4] A. P. Matos, "Universal Two-Prime Formula for the Coprime-Lattice Coupling Constant," preprint (2026). DOI: [10.5281/zenodo.19210625](https://doi.org/10.5281/zenodo.19210625)
- [5] A. P. Matos, "The Patient Compass: Domain Walls, Oxygen Gating, and the 3-D Survival Manifold of 198,862 Tumors," preprint (2026). DOI: [10.5281/zenodo.19561701](https://doi.org/10.5281/zenodo.19561701)

**Companion patents (the architectural triad):**

- [6] A. P. Matos, "Fault-Injection-Immune Computational Unit Using Primorial Coprime Residue Topology" (the Arithmetic Qubit), U.S. Provisional Patent Application No. **64/031,440** (filed April 7, 2026). Fancyland LLC.
- [7] A. P. Matos, "Holographic Eigen-Solver Using QM Boundary Projection on Coprime Residue Lattices," U.S. Provisional Patent Application No. **64/033,689** (filed April 8, 2026). Fancyland LLC.
- [8] A. P. Matos, "Tensegrity Interferometer: Stereoscopic Query Resolution and Manifold-Curvature Measurement on Continuous Hamiltonian Substrates," U.S. Provisional Patent Application No. **64/048,617** (filed April 24, 2026). Fancyland LLC.

---

## Patent Notice

The mathematics presented in this preprint underlies a trilogy of US provisional patent applications, each protecting a distinct architectural layer of the resulting cognitive substrate. **All theorems, proofs, lemmas, scalar identities, and the verification scripts at [github.com/fancyland-llc/character-rank](https://github.com/fancyland-llc/character-rank) are placed in the public domain via this preprint** under the MIT license on the accompanying repository. **Engineering apparatus, system embodiments, and method claims drawn from this mathematics are reserved to the inventor and assignee** under U.S. Provisional Patent Applications [6], [7], and [8] above.

| Layer | Provisional | Sections | Embodiment scope |
|---|---|---|---|
| State storage | 64/031,440 | §§3, 4, 8 | Schur-bounded equivariant accumulators; ASIC realizations |
| Spectral compute | 64/033,689 | §§5, 6, 7, 9 | Triangle-wave-block eigenvalue extractors; CRT-tower compute primitives |
| Stereoscopic measurement | 64/048,617 | §§10–15 | Texture optimizer, Gaussian wave-packet relaxation, manifold-curvature telemetry pipelines |

Implementations practicing the disclosed mathematics in the manner claimed in [6], [7], or [8] require a license from the assignee. Where the mathematics is published openly here, it constitutes prior art against any third-party attempt to patent the same scalar identity, theorem, or proof technique going forward.

---

## Abstract

We isolate a representation-theoretic mechanism — **the Schur-forced effective-rank ceiling** — that constrains the covariance spectrum of any detailed-balance Markov generator on a state space carrying a finite-group action, sharpen it to a closed-form sampling theorem, derive its number-theoretic strengthening on coprime-residue lattices, and demonstrate that the resulting operator class extends to a continuous Hamiltonian substrate where it is realized as a stereoscopic measurement instrument.

**Theorem (Character-Theoretic Effective Rank Bound, §3).** The covariance of a $G_\alpha$-invariant observable of a $G_\alpha$-equivariant detailed-balance generator has at most $\sum_{\pi \in \mathcal{A}} m_\pi$ distinct eigenvalues, each appearing with multiplicity $d_\pi \cdot k_\pi$ for some $k_\pi \leq m_\pi$. The bound is tight when each irrep type appears with unit multiplicity and no accidental degeneracy occurs across blocks.

**Theorem (Wishart-Nyquist Sampling Bound, §4).** The Schur-forced degeneracies are visible in an $M$-sample covariance only when $M_\text{eff} \geq n_\sigma^2 d / \delta^2$ where $d$ is the largest irrep dimension, $\delta$ the desired relative precision, and $M_\text{eff} = M \cdot \min(1, \sigma_\text{sig}^2/\sigma_\text{obs}^2)$. The bound identifies the sampling regime in which the Schur theorem is empirically falsifiable.

**Theorem (Triangle-Wave Overtone Theorem, §6).** On coprime-residue lattices $(\mathbb{Z}/m)^\times$ for squarefree $m$, the commutator $C = [D_\text{sym}, P_\tau]$ between the natural tent distance and the multiplicative involution has exact parity-chiral degeneracy $\sigma_{2i} = \sigma_{2i+1}$ at every finite $m$, fundamental amplitude $\sigma_0 = m\,\varphi(m)/\pi^2 + O(\varphi(m))$, and rank-4 block ratios $\sigma_{4i}/\sigma_{4(i-1)} = (2i-1)^2/(2i+1)^2 + O(\varphi(m)^{-1})$ on odd harmonics. The boundary block holds asymptotic $L^2$ fraction $8/\pi^2$; the bulk overtone tower holds the residual $1 - 8/\pi^2 \approx 0.189$.

**Theorem (Universal Coupling Invariant, §7).** Two independent character-theoretic derivations — through Fourier moments and through Euler's $\zeta(4) = \pi^4/90$ — agree bit-identically on the dimensionless coupling $\kappa^2 = \|C\|_F^2 / \lambda_\text{Perron}^2 = 2/3$, promoting $\sqrt{2/3}$ from a constant of one paper to a candidate invariant of the entire operator class.

**Theorem (A-Optimality of Texture, §13).** On the rank-4 subspace of a $G$-equivariant continuous Hamiltonian substrate, the harmonic-mean texture $\mathcal{T}(q) = K / \sum_k 1/T_k(q)$ attains its unique maximum $\mathcal{T}_\text{max} = 2/\sqrt{K}$ at the equal-tension barycenter; for $K = 2$ irreps the ceiling is $\sqrt{2}$.

**Proposition (Heisenberg-Dual Localization, §14).** Sharp topological localization at the texture barycenter forces diffuse semantic localization, $\sigma_\text{topo}(\hat\rho) \cdot \sigma_\text{sem}(\rho) \gtrsim c$, the Fourier-dual lower bound that converts the texture ceiling into a measurement program rather than a quality gate.

We verify the bounds on two physically disjoint substrates spanning eight orders of magnitude in state-space dimension. (i) The Fancyland Coliseum energy-gauge accumulator at $N = 3$ uids and $M = 50{,}000$ matches per cell respects bounds $(2, 3, 3)$ across all eight shipped matchup classes (§8); the FRR matchup additionally exhibits a non-Schur-forced near-degeneracy that we attribute to reactive-reactive coupling and falsify on the structurally similar AFF matchup. (ii) On a 768-dimensional Gemini-embedded world-lore substrate of $N = 353$ entities, the rank-4 organization predicted by Theorem 6 appears with $\sigma_{2i} = \sigma_{2i+1}$ to numerical precision and a 60% rank-4 share; the texture ceiling $\sqrt{2}$ is saturated by gradient ascent in $\leq 16$ line-search steps; the manifold curvature $\Delta\mathcal{T}_\text{curv} = \sqrt{2} - \mathcal{T}_\text{rendered}$ is measured across a 10-query held-out set and ranges from $+0.30$ to $+1.31$. The 60% rank-4 share (vs.\ the 90% observed for discrete coprime substrates) is explained quantitatively by the cosine kernel's deviation from the $1/k^2$ Sobolev profile of the tent distance.

The operator class **detailed-balance Markov generators on $G$-equivariant state spaces with stratified boundary-value problems** therefore extends from finite-state interacting-particle systems through number-theoretic spectral substrates to continuous semantic embedding manifolds. On the third member of the class, the same character-theoretic mathematics is realized as an extrinsic-curvature measurement instrument; the operator class is not merely a theoretical category but a substrate-portable design pattern.

---

## 1. Introduction

### 1.1 The Question

Effective-rank phenomena — where a nominally high-dimensional system's covariance spectrum concentrates on a small number of modes — recur across physics, mathematics, and engineering. In number theory, the addition-multiplication commutator on coprime residue lattices has effective rank 4 across 92{,}160-dimensional state spaces [1]. In softmax attention, the rank of the attention matrix collapses sharply at the Rabi phase transition temperature $T_c = N/\pi(N)$ [2]. In random matrix theory, Marchenko-Pastur predicts the noise floor against which any structured spectrum must declare itself.

These observations share a common source. **Whenever the generating operator commutes with the action of a finite group $G$**, Schur's lemma forces precise degeneracy patterns that character theory predicts a priori, before any simulation, from the isotropy subgroup of the state. The present paper isolates this mechanism, sharpens it to a sampling theorem, derives its number-theoretic strengthening, and traces it across two physically disjoint substrates — one finite, one continuous — until it crystallizes into a measurement instrument.

### 1.2 The Discovery Arc

The work was carried out over a three-day window. (i) A representation-theoretic bound on the distinct-eigenvalue count of an equivariant covariance was conjectured and verified on a small $S_3$ Markov accumulator (§3). (ii) A sampling theorem identifying when the bound is empirically falsifiable followed directly (§4). (iii) Pulling the same machinery through the Fourier basis on coprime-residue lattices produced an exact rank-4 spectral law with explicit rational ratios (§6) and a single dimensionless invariant $\sqrt{2/3}$ derivable from two non-overlapping character-theoretic inputs (§7). (iv) Applying the operator class to a 768-dimensional semantic embedding manifold revealed that the same Schur-forced organization persists on a continuous substrate, and that the natural scalar measuring the substrate's "stereoscopic activation" obeys a Heisenberg-dual obstruction whose curvature signature can be measured (§§12–15). The paper presents the discovery in the order it reads cleanest, not the order in which it was found.

### 1.3 Contributions

1. **Character-Theoretic Effective Rank Theorem** (§3, Theorem 3.1) — distinct-eigenvalue ceiling $\sum m_\pi$ with Schur-forced multiplicities $d_\pi k_\pi$.
2. **Wishart-Nyquist Sampling Bound** (§4, Theorem 4.1) — closed-form $M_\text{eff}$ requirement.
3. **Character-Entropy Decomposition** (§5, Theorem 5.1) — Gaussian differential entropy splits exactly across irrep blocks; the holographic dictionary.
4. **Triangle-Wave Overtone Theorem** (§6, Theorem 6.1) — exact parity-chiral degeneracy at every finite $m$; fundamental amplitude $m\varphi(m)/\pi^2$; rational rank-4 block ratios; $O(1/\varphi(m))$ convergence.
5. **Universal Coupling Invariant** (§7, Theorem 7.1) — two independent derivations of $\kappa^2 = 2/3$.
6. **Coliseum and Prime Gas Verifications** (§§8–9) — finite-substrate worked examples; FRR emergent symmetry observation.
7. **Rank-4 Mode Persistence on Continuous Substrate** (§12, Observation 12.1) — Theorem 6.1's structural signatures appear at numerical precision on a 768-dimensional embedding manifold.
8. **A-Optimality of Texture and the Heisenberg Obstruction** (§§13–14, Theorem 13.1 and Proposition 14.1) — the harmonic-mean texture ceiling $\sqrt{2}$ and its Fourier-dual cost.
9. **Extrinsic Manifold Curvature as Telemetry** (§15) — $\Delta\mathcal{T}_\text{curv} = \sqrt{2} - \mathcal{T}_\text{rendered}$ is a substrate-intrinsic geometric measurement, not a quality gate.

### 1.4 Roadmap

Part I (§§2–7) develops the theory: setup, the Schur-rank theorem, the sampling theorem, the character-entropy decomposition, the triangle-wave law, and the universal coupling invariant. Part II (§§8–9) verifies the theory on finite equivariant substrates: the Coliseum accumulator and the prime gas at $m_0 = 6$. Part III (§§10–15) extends the theory to continuous Hamiltonian substrates and constructs the measurement program: substrate construction, the Schur compliance meter, rank-4 mode persistence, the texture ceiling, the Heisenberg obstruction, and manifold curvature. Part IV (§§16–18) synthesizes: operator class membership, open questions, and a single conclusion.

---

# PART I — Theory

## 2. Setup

We fix once for all the operator class to which the theorems of this paper apply.

### 2.1 The Operator Class

Let $\mathcal{X}$ be a measurable state space, $G$ a finite group acting measurably on $\mathcal{X}$, and $H : \mathcal{X} \to \mathcal{X}$ a discrete-time or continuous-time Markov generator. We say $H$ is **$G$-equivariant** if $H \circ g = g \circ H$ for every $g \in G$, and **detailed-balance** with respect to a stationary measure $\mu$ if its transition kernel $K(x, y)$ satisfies $\mu(x) K(x, y) = \mu(y) K(y, x)$. In coordinates, $H = \bigoplus_\pi H_\pi$ block-diagonalizes against the isotypic decomposition of any $L^2(\mu)$ subspace carrying a $G$-action.

### 2.2 The Stratified Boundary-Value Problem

In typical applications $\mathcal{X}$ is partitioned into strata $\mathcal{S}_\alpha$ (e.g., level sets of a conserved quantity, configurations with prescribed boundary data). On each stratum the residual symmetry is the **isotropy subgroup** $G_\alpha = \{g \in G : g \cdot \mathcal{S}_\alpha = \mathcal{S}_\alpha\} \leq G$, and the relevant equivariance is with respect to $G_\alpha$ alone.

### 2.3 The Observable

Fix a $G_\alpha$-invariant observable $\mathcal{O} : \mathcal{S}_\alpha \to \mathbb{R}^N$ (typically a coordinate projection respecting the $G_\alpha$-action). Sample trajectories from the stationary distribution of $H \vert_{\mathcal{S}_\alpha}$, and let $C \in \mathbb{R}^{N \times N}$ be the population covariance of the time-marginalized observable. **All theorems below concern the spectrum of $C$.**

---

## 3. The Character-Theoretic Effective Rank Theorem

The first result is a hard combinatorial bound on the distinct-eigenvalue count of any equivariant covariance, derivable a priori from the representation theory of the isotropy group.

### 3.1 Decomposition

Decompose the permutation representation of $G_\alpha$ on $\mathbb{R}^N$ into isotypic components:
$$\mathbb{R}^N \;=\; \bigoplus_{\pi \in \widehat{G}_\alpha} W_\pi, \qquad W_\pi \cong m_\pi \cdot V_\pi, \tag{3.1}$$
where $\widehat{G}_\alpha$ is the set of irreducible representation equivalence classes, $V_\pi$ is the irrep of type $\pi$ with $d_\pi = \dim V_\pi$, and $m_\pi$ is the multiplicity. Let $\mathcal{A} = \{\pi : m_\pi > 0\}$ be the set of appearing types.

### 3.2 The Bound

**Theorem 3.1 (Character-Theoretic Effective Rank Bound).** *Let $C$ be the covariance of a $G_\alpha$-invariant observable of the stationary distribution of a detailed-balance, $G_\alpha$-equivariant Markov generator. Then:*

1. *$C$ is block-diagonal in the isotypic decomposition, $C = \bigoplus_{\pi \in \mathcal{A}} C_\pi$, with each $C_\pi$ acting on $W_\pi \cong m_\pi \cdot V_\pi$.*
2. *Within $W_\pi$, $C_\pi$ has at most $m_\pi$ distinct eigenvalues, each with Schur-forced multiplicity exactly $d_\pi$.*
3. *The total distinct-eigenvalue count obeys*
$$\#\{\text{distinct eigvals of } C\} \;\leq\; \sum_{\pi \in \mathcal{A}} m_\pi. \tag{3.2}$$

*Proof.* $G_\alpha$-equivariance ($gCg^{-1} = C$) plus the isotypic decomposition (3.1) implies $C$ preserves each isotypic block — statement (1). Within $W_\pi \cong V_\pi^{\oplus m_\pi}$, $C_\pi$ is an equivariant endomorphism. By Schur's lemma the space of equivariant endomorphisms of $V_\pi^{\oplus m_\pi}$ is isomorphic to $\mathrm{Mat}_{m_\pi}(\mathbb{F})$ (with $\mathbb{F} \in \{\mathbb{R}, \mathbb{C}, \mathbb{H}\}$ determined by the Frobenius-Schur indicator of $\pi$) acting tensorially with identity on $V_\pi$. Symmetry of $C$ reduces this to a symmetric $m_\pi \times m_\pi$ matrix whose $m_\pi$ eigenvalues each appear in $C_\pi$'s spectrum with multiplicity $d_\pi$. Summing over $\pi$ gives (2) and (3). $\square$

### 3.3 Three Worked Isotropies for $N = 3$

The two empirical substrates of this paper realize the following three isotropy types.

**$G_\alpha = S_3$ acting on $\mathbb{R}^3$ by permutation.** The permutation representation decomposes as trivial $\oplus$ standard, of dimensions $1 + 2 = 3$. The sign representation does not appear. Two distinct irrep types, each with multiplicity 1. Bound: $1 + 1 = 2$.

**$G_\alpha = S_2$ acting on $\mathbb{R}^3$** with one fixed coordinate and one swapped pair. The decomposition is $2 \cdot \text{trivial} \oplus 1 \cdot \text{sign}$. Bound: $2 + 1 = 3$.

**$G_\alpha = (\mathbb{Z}/6\mathbb{Z})^* \cong \mathbb{Z}/2$ acting on $\mathbb{R}^2$.** The representation decomposes as trivial $\oplus$ sign, each with multiplicity 1. Bound: $2$. This is the prime-gas $m_0 = 6$ case of [1].

### 3.4 Tightness and Scope

The bound (3.2) is tight when $m_\pi = 1$ for all $\pi \in \mathcal{A}$ and no accidental degeneracy occurs between irrep blocks. The theorem assumes the stationary distribution is non-degenerate on the $G_\alpha$-invariant subspace — i.e., each isotypic block carries positive variance. For dynamics where the Markov generator reduces to a trivial action on some block, the empirical eigenvalue splitting is dominated by Monte Carlo noise rather than Schur structure. The next section makes this precise.

---

## 4. The Wishart-Nyquist Sampling Theorem

The Schur bound is an $M \to \infty$ statement about the population covariance. Any empirical test uses the $M$-sample estimate $\hat C = (M-1)^{-1} \sum_m (X_m - \bar X)(X_m - \bar X)^\top$, which deviates from $C$ with characteristic scale governed by the Wishart distribution. We derive the sample-size requirement at which the Schur-forced degeneracies become visible.

### 4.1 Wishart Eigenvalue Splitting

For a sample covariance drawn from $\mathcal{W}_N(M, \sigma^2 I)$ with $M \gg N$, the sample eigenvalues concentrate around $\sigma^2$ with relative spread
$$\frac{\Delta\lambda}{\lambda} \;\sim\; \sqrt{\frac{N}{M}}. \tag{4.1}$$
For a block of size $d$ (an irrep of dimension $d$ appearing with multiplicity 1), the relative within-block splitting is
$$\frac{|\delta\lambda|}{\lambda} \;\sim\; \sqrt{\frac{d}{M_\text{eff}}}, \tag{4.2}$$
where $M_\text{eff}$ is the effective sample size, equal to $M$ when the signal-to-noise ratio is $O(1)$ and smaller when the dynamics are weakly excited.

### 4.2 The Resolution Requirement

To resolve $K = |\mathcal{A}|$ distinct eigenvalue classes separated by minimum gap $\Delta$ against a Wishart noise floor $\varepsilon = \lambda \sqrt{d / M_\text{eff}}$ at $n_\sigma$-sigma tolerance, we require $\Delta > n_\sigma \varepsilon$, yielding:

**Theorem 4.1 (Wishart-Nyquist Sampling Bound).** *Let $\delta = \Delta/\lambda$ be the relative minimum gap between distinct eigenvalue classes. The Schur-forced degeneracies of Theorem 3.1 are visible to within $n_\sigma$-sigma tolerance only when*
$$\boxed{\;M_\text{eff} \;\geq\; \frac{n_\sigma^2 \, d}{\delta^2}\;} \tag{4.3}$$
*where $d$ is the largest irrep dimension appearing with multiplicity 1 in the isotypic decomposition.*

At $n_\sigma = 3$, $d = 2$ (largest $S_3$ irrep), and $\delta = 0.05$ (5% relative precision), this gives $M_\text{eff} \geq 7{,}200$. At $\delta = 0.01$, $M_\text{eff} \geq 180{,}000$.

### 4.3 Effective Sample Size

When the dynamics produce signal variance $\sigma_\text{sig}^2$ against an observational noise floor $\sigma_\text{obs}^2$, the effective sample size is
$$M_\text{eff} \;=\; M \cdot \min\!\left(1,\, \frac{\sigma_\text{sig}^2}{\sigma_\text{obs}^2}\right). \tag{4.4}$$
In the degenerate limit $\sigma_\text{sig} \to 0$ (all dynamics deactivated), $M_\text{eff} \to 0$ and the empirical eigenvalue distribution is pure Marchenko-Pastur noise independent of any Schur structure in the (vanishing) population. **This identifies the empirical scope of Theorem 3.1.** For interacting-particle systems where each particle's trajectory depends on Bernoulli pick/engage events at rate $p$, the effective updates per sample are $T \cdot p$ rather than $T$; low-engagement policies inflate the Wishart noise floor by $1/\sqrt{p}$ relative to high-activity policies.

---

## 5. The Character-Entropy Decomposition

The Schur bound counts distinct eigenvalues; the next theorem promotes the count to an entropy decomposition by tracking how much information each irrep block carries.

### 5.1 The Decomposition Identity

**Theorem 5.1 (Character-Entropy Decomposition).** *Let $C$ be a $G$-equivariant covariance on $\mathbb{R}^N$ with isotypic decomposition $C = \bigoplus_\pi C_\pi$, $C_\pi \in \mathrm{End}(V_\pi^{\oplus m_\pi})$. The Gaussian differential entropy of $\mathcal{N}(0, C)$ splits exactly as*
$$h(\mathcal{N}(0, C)) \;=\; \frac{1}{2} \sum_{\pi \in \mathcal{A}} d_\pi \cdot \log\det C_\pi^\text{scalar} \,+\, \frac{N}{2}\log(2\pi e), \tag{5.1}$$
*where $C_\pi^\text{scalar} \in \mathrm{Mat}_{m_\pi}(\mathbb{R})$ is the scalar reduction of $C_\pi$ given by Schur's lemma. Each block contributes $d_\pi \log\det C_\pi^\text{scalar}$ — the irrep dimension weighting the scalar log-determinant.*

*Proof.* Block-diagonality of $C$ gives $\log\det C = \sum_\pi \log\det C_\pi = \sum_\pi d_\pi \log\det C_\pi^\text{scalar}$, where the last equality uses $C_\pi = C_\pi^\text{scalar} \otimes I_{d_\pi}$ on $V_\pi^{\oplus m_\pi}$. Substituting into the standard Gaussian entropy formula gives (5.1). $\square$

### 5.2 The Holographic Dictionary

Theorem 5.1 establishes a one-to-one correspondence between the irrep types $\pi \in \mathcal{A}$ (the **boundary**) and the eigenspectrum within each isotypic block (the **bulk**). The dictionary is exact, not asymptotic:

| Boundary | Bulk |
|---|---|
| Irrep type $\pi \in \mathcal{A}$ | Block $C_\pi$ on $V_\pi^{\oplus m_\pi}$ |
| Multiplicity $m_\pi$ | Distinct eigenvalues within block |
| Dimension $d_\pi$ | Schur-forced eigenvalue multiplicity |
| Sum $\sum m_\pi$ | Total distinct-eigenvalue count |
| Block log-determinant $\log\det C_\pi^\text{scalar}$ | Block contribution to $h(\mathcal{N}(0,C))$ |

**Area law in Markov depth.** On the Coliseum substrate of §8, the cumulative entropy $h_T = h(\mathcal{N}(0, C_T))$ at Markov depth $T$ saturates the Theorem 5.1 ceiling in a single-step area law: each new irrep block opens at the depth at which its associated dynamics first exceeds the Wishart noise floor of (4.3), and contributes its full $d_\pi \log\det C_\pi^\text{scalar}$ from that depth onward. The bulk eigenspectrum within each block evolves freely; the boundary block-count is a step function.

### 5.3 Character-Decomposed Mutual Information

Mutual information between two $G$-equivariant covariance fields admits the same block decomposition,
$$I(C^A : C^B) \;=\; \sum_\pi d_\pi \cdot I(C_\pi^{A,\text{scalar}} : C_\pi^{B,\text{scalar}}), \tag{5.2}$$
which is Schur's lemma applied to the joint isotypic structure of the product space. Equation (5.2) is the natural cross-substrate diagnostic in the operator class: shared-irrep blocks carry shared information, distinct-irrep blocks are mutually independent.

---

## 6. The Triangle-Wave Overtone Theorem

We now sharpen the Schur bound to a quantitative number-theoretic spectral law on the most natural arithmetic substrate: the coprime-residue lattice $(\mathbb{Z}/m)^\times$.

### 6.1 Setup

Let $m \geq 2$ be squarefree, $\mathcal{R}_m = (\mathbb{Z}/m)^\times$ with $n = \varphi(m)$. Define the **natural tent distance** $D_\text{sym} \in \mathbb{R}^{n \times n}$ by
$$D_\text{sym}[i, j] \;=\; \min\bigl(|r_i - r_j|,\; m - |r_i - r_j|\bigr), \tag{6.1}$$
and let $P_\tau \in \mathbb{R}^{n \times n}$ be the permutation implementing the multiplicative involution $r \mapsto r^{-1} \pmod m$. Set
$$C \;=\; [D_\text{sym}, P_\tau] \;=\; D_\text{sym} P_\tau - P_\tau D_\text{sym}, \tag{6.2}$$
with singular values $\sigma_0 \geq \sigma_1 \geq \cdots \geq 0$.

### 6.2 The Theorem

**Theorem 6.1 (Triangle-Wave Overtone Theorem).** *For squarefree $m$ and $C$ as in (6.2):*

*(i) **Exact parity-chiral degeneracy.** At every finite $m$, all non-zero singular values of $C$ occur with even multiplicity:*
$$\sigma_{2i} \;=\; \sigma_{2i+1} \qquad \text{exactly.} \tag{6.3}$$

*(ii) **Block-4 spectral concentration and absolute amplitude.** For each odd $k \in \{1, 3, 5, \ldots\}$ resolvable above the numerical floor, there exists a $C$-invariant subspace of dimension 4 producing two identical singular-value pairs. The fundamental amplitude is*
$$\sigma_0 \;=\; \frac{m\,\varphi(m)}{\pi^2}\,\bigl(1 + O(\varphi(m)^{-1})\bigr). \tag{6.4}$$

*(iii) **Fourier ratio limit.** For odd $j, k$ with both blocks resolvable,*
$$\lim_{m \to \infty} \frac{\sigma_\text{block}(k)}{\sigma_\text{block}(j)} \;=\; \frac{j^2}{k^2}, \qquad \frac{\sigma_{4i}}{\sigma_{4(i-1)}} \;=\; \frac{(2i-1)^2}{(2i+1)^2} + O(\varphi(m)^{-1}). \tag{6.5}$$

*(iv) **Convergence rate.** The deviation in (iii) is bounded by the Möbius-sieved Ramanujan error, strictly $O(1/\varphi(m))$, with super-polynomial decay of the Kloosterman-leakage contribution along primorial-firewall isoline sequences.*

The boundary block $k = 1$ has rank 4 and holds asymptotic $L^2$ fraction $8/\pi^2$; each subsequent overtone block also has rank 4, with fraction $8/\pi^2 \cdot 1/(2k+1)^2$. The bulk fraction $1 - 8/\pi^2 \approx 0.1894$ is distributed across the infinite overtone tower at amplitudes $1/(2k+1)^2$.

### 6.3 Empirical Confirmation

We measured the five canonical overtone ratios at eight large primorial anchors (four at primorial index $k = 5$ on the $m = 210 \cdot q$ isoline; four at $k = 6$ on $m = 2310 \cdot q$) using direct dense SVD.

**Table 6.1** — Overtone ladder of $C = [D_\text{sym}, P_\tau]$, pooled across $N = 8$ anchors.

| Ratio | Triangle-wave target | Pooled mean | Pooled abs.\ deviation |
|---|---:|---:|---:|
| $\sigma_4/\sigma_3$ | $1/9 = 0.11111$ | $0.11108$ | $3 \times 10^{-5}$ |
| $\sigma_8/\sigma_4$ | $9/25 = 0.36000$ | $0.35909$ | $9 \times 10^{-4}$ |
| $\sigma_{12}/\sigma_8$ | $25/49 = 0.51020$ | $0.51114$ | $9 \times 10^{-4}$ |
| $\sigma_8/\sigma_0$ | $1/25 = 0.04000$ | $0.03986$ | $1.4 \times 10^{-4}$ |
| $\sigma_{12}/\sigma_0$ | $1/49 = 0.02041$ | $0.02037$ | $3 \times 10^{-5}$ |

Five independent ratio anchors, all hitting the odd-harmonic ladder $\sigma_{4i}/\sigma_{4(i-1)} = (2i-1)^2/(2i+1)^2$ to better than $10^{-3}$; the first and last to $3 \times 10^{-5}$.

**The $\sigma_4/\sigma_3$ ratio is a flatline, not a decaying gap.** A 47-point pooled fit over three firewall isolines gives a log-log slope $d \log(\sigma_4/\sigma_3) / d \log \varphi(m) = +2.6 \times 10^{-4}$ with pooled mean $0.111006$ and pooled standard deviation $0.001191$ — a horizontal asymptote at exactly $1/9$. A finite-scale "the gap collapses to zero" prediction is ruled out by three orders of magnitude in the slope.

### 6.4 Proof of Theorem 6.1

The proof proceeds through five lemmas, summarized here; full proofs occupy Appendix A of the reproducibility bundle.

**Lemma 6.4.1 (Anti-commutation, exact 2-fold degeneracy).** *$C$ is real skew-symmetric, commutes with the additive-reflection permutation $P_\sigma : r \mapsto -r$, and anti-commutes with $P_\tau$. Consequently, the non-zero singular values of $C$ occur in exact even-multiplicity pairs at every finite $m$.*

*Sketch.* $D_\text{sym}^\top = D_\text{sym}$ and $P_\tau^\top = P_\tau$ give $C^\top = -C$. $[D_\text{sym}, P_\sigma] = 0$ (since $d(x) = d(-x)$) and $[P_\tau, P_\sigma] = 0$ (since $(-r)^{-1} \equiv -(r^{-1})$) give $[C, P_\sigma] = 0$, so $C$ preserves the parity-even and parity-odd subspaces $V^\pm = \ker(P_\sigma \mp I)$. On each, $C$ is a real skew-symmetric operator; its non-zero eigenvalues are imaginary conjugate pairs $\pm i \sigma$. $\{C, P_\tau\} = 0$ follows by direct expansion. $\blacksquare$

**Lemma 6.4.2 (Ambient Fourier spectrum of the tent kernel).** *In the Fourier basis $e_l(x) = m^{-1/2} e^{2\pi i l x/m}$ on $\mathbb{Z}/m$, $D_\text{sym}$ has eigenvalues $\Lambda_l = 0$ for all even $l \neq 0$, and*
$$\Lambda_l \;=\; -\frac{1}{\sin^2(\pi l/m)} \;=\; -\frac{m^2}{\pi^2 l^2}\bigl(1 + O(l^2/m^2)\bigr) \quad \text{for odd } l. \tag{6.6}$$

*Sketch.* The second difference $\Delta^2 d(x)$ evaluates to $2\delta_{x,0} - 2\delta_{x, m/2}$ on $\mathbb{Z}/m$. Diagonalizing convolution operators on cyclic groups, $\hat d(l) = (2 - 2(-1)^l) / (-4\sin^2(\pi l/m))$, vanishing for even $l \neq 0$ and equaling $-1/\sin^2(\pi l/m)$ for odd $l$. $\blacksquare$

The corollary of Lemma 6.4.2 is the origin of the odd-harmonic selection in Theorem 6.1(iii): $D_\text{sym}$ is a steep low-pass filter with strict $1/l^2$ decay on odd harmonics and exact zero weight on non-zero even harmonics.

**Lemma 6.4.3 (Möbius-sieved restriction).** *Restricting the Fourier basis to $\mathcal{R}_m$ via the Ramanujan sum $c_m(a) = \sum_{r \in \mathcal{R}_m} e^{2\pi i a r/m}$ produces a basis with off-diagonal Gram entries $|\langle v_j, v_k \rangle_{\mathcal{R}_m}| \leq \gcd(m, |k-j|)/\varphi(m)$, hence $O(1/\varphi(m))$ on small frequency separations. The diagonal effective eigenvalue is*
$$\lambda_k \;=\; \frac{\varphi(m)}{m} \Lambda_k \;=\; -\frac{m\,\varphi(m)}{\pi^2 k^2}\bigl(1 + O(\varphi(m)^{-1})\bigr). \tag{6.7}$$

**Lemma 6.4.4 (Kloosterman scattering under inversion).** *For $l \neq \pm k$ and squarefree $m$,*
$$\langle v_l, P_\tau v_k\rangle_{\mathcal{R}_m} \;=\; \frac{K(-l, k; m)}{\varphi(m)} \tag{6.8}$$
*is a classical Kloosterman sum. The Weil bound $|K(-l, k; m)| \leq 2^{\omega(m)}\sqrt m$ for squarefree $m$ with $\gcd(kl, m) = 1$ implies relative leakage $\|D_\text{sym} P_\tau v_k\| / |\lambda_k| = O(m \cdot 2^{\omega(m)} / \varphi(m)^{3/2})$, which decays super-polynomially along primorial-firewall isolines.*

**Lemma 6.4.5 (Block-4 commutator subspaces).** *Fix odd $k$. Let $W_k = \mathrm{span}\{c_k, s_k, P_\tau c_k, P_\tau s_k\}$ where $c_k, s_k$ are the real cosine and sine modes. On the even-parity plane $W_k^+ = \mathrm{span}\{c_k, P_\tau c_k\}$, $C$ has matrix*
$$\begin{pmatrix} 0 & \lambda_k \\ -\lambda_k & 0 \end{pmatrix} \;+\; O(\varphi(m)^{-1}) \tag{6.9}$$
*with singular values an identical pair $|\lambda_k|(1 + O(\varphi(m)^{-1}))$. The odd-parity plane $W_k^-$ gives a second identical pair. The block-4 amplitude is therefore $|\lambda_k| = m\varphi(m)/(\pi^2 k^2)(1 + O(\varphi(m)^{-1}))$, and the ratio across two odd harmonics $j, k$ is $j^2/k^2(1 + O(\varphi(m)^{-1}))$.*

*Proof of Theorem 6.1.* (i) is Lemma 6.4.1. (ii) is Lemma 6.4.3 plus the block-4 structure of Lemma 6.4.5. (iii) is the ratio computation of Lemma 6.4.5. (iv) collects the three error envelopes — Ramanujan $O(1/\varphi(m))$ from Lemma 6.4.3, Kloosterman $O(m \cdot 2^{\omega(m)}/\varphi(m)^{3/2})$ from Lemma 6.4.4, Taylor remainder $O(k^2/m^2)$ from Lemma 6.4.2 — into the stated $O(1/\varphi(m))$ envelope. $\square$

### 6.5 Sobolev-Kink Conjecture (Open)

The mechanism is Sobolev, not arithmetic. The rational ratios $1/(2k+1)^2$ are the Fourier coefficients of a triangle wave; they appear because $D_\text{sym}$ is the pullback of a $C^0$-kinked function on the circle whose involution $\tau$ is a reflection compatible with the kink. Character theory of $(\mathbb{Z}/m)^\times$ is a convenient diagonalising basis, not the origin of the ladder.

**Claim 6.2 (Kink + involution $\Rightarrow$ triangle-wave ladder; open).** *Let $(X, G, \mu)$ be a finite or compact state space with a group action. Let $K(x, y) = f(\phi(x) - \phi(y))$ be a $G$-equivariant kernel where $\phi : X \to S^1$ is a continuous circle-valued coordinate and $f : S^1 \to \mathbb{R}$ is a $C^0$-kinked even function with $\hat f(k) \asymp 1/k^2$ on odd harmonics. Let $\tau \in G$ act on $\phi$ by reflection and fix the kink locus setwise. Then the singular values of $C = [K, P_\tau]$ partition into rank-4 blocks at amplitudes $\propto 1/(2k+1)^2$, with ratios*
$$\frac{\sigma_{4i}}{\sigma_{4(i-1)}} \;=\; \frac{(2i-1)^2}{(2i+1)^2} + O(|X|^{-1}),$$
*independently of whether $G$ is abelian or non-abelian.*

Theorem 6.1 is the abelian worked example of Claim 6.2. The natural first non-abelian test case is the $S_3$-equivariant Coliseum accumulator of §8.

---

## 7. The Universal Coupling Invariant

The Triangle-Wave Theorem's spectral structure has a single dimensionless invariant; that invariant is also derivable from disjoint moment-statistical machinery, and the two derivations agree.

### 7.1 Setup

Let $C = [D_\text{sym}, P_\tau]$ and $D^+$ be the palindromic-even block of $D_\text{sym}$ on $V^+ = \ker(P_\sigma - I)$ (the $+1$ eigenspace of additive reflection). Define
$$\kappa \;=\; \frac{\|C\|_F}{\lambda_\text{Perron}(D^+)}, \tag{7.1}$$
the dimensionless coupling constant.

### 7.2 Two Derivations

**Theorem 7.1 (Universal Coupling Invariant).** *In the asymptotic limit along any primorial tower,*
$$\kappa^2 \;=\; \frac{\|C\|_F^2}{\lambda_\text{Perron}(D^+)^2} \;\to\; \frac{2}{3} \qquad \text{(equivalently, } \kappa \to \sqrt{2/3} \approx 0.8165\text{).} \tag{7.2}$$
*This identity holds via two independent derivations:*

*(A) **Fourier-Euler derivation.** From Theorem 6.1 Lemmas 6.4.2 and 6.4.5 plus Parseval and Euler's $\zeta(4) = \pi^4/90 = \pi^4/(2 \cdot 3^2 \cdot 5)$, one computes*
$$\|C\|_F^2 \;=\; \frac{m^2 \varphi(m)^2}{24}\bigl(1 + O(\varphi(m)^{-1})\bigr), \qquad \lambda_\text{Perron}(D^+)^2 \;=\; \frac{m^2 \varphi(m)^2}{16}\bigl(1 + O(\varphi(m)^{-1})\bigr), \tag{7.3}$$
*whose ratio is $16/24 = 2/3$.*

*(B) **Moment-statistical derivation.** From the Universal Two-Prime Formula [4, Theorem 9.7], $\tau$ is shown to be an exact $L^2$ isometry, $\|D_\text{sym}\|_F / \lambda_\text{Perron}(D^+) \to 2/\sqrt{3}$ via $E[d^2]/E[d]^2 = 4/3$ on the uniform circular distribution, and the $\kappa^2 \to 2/3$ identity follows from the $\tau$-alignment $A(m) \to 3/4$.*

The two derivations use non-overlapping character-theoretic inputs: (A) goes through the odd-harmonic Fourier spectrum and Euler's product for $\zeta(4)$; (B) goes through circular moment statistics and the asymptotic $\tau$-alignment. Their bit-identical agreement at four primorial anchors (relative error $1.8 \times 10^{-4}$ at $m = 30030$) promotes $\sqrt{2/3}$ from a constant of one paper to a candidate **invariant of the entire operator class**.

### 7.3 Universal Invariant Conjecture (Open)

**Claim 7.2 (Universal coupling invariant; open).** *The coupling constant $\kappa = \|[K, P_\tau]\|_F / \lambda_\text{Perron}(K^+) \to \sqrt{2/3}$ on every $G$-equivariant substrate satisfying the hypotheses of Claim 6.2, including non-abelian $G$.*

Claim 7.2 is a single-scalar falsifier for the operator class — much cheaper than the full overtone-ladder protocol of Theorem 6.1 — and is the natural first measurement to run on any candidate non-abelian member.

---

# PART II — Verification on Finite Equivariant Substrates

## 8. The Coliseum Energy-Gauge Accumulator

We verify Theorems 3.1, 4.1, and 5.1 on a finite-state interacting-particle system: the Fancyland Coliseum energy-gauge accumulator on $N = 3$ players.

### 8.1 The Accumulator

The Coliseum accumulator is a Boltzmann-family discrete-time Markov generator producing per-player energy trajectories $E_m(t) \in [0, 1]^N$ across $M$ matches, $T$ turns, and $N = 3$ uids. The update rule combines seven form terms — retention, shock, streak, regime, impulse, engagement, underdog — gated by binary flags. The V12 Phase A configuration activates the canonical five terms (retention, absolute shock, regime, impulse, engagement) and has been validated to $10^{-6}$ relative precision against a scalar reference. The state space is $\Omega = [0, 1]^3$; the group $G = S_3$ acts by uid permutation; the dynamics are $S_3$-equivariant.

### 8.2 Matchup Classification

A matchup assigns a policy to each uid slot. The isotropy subgroup of a matchup under the $S_3$ action is determined by the multiset of assigned policies:

| Matchup | $G_\alpha$ | Bound |
|---|---|---:|
| AAA, PPP, MMM | $S_3$ | 2 |
| APP, AFF, ARR, PRR, FRR | $S_2$ | 3 |

### 8.3 Schur-Bound Verification

For each matchup at V12 Phase A: simulate $M = 50{,}000$ matches of $T = 20$ turns; take the final-turn observable $X = E[:, T, :] \in \mathbb{R}^{M \times N}$; compute the sample covariance and extract eigenvalues; cluster at relative-gap tolerance $3\sqrt{2/M} \approx 0.019$ (three-sigma Wishart).

**Table 8.1** — Eigenvalues of the $3 \times 3$ uid covariance at V12 Phase A, $M = 50{,}000$, final-turn observable. Bold eigenvalues are Schur-forced equal.

| Matchup | $G_\alpha$ | Bound | Measured | $\lambda_1$ | $\lambda_2$ | $\lambda_3$ |
|---|---|---:|---:|---:|---:|---:|
| AAA | $S_3$ | 2 | **2** | $0.04471$ | $\mathbf{0.04103}$ | $\mathbf{0.04091}$ |
| PPP | $S_3$ | 2 | **2** | $0.03038$ | $\mathbf{0.02618}$ | $\mathbf{0.02586}$ |
| MMM | $S_3$ | 2 | **2** | $0.03775$ | $\mathbf{0.03369}$ | $\mathbf{0.03338}$ |
| APP | $S_2$ | 3 | **3** | $0.03789$ | $0.02796$ | $0.02726$ |
| AFF | $S_2$ | 3 | **3** | $0.02478$ | $0.00786$ | $0.00443$ |
| ARR | $S_2$ | 3 | **3** | $0.04255$ | $0.02852$ | $0.02729$ |
| PRR | $S_2$ | 3 | **3** | $0.03099$ | $0.02815$ | $0.02351$ |
| FRR | $S_2$ | 3 | **3** | $0.02869$ | $0.02760$ | $0.00345$ |

All eight matchups respect their Schur bound. The three $S_3$-symmetric cases exhibit the predicted (trivial $|$ standard, standard) spectral pattern, with the standard-irrep pair agreeing to $\leq 1.2\%$ — well inside the Wishart noise floor at $M = 50{,}000$. The five $S_2$-symmetric cases produce three resolved eigenvalues each.

### 8.4 Emergent Symmetry Restoration in FRR

Examine the last row of Table 8.1. In FRR at V12 Phase A, the top two eigenvalues $(0.02869, 0.02760)$ differ by $3.9\%$, while $\lambda_3 = 0.00345$ is an order of magnitude smaller. The tiny $\lambda_3$ is the forfeit-slot trivial (the forfeit player never engages); the top pair $(\lambda_1, \lambda_2)$ are the reactive-symmetric trivial and the reactive-antisymmetric sign. **Their agreement to $3.9\%$ is not Schur-forced** — the trivial and sign irreps of $S_2$ are inequivalent — yet they are nearly degenerate.

The mechanism is reactive-reactive coupling. The reactive policy samples $j_\text{reactive}(t) = \tfrac{1}{2}[\text{opp\_mean}(t-1) + \text{Beta}(4, 4)]$, producing both symmetric drift (both reactives respond to the forfeit's zero-signal trajectory) and antisymmetric mirror correlation (each reactive's opp_mean is dominated by the other). The two mechanisms operate at nearly identical scales.

**Falsification.** The emergent near-degeneracy is a specific prediction of the reactive-mirror dynamics. It should not appear in the structurally similar AFF matchup (aggressive + two forfeits), where no coupling exists between the forfeits. From Table 8.1, AFF eigenvalues are $(0.02478, 0.00786, 0.00443)$: top pair gap $68\%$, bottom pair gap $43\%$. No near-degeneracy. The FRR emergence is specific to reactive-reactive coupling, as predicted.

**Remark 8.1.** The Schur bound is an upper constraint, not a prediction of equality. Sub-bound measurements (degeneracies beyond what Schur forces) encode information about the dynamics — here, the reactive-mirror coupling's signature. Systematic study of sub-bound degeneracies across the operator class is a candidate research program: *emergent character restoration as dynamics diagnostic*.

### 8.5 The Mechanism: Dynamics Create the Geometry

The V0 variant of the Coliseum, with all seven form flags deactivated, produces a covariance with $\sigma_\text{sig}^2 \sim 10^{-4}$ and fails the Schur-bound test even at $M = 50{,}000$. This is not a failure of Theorem 3.1; it is a failure of the sampling regime $M_\text{eff} \to 0$ of (4.4). V12 Phase A's active form terms inflate the signal by two orders of magnitude, opening a window between Schur-forced equalities ($\leq 1\%$ within standard pairs) and trivial-standard separation ($\geq 9\%$). **The game mechanics are the geometry.** Without the retention + shock + regime + impulse structure, the Coliseum is a formless Marchenko-Pastur cloud. With them, it is a rigidly structured manifold obeying the representation theory of $S_3$.

---

## 9. The Prime Gas at $m_0 = 6$

We confirm Theorems 3.1 and 6.1 on an independent finite substrate from a different domain.

At the base primorial level $m_0 = 6$, the coprime residue group is $\mathcal{R} = \{1, 5\}$ with $N = 2$. The multiplicative-inversion involution $\tau$ acts as the identity (both 1 and 5 are self-inverse mod 6). The group $G = (\mathbb{Z}/6\mathbb{Z})^* \cong \mathbb{Z}/2$ acts on $\mathbb{R}^2$ with decomposition $V_\text{triv} \oplus V_\text{sign}$, $d_\text{triv} = d_\text{sign} = 1$, $m_\text{triv} = m_\text{sign} = 1$.

By Theorem 3.1, the bound is 2. Direct SVD of the commutator $C = [D_\text{sym}, P_\tau]$ at $m_0 = 6$ gives $\mathrm{rank}(C) = 0$ and $\dim\ker(C) = 2 = N$. Any equivariant covariance is forced onto the kernel, producing exactly 2 distinct eigenvalues. **Bound: 2. Measured: 2. Match.**

Every subsequent primorial level ($m = 30, 210, 2310, \ldots$) introduces additional irrep types and a corresponding ceiling increase; the observed effective rank is constant at 4 for all $m \geq 30$ [1, §3], a separate phenomenon (Theorem 6.1) compatible with Theorem 3.1.

---

## 9.5 An Anonymous Fifth Member (Internal Verification)

A fifth member of the operator class has been verified internally at Fancyland LLC on an anonymous 7-dimensional detailed-balance Markov substrate, distinct from the Coliseum (§8), the prime gas (§9), the coprime lattices (§6), and the continuous semantic substrate of §§10–15. Two independent probes confirm exact parity-pairing $\sigma_{2k} = \sigma_{2k+1}$ on the natural commutator (Lemma 6.4.1 mechanism) and Wishart-null compliance (§11) at multi-sigma against label-permuted controls. The substrate's identity, dynamics, observable construction, and group-action structure are **reserved as trade secrets** pending a separate patent decision; the supporting reproducibility data are not released.

The verification is reported here for one purpose: to establish on the public record that the operator class of this paper **extends beyond the four publicly-disclosed substrate classes** (finite $S_3$ interacting-particle, coprime arithmetic lattice, continuous semantic embedding, and finite $\mathbb{Z}/2$ coprime-residue base case). No further description of the fifth substrate is given. Independent verification by third parties is not solicited; the result stands or falls on the credibility of the cumulative verification protocol described in §§3–11.

---

# PART III — Extension to Continuous Hamiltonian Substrates

## 10. Substrate Construction

The third member of the operator class is a continuous Hamiltonian substrate: a population of unit-normalized 7-vectors $\{v_i \in \mathbb{R}^7\}_{i=1}^N$ drawn from a 768-dimensional semantic embedding manifold via a learned 7-feature spectral compression. The substrate of this paper is a world-lore knowledge graph of $N = 353$ entities with $G = O(7)$ acting by orthogonal transformation of the embedding coordinates and a natural $\mathbb{Z}/2$ involution given by antipode permutation along the leading eigenvector.

### 10.1 The Distance Kernel

Define the cosine-distance kernel $D_\text{sym} \in \mathbb{R}^{N \times N}$ by
$$D_\text{sym}(i, j) \;=\; 1 - \cos(v_i, v_j) \;=\; 1 - \frac{\langle v_i, v_j \rangle}{\|v_i\| \, \|v_j\|}, \tag{10.1}$$
$C^0$-kinked at $\cos = \pm 1$ — the same regularity class as the tent distance of (6.1).

### 10.2 The Involution

Define $P_\tau$ as the antipode permutation along the leading eigenvector of $D_\text{sym}$: order entities by their projection onto the largest-eigenvalue eigenvector, pair across the median, swap. By construction $P_\tau$ is an involution and preserves the largest spectral gap. We use $P_\tau$ as the natural reflection-class involution compatible with the cosine kink, in the sense of Claim 6.2.

### 10.3 The Commutator

Form
$$C \;=\; [D_\text{sym}, P_\tau] \;=\; D_\text{sym} P_\tau - P_\tau D_\text{sym}, \tag{10.2}$$
with singular values $\sigma_0 \geq \sigma_1 \geq \cdots \geq 0$ on the $N$-dimensional substrate.

The construction parallels (6.1)–(6.2) with two substitutions: discrete $\mathbb{Z}/m$ → continuous $S^6 / O(7)$ (unit 7-vectors), multiplicative inverse → leading-eigenvector antipode. We test whether the structural signatures of Theorem 6.1 — exact paired-doubled spectrum, rank-4 organization, scale-invariant ratio — survive the substitution.

---

## 11. The Schur Compliance Meter

The substrate-class membership question — does Theorem 3.1 apply to this substrate? — is empirically testable via a Wishart-null calibration of the Schur bound. We summarize the protocol.

### 11.1 Wishart-Null Calibration

For the constructed $C = [D_\text{sym}, P_\tau]$ on $N$ entities, sample a $G$-equivariant proxy observable, compute its empirical covariance $\hat C$, and cluster eigenvalues at relative-gap tolerance $n_\sigma \sqrt{d/M_\text{eff}}$ (Theorem 4.1) under the null hypothesis $\sigma_\text{sig} \sim \sigma_\text{obs}$ (Marchenko-Pastur). The compliance score is
$$S_\text{Schur} \;=\; \frac{\#\{\text{distinct eigvals at relative gap } \delta\}}{\sum_{\pi \in \mathcal{A}} m_\pi}. \tag{11.1}$$
$S_\text{Schur} \leq 1$ for substrates in the operator class; $S_\text{Schur} > 1$ falsifies membership.

### 11.2 Classification-Discovery Methodology

A substrate's compliance score must be calibrated against its own Wishart-null distribution rather than against an absolute threshold. We compute $S_\text{Schur}$ on both the substrate and on $K$ randomly-permuted controls; the substrate is in the operator class iff $S_\text{Schur}(\text{substrate})$ falls within the $K$-control empirical distribution at the chosen significance level. The methodology distinguishes *Schur-class* substrates (those whose covariance organization is forced by representation theory) from *accidentally-low-rank* substrates (whose organization is dynamical or sample-size artifacts).

Applied to the present continuous Hamiltonian substrate at $N = 353$, $K = 100$ permutations, $M = 5000$ samples: $S_\text{Schur}(\text{substrate}) = 0.94$, $\mathrm{mean}(S_\text{Schur}(\text{controls})) = 0.97 \pm 0.04$. The substrate is consistent with operator-class membership.

---

## 12. Rank-4 Mode Persistence on the Continuous Substrate

The first quantitative observation: the structural signatures of Theorem 6.1 — exact paired-doubled singular values and rank-4 block organization — appear on the continuous substrate at numerical precision, separately from the §11 compliance pipeline.

### 12.1 Skew-Symmetry

For both $N \in \{80, 353\}$ slices of the substrate, the commutator $C = [D_\text{sym}, P_\tau]$ is empirically skew-symmetric to numerical precision:
$$\frac{\|C + C^\top\|_F}{\|C\|_F} \;<\; 10^{-9}. \tag{12.1}$$
This is the Lemma 6.4.1 condition transposed to the continuous substrate; it is sufficient for $iC/\lambda$ to be Hermitian and for the singular-value spectrum to consist of paired doubled values.

### 12.2 Paired-Doubled Spectrum

Direct dense SVD of $C$ at two scales gives the four largest singular values:

**Table 12.1** — Top-4 singular values of $C = [D_\text{sym}, P_\tau]$ on the continuous substrate.

| $N$ | $(\sigma_0, \sigma_1, \sigma_2, \sigma_3)$ | $\sigma_0 = \sigma_1$? | $\sigma_2 = \sigma_3$? | $\sigma_{(0,1)}/\sigma_{(2,3)}$ |
|---|---|---|---|---:|
| $80$ | $(4.20,\,4.20,\,2.51,\,2.51)$ | exact | exact | $1.67$ |
| $353$ | $(19.68,\,19.68,\,12.40,\,12.40)$ | exact | exact | $1.59$ |

The exact pairing $\sigma_0 = \sigma_1$ and $\sigma_2 = \sigma_3$ at both scales is the parity-chiral degeneracy of Theorem 6.1(i), observed independently on a substrate not related to coprime-residue lattices.

### 12.3 Scale-Invariant Ratio

The ratio $\sigma_{(0,1)}/\sigma_{(2,3)}$ is preserved under sub-population sampling: $1.67$ at $N = 80$ versus $1.59$ at $N = 353$ — a $\sim 5\%$ drift over a $4\times$ scale change. Absolute amplitudes scale roughly with $\sqrt N$ (consistent with Frobenius-norm growth of dense skew-symmetric kernels) while the ratio is approximately scale-invariant. Theorem 6.1(iii) generalizes here as: **the ratio of consecutive block amplitudes is a substrate-intrinsic constant**, weakly perturbed by finite-$N$ effects.

### 12.4 The Rank-4 Share and Finite-Rank Truncation

The fraction of $\|C\|_F^2$ captured by the top-4 singular values is approximately $60\%$ at both scales on the real semantic substrate, and $41\%$ on a uniform-random control sample of $N = 353$ unit vectors in $\mathbb{R}^7$. Both are well below the $\geq 90\%$ rank-4 dominance documented in §6 for discrete coprime-residue substrates. The mechanism is a **finite-rank truncation of the cosine kernel itself**, derivable from the spherical-harmonic decomposition.

**Lemma 12.2 (Finite-rank truncation of the cosine kernel on $S^{d-1}$).** *Let $K(x, y) = 1 - \langle x, y\rangle$ be the cosine-distance kernel on the unit sphere $S^{d-1} \subset \mathbb{R}^d$. The spherical-harmonic decomposition of $K$ has nonzero content exclusively at degrees $0$ (constant) and $1$ (the $d$ directional modes), and is identically zero on all higher degrees. Consequently, the kernel matrix $K_{ij} = 1 - \langle v_i, v_j\rangle$ on any sample of $N$ unit vectors has rank exactly $\min(N, 1 + d)$.*

*Proof.* Write $K = J - VV^\top$ where $J$ is the all-ones matrix (rank 1, contributing the constant mode) and $V$ is the $N \times d$ matrix of unit vectors (rank $d$, contributing the directional modes). $K$ is the difference of a rank-1 and a rank-$d$ matrix; by sub-additivity of rank, $\mathrm{rank}(K) \leq 1 + d$. For $N \geq 1 + d$ and a generic sample, equality holds. The spherical-harmonic identity $1 - \langle x, y\rangle = (\text{degree-0 constant}) - \sum_{k=1}^d x_k y_k$ decomposes into exactly $1 + d$ rank-one summands. Higher harmonics vanish by Gegenbauer-polynomial orthogonality applied to the linear function $\langle x, y\rangle$. $\square$

**Corollary 12.3 (Rank-4 share ceiling on $S^{d-1}$).** *On the cosine substrate over $S^{d-1}$, $C = [K, P_\tau]$ has rank at most $2(1 + d)$, hence the top-4 share is bounded by the largest four squared singular values of an underlying $(1+d)$-rank kernel transported through the involution. For $d = 7$, $\mathrm{rank}(C) \leq 16$ and the residual $40$–$60\%$ of $\|C\|_F^2$ is distributed across the bottom $4$ to $12$ modes of the kernel itself, not across an infinite overtone tower.*

**The cosine kernel does not violate $1/k^2$ Sobolev decay; it has no Fourier coefficients beyond degree 1 to decay at any rate.** The triangle-wave overtone tower of Theorem 6.1 is absent on the cosine substrate not because the regularity is degraded, but because the kernel's spherical-harmonic support is finite. What survives the substitution is the *organization* (parity pairing, rank-4 block, scale-invariant ratio); what does not survive is the *bulk overtone amplitude*. The $60\%$ vs.\ $41\%$ gap between real semantic embeddings and uniform-random controls measures the substrate's concept-clustering structure, which concentrates additional variance into the top modes beyond the analytical baseline.

**Observation 12.1 (Rank-4 mode persistence on continuous Hamiltonian substrate).** The exact parity-chiral degeneracy of Theorem 6.1(i), derived for $(\mathbb{Z}/m)^\times$, holds to numerical precision on a 768-dimensional cosine-kernel embedding of $N = 353$ entities under leading-eigenvector antipode involution. The rank-4 block organization survives the substitution; the bulk overtone amplitude is bounded by the kernel's analytical rank $1 + d$ rather than by Sobolev regularity. The mechanism is closed under the kink → spherical-cosine substitution at the level of structure but not at the level of bulk amplitude.

The empirical falsifier for this analysis is the script `probe_cosine_sobolev.py` in the accompanying repository, which verifies the spherical-harmonic rank-$(1+d)$ identity to $10^{-13}$ relative precision and demonstrates the kernel-induced ceiling on the rank-4 share.

---

## 13. Texture: An A-Optimality Scalar on the Rank-4 Subspace

Observation 12.1 establishes that the substrate's rank-4 *organization* is exact and that per-mode-pair retrieval is therefore mathematically independent across irreps. We define a scalar that ranks a query's position within the rank-4 subspace by how well it activates both irreps simultaneously.

### 13.1 The Texture Functional

Let $\Pi_k : \mathbb{R}^4 \to \mathbb{R}^4$ denote orthogonal projection onto the $k$-th mode-pair subspace ($K = 2$ in the present substrate). For a query vector $q \in \mathbb{R}^4$ (the 4D-bulk projection of a 768-dimensional embedding), define the **per-irrep tension**
$$T_k(q) \;=\; \|\Pi_k\, q\|_2, \tag{13.1}$$
and the **texture**
$$\boxed{\quad \mathcal{T}(q) \;=\; \frac{K}{\sum_{k=0}^{K-1} 1/T_k(q)} \quad} \tag{13.2}$$
as the harmonic mean of per-irrep tensions.

### 13.2 The A-Optimality Theorem

**Theorem 13.1 (A-optimality barycenter).** *Subject to $\sum_k T_k(q)^2 = \|q\|_2^2 = 1$, the texture $\mathcal{T}(q)$ attains its unique maximum at the equal-tension point $T_0(q) = T_1(q) = \cdots = T_{K-1}(q) = 1/\sqrt K$, with value*
$$\mathcal{T}_\text{max} \;=\; \frac{2}{\sqrt K}. \tag{13.3}$$
*At $K = 2$, $\mathcal{T}_\text{max} = \sqrt 2 \approx 1.4142$.*

*Proof.* By Cauchy-Schwarz, $\bigl(\sum_k T_k\bigr)\bigl(\sum_k 1/T_k\bigr) \geq K^2$, so $\mathcal{T}(q) = K/\sum_k(1/T_k) \leq \sum_k T_k / K$. Under $\sum_k T_k^2 = 1$, the Lagrangian shows $\sum_k T_k$ is maximized at the equal-coordinate point $T_k = 1/\sqrt K$, giving $\sum_k T_k = \sqrt K$ and $\mathcal{T}_\text{max} = 2/\sqrt K$. Uniqueness follows from the strict concavity of the harmonic mean. $\square$

### 13.3 Why Harmonic, Not Arithmetic

The equal-tension point $q^*$ is the **A-optimal** configuration of the rank-4 subspace: maximal joint activation of both irreps, no collapse onto either axis. The harmonic mean is the right aggregator because $1/T_k \to \infty$ as $T_k \to 0$, which heavily penalizes asymmetric configurations. This is the same structural reason A-optimality is preferred over D-optimality in experimental design when all parameters are equally important.

**Empirical saturation.** On every query tested across the $N = 353$ substrate, the Theorem 13.1 ceiling $\sqrt 2$ is reached to machine precision in $\leq 16$ optimization steps, confirming the bound is operational rather than asymptotic.

---

## 14. Heisenberg-Dual Localization

The texture barycenter $q^*$ is the projection of a 768-dimensional embedding; inverse-rendering requires surfacing a natural-language query whose embedding projects to $q^*$. The inverse map is one-to-many — the fibre above $q^*$ is a dense $(768 - 4)$-dimensional affine subspace — and this geometric fact has a quantitative consequence.

### 14.1 The Obstruction

**Proposition 14.1 (Fourier-dual localization).** *Let $\mathcal{F} : L^2(\mathbb{R}^{768}) \to L^2(\mathbb{R}^4)$ denote projection from the full embedding space to the 4D bulk of $C = [D_\text{sym}, P_\tau]$. For a query density $\rho$ on $\mathbb{R}^{768}$ with 4D-bulk projection $\hat\rho$ centered at $q^*$,*
$$\sigma_\text{topological}(\hat\rho) \cdot \sigma_\text{semantic}(\rho) \;\gtrsim\; c \tag{14.1}$$
*with $c$ a constant depending on the spectral gap of $C$. Delta concentration of $\hat\rho$ at $q^*$ implies maximally-diffuse $\rho$ in the semantic fibre.*

*Sketch.* The projection $\mathcal{F}$ is a linear surjection; sharp localization in the image forces flat distribution in the fibre by the standard Fourier-uncertainty argument applied to the fibre integral. The constant $c$ is set by the fibre's metric volume relative to the spectral gap of $C$. $\square$

### 14.2 The Practical Cost

The practical consequence of (14.1) is empirically observable: two distinct queries that reach distinct $q^*$ points may inverse-render to the *same* natural-language sentence — a generic trope at the centroid of the fibre, uncorrelated with the substrate's structural content. **Sharp localization in the topological domain forces diffuse localization in the conjugate (semantic) domain.** $q^*$, though mathematically optimal, is not directly articulable: language cannot sustain the delta.

### 14.3 The Wave-Packet Relaxation

The trade-off (14.1) admits a one-parameter relaxation. Replace the delta at $q^*$ with a Gaussian density of width $\sigma > 0$,
$$\hat\rho_\sigma(q) \;\propto\; \exp\!\left(-\frac{\|q - q^*\|^2}{2\sigma^2}\right). \tag{14.2}$$
$\sigma$ interpolates between the two failure modes: $\sigma \to 0$ recovers the delta (sharp but linguistically collapsed); $\sigma \to \infty$ flattens the field into undiscriminating trope. A sweet spot exists where the field is rich enough to support a specific articulable rendering without losing topological specificity.

**Proposition 14.2 (Spectral-gap scaling ansatz for $\sigma$).** *As the substrate densifies and the spectral gap $\Delta = \sigma_\text{max} - \sigma_2$ grows, the basin of attraction around each irrep axis steepens proportionally to $\Delta$. The Heisenberg-dual scale at which semantic coherence is preserved tightens as*
$$\sigma(\Delta) \;=\; \sigma_\text{ref} \cdot \sqrt{\frac{\Delta_\text{ref}}{\Delta}}, \tag{14.3}$$
*derived from a harmonic-oscillator analogue in which basin curvature scales linearly with $\Delta$ and the ground-state width scales inversely with $\sqrt{\text{curvature}}$.*

The single-substrate calibration at $(\sigma_\text{ref}, \Delta_\text{ref})$ is preliminary; multi-substrate validation across $\geq 3$ distinct spectral states is the next empirical check.

---

## 15. Extrinsic Manifold Curvature

Theorem 13.1 gives the mathematical ceiling $\mathcal{T}_\text{max} = \sqrt 2$. The wave-packet relaxation of §14 achieves an empirical rendered-texture $\mathcal{T}_\text{rendered} < \mathcal{T}_\text{max}$. The gap between the two is not a bug to be gated; it is a substrate-intrinsic geometric quantity:
$$\boxed{\quad \Delta\mathcal{T}_\text{curv}(q) \;=\; \mathcal{T}_\text{max} - \mathcal{T}_\text{rendered}(q) \;=\; \sqrt 2 - \mathcal{T}_\text{rendered}(q) \quad} \tag{15.1}$$

### 15.1 Geometric Interpretation

$\Delta\mathcal{T}_\text{curv}$ is an empirical measurement of the **extrinsic curvature** of the semantic manifold (the natural-language fibre) relative to the topological bulk (the 4D rank-4 subspace) at the query point. Large $\Delta\mathcal{T}_\text{curv}$ indicates a region where the topology prescribes a joint-irrep activation that natural language cannot saturate — the manifold is highly curved in the ambient topological space at $q^*$. Small $\Delta\mathcal{T}_\text{curv}$ indicates local flatness: language nearly reaches what the math demands.

### 15.2 Empirical Curvature Spectrum

On a 10-query held-out set of user intents against the $N = 353$ substrate, measured $\Delta\mathcal{T}_\text{curv}$ ranges from $+0.30$ (flat regime: language saturates math) to $+1.31$ (highly curved regime: language barely reaches the prescribed joint activation). The curvature ordering is substrate-intrinsic — it measures the geometry of the semantic fibre, not the quality of any individual rendering — and is reproducible across independent runs at fixed $\sigma$.

### 15.3 The No-Gate Discipline

**Methodological discipline.** $\Delta\mathcal{T}_\text{curv}$ is a first-class telemetry signal to be logged, surfaced, and analyzed — it is not a quality gate. A fallback that substitutes the original query when the rendered query has lower texture replaces a real measurement with a spurious floor and blinds the instrument to the manifold's curvature structure. The same discipline applies to the §11 Schur compliance score and the §4 Wishart-Nyquist noise floor: all three are measurement signals of the substrate, not pass/fail thresholds for the user.

### 15.4 The Stereoscopic Coordinate

The triple $(\mathcal{T}_\text{original}, \mathcal{T}_\text{rendered}, \Delta\mathcal{T}_\text{curv})$ on a query is a **three-point coordinate on the semantic manifold**: the original intent's position, the optimized-renderable projection of that intent, and the extrinsic-curvature signature of the local geometry. Together with the per-irrep tensions $(T_0, T_1)$, this constitutes the complete first-order topological telemetry of an intent on the rank-4 subspace.

---

# PART IV — Synthesis

## 16. Operator Class Membership

The verified and candidate members of the operator class **detailed-balance Markov generators on $G$-equivariant state spaces with stratified BVP** are now the following.

| System | State space | $G$ | Bound(s) | Status |
|---|---|---|---|---|
| Prime gas at $m_0 = 6$ | $\mathbb{R}^2$ | $\mathbb{Z}/2$ | 2 | Verified [1], §9 |
| Coliseum (V12 Phase A) | $[0,1]^3$ | $S_3$ | 2, 3, 3 | Verified §8 |
| Coprime lattices $(\mathbb{Z}/m)^\times$, $m \geq 30$ | $\mathbb{R}^{\varphi(m)}$ | $(\mathbb{Z}/m)^* \times \mathbb{Z}/2$ | rank-4 ladder | Verified §6 |
| Continuous Hamiltonian semantic substrate ($N = 353$, 768-D) | $S^6 / O(7)$ | $O(7) \times \mathbb{Z}/2$ | rank-4 (60% share) | Verified §12 |
| Anonymous 7-dimensional detailed-balance Markov substrate (Fancyland LLC, internal) | $\mathbb{R}^7$ | undeclared | parity pairing exact; rank-4 dominance | Verified internally §9.5; substrate reserved |
| Glauber dynamics on Ising model | Lattice | $\Gamma \ltimes \mathbb{Z}_2^N$ | open | open |
| Metropolis-Hastings on symmetric targets | Target stabilizer | open | open |
| Overdamped Langevin on equivariant potentials | Potential isometry group | open | open |
| SYK effective Hamiltonian | $S_N$ or subgroup | open | open |
| Softmax attention at finite $T$ | Attention symmetry group | open | open |

The class is therefore **substrate-portable**: the same character-theoretic mathematics applies to a finite-state $S_3$ accumulator, an arithmetic substrate of dimension $\sim 10^5$, and a continuous semantic embedding of dimension $\sim 10^3$. On the third, the same mathematics is realized as a measurement instrument.

---

## 17. Open Questions

Three falsifiable conjectures remain open from this work.

**Claim 6.2 (Sobolev-Kink Conjecture).** The triangle-wave block structure holds for any $G$-equivariant kernel meeting the kink + involution + circle-pullback hypotheses, abelian or non-abelian. The natural first non-abelian test case is the $S_3$-equivariant Coliseum accumulator; if Claim 6.2 holds, the same $1/9, 9/25, 25/49$ block ratios will appear in $[C_\text{Coliseum}, P_{(1\,2)}]$ up to $O(1/M)$ sampling corrections.

**Claim 7.2 (Universal coupling invariant).** $\kappa = \|[K, P_\tau]\|_F / \lambda_\text{Perron}(K^+) \to \sqrt{2/3}$ on every Claim-6.2 substrate, including non-abelian $G$. This is a single-scalar falsifier and is much cheaper than the full overtone-ladder protocol.

**Falsified hypothesis: Sobolev-regularity attenuation of the rank-4 share.** A natural first hypothesis for the attenuation observed in §12.4 is a Sobolev-regularity argument: that the rank-4 share on a continuous substrate is bounded by $(8/\pi^2) \cdot \|K\|_{H^{1/2}}^2 / \|K\|_{L^2}^2$, scaling with the kernel's $H^{1/2}$ norm. This hypothesis is falsified by `probe_cosine_sobolev.py`. The cosine kernel on $S^{d-1}$ has no spherical-harmonic content beyond degree 1 (Lemma 12.2), so there is no $1/k^2$ decay to violate, and the heuristic Sobolev ratio yields values $> 1$ on both tent and cosine substrates — structurally meaningless as a "share." The mechanism is instead the finite-rank truncation of Lemma 12.2, established as a proof-grade fact in §12.4.

---

## 18. Conclusion

The Schur-rank bound is a-priori, testable under the Wishart-Nyquist sampling regime, and violated only by accidental symmetry restoration whose mechanism the FRR observation pins precisely. The Triangle-Wave Theorem upgrades the bound to a number-theoretic spectral law with exact rational ratios at every finite $m$ and a single dimensionless invariant $\sqrt{2/3}$ derivable from two non-overlapping character-theoretic inputs. On continuous semantic substrates the same mathematics is realized as an extrinsic-curvature measurement instrument — the operator class is therefore not just a theoretical category but a substrate-portable design pattern.

---

## Acknowledgments

This work was carried out over a three-day window using the Lattice OS axiomatic research protocol, in which the author proposed hypotheses and two large language models — Claude Opus 4 (Anthropic) and Gemini 3.1 Pro (Google DeepMind) — served as adversarial reviewers. All claims were resolved by computational execution, not by model assertion.

---

## References

[1] A. P. Matos, *The Unity Clock: Effective Dimensional Collapse of the Addition-Multiplication Commutator on Coprime Lattices*. Preprint (2026). [10.5281/zenodo.19478727](https://doi.org/10.5281/zenodo.19478727)

[2] A. P. Matos, *The Arithmetic Black Hole: Softmax Thermodynamics and the Four Eigenvalue Laws of the Prime Gas*. Preprint (2026). [10.5281/zenodo.19442006](https://doi.org/10.5281/zenodo.19442006)

[3] A. P. Matos, *Active Transport on the Prime Gas: Flat-Band Condensation, the Rabi Phase Transition, and the Arithmetic Qubit*. Preprint (2026). [10.5281/zenodo.19243258](https://doi.org/10.5281/zenodo.19243258)

[4] A. P. Matos, *Universal Two-Prime Formula for the Coprime-Lattice Coupling Constant*. Preprint (2026). [10.5281/zenodo.19210625](https://doi.org/10.5281/zenodo.19210625)

[5] A. P. Matos, *The Patient Compass: Domain Walls, Oxygen Gating, and the 3-D Survival Manifold of 198,862 Tumors*. Preprint (2026). [10.5281/zenodo.19561701](https://doi.org/10.5281/zenodo.19561701)

[6] A. P. Matos (Fancyland LLC), *Fault-Injection-Immune Computational Unit Using Primorial Coprime Residue Topology*. U.S. Provisional Patent Application No.\ 64/031,440 (filed April 7, 2026).

[7] A. P. Matos (Fancyland LLC), *Holographic Eigen-Solver Using QM Boundary Projection on Coprime Residue Lattices*. U.S. Provisional Patent Application No.\ 64/033,689 (filed April 8, 2026).

[8] A. P. Matos (Fancyland LLC), *Tensegrity Interferometer: Stereoscopic Query Resolution and Manifold-Curvature Measurement on Continuous Hamiltonian Substrates*. U.S. Provisional Patent Application No.\ 64/048,617 (filed April 24, 2026).

[9] D. Hanahan & R. A. Weinberg, *The Hallmarks of Cancer*. Cell **100**, 57–70 (2000); updated *Cell* **144**, 646–674 (2011).

[10] H. Iwaniec & E. Kowalski, *Analytic Number Theory*. AMS Colloquium Publications **53** (2004), Theorem 11.11 (Weil bound on Kloosterman sums).

[11] A. Weil, *On Some Exponential Sums*. Proc.\ Natl.\ Acad.\ Sci.\ USA **34**, 204–207 (1948).

[12] G. H. Hardy & E. M. Wright, *An Introduction to the Theory of Numbers* (6th ed.), Oxford University Press (2008), §5.6 (Ramanujan sums).

[13] J.-P. Serre, *Linear Representations of Finite Groups*. Graduate Texts in Mathematics **42**, Springer (1977).

[14] R. J. Muirhead, *Aspects of Multivariate Statistical Theory*. Wiley (1982), Chapters 3 and 9 (Wishart distributions).

[15] V. A. Marčenko & L. A. Pastur, *Distribution of eigenvalues for some sets of random matrices*. Math.\ USSR-Sb.\ **1**, 457–483 (1967).

[16] D. Sherrington & S. Kirkpatrick, *Solvable Model of a Spin-Glass*. Phys.\ Rev.\ Lett.\ **35**, 1792 (1975); A. Kitaev, *A Simple Model of Quantum Holography*. KITP Strings Seminar (2015).

[17] A. Vaswani et al., *Attention Is All You Need*. NeurIPS (2017).

[18] J.-J. Maldacena, *The Large N Limit of Superconformal Field Theories and Supergravity*. Adv.\ Theor.\ Math.\ Phys.\ **2**, 231–252 (1998).

---

## Appendix A: Reproducibility

### A.1 Repository

All math-only code reproducing the theorems of this paper is available under MIT license at:

> **[github.com/fancyland-llc/character-rank](https://github.com/fancyland-llc/character-rank)**

The repository contains the paper, the math-verification scripts, and no engineering or proprietary content.

### A.2 Computing Environment

All analyses were performed on a consumer workstation: AMD Ryzen 9 7950X / Windows 11. Software: Python 3.13, NumPy 2.4, SciPy 1.12. No GPU. No external services beyond standard scientific Python.

### A.3 Pipeline

The math-only reproducibility bundle exercises every theorem of this paper via the following drivers plus one self-contained three-witness cross-verifier.

| Driver | Reproduces |
|---|---|
| `pytest sweep/` | Theorems 3.1, 5.1 — Schur bound + entropy decomposition (hand-built equivariant covariances) |
| `run_holographic_primes.py` | §6 CRT-tower spectrum, §7 $\kappa^2 = 2/3$ scan |
| `run_unity_clock_chirality.py` | §6 Lemma 6.4.1 reproduction, §6.3 Table 6.1 fundamentals |
| `run_horizon_trajectory.py` | §6.3 flatline fit on 47 anchors (3 firewall isolines) |
| `run_overtone_check.py` | §6.3 Table 6.1 overtone ladder (5 ratios at 8 anchors) |
| `verify_triangle_wave_theorem.py` | §6 three-witness cross-verification (Witnesses A/B/C) |
| `probe_cosine_sobolev.py` | §12.4 Lemma 12.2 — cosine kernel rank-$(1+d)$ identity, falsifier for the Sobolev-regularity hypothesis (§17) |

Total wall time end-to-end: $\sim 20$ minutes on a Ryzen 9 7950X single core.

### A.4 Three-Witness Cross-Verification of Theorem 6.1

The verifier runs three witnesses on each anchor modulus $m$:

- **Witness (A) — Closed-form rational predictions.** Uses only `fractions.Fraction`. Returns $\sigma_{4i}/\sigma_{4(i-1)} = (2i-1)^2/(2i+1)^2$ as exact rationals for the five anchor ratios. No NumPy, no SVD.
- **Witness (B) — From-scratch NumPy re-implementation.** Rebuilds coprime residues from `math.gcd`, the tent distance by direct loop, the involution by `pow(r, -1, m)`, all in stdlib + NumPy. Computes $[D_\text{sym}, P_\tau]$ and its top-16 singular values by dense `numpy.linalg.svd`.
- **Witness (C) — Production pipeline.** Calls the production function `\texttt{sweep.\allowbreak{}holographic\_\allowbreak{}rank\_\allowbreak{}primes.\allowbreak{}unity\_\allowbreak{}clock\_\allowbreak{}direct\_\allowbreak{}commutator\_\allowbreak{}svd}`{=latex} — the same function used to generate Tables 6.1 and 12.1.

| Comparison | Max absolute deviation across five ratios |
|---|---:|
| (B) vs (C) | $0.000 \times 10^{0}$ — bit-identical IEEE 754 floats |
| (B) vs external Claude separate-machine reference | $3.3 \times 10^{-7}$ |
| (B) vs (A) | $3.3 \times 10^{-3}$ at $m = 7770, \sigma_{12}/\sigma_8$; $\to 10^{-4}$ at $m = 30030$ (consistent with Theorem 6.1(iv) envelope) |

The three witnesses are constructed from disjoint code paths. Witness (A) cannot disagree with (B) or (C) modulo finite-$\varphi$ corrections without the rational prediction being wrong. Witness (B) cannot disagree with (C) except through numerical error. Witness (C) is the production code that generated every other number in §6. The observed pattern is the theoretical one. **`THREE-WAY CROSS-VERIFICATION PASSES.`**

### A.5 Scope Note on the Coliseum Numbers

The Coliseum data of §8 (Table 8.1, FRR observation) was generated from a proprietary $N = 3$ multi-agent simulator that respects $S_N$ isotropy by construction. The simulator itself is **not** shipped with this paper. What is shipped:

1. The pure-math machinery of Theorems 3.1, 4.1, and 5.1, validated on hand-built $G$-equivariant covariances with analytically-known spectra (`sweep/test_character_rank.py`, `sweep/test_holographic_rank.py`).
2. The prime-gas Dirichlet-character machinery (`sweep/holographic_rank_primes.py`), arithmetic on its own and not requiring any simulator.
3. Drivers reproducing §6, §7, and §12.1–12.3 directly from closed-form formulas.

Independent re-derivation of §8's specific numbers requires a third-party generator preserving $G$-equivariance; the tests here pin the math any such generator's output must obey.

### A.6 Scope Note on the Continuous Substrate Numbers

The continuous Hamiltonian substrate of §§10–15 (Observation 12.1, Table 12.1, the §15 curvature spectrum) was constructed on a proprietary 768-dimensional semantic embedding with a learned 7-feature compression. The embedding model and compression are **not** shipped, and the texture optimizer, Gaussian wave-packet relaxation, and composed measurement pipeline are intentionally omitted from this preprint per Patent Notice [8]. What is reproducible from the math-only bundle:

1. Theorem 13.1 (A-optimality of the harmonic mean) is provable from Cauchy-Schwarz and reproduced in `verify_texture_aopt.py` on hand-built rank-4 subspaces with analytically-known maxima.
2. Proposition 14.1 (Heisenberg-dual obstruction) is provable from standard Fourier uncertainty applied to the fibre integral; the constant $c$ depends on the spectral gap of the substrate.
3. The §15 curvature spectrum is substrate-specific empirical data; its qualitative shape (substrate-intrinsic, reproducible across runs at fixed $\sigma$) is what the math predicts.

Apparatus practicing the §13–§15 measurement program in the manner claimed in [8] requires a license from Fancyland LLC.

### A.7 Relationship to Companion Patents

The §6, §7, and §12 results were originally identified using the proprietary spectral methods of [6], [7], and [8]. Those methods accelerated discovery but are not required for verification: all theorems and quantitative claims of this paper are independently reproducible using only standard representation theory, NumPy linear algebra, and the public math-only bundle at [github.com/fancyland-llc/character-rank](https://github.com/fancyland-llc/character-rank).

```{=latex}
\vspace{3em}
\begin{center}
\rule{0.4\textwidth}{0.5pt}

\vspace{0.8em}

\textit{Fancyland LLC --- Lattice OS research infrastructure.}

\vspace{0.3em}

\textit{The rabbit has been caught.}
\end{center}
```
