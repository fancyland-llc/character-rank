"""Unity Clock chiral decomposition — with the rigorous scope distinction.

The chiral identity

    C u_χ = (λ(χ̄) − λ(χ)) u_{χ̄} = −2i Im(λ(χ)) u_{χ̄}

holds *only* for G-equivariant kernels, i.e., kernels K(r,s) = k(s r^{-1})
invariant under (ℤ/m)^× acting by multiplication.

  • D^CRT (CRT-local, direct-sum per-prime circular distance) IS
    G-equivariant: Dirichlet characters diagonalize it exactly and the
    chiral formula applies rigorously.
  • D_sym (natural tent distance on coprime residues) is NOT
    G-equivariant (it's additively invariant on ℤ/m, not multiplicatively
    invariant on (ℤ/m)^×). Its eigenvectors are Möbius-sieved ambient
    Fourier modes (Unity Clock §3.4), not Dirichlet characters.

So we run two tracks:

  (A) CRT-local chiral decomposition — character basis, closed form,
      should confirm σ_C(D^CRT) = 2|Im λ(χ)| with exact 2-degeneracy.

  (B) Direct D_sym commutator SVD — matrix-level O(n³), reproduces
      Unity Clock §3.1-§3.3 σ_0/σ_2 → 1 convergence and 4-mode
      dominance of C = [D_sym, P_τ].

The "chiral imaginary wavefunction" is real on both tracks, but in
different bases: Dirichlet characters for D^CRT (exact factorization),
and Fourier tent-harmonics for D_sym (which the UC paper already
formalized via σ₀/σ₂ → 1).
"""
from __future__ import annotations

import sys
import time

import numpy as np

from sweep.holographic_rank_primes import (
    analytic_tc_vs_time,
    unity_clock_chiral_decomposition,
    unity_clock_direct_commutator_svd,
)


PRIMORIALS = [6, 30, 210, 2310, 30030]


def fmt_e(v: float, width: int = 11) -> str:
    if v != v:
        return " " * (width - 3) + "nan"
    if abs(v) < 1e-3 or abs(v) > 1e6:
        return f"{v:>{width}.3e}"
    return f"{v:>{width}.3f}"


def scan_peak(m: int, kernel: str, lo: float, hi: float,
              n: int = 80) -> tuple[float, float]:
    schedule = np.logspace(np.log10(lo), np.log10(hi), n).tolist()
    rows = analytic_tc_vs_time(m, schedule, kernel=kernel)
    tc = np.array([r["TC"] for r in rows])
    bt = np.array([r["beta_t"] for r in rows])
    i = int(np.argmax(tc))
    return float(bt[i]), float(tc[i])


def main() -> None:
    print()
    print("═" * 82)
    print(" UNITY CLOCK / CHIRAL-IMAGINARY-WAVEFUNCTION — FULL TOWER")
    print("═" * 82)
    print()
    print(" The chiral identity  C u_χ = −2i Im(λ(χ)) u_{χ̄}  holds rigorously only")
    print(" on G-equivariant kernels. This splits the analysis into two tracks:")
    print()
    print("   (A) D^CRT-local   — G-equivariant by construction,")
    print("                        Dirichlet characters diagonalize exactly.")
    print("   (B) D_sym natural — NOT G-equivariant (additively invariant, not")
    print("                        multiplicatively). Eigenvectors are Möbius-")
    print("                        sieved ambient Fourier modes (UC §3.4).")
    print()

    t0 = time.time()

    # ── TRACK A: CRT-local chiral decomposition ────────────────────────────
    print("┏" + "━" * 80 + "┓")
    print("┃ (A) CRT-LOCAL:  σ_C(χ) = 2|Im λ(χ)| on Dirichlet characters" +
          " " * 19 + "┃")
    print("┗" + "━" * 80 + "┛")
    print()
    hdr = (
        f"  {'m':>6} {'φ(m)':>6}  {'#real χ':>8}  {'#chiral prs':>12}  "
        f"{'max σ_C':>11}  {'note':>28}"
    )
    print(hdr)
    print("  " + "─" * (len(hdr) - 2))
    print("  (Expected max σ_C ≈ 0: CRT-local dynamics don't couple primes.)")
    for m in PRIMORIALS:
        sc = unity_clock_chiral_decomposition(m, kernel="crt")
        n_real = int(sc.is_real_char.sum())
        n_chi_prs = (sc.phi_m - n_real) // 2
        max_sigma = float(sc.sigma_C[0]) if len(sc.sigma_C) > 0 else 0.0
        note = "chirally zero ✓" if max_sigma < 1e-9 else "↳ per-prime chirality only"
        print(
            f"  {m:>6d} {sc.phi_m:>6d}  {n_real:>8d}  {n_chi_prs:>12d}  "
            f"{fmt_e(max_sigma)}  {note:>28}"
        )
    print()

    # ── TRACK B: Direct D_sym commutator SVD (matches UC paper) ────────────
    print("┏" + "━" * 80 + "┓")
    print("┃ (B) NATURAL D_sym:  direct SVD of C = [D_sym, P_τ]  (UC §3.1-§3.3)" +
          " " * 13 + "┃")
    print("┗" + "━" * 80 + "┛")
    print()
    hdr = (
        f"  {'m':>6}  {'σ_0':>11}  {'σ_2/σ_0':>8}  {'top4 %':>7}  "
        f"{'rank(C)':>8}  {'dim ker':>8}  {'eff_90':>7}"
    )
    print(hdr)
    print("  " + "─" * (len(hdr) - 2))
    directs = {}
    for m in PRIMORIALS:
        sd = unity_clock_direct_commutator_svd(m)
        directs[m] = sd
        print(
            f"  {m:>6d}  {fmt_e(sd.C_sigma[0] if len(sd.C_sigma) > 0 else 0.0)}  "
            f"{sd.C_sigma_ratio_pair_02:>8.4f}  "
            f"{100.0*sd.C_top4_frac:>6.2f}%  "
            f"{sd.C_rank:>8d}  {sd.dim_kernel:>8d}  {sd.C_effective_rank_90:>7d}"
        )
    print()
    print("  UC §3.1: dim ker(C) = 2^{k-1} for primorial with k primes.")
    print("  UC §3.3: σ_2/σ_0 → 1 as m → ∞  (top 4 modes degenerate into 2 pairs).")
    print()

    # ── Compare to UC paper's reported σ_0 ─────────────────────────────────
    uc_reported = {30: None, 210: None, 2310: None,
                   30030: 17526698.0, 510510: 4.767e9}
    print("  Cross-check against UC paper §3.3 reported σ_0:")
    for m in PRIMORIALS:
        if uc_reported.get(m) is not None:
            mine = float(directs[m].C_sigma[0])
            print(f"    m={m}: direct σ_0 = {mine:.3e};  UC paper {uc_reported[m]:.3e};"
                  f"  ratio {mine/uc_reported[m]:.4f}")
    print()

    # ── Show floor-by-floor Unity Clock 4-mode spectrum from direct ───────
    print("┏" + "━" * 80 + "┓")
    print("┃ UNITY-CLOCK 4-MODE STRUCTURE (direct SVD of C = [D_sym, P_τ])" +
          " " * 17 + "┃")
    print("┗" + "━" * 80 + "┛")
    print()
    hdr = f"  {'m':>6}   {'σ_0':>11} {'σ_1':>11} {'σ_2':>11} {'σ_3':>11}   {'top4 frac':>10}"
    print(hdr)
    print("  " + "─" * (len(hdr) - 2))
    for m in PRIMORIALS:
        sd = directs[m]
        s = sd.C_sigma
        pad = list(s[:4]) + [0.0] * (4 - min(len(s), 4))
        print(
            f"  {m:>6d}   "
            f"{fmt_e(float(pad[0]))} {fmt_e(float(pad[1]))} "
            f"{fmt_e(float(pad[2]))} {fmt_e(float(pad[3]))}   "
            f"{100.0*sd.C_top4_frac:>9.3f}%"
        )
    print()

    # ── Peak-βt of analytic_tc (character Fourier dynamic) across tower ────
    print("┏" + "━" * 80 + "┓")
    print("┃ CHARACTER-MARGINAL CROSS-CORRELATION: peak TC vs βt across tower" +
          " " * 14 + "┃")
    print("┗" + "━" * 80 + "┛")
    print()
    print("  (This is the analytic_tc_vs_time dynamic with kernel='natural'.")
    print("  Its 'eigenvalues' are the Dirichlet-character Fourier coefficients")
    print("  of k(g) = K(1,g), not eigenvalues of D_sym. It's the cross-prime")
    print("  character-marginal correlation observable — well-defined, but")
    print("  NOT the physical D_sym diffusion.)")
    print()
    hdr = f"  {'m':>6}  {'meas βt*':>12}  {'peak TC':>8}  {'σ_0(C)':>12}  {'σ_0/m':>10}"
    print(hdr)
    print("  " + "─" * (len(hdr) - 2))
    for m in PRIMORIALS:
        sd = directs[m]
        meas_bt, peak_tc = scan_peak(m, "natural", 1e-2, 1e7, n=200)
        ratio = float(sd.C_sigma[0]) / m if m > 0 else 0.0
        print(
            f"  {m:>6d}  {fmt_e(meas_bt)}  {peak_tc:>8.4f}  "
            f"{fmt_e(float(sd.C_sigma[0]) if len(sd.C_sigma) > 0 else 0.0)}  "
            f"{fmt_e(ratio)}"
        )
    print()

    print("═" * 82)
    print(f"  wall time: {time.time() - t0:.1f}s")
    print("═" * 82)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.exit(130)
