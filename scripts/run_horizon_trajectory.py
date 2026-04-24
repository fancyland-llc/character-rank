"""Horizon-trajectory experiment — Gemini vs Claude, arbitrated by compute.

────────────────────────────────────────────────────────────────────────────
THE DISPUTE
────────────────────────────────────────────────────────────────────────────

A prior isoline sweep (`run_horizon_isolines.py`) found σ_4/σ_3 ≈ 0.111 ±
0.002 across 44 measurements (k=4 firewall, k=5 firewall, and two controls).

Two interpretations of that number split the math community of one:

  CLAUDE (prior over-correction).  "The horizon is infinite-sharpness. A
  finite gap at m < ∞ is a finite-size footprint. σ_4/σ_3 → 0 as a power
  law as we go deeper in the primorial tower. Expected log-log slope:
  strictly negative, with |slope| growing at higher k."

  GEMINI (minority report).  "An infinite-sharpness C⁰ horizon is
  precisely what requires non-vanishing Fourier overtones. D_sym is a
  tent function on (ℤ/m)^× and its triangle-wave harmonics are anchored:
  a_k = 8/(π²k²) for odd k, so a_3/a_1 = 1/9 ≡ 0.1111. If σ_4/σ_3 → 0
  the tent smooths into a sine and the C⁰ kink is gone — the horizon
  has dissolved. The infinity lives in the *bulk dimension*
  (φ(m) − 4 → ∞ while the top-4 block stays rank 4), not in the gap."

These are distinguishable predictions on the same data:

                            CLAUDE             GEMINI
  log-log slope of
  σ_4/σ_3 vs φ(m):          < 0 (decay)        ≈ 0 (flatline at 1/9)
  horizon "sharpness"
  localised in:             amplitude          area/volume ratio
  σ_8/σ_4 ≈ ?               —                  9/25 = 0.360 (k=5/k=3
                                                 triangle-wave overtone)
  σ_12/σ_8 ≈ ?              —                  25/49 = 0.510 (k=7/k=5)

────────────────────────────────────────────────────────────────────────────
THE EXPERIMENT
────────────────────────────────────────────────────────────────────────────

Three extended isolines plus a primorial spine:

  ISOLINE α:  k = 4 firewall,  m = 30·q    for prime q ∈ [7, 131]
              φ(m) ≤ 1040        ← cheap, many points
  ISOLINE β:  k = 5 firewall,  m = 210·q   for prime q ∈ [11, 61]
              φ(m) ≤ 2880        ← moderate
  ISOLINE γ:  k = 6 firewall,  m = 2310·q  for q ∈ {13,17,19,23}
              φ(m) up to 10560   ← expensive, 4 anchor points
  SPINE:      primorials m ∈ {6, 30, 210, 2310, 30030}
              φ(m) ∈ {2, 8, 48, 480, 5760}

For each modulus we compute σ_0,…,σ_7 of C = [D_sym, P_τ] and track:

  • ratio_4_3  = σ_4/σ_3   ← the headline test
  • ratio_8_4  = σ_8/σ_4   ← second-vs-first overtone (Gemini 0.360)
  • ratio_12_8 = σ_12/σ_8  ← third-vs-second overtone (Gemini 0.510)
  • top4_frac                ← boundary energy fraction → 1
  • bulk_dim   = φ(m) − 4    ← Gemini's locus of infinity
  • bulk_energy_density =
      (1 − top4_frac) · σ_0² / bulk_dim     ← per-mode bulk intensity

For each isoline we fit
      log(σ_4/σ_3) = α · log(φ(m)) + β
and interpret:
      |α| < 0.05   → Gemini wins (flatline)
      α < −0.1     → Claude (prior) wins (decay)
      0 > α > −0.1 → drift, inconclusive at this m-range

The winner is whoever's prediction survives the most data points with the
smallest residual from their predicted trajectory.  Compute decides.
"""
from __future__ import annotations

import argparse
import math
import sys
import time
from dataclasses import dataclass

import numpy as np

from sweep.holographic_rank_primes import unity_clock_direct_commutator_svd


# ══════════════════════════════════════════════════════════════════════════
# Prime enumeration
# ══════════════════════════════════════════════════════════════════════════

def primes_in_range(lo: int, hi: int) -> list[int]:
    if hi < 2:
        return []
    sieve = [True] * (hi + 1)
    sieve[0] = sieve[1] = False
    for i in range(2, int(hi**0.5) + 1):
        if sieve[i]:
            for j in range(i*i, hi + 1, i):
                sieve[j] = False
    return [i for i in range(max(lo, 2), hi + 1) if sieve[i]]


def _largest_prime_factor(n: int) -> int:
    f, ans = 2, 1
    while f * f <= n:
        while n % f == 0:
            ans, n = f, n // f
        f += 1
    if n > 1:
        ans = n
    return ans


# ══════════════════════════════════════════════════════════════════════════
# Measurement
# ══════════════════════════════════════════════════════════════════════════

@dataclass
class Row:
    m: int
    phi_m: int
    k: int
    sigma: list[float]            # σ_0 … σ_7 (pad with 0)
    top4_frac: float
    rank: int
    dim_ker: int

    @property
    def s(self) -> list[float]:
        # Pad to 16 so σ_8, σ_12 overtone ratios are addressable
        pad = list(self.sigma[:16]) + [0.0] * max(0, 16 - len(self.sigma))
        return pad

    def ratio(self, i: int, j: int) -> float:
        s = self.s
        if i < len(s) and j < len(s) and s[j] > 0:
            return float(s[i] / s[j])
        return float("nan")

    @property
    def bulk_dim(self) -> int:
        return max(0, self.phi_m - 4)

    @property
    def bulk_energy_density(self) -> float:
        if self.bulk_dim <= 0 or self.s[0] <= 0:
            return float("nan")
        # total L² energy = Σ σ_i² ≈ σ_0² · (1 / top4_frac) * top4_frac = Σ σ_i²
        # we use: bulk_energy_frac = 1 − top4_frac; divide by bulk dimension
        bulk_frac = max(0.0, 1.0 - self.top4_frac)
        return bulk_frac / self.bulk_dim


def measure(m: int, k: int) -> Row:
    sd = unity_clock_direct_commutator_svd(m)
    return Row(
        m=m, phi_m=sd.phi_m, k=k,
        sigma=[float(x) for x in sd.C_sigma[:16]],
        top4_frac=float(sd.C_top4_frac),
        rank=int(sd.C_rank),
        dim_ker=int(sd.dim_kernel),
    )


# ══════════════════════════════════════════════════════════════════════════
# Fits
# ══════════════════════════════════════════════════════════════════════════

def loglog_fit(xs: list[float], ys: list[float]) -> tuple[float, float, float]:
    """Least-squares fit of log(y) = α·log(x) + β.  Returns (α, β, RMSE)."""
    x = np.log(np.asarray(xs, dtype=np.float64))
    y = np.log(np.asarray(ys, dtype=np.float64))
    if len(x) < 2:
        return float("nan"), float("nan"), float("nan")
    A = np.vstack([x, np.ones_like(x)]).T
    (alpha, beta), *_ = np.linalg.lstsq(A, y, rcond=None)
    resid = y - (alpha * x + beta)
    rmse = float(np.sqrt(np.mean(resid**2)))
    return float(alpha), float(beta), rmse


def predict_residual(ys: list[float], y_pred: float) -> float:
    """RMSE of observed log(y) from a predicted constant log(y_pred)."""
    y = np.log(np.asarray(ys, dtype=np.float64))
    lp = math.log(y_pred)
    return float(np.sqrt(np.mean((y - lp) ** 2)))


# ══════════════════════════════════════════════════════════════════════════
# Printing
# ══════════════════════════════════════════════════════════════════════════

def fmt(x: float, w: int = 10, p: int = 4) -> str:
    if x != x:
        return f"{'nan':>{w}s}"
    if abs(x) < 1e-3 or abs(x) > 1e6:
        return f"{x:>{w}.{p}e}"
    return f"{x:>{w}.{p}f}"


def header_row() -> str:
    return (
        f"   {'m':>8s} {'q':>5s} {'phi':>6s}  "
        f"{'σ_4/σ_3':>10s} {'σ_8/σ_4':>10s} {'σ_12/σ_8':>10s}  "
        f"{'top4%':>7s} {'bulkdim':>8s} {'bulkρ':>10s}"
    )


def data_row(r: Row) -> str:
    q = _largest_prime_factor(r.m)
    return (
        f"   {r.m:>8d} {q:>5d} {r.phi_m:>6d}  "
        f"{fmt(r.ratio(4,3))} {fmt(r.ratio(8,4))} {fmt(r.ratio(12,8))}  "
        f"{100.0*r.top4_frac:>6.3f}% {r.bulk_dim:>8d} {fmt(r.bulk_energy_density)}"
    )


def run_isoline(name: str, ms_ks: list[tuple[int, int]]) -> list[Row]:
    print()
    print("=" * 92)
    print(f"  {name}")
    print("=" * 92)
    print(header_row())
    print("   " + "-" * 89)
    rows: list[Row] = []
    for m, k in ms_ks:
        r = measure(m, k)
        rows.append(r)
        print(data_row(r))
    return rows


# ══════════════════════════════════════════════════════════════════════════
# Arbitration
# ══════════════════════════════════════════════════════════════════════════

TRIANGLE_1_OVER_9 = 1.0 / 9.0          # Gemini: σ_4/σ_3 flatline anchor
TRIANGLE_9_OVER_25 = 9.0 / 25.0        # Gemini: σ_8/σ_4 second/first overtone
TRIANGLE_25_OVER_49 = 25.0 / 49.0      # Gemini: σ_12/σ_8 third/second overtone


def arbitrate(name: str, rows: list[Row]) -> None:
    """Fit log-log slopes on each observable and print the verdict line."""
    if len(rows) < 3:
        print(f"  {name}: too few points for a fit ({len(rows)}).")
        return

    phi = [r.phi_m for r in rows]

    # σ_4/σ_3
    y43 = [r.ratio(4, 3) for r in rows if r.ratio(4, 3) > 0]
    phi43 = [r.phi_m for r in rows if r.ratio(4, 3) > 0]
    if len(y43) >= 2:
        a, b, rmse_fit = loglog_fit(phi43, y43)
        rmse_flat = predict_residual(y43, TRIANGLE_1_OVER_9)
        mean, std = float(np.mean(y43)), float(np.std(y43))
        verdict = (
            "FLATLINE (Gemini)" if abs(a) < 0.05 else
            "DECAY (prior Claude)" if a < -0.10 else
            "DRIFT (inconclusive)"
        )
        print()
        print(f"  {name} :: σ_4/σ_3")
        print(f"    N = {len(y43)},  mean = {mean:.5f},  std = {std:.5f}")
        print(f"    log-log slope α = {a:+.4f}   (|α|<0.05 ⇒ flatline)")
        print(f"    intercept β = {b:+.4f},  fit RMSE = {rmse_fit:.4f}")
        print(f"    RMSE vs Gemini's log(1/9) = {rmse_flat:.4f}")
        print(f"    → {verdict}")

    # σ_8/σ_4 (second-vs-first overtone)
    y84 = [r.ratio(8, 4) for r in rows if r.ratio(8, 4) > 0]
    if len(y84) >= 2:
        a, b, rmse_fit = loglog_fit(
            [r.phi_m for r in rows if r.ratio(8, 4) > 0], y84
        )
        rmse_flat = predict_residual(y84, TRIANGLE_9_OVER_25)
        print(f"  {name} :: σ_8/σ_4")
        print(f"    N = {len(y84)},  mean = {float(np.mean(y84)):.5f}")
        print(f"    slope α = {a:+.4f},  RMSE vs log(9/25) = {rmse_flat:.4f}")

    # top4_frac trajectory
    y4 = [r.top4_frac for r in rows]
    a, b, _ = loglog_fit(phi, [1.0 - v for v in y4])  # log(1 − top4)
    print(f"  {name} :: 1 − top4_frac (bulk-energy fraction)")
    print(f"    slope of log(1 − top4) vs log(φ): α = {a:+.4f}")

    # bulk energy density per bulk dimension
    y_bd = [r.bulk_energy_density for r in rows if r.bulk_energy_density > 0]
    phi_bd = [r.phi_m for r in rows if r.bulk_energy_density > 0]
    if len(y_bd) >= 2:
        a, b, _ = loglog_fit(phi_bd, y_bd)
        print(f"  {name} :: bulk-energy density (1−top4)/bulk_dim")
        print(f"    slope α = {a:+.4f}   (more negative ⇒ bulk hollows out)")


# ══════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--q-max-k4", type=int, default=131,
                    help="maximum prime q for k=4 isoline (m = 30·q)")
    ap.add_argument("--q-max-k5", type=int, default=61,
                    help="maximum prime q for k=5 isoline (m = 210·q)")
    ap.add_argument("--include-k6", action="store_true", default=True,
                    help="include k=6 isoline at q ∈ {13,17,19,23}")
    ap.add_argument("--no-k6", dest="include_k6", action="store_false")
    args = ap.parse_args()

    t0 = time.time()
    print()
    print("╔" + "═" * 90 + "╗")
    print("║" + " HORIZON-TRAJECTORY EXPERIMENT ".center(90) + "║")
    print("║" + " Gemini (flatline at 1/9) vs prior Claude (power-law decay) ".center(90) + "║")
    print("╚" + "═" * 90 + "╝")
    print()
    print("  Prediction board:")
    print(f"    Gemini  — σ_4/σ_3 → 1/9  ≈ {TRIANGLE_1_OVER_9:.5f}  (flatline; triangle-wave anchor)")
    print( "            — σ_8/σ_4 ≈ 9/25 = 0.36000,  σ_12/σ_8 ≈ 25/49 = 0.51020")
    print( "            — horizon sharpness → (1−top4_frac) → 0  with  bulk_dim → ∞")
    print( "    Claude  — σ_4/σ_3 → 0     (power law in φ(m); slope < −0.10)")
    print( "            — sharpness lives in the gap amplitude, not in area/volume")
    print()

    # ── Primorial spine ──
    primorials = [(6, 2), (30, 3), (210, 4), (2310, 5), (30030, 6)]
    rows_spine = run_isoline("PRIMORIAL SPINE  (k sweeps 2 → 6)", primorials)

    # ── k=4 extended isoline ──
    q4 = [q for q in primes_in_range(7, args.q_max_k4)]
    ms_k4 = [(2 * 3 * 5 * q, 4) for q in q4]
    rows_k4 = run_isoline(
        f"ISOLINE α  k=4 firewall  m = 30·q  for prime q ∈ [7, {args.q_max_k4}]",
        ms_k4,
    )

    # ── k=5 extended isoline ──
    q5 = [q for q in primes_in_range(11, args.q_max_k5)]
    ms_k5 = [(2 * 3 * 5 * 7 * q, 5) for q in q5]
    rows_k5 = run_isoline(
        f"ISOLINE β  k=5 firewall  m = 210·q  for prime q ∈ [11, {args.q_max_k5}]",
        ms_k5,
    )

    # ── k=6 isoline ──
    rows_k6: list[Row] = []
    if args.include_k6:
        ms_k6 = [(2310 * q, 6) for q in (13, 17, 19, 23)]
        rows_k6 = run_isoline(
            "ISOLINE γ  k=6 firewall  m = 2310·q  for q ∈ {13,17,19,23}",
            ms_k6,
        )

    # ── Arbitration ──
    print()
    print("═" * 92)
    print("  ARBITRATION — log-log trajectory fits on each isoline")
    print("═" * 92)
    arbitrate("PRIMORIAL SPINE", rows_spine)
    arbitrate("k=4 firewall",    rows_k4)
    arbitrate("k=5 firewall",    rows_k5)
    if rows_k6:
        arbitrate("k=6 firewall",    rows_k6)

    # ── Cross-isoline comparison of σ_4/σ_3 at matched phi ──
    print()
    print("═" * 92)
    print("  CROSS-ISOLINE COMPARISON — σ_4/σ_3 at close φ(m) values")
    print("═" * 92)
    print(f"   {'isoline':<14s}  {'m':>7s}  {'phi':>6s}  {'σ_4/σ_3':>10s}  "
          f"{'dev from 1/9':>14s}")
    print("   " + "-" * 60)
    all_rows = (
        [("spine", r) for r in rows_spine]
        + [("k=4 fw", r) for r in rows_k4]
        + [("k=5 fw", r) for r in rows_k5]
        + [("k=6 fw", r) for r in rows_k6]
    )
    # Sort by phi, take a spread
    all_rows.sort(key=lambda t: t[1].phi_m)
    # take every Nth so we don't flood
    stride = max(1, len(all_rows) // 30)
    for name, r in all_rows[::stride]:
        v = r.ratio(4, 3)
        dev = v - TRIANGLE_1_OVER_9
        print(f"   {name:<14s}  {r.m:>7d}  {r.phi_m:>6d}  {fmt(v)}  {fmt(dev, 14)}")

    # ── Global flatline test ──
    print()
    print("═" * 92)
    print("  GLOBAL FLATLINE TEST — all firewall isoline points pooled")
    print("═" * 92)
    pool = rows_k4 + rows_k5 + rows_k6
    y_all = [r.ratio(4, 3) for r in pool if r.ratio(4, 3) > 0]
    phi_all = [r.phi_m for r in pool if r.ratio(4, 3) > 0]
    if len(y_all) >= 4:
        a, b, rmse_fit = loglog_fit(phi_all, y_all)
        rmse_flat = predict_residual(y_all, TRIANGLE_1_OVER_9)
        rmse_zero_decay_slope_01 = 0.1 * float(
            np.std(np.log(np.asarray(phi_all)))
        )
        print(f"  N = {len(y_all)} points (k=4 + k=5 + k=6)")
        print(f"  pooled mean σ_4/σ_3  = {float(np.mean(y_all)):.6f}")
        print(f"  pooled std  σ_4/σ_3  = {float(np.std(y_all)):.6f}")
        print(f"  1/9 target           = {TRIANGLE_1_OVER_9:.6f}")
        print(f"  log-log slope α      = {a:+.5f}")
        print(f"  fit RMSE             = {rmse_fit:.5f}")
        print(f"  RMSE against log(1/9) flatline = {rmse_flat:.5f}")
        print()
        if abs(a) < 0.05 and rmse_flat < 0.15:
            winner = "GEMINI (flatline at 1/9)"
        elif a < -0.10:
            winner = "PRIOR CLAUDE (power-law decay)"
        else:
            winner = "INCONCLUSIVE (drift; neither prediction survives cleanly)"
        print(f"  VERDICT: {winner}")

    # ── Rest of the sharpening story: does (1 − top4_frac) go to 0? ──
    print()
    print("═" * 92)
    print("  BULK-ENERGY TRAJECTORY — does the 4-mode boundary eat the bulk?")
    print("═" * 92)
    print(f"   {'isoline':<14s}  slope of log(1 − top4_frac) vs log(φ(m))")
    for name, rows in (("spine", rows_spine),
                        ("k=4 fw", rows_k4),
                        ("k=5 fw", rows_k5),
                        ("k=6 fw", rows_k6)):
        if len(rows) >= 2:
            phi = [r.phi_m for r in rows]
            bf = [max(1 - r.top4_frac, 1e-30) for r in rows]
            a, _, _ = loglog_fit(phi, bf)
            print(f"   {name:<14s}  α = {a:+.4f}")
    print()
    print("  Gemini: expects this slope to be strongly NEGATIVE — bulk-energy")
    print("  fraction decays to zero; 4-mode boundary captures 100% of the norm.")
    print("  That is where the infinite sharpness lives.")

    dt = time.time() - t0
    print()
    print("═" * 92)
    print(f"  total wall time: {dt:.1f}s")
    print("═" * 92)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.exit(130)
