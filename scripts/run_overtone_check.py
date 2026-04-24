"""Quick overtone-ladder check on Gemini's triangle-wave anchors.

Triangle wave Fourier amplitudes are a_k = 8/(π²k²) for odd k.
If the 4-mode SVD blocks of C = [D_sym, P_τ] correspond to k=1, k=3,
k=5, k=7 harmonics (four modes per block, from two-parity × two
reflection-symmetry), then:

    block k=1:  σ_0, σ_1, σ_2, σ_3   amplitude ∝ 1
    block k=3:  σ_4, σ_5, σ_6, σ_7   amplitude ∝ 1/9
    block k=5:  σ_8, σ_9, σ_10, σ_11 amplitude ∝ 1/25
    block k=7:  σ_12, σ_13, σ_14, σ_15 amplitude ∝ 1/49

Predicted block-to-block ratios (all four modes land at roughly the same
value within a block):

    σ_3:σ_0 ratio is internal to the k=1 block; → 1 asymptotically (UC §3.3)
    σ_4/σ_3 → 1/9  = 0.11111         (first overtone / fundamental)
    σ_8/σ_4 → 9/25 = 0.36000         (second overtone / first)
    σ_8/σ_7 → 1/25 / 1/9 = 9/25 = 0.36000   (same)
    σ_12/σ_8 → 25/49 = 0.51020       (third overtone / second)
    σ_12/σ_11 → 25/49 = 0.51020      (same)

Also the intra-block gaps should all be ≈ 1 (within block, σ_{4k+3}/σ_{4k}).

This only needs the large-φ anchor points (k=5, k=6 isoline tail).
"""
from __future__ import annotations

import time

import numpy as np

from sweep.holographic_rank_primes import unity_clock_direct_commutator_svd


# Gemini's triangle-wave targets
T_1_9   = 1.0/9.0      # σ_4/σ_3
T_9_25  = 9.0/25.0     # σ_8/σ_4, σ_8/σ_7
T_1_25  = 1.0/25.0     # σ_8/σ_0
T_25_49 = 25.0/49.0    # σ_12/σ_8, σ_12/σ_11
T_1_49  = 1.0/49.0     # σ_12/σ_0


def main() -> None:
    # Use the existing k=5 tail + k=6 points (all have φ ≥ 1000, so
    # σ_12 is populated comfortably)
    modali = [
        ("k=5", 210 * 29),    # m = 6090,   φ = 1344
        ("k=5", 210 * 37),    # m = 7770,   φ = 1728
        ("k=5", 210 * 47),    # m = 9870,   φ = 2208
        ("k=5", 210 * 61),    # m = 12810,  φ = 2880
        ("k=6", 2310 * 13),   # m = 30030,  φ = 5760
        ("k=6", 2310 * 17),   # m = 39270,  φ = 7680
        ("k=6", 2310 * 19),   # m = 43890,  φ = 8640
        ("k=6", 2310 * 23),   # m = 53130,  φ = 10560
    ]

    t0 = time.time()
    print()
    print("=" * 96)
    print("  OVERTONE LADDER — three triangle-wave anchors on a single kernel")
    print("=" * 96)
    print()
    print("  Gemini's predicted ratios from triangle-wave Fourier coefficients")
    print(f"     σ_4/σ_3    →  1/9   = {T_1_9:.5f}     (block k=3 / block k=1)")
    print(f"     σ_8/σ_4    →  9/25  = {T_9_25:.5f}    (block k=5 / block k=3)")
    print(f"     σ_12/σ_8   →  25/49 = {T_25_49:.5f}   (block k=7 / block k=5)")
    print()
    print(f"  Cross-tier:")
    print(f"     σ_8/σ_0    →  1/25  = {T_1_25:.5f}")
    print(f"     σ_12/σ_0   →  1/49  = {T_1_49:.5f}")
    print()
    print(f"  {'kind':>6s} {'m':>8s} {'phi':>6s}    "
          f"{'σ_4/σ_3':>10s} {'σ_8/σ_4':>10s} {'σ_12/σ_8':>10s}  "
          f"{'σ_8/σ_0':>10s} {'σ_12/σ_0':>10s}")
    print("  " + "-" * 94)

    r43, r84, r128 = [], [], []
    r80, r120 = [], []
    for kind, m in modali:
        sd = unity_clock_direct_commutator_svd(m)
        s = list(sd.C_sigma)
        if len(s) < 16:
            s = s + [0.0] * (16 - len(s))
        sigma_0, sigma_3, sigma_4 = s[0], s[3], s[4]
        sigma_7, sigma_8 = s[7], s[8]
        sigma_11, sigma_12 = s[11], s[12]

        def rr(a, b):
            return (a / b) if b > 0 else float("nan")

        v43   = rr(sigma_4,  sigma_3)
        v84   = rr(sigma_8,  sigma_4)
        v128  = rr(sigma_12, sigma_8)
        v80   = rr(sigma_8,  sigma_0)
        v120  = rr(sigma_12, sigma_0)
        r43.append(v43); r84.append(v84); r128.append(v128)
        r80.append(v80); r120.append(v120)
        print(
            f"  {kind:>6s} {m:>8d} {sd.phi_m:>6d}    "
            f"{v43:>10.5f} {v84:>10.5f} {v128:>10.5f}  "
            f"{v80:>10.5f} {v120:>10.5f}"
        )

    print()
    print("  MEAN / DEVIATION FROM GEMINI TARGETS")
    print(f"     σ_4/σ_3    mean = {np.mean(r43):.5f}   target 1/9   = {T_1_9:.5f}   "
          f"abs dev = {abs(np.mean(r43) - T_1_9):.5f}")
    print(f"     σ_8/σ_4    mean = {np.mean(r84):.5f}   target 9/25  = {T_9_25:.5f}   "
          f"abs dev = {abs(np.mean(r84) - T_9_25):.5f}")
    print(f"     σ_12/σ_8   mean = {np.mean(r128):.5f}   target 25/49 = {T_25_49:.5f}   "
          f"abs dev = {abs(np.mean(r128) - T_25_49):.5f}")
    print(f"     σ_8/σ_0    mean = {np.mean(r80):.5f}   target 1/25  = {T_1_25:.5f}   "
          f"abs dev = {abs(np.mean(r80) - T_1_25):.5f}")
    print(f"     σ_12/σ_0   mean = {np.mean(r120):.5f}   target 1/49  = {T_1_49:.5f}   "
          f"abs dev = {abs(np.mean(r120) - T_1_49):.5f}")
    print()
    print(f"  wall time: {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
