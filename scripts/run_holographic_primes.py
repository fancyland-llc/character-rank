"""Bonkers-scale holographic experiment on the prime gas.

This is the Shadow-3 dictionary-row-2 test: apply Theorem 11.1 (character-
entropy decomposition) to the prime gas at a ladder of primorial bases,
subgroup refinements, evolution depths, Haar baselines, CRT factorizations,
and a functorial comparison to the Coliseum at m = 6.

Design goal: not "does the decomposition hold?" — it *has* to hold by
construction whenever the ensemble is G-equivariant. The interesting
data is WHERE the decomposition diverges from the naive expectation,
what the subgroup refinement tower reveals, and whether the dictionary
entry at m = 6 matches the Coliseum APP signature. These are the places
the next rabbit hides.

Output: a tower of tables, each emphasising a different invariant.
Grep the output for 'RABBIT:' lines — those flag observations that
open the solution space to higher-order systems.
"""
from __future__ import annotations

import argparse
import time
from math import log

import numpy as np

from sweep.holographic_rank import (
    character_entropy_decomposition,
)
from sweep.holographic_rank_primes import (
    abelian_character_decomposition,
    analytic_tc_vs_time,
    boundary_covariance,
    build_crt_local_distance_kernel,
    build_distance_kernel_boundary,
    coprime_residues,
    crt_character_marginal_check,
    crt_factorization_check,
    dirichlet_character_table,
    enumerate_subgroups,
    group_order,
    haar_random_entropy_baseline,
    page_curve_on_residues,
    prime_gas_ensemble,
    subgroup_refinement,
    von_neumann_entropy_complex,
)


PRIMORIALS = (6, 30, 210, 2310)
T_LADDER = (0.0, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0)


# ──────────────────────────────────────────────────────────────────────────
# Pretty printers
# ──────────────────────────────────────────────────────────────────────────

def print_header(title: str) -> None:
    print()
    print("=" * 78)
    print(f"  {title}")
    print("=" * 78)


def print_sub(title: str) -> None:
    print()
    print(f"-- {title} --")


def rabbit(msg: str) -> None:
    """Flag an observation worth hunting — a rabbit hole seed."""
    print(f"   RABBIT: {msg}")


# ──────────────────────────────────────────────────────────────────────────
# Run sections
# ──────────────────────────────────────────────────────────────────────────

def section_1_character_tables(primorials: tuple[int, ...]) -> None:
    print_header("1. Dirichlet character tables — the dictionary's Fourier basis")
    for m in primorials:
        G = coprime_residues(m)
        n = len(G)
        if n > 32:
            print(f"  m = {m}:  |G| = {n}  (character table suppressed, too big)")
            continue
        Lam = dirichlet_character_table(m)
        # Schur orthogonality check
        ortho = Lam @ Lam.conj().T / n
        err = float(np.max(np.abs(ortho - np.eye(n))))
        print(f"  m = {m}:  |G| = φ(m) = {n}   Schur-orthogonality err = {err:.2e}")
        if n <= 8:
            print(f"           coprimes = {G}")
            # Print character magnitudes rounded
            for r in range(n):
                row = "  ".join(f"{c: 6.3f}" for c in np.real(Lam[r]))
                print(f"           χ_{r}  = [{row}]")
        if err > 1e-6:
            rabbit(
                f"character table at m={m} has Schur-orthogonality error "
                f"{err:.2e} — investigate numerical conditioning of the "
                f"eigendecomposition"
            )


def section_2_character_entropy_primorial_ladder(
    primorials: tuple[int, ...], M: int, t_schedule: tuple[float, ...],
) -> list[dict]:
    print_header(
        "2. Character-entropy decomposition across the primorial ladder"
    )
    rows = []
    for m in primorials:
        n = group_order(m)
        print_sub(f"m = {m}   |G| = {n}   log|G| = {log(n):.4f}")
        print("   t        S(ρ)      H(p)      leakage   top_char_p   entropy_ratio")
        print("   ----    -------   -------   --------   ----------   -------------")
        psi = prime_gas_ensemble(m, M=M, t_schedule=t_schedule)
        for ti, t in enumerate(t_schedule):
            C = boundary_covariance(psi, time_idx=ti)
            full = abelian_character_decomposition(C, m)
            ratio = full.shannon / log(n) if n > 1 else 0.0
            top_p = float(np.max(full.p_char))
            print(
                f"   {t:5.2f}    {full.direct_entropy:7.4f}   "
                f"{full.shannon:7.4f}   {full.nonabelian_leakage:8.1e}   "
                f"{top_p:10.4f}   {ratio:13.4f}"
            )
            rows.append({
                "m": m, "n": n, "t": t,
                "S_direct": full.direct_entropy,
                "S_shannon": full.shannon,
                "leakage": full.nonabelian_leakage,
                "top_char_p": top_p,
                "ratio_to_logn": ratio,
            })
        last = rows[-1]
        prev = rows[-2]
        rel = abs(last["S_direct"] - prev["S_direct"]) / max(1e-9, prev["S_direct"])
        print(
            f"   [m={m}] last-doubling change in S(ρ): {100.0 * rel:.2f}%  "
            f"({'SATURATES' if rel < 0.02 else 'still moving'})"
        )
        if last["ratio_to_logn"] > 0.999:
            rabbit(
                f"at m={m}, ratio_to_logn = {last['ratio_to_logn']:.6f} — "
                f"the prime-gas boundary covariance achieves MAXIMUM possible "
                f"entropy log|G|. Dynamics are indistinguishable from Haar "
                f"random at the abelian-character resolution."
            )
        if last["ratio_to_logn"] < 0.50:
            rabbit(
                f"at m={m}, ratio_to_logn = {last['ratio_to_logn']:.4f} — "
                f"well below the log|G| ceiling. Strong character-localization."
            )
    return rows


def section_3_subgroup_refinement_tower(m: int, M: int) -> None:
    print_header(
        f"3. Subgroup refinement tower at m = {m} — "
        "where the non-trivial Schur bonuses live"
    )
    n = group_order(m)
    if n > 24:
        print(f"  skipped: |G| = {n} too many subgroups to enumerate")
        return
    psi = prime_gas_ensemble(m, M=M, t_schedule=(2.0,))
    C = boundary_covariance(psi, time_idx=0)
    subs = enumerate_subgroups(m)
    print("   |H|   [G:H]    H(p)       internal    total        direct      recon")
    print("   ---   -----    -------    --------    -------      -------     --------")
    S_refined = []
    for H in subs:
        ref = subgroup_refinement(C, m, H)
        S_refined.append((ref, H))
        print(
            f"   {len(H):3d}   {ref.subgroup_index:5d}    "
            f"{ref.shannon:7.4f}    {ref.internal:8.4f}    "
            f"{ref.total:7.4f}     {ref.direct:7.4f}     "
            f"{ref.reconstruction_err:.2e}"
        )
    # Refinement bonus: (internal entropy at H) - 0 (at G)
    # Reconstruction — largest H reaches direct entropy; smaller H misses.
    # INTERESTING: what subgroup index MAXIMIZES internal entropy / minimizes
    # reconstruction error? Those identify the "natural" effective symmetry.
    best_internal = max(S_refined, key=lambda r: r[0].internal)
    best_recon = min(S_refined, key=lambda r: r[0].reconstruction_err)
    print()
    print(
        f"   most internal at  |H| = {len(best_internal[1])}:  "
        f"internal = {best_internal[0].internal:.4f}"
    )
    print(
        f"   best reconstruction at  |H| = {len(best_recon[1])}:  "
        f"err = {best_recon[0].reconstruction_err:.2e}"
    )
    if best_internal[0].internal > 0.1 * log(n):
        rabbit(
            f"at m={m}, some subgroup of index {best_internal[0].subgroup_index} "
            f"carries internal entropy {best_internal[0].internal:.4f} — "
            f"a measurable 'missed-symmetry' signature"
        )


def section_4_haar_baseline(primorials: tuple[int, ...], M: int) -> None:
    print_header("4. Haar-random baseline — how ordered is the prime gas?")
    print("   m      φ(m)    Haar_mean ± std    log|G|     ratio_to_sat")
    print("   ---    ----    ---------------    -------    ------------")
    for m in primorials:
        n = group_order(m)
        if n > 50:
            continue
        hb = haar_random_entropy_baseline(n, M=M)
        print(
            f"   {m:3d}    {n:4d}    "
            f"{hb['mean_entropy']:7.4f} ± {hb['std_entropy']:.4f}   "
            f"{hb['saturation_log_n']:7.4f}    {hb['ratio_to_saturation']:12.4f}"
        )
        if hb['ratio_to_saturation'] < 0.95:
            rabbit(
                f"Haar baseline at n={n} sits at {hb['ratio_to_saturation']:.4f} "
                f"of log|G| saturation — finite-M Wishart curvature, not a "
                f"violation. Use this as the noise floor for prime-gas comparison."
            )


def section_5_crt_factorization(
    primorials: tuple[int, ...], M: int,
) -> list[dict]:
    print_header(
        "5. CRT factorization — does entropy split across prime factors?"
    )
    results = []
    for m in primorials:
        n = group_order(m)
        if n > 24:
            continue
        psi = prime_gas_ensemble(m, M=M, t_schedule=(2.0,))
        C = boundary_covariance(psi, time_idx=0)
        crt = crt_factorization_check(m, C)
        print_sub(f"m = {m}  factors = {crt['factors']}")
        print(f"   full Shannon:        {crt['full_shannon']:.4f}")
        for p, s in crt['per_factor_shannon'].items():
            print(f"   factor p = {p:3d}:       {s:.4f}")
        print(f"   sum of factors:      {crt['sum_of_factor_shannons']:.4f}")
        print(f"   residual:            {crt['factorization_residual']:+.4f}  "
              f"({100 * crt['factorization_residual'] / max(1e-9, crt['full_shannon']):.1f}%)")
        results.append(crt)
        if abs(crt['factorization_residual']) < 0.02 * crt['full_shannon']:
            rabbit(
                f"at m={m}: entropy factorizes across prime factors to "
                f"within 2% — the dynamics respect CRT. This is a "
                f"functorial-factorization theorem candidate."
            )
        elif abs(crt['factorization_residual']) > 0.1 * crt['full_shannon']:
            rabbit(
                f"at m={m}: entropy DOES NOT factorize across CRT — residual "
                f"is {100 * crt['factorization_residual'] / crt['full_shannon']:.1f}%. "
                f"The CRT isomorphism at the group level is NOT preserved by "
                f"the Markov generator. Localize which prime breaks it."
            )
    return results


def section_6_page_curve(m: int, M: int) -> list[dict]:
    print_header(f"6. Page curve on residues at m = {m}")
    psi = prime_gas_ensemble(m, M=M, t_schedule=(4.0,))
    C = boundary_covariance(psi, time_idx=0)
    n = C.shape[0]
    pc = page_curve_on_residues(C, m)
    print("   |A|    mean S(ρ_A)    std          ratio_to_log|A|")
    print("   ---    -----------    ---------    --------------")
    for row in pc:
        print(
            f"   {row['subset_size']:3d}    "
            f"{row['mean_entropy']:11.4f}    {row['std_entropy']:.4f}       "
            f"{row['ratio_to_log_k']:14.4f}"
        )
    # Check: does the max sit at |A| = n/2 (Page peak) or elsewhere?
    peak = max(pc, key=lambda r: r['mean_entropy'])
    if peak['subset_size'] == n // 2:
        rabbit(
            f"at m={m}, Page curve peaks at |A| = n/2 = {n // 2} — canonical "
            f"Page behavior (maximum entanglement at half-system)"
        )
    else:
        rabbit(
            f"at m={m}, Page curve peaks at |A| = {peak['subset_size']} "
            f"(not n/2 = {n // 2}) — asymmetric subset dependence, "
            f"possible evidence of preferred residues"
        )
    return pc


def section_7_m6_prime_gas_signature(M: int) -> None:
    """Print the prime-gas (ℤ/6)^× character signature at t = 4.

    (Earlier revisions compared this signature against an external S_2
    simulator; the proprietary cross-check has been removed. The pure
    prime-gas readout is kept as-is — the character weights on the
    trivial and sign representations, plus a noise-floor leakage check.)
    """
    print_header(
        "7. Prime-gas (ℤ/6)^× character signature at t = 4"
    )
    psi = prime_gas_ensemble(6, M=M, t_schedule=(4.0,))
    C_primes = boundary_covariance(psi, time_idx=0)
    full = abelian_character_decomposition(C_primes, 6)
    # For (ℤ/6)^× ≅ ℤ/2 there are exactly two characters: trivial and sign.
    p_triv = float(full.p_char[0])
    p_sign = float(full.p_char[1]) if full.p_char.size > 1 else 0.0
    print(f"   p_triv             = {p_triv:.4f}")
    print(f"   p_sign             = {p_sign:.4f}")
    print(f"   direct entropy     = {full.direct_entropy:.4f}")
    print(f"   Shannon (chars)    = {full.shannon:.4f}")
    print(f"   leakage            = {full.nonabelian_leakage:.2e}")


def section_8_alpha_sweep_rabi_signature(M: int) -> None:
    print_header(
        "8. Rabi phase transition: character entropy vs α (shear parameter)"
    )
    m = 30
    n = group_order(m)
    alpha_ladder = (0.0, 0.25, 0.5, 0.75, 1.0, 1.239, 1.5, 2.0)
    print(f"   m = {m}  n = {n}  α_c ≈ √(135/88) ≈ 1.239")
    print("   α        S(ρ)       H(p)       leakage       top_p        localization")
    print("   ----    -------    -------    --------    ----------    ------------")
    top_p_vs_alpha = []
    for a in alpha_ladder:
        psi = prime_gas_ensemble(m, M=M, t_schedule=(4.0,), alpha=a)
        C = boundary_covariance(psi, time_idx=0)
        full = abelian_character_decomposition(C, m)
        top_p = float(np.max(full.p_char))
        # "Localization": how concentrated is p on its top character?
        # Inverse participation ratio on characters.
        ipr = float(np.sum(full.p_char ** 2))
        top_p_vs_alpha.append((a, top_p, full.shannon))
        print(
            f"   {a:5.3f}   {full.direct_entropy:7.4f}    "
            f"{full.shannon:7.4f}    {full.nonabelian_leakage:.2e}    "
            f"{top_p:10.4f}    {1.0/ipr:12.4f}"
        )
    # Sharp transition in top_p or Shannon?
    top_ps = np.array([t[1] for t in top_p_vs_alpha])
    derivs = np.abs(np.diff(top_ps))
    max_deriv_idx = int(np.argmax(derivs))
    rabbit(
        f"largest |Δtop_p| between α = {top_p_vs_alpha[max_deriv_idx][0]:.3f} "
        f"and α = {top_p_vs_alpha[max_deriv_idx + 1][0]:.3f}: "
        f"|Δ| = {derivs[max_deriv_idx]:.4f}. "
        f"If this straddles α_c = 1.239, the Rabi transition imprints on the "
        f"character-entropy decomposition."
    )


def section_9_spectral_fingerprints(primorials: tuple[int, ...], M: int) -> None:
    """Exotic: the distance matrix D_sym on the boundary IS G-equivariant by
    construction. Its eigenvalues are labeled by Dirichlet characters. What
    do those eigenvalues look like across the primorial ladder? This is the
    'character spectrum' of the prime gas at α=0.
    """
    print_header("9. Character spectrum of D_sym (the G-equivariant generator)")
    for m in primorials:
        n = group_order(m)
        if n > 24:
            continue
        D = build_distance_kernel_boundary(m)
        evals = np.sort(np.linalg.eigvalsh(D))
        print_sub(f"m = {m}  n = {n}")
        print(f"   eigenvalues of D_sym (sorted): {np.round(evals, 4).tolist()}")
        # Perron gap: λ_top - λ_second (in absolute magnitude)
        top_two = np.sort(np.abs(evals))[::-1][:2]
        gap = (top_two[0] - top_two[1]) / top_two[0]
        print(f"   relative Perron gap = {gap:.4f}")
        # Degeneracies — look for repeated eigenvalues (character
        # multiplicities in the actual spectrum)
        dedup = sorted(set(np.round(evals, 6)))
        mults = {v: int(np.sum(np.isclose(evals, v, atol=1e-5))) for v in dedup}
        degen = {v: mu for v, mu in mults.items() if mu > 1}
        print(f"   distinct eigvals = {len(dedup)}  (bound: n = {n})")
        if degen:
            print(f"   degeneracies     = {degen}")
            rabbit(
                f"at m={m}, D_sym has true eigenvalue degeneracies "
                f"{degen} — pairs of characters share an eigenvalue! "
                f"This is extra symmetry beyond (ℤ/m)^×; probably Galois/"
                f"complex-conjugation pairing χ ↔ χ̄."
            )


def section_10_dictionary_summary_table(primorials: tuple[int, ...], M: int) -> None:
    print_header("10. Dictionary summary — prime-gas rows")
    print("   system          G              |G|    bound   meas    ratio_to_logN")
    print("   -------------   ------------   ----   -----   -----   -------------")
    # Prime-gas rows
    for m in primorials:
        n = group_order(m)
        if n > 32:
            continue
        psi = prime_gas_ensemble(m, M=min(M, 20000), t_schedule=(4.0,))
        C = boundary_covariance(psi, time_idx=0)
        full = abelian_character_decomposition(C, m)
        G_name = f"(Z/{m})*"
        ratio = full.shannon / log(n) if n > 1 else 0.0
        print(f"   prime-gas m={m:4d} {G_name:12s}   {n:4d}   "
              f"   {n}    {full.shannon:.2f}     {ratio:.4f}")
    print()
    print("   (bound here is the Theorem-3.1 Schur count = # irrep types; "
          "for abelian G it equals |G|.)")


# ──────────────────────────────────────────────────────────────────────────
# §11. CRT-split comparison — honest version
#
# CORRECTION NOTE. A prior version of this script committed predictions
# against `crt_factorization_check.factorization_residual`. That quantity
# is NOT a valid CRT-factorization test — it returns −2.079 at uniform
# initialization, where a correct test must return 0. It was measuring
# a different thing (H_full minus a sum of H_p-subgroup-refinement
# shannons), not the cross-prime mutual information.
#
# The correct observable is the total correlation on the Dirichlet joint:
#     TC(m, t)  ≡  Σ_p H(p_p) − H(p_full)   ≥ 0
# where p(χ) is the Born weight on G-characters χ = (χ_{p_1},...,χ_{p_k}),
# and p_p is its marginal on the p-th factor (ℤ/p)^×. TC = 0 iff the
# character distribution factorizes across primes.
#
# PREDICTIONS COMMITTED IN ADVANCE (take 2, against the correct observable):
#
#   (P1) With the CRT-local generator D^CRT, TC(m, t) = 0 at all m, all t,
#        up to Monte-Carlo noise. Theoretical proof: under D^CRT the
#        relaxation amplitude on χ = Π χ_p is exp(−2βt Σ_p λ_p(χ_p)/λ_P),
#        a product over p ⇒ p(χ) = Π_p p_p(χ_p) exactly.
#        Tolerance: |TC| < 0.05 nats.
#
#   (P2) With natural D_sym on m = 30, TC(m, t) > 0 and grows in t. The
#        magnitude at t = 16 is TBD by measurement — the previous run's
#        −1.94 was the wrong observable and is not a valid prediction.
#        We measure TC honestly and let the number speak.
#
#   (P3) TC at m = 30 scales subadditively in t (bounded by log φ(m) minus
#        H_full at saturation).
#
# Outcome: if P1 passes, the CRT-local construction is verified and TC is
# a clean observable. The magnitude of TC under natural D_sym is the
# *correct* Arithmetic-RT-analog candidate — not the old residual.
# ──────────────────────────────────────────────────────────────────────────

def section_11_crt_split_comparison(
    primorials: tuple[int, ...], *, M: int, t_schedule: tuple[float, ...],
) -> None:
    print()
    print("=" * 78)
    print("  11. CRT-split comparison — natural D_sym vs CRT-local D^CRT")
    print("=" * 78)
    print()
    print("  Observable: total correlation TC = Σ H(p_p) − H(p_full)  on Dirichlet joint")
    print("  P1: CRT-local TC ≈ 0 at all m, all t   (|TC| < 0.05 nats)")
    print("  P2: Natural   TC > 0, grows with t       (measure, do not predict)")
    print("  P3: TC_nat(m=30, t=16) is the Arithmetic-RT candidate")
    print()

    summary_rows = []

    for m in primorials:
        if group_order(m) > 64:
            print(f"  -- m = {m}: group too large; skipping --")
            continue
        print(f"-- m = {m}   φ(m) = {group_order(m)} --")

        D_nat = build_distance_kernel_boundary(m)
        D_crt = build_crt_local_distance_kernel(m)
        fro_nat = float(np.linalg.norm(D_nat))
        fro_crt = float(np.linalg.norm(D_crt))
        fro_diff = float(np.linalg.norm(D_nat - D_crt))
        print(f"   ‖D_nat‖_F = {fro_nat:.3f}   ‖D_CRT‖_F = {fro_crt:.3f}   "
              f"‖D_nat − D_CRT‖_F = {fro_diff:.3f}")
        print()
        hdr = (f"   t       H_full_nat  H_full_crt   TC_nat    TC_crt    "
               f"ΔTC     per-prime marginals (nat | crt)")
        print(hdr)
        print("   " + "-" * (len(hdr) - 3))

        psi_nat = prime_gas_ensemble(
            m, M=M, t_schedule=t_schedule,
            rng=np.random.default_rng(hash(("nat", m)) & 0xFFFFFFFF),
        )
        psi_crt = prime_gas_ensemble(
            m, M=M, t_schedule=t_schedule, D_override=D_crt,
            rng=np.random.default_rng(hash(("crt", m)) & 0xFFFFFFFF),
        )

        for ti, t in enumerate(t_schedule):
            C_nat = boundary_covariance(psi_nat, time_idx=ti)
            C_crt = boundary_covariance(psi_crt, time_idx=ti)

            tc_nat = crt_character_marginal_check(m, C_nat)
            tc_crt = crt_character_marginal_check(m, C_crt)

            H_nat = tc_nat["H_full"]
            H_crt = tc_crt["H_full"]
            TC_nat = tc_nat["total_correlation"]
            TC_crt = tc_crt["total_correlation"]
            dTC = TC_nat - TC_crt

            marg_nat_s = ",".join(f"{tc_nat['H_marginals'][p]:.2f}" for p in tc_nat["factors"])
            marg_crt_s = ",".join(f"{tc_crt['H_marginals'][p]:.2f}" for p in tc_crt["factors"])

            print(
                f"   {t:5.2f}   {H_nat:8.4f}   {H_crt:8.4f}   "
                f"{TC_nat:+7.4f}   {TC_crt:+7.4f}   {dTC:+7.4f}   "
                f"[{marg_nat_s}] | [{marg_crt_s}]"
            )
            summary_rows.append({
                "m": m, "t": t,
                "H_nat": H_nat, "H_crt": H_crt,
                "TC_nat": TC_nat, "TC_crt": TC_crt,
                "dTC": dTC,
            })

        # Per-m verdict at the tail
        tail = summary_rows[-1]
        p1_pass = abs(tail["TC_crt"]) < 0.05
        print(f"   [m={m}] tail TC_crt = {tail['TC_crt']:+.4f} nats  "
              f"({'PASS' if p1_pass else 'FAIL'} P1, threshold 0.05)")
        print(f"   [m={m}] tail TC_nat = {tail['TC_nat']:+.4f} nats  "
              f"(natural cross-prime mutual information at t={tail['t']})")
        print(f"   [m={m}] ΔTC       = {tail['dTC']:+.4f} nats  "
              f"(how much extra coupling the natural generator carries)")
        print()

    # Cross-m summary
    print("   === TC at the ladder tail across primorials ===")
    print("   m       φ(m)   H_full_nat   H_full_crt   TC_nat      TC_crt      ΔTC")
    print("   ---     ----   ----------   ----------   --------    --------    --------")
    last_by_m: dict = {}
    for row in summary_rows:
        last_by_m[row["m"]] = row
    for m, row in last_by_m.items():
        print(f"   {m:5d}   {group_order(m):4d}   "
              f"{row['H_nat']:9.4f}    {row['H_crt']:9.4f}    "
              f"{row['TC_nat']:+7.4f}    {row['TC_crt']:+7.4f}    {row['dTC']:+7.4f}")
    print()
    print("   Interpretation:")
    print("     TC_crt ≈ 0 at all m  →  P1 holds, CRT-local factorization verified")
    print("     TC_nat > 0          →  natural generator carries cross-prime MI")
    print("     ΔTC = TC_nat − TC_crt is the cleanest candidate for the")
    print("     Arithmetic-Ryu-Takayanagi deficit: bulk information the dynamics")
    print("     forces across CRT boundaries that the group structure alone forbids.")


# ──────────────────────────────────────────────────────────────────────────
# §12. Asymptotic scaling — TC_nat across the full primorial well
#
# The m ∈ {6, 30, 210} data from §11 showed non-monotone TC_nat:
# 0.00 → 0.30 → 0.03. Two competing hypotheses:
#
#   H0 (Claude.AI "dynamical concentration"):
#       TC_nat → 0 as m → ∞. The cross-prime coupling is there in the
#       generator but invisible at fixed t because p(χ) concentrates on
#       the trivial character. At matched-entropy t* (where H_full ≈ const
#       across m), TC is non-decreasing in m.
#
#   H1 (Tony's "3-adic polarity resonance"):
#       TC_nat does NOT decay. At the primorial m = 30030 the 3-adic
#       valuation hierarchy flips sign (reported by v3_character_hunt),
#       and TC should hit a spike or phase-shift, not a decay.
#
#   H2 (null "noise"):
#       m=30's bump is a small-m artifact; the curve is essentially flat
#       and near zero everywhere.
#
# PREDICTIONS COMMITTED IN ADVANCE for this driver:
#
#   (D1) Peak TC over t (sweeping βt ∈ [0, 64]) at each m.
#        If H0:  peak TC stays roughly constant or grows slowly; tail
#                TC at large βt → 0.
#        If H1:  peak TC spikes at m = 30030 vs m = 2310.
#        If H2:  peak TC ≤ 0.05 everywhere.
#
#   (D2) Matched-entropy TC at H_full ≈ 1.2 nats for each m ≥ 30. We
#        find the βt where the natural kernel's full Shannon crosses
#        1.2 nats and report TC there. Call this TC*(m).
#        If H0:  TC*(m) is non-decreasing in m.
#        If H1:  TC*(m) spikes at 30030.
#        If H2:  TC*(m) ≤ 0.05 everywhere.
#
# Method: analytic, Monte-Carlo-free. For primorial m, λ(χ) is computed
# via the CRT-tensor contraction in `primorial_kernel_eigenvalue_tensor`.
# Born weights p(χ) = exp(-2βt λ(χ)/λ_P) / Z(t) exactly. No ensembles,
# no Wishart noise. O(n²) at worst; m = 30030 (n = 5760) runs in seconds.
# ──────────────────────────────────────────────────────────────────────────

def section_12_asymptotic_scaling(
    primorials: tuple[int, ...],
    beta_t_schedule: tuple[float, ...],
    matched_entropy_target: float = 1.2,
) -> None:
    print()
    print("=" * 78)
    print("  12. Asymptotic TC scaling — does the bulk topology persist at m=30030?")
    print("=" * 78)
    print()
    print("  H0: dynamical concentration (TC → 0 at fixed t as m → ∞)")
    print("  H1: 3-adic polarity resonance (TC spikes at m = 30030)")
    print("  H2: null (TC ≤ 0.05 everywhere)")
    print()
    print("  Observable: TC on Dirichlet joint, computed analytically "
          "(no MC noise).")
    print("  Predictions committed in code comments above.")
    print()

    summary = []
    for m in primorials:
        n = group_order(m)
        print(f"-- m = {m}   φ(m) = {n} --")
        # Natural kernel sweep
        rows_nat = analytic_tc_vs_time(m, beta_t_schedule, kernel="natural")
        rows_crt = analytic_tc_vs_time(m, beta_t_schedule, kernel="crt")

        peak_tc_nat = max(r["TC"] for r in rows_nat)
        peak_bt = rows_nat[int(np.argmax([r["TC"] for r in rows_nat]))]["beta_t"]
        peak_tc_crt = max(abs(r["TC"]) for r in rows_crt)
        tail_tc_nat = rows_nat[-1]["TC"]
        tail_H_nat = rows_nat[-1]["H_full"]
        tail_top = rows_nat[-1]["top_char_prob"]

        print(f"   peak TC_nat over βt ladder      = {peak_tc_nat:+.4f} nats"
              f"   (at βt = {peak_bt:.2f})")
        print(f"   peak |TC_crt| over βt ladder    = {peak_tc_crt:+.4e} nats"
              f"   (should be ~0)")
        print(f"   TC_nat at βt = {beta_t_schedule[-1]:.1f}           = "
              f"{tail_tc_nat:+.4f} nats   (H_full = {tail_H_nat:.3f}, "
              f"top p(χ) = {tail_top:.3f})")

        # Matched-entropy t* — find βt where H_full is closest to target
        H_gap = [abs(r["H_full"] - matched_entropy_target) for r in rows_nat]
        ti_star = int(np.argmin(H_gap))
        r_star = rows_nat[ti_star]
        tc_star = r_star["TC"]
        print(f"   matched entropy H_full ≈ {matched_entropy_target:.2f}: "
              f"βt = {r_star['beta_t']:.2f}, H_full = {r_star['H_full']:.3f}, "
              f"TC* = {tc_star:+.4f} nats")

        # Per-prime marginal Shannon at peak TC
        peak_row = rows_nat[int(np.argmax([r["TC"] for r in rows_nat]))]
        marg_str = ", ".join(
            f"p={p}:{h:.3f}" for p, h in peak_row["H_marginals"].items()
        )
        print(f"   per-prime marginals at peak βt:  {marg_str}")

        summary.append({
            "m": m, "n": n,
            "peak_tc_nat": peak_tc_nat, "peak_bt": peak_bt,
            "peak_tc_crt": peak_tc_crt,
            "tail_tc_nat": tail_tc_nat,
            "tc_star": tc_star,
            "bt_star": r_star["beta_t"],
            "H_star": r_star["H_full"],
        })
        print()

    # Cross-m summary
    print("   === Cross-m summary ===")
    print(f"   {'m':>6}  {'φ(m)':>5}  {'peak_TC':>9}  "
          f"{'peak_βt':>8}  {'tail_TC':>9}  {'TC*(H=1.2)':>11}  "
          f"{'peak|TC_crt|':>12}")
    for s in summary:
        print(f"   {s['m']:>6d}  {s['n']:>5d}  "
              f"{s['peak_tc_nat']:+8.4f}  {s['peak_bt']:>8.2f}  "
              f"{s['tail_tc_nat']:+8.4f}  {s['tc_star']:+10.4f}   "
              f"{s['peak_tc_crt']:.2e}")
    print()

    # Verdict
    peak_seq = [s["peak_tc_nat"] for s in summary]
    star_seq = [s["tc_star"] for s in summary]
    is_monotone_decay = all(peak_seq[i] >= peak_seq[i + 1]
                             for i in range(len(peak_seq) - 1))
    is_flat = max(peak_seq) < 0.05
    has_spike_at_30030 = False
    if len(summary) >= 2 and summary[-1]["m"] == 30030:
        has_spike_at_30030 = (
            summary[-1]["peak_tc_nat"] > 1.1 * summary[-2]["peak_tc_nat"]
            and summary[-1]["peak_tc_nat"] > 0.05
        )
    print("   Verdict on competing hypotheses:")
    print(f"     monotone-decay in peak TC?  {is_monotone_decay}")
    print(f"     flat (all peaks < 0.05)?    {is_flat}")
    print(f"     spike at m = 30030?         {has_spike_at_30030}")
    if is_flat:
        print("     ⇒ H2 (null) consistent with data")
    elif has_spike_at_30030:
        print("     ⇒ H1 (3-adic polarity resonance) supported")
    elif is_monotone_decay and summary[-1]["peak_tc_nat"] < 0.05:
        print("     ⇒ H0 (dynamical concentration) supported")
    else:
        print("     ⇒ none of H0/H1/H2 cleanly — non-monotone with no spike")


# ──────────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--M", default=5000, type=int)
    ap.add_argument("--M-haar", default=3000, type=int)
    ap.add_argument(
        "--primorials", default=",".join(str(m) for m in PRIMORIALS),
        help="comma-separated list of m values",
    )
    ap.add_argument(
        "--big-primorials", default="6,30,210,2310,30030",
        type=lambda s: [int(x) for x in s.split(",")],
        help="primorials for §12 analytic scaling sweep",
    )
    ap.add_argument(
        "--big-bt-ladder",
        default="0.01,0.03,0.1,0.3,1.0,3.0,10.0,30.0,100.0,300.0",
        help="βt values for analytic §12 sweep",
    )
    ap.add_argument(
        "--matched-entropy", default=1.2, type=float,
        help="target H_full in nats for matched-entropy TC comparison",
    )
    args = ap.parse_args()
    primorials = tuple(int(x) for x in args.primorials.split(","))

    t0 = time.time()
    print()
    print("╔══════════════════════════════════════════════════════════════════════════╗")
    print("║   HOLOGRAPHIC PRIME GAS — dictionary row 2 test                          ║")
    print("║   Character-entropy decomposition on (Z/m)^×                             ║")
    print("╚══════════════════════════════════════════════════════════════════════════╝")
    print(f"   M = {args.M}   primorials = {primorials}")

    section_1_character_tables(primorials)
    section_2_character_entropy_primorial_ladder(
        primorials, M=args.M, t_schedule=T_LADDER,
    )
    for m in primorials:
        if group_order(m) <= 24:
            section_3_subgroup_refinement_tower(m, M=args.M)
    section_4_haar_baseline(primorials, M=args.M_haar)
    section_5_crt_factorization(primorials, M=args.M)
    for m in primorials:
        if group_order(m) <= 20:
            section_6_page_curve(m, M=args.M)
    section_7_m6_prime_gas_signature(M=args.M)
    section_8_alpha_sweep_rabi_signature(M=args.M)
    section_9_spectral_fingerprints(primorials, M=args.M)
    section_10_dictionary_summary_table(primorials, M=args.M)
    # §11: MC-based comparison — restricted to m where the group is small enough
    crt_primorials = tuple(m for m in primorials if group_order(m) <= 48)
    section_11_crt_split_comparison(
        crt_primorials, M=args.M, t_schedule=T_LADDER,
    )
    # §12: analytic asymptotic scaling — pushes to m=30030
    big_primorials = tuple(args.big_primorials)
    big_bt_ladder = tuple(float(x) for x in args.big_bt_ladder.split(","))
    section_12_asymptotic_scaling(
        big_primorials, big_bt_ladder,
        matched_entropy_target=args.matched_entropy,
    )

    dt = time.time() - t0
    print()
    print("=" * 78)
    print(f"  total wall time: {dt:.1f}s")
    print("  grep 'RABBIT:' in this output for the leads worth chasing.")
    print("=" * 78)


if __name__ == "__main__":
    main()
