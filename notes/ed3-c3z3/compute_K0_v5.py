#!/usr/bin/env python3
"""
K_0 computation for C^3/Z_3 orientifold - Version 5.
Closed string channel computation with explicit tadpole cancellation.

Key insight: Transform to closed string channel where:
- Annulus: t_open -> l = 1/(2t) closed
- Mobius: t_open -> l = 1/(8t) closed (standard convention)

In closed string channel, tadpole cancellation is manifest as the
cancellation of the l -> infinity (massless) divergence.
"""

import numpy as np
from scipy import integrate
import warnings
warnings.filterwarnings('ignore')

def eta(tau, n_max=100):
    """Dedekind eta function: eta(tau) = q^{1/24} prod_{n>=1}(1-q^n)"""
    q = np.exp(2j * np.pi * tau)
    if np.abs(q) >= 1 - 1e-10:
        return np.nan
    result = q**(1/24)
    for n in range(1, n_max + 1):
        result *= (1 - q**n)
    return result

def vartheta_00(tau, n_max=100):
    """theta_00(tau) = sum_n q^{n^2/2}"""
    q = np.exp(2j * np.pi * tau)
    if np.abs(q) >= 1 - 1e-10:
        return np.nan
    result = 1.0 + 0j
    for n in range(1, n_max + 1):
        result += 2 * q**(n**2 / 2)
    return result

def vartheta_01(tau, n_max=100):
    """theta_01(tau) = sum_n (-1)^n q^{n^2/2}"""
    q = np.exp(2j * np.pi * tau)
    if np.abs(q) >= 1 - 1e-10:
        return np.nan
    result = 1.0 + 0j
    for n in range(1, n_max + 1):
        result += 2 * ((-1)**n) * q**(n**2 / 2)
    return result

def vartheta_10(tau, n_max=100):
    """theta_10(tau) = sum_n q^{(n+1/2)^2/2}"""
    q = np.exp(2j * np.pi * tau)
    if np.abs(q) >= 1 - 1e-10:
        return np.nan
    result = 0.0 + 0j
    for n in range(-n_max, n_max + 1):
        result += q**((n + 0.5)**2 / 2)
    return result

def Z_A_closed(l, N=24, n_max=100):
    """
    Annulus amplitude in closed string channel.

    The closed string channel has tau = il (purely imaginary).
    The amplitude includes the modular transformation factor.

    For D7-D(-1) strings on C^3/Z_3:
    Z_A^{closed}(l) = N * [theta functions at tau=il] * l^{-2}

    The l^{-2} comes from the transformation of the measure.
    """
    if l < 1e-12:
        return 0.0

    tau = 1j * l
    q = np.exp(-2 * np.pi * l)

    th00 = np.real(vartheta_00(tau, n_max))
    th01 = np.real(vartheta_01(tau, n_max))
    th10 = np.real(vartheta_10(tau, n_max))
    eta_val = np.real(eta(tau, n_max))

    if np.isnan(th00) or np.isnan(eta_val) or eta_val == 0:
        return np.nan

    # For 4 transverse directions (D7 to D(-1)):
    # Z_A ~ N * (th00^4 - th01^4 - th10^4) / eta^{12}
    # But with GSO and orientifolding, the structure changes.

    # Simple model: Z_A ~ N * (massive modes only)
    # At large l (IR): the massless modes give divergence ~ l

    # For massive modes only (subtracting tachyon and massless):
    # Z_A^{massive} = Z_A - (massless contribution)

    # The massless contribution at large l goes as:
    # Z_massless ~ N * (constant from oscillator zero modes)

    # For now, use a simplified form based on the spectrum:
    # 4 bosons (transverse) + 2 fermions from broken SUSY

    # Using the SO(2,2) character formula:
    # Z = (1/eta^4) * [massless + massive]

    # For D7-D(-1) strings in the DN sector:
    # The partition function in closed channel is:
    Z_bosons = (th00**4) / (eta_val**8)  # 4 transverse directions
    Z_fermions = (th01**4 - th10**4) / (eta_val**8)  # R-NS sectors

    # Total with N D7-branes:
    Z_total = N * np.real(Z_bosons + Z_fermions) / 2

    return Z_total

def Z_M_closed(l, n_max=100):
    """
    Mobius strip amplitude in closed string channel.

    For the O7-plane contribution:
    Z_M^{closed}(l) = +-32 * [hatted theta functions]

    The sign and coefficient come from the O-plane charge.
    """
    if l < 1e-12:
        return 0.0

    # Mobius strip modulus: hat_tau = i*l + 1/2
    hat_tau = 1j * l + 0.5

    th00 = vartheta_00(hat_tau, n_max)
    th01 = vartheta_01(hat_tau, n_max)
    th10 = vartheta_10(hat_tau, n_max)
    eta_val = eta(hat_tau, n_max)

    if np.isnan(np.abs(th00)) or np.isnan(np.abs(eta_val)) or np.abs(eta_val) < 1e-100:
        return np.nan

    # O7-plane contribution (with charge Q_O7 = -8 in D7-brane units):
    Z_bosons = (th00**4) / (eta_val**8)
    Z_fermions = (th01**4 - th10**4) / (eta_val**8)

    # The prefactor is related to O-plane charge
    # For O7-plane: Q_O7 = -8, contributing -8 to the tadpole
    Z_total = -8 * np.real(Z_bosons + Z_fermions) / 2

    return Z_total

def test_tadpole_cancellation():
    """
    Test that the massless closed string tadpole cancels.

    At large l (IR), the amplitude should behave as:
    A + M ~ (N - 8*k) * l + O(1)

    where k depends on normalization. For N = 24, we need k = 3.
    """
    print("Testing tadpole cancellation in closed string channel:")
    print("-"*60)

    for l in [0.5, 1.0, 2.0, 5.0, 10.0]:
        z_a = Z_A_closed(l)
        z_m = Z_M_closed(l)
        total = z_a + z_m

        # The massless tadpole is proportional to l at large l
        tadpole_coeff = total / l if l > 0 else 0

        print(f"l = {l:5.1f}: Z_A = {z_a:12.4f}, Z_M = {z_m:12.4f}, "
              f"total = {total:12.4f}, total/l = {tadpole_coeff:10.4f}")

def Z_total_openstring(t, N=24, n_max=100):
    """
    Total partition function in open string channel.
    Recomputed carefully following the conventions.
    """
    if t < 1e-12:
        return 0.0

    tau = 1j * t
    tau_2 = 2j * t  # Doubled modulus for DN sector
    q = np.exp(-2 * np.pi * t)

    # Annulus DN sector
    th10 = np.real(vartheta_10(tau_2, n_max))
    th00 = np.real(vartheta_00(tau_2, n_max))

    if np.isnan(th10) or np.isnan(th00):
        return np.nan

    # From the vacuum representation spectrum:
    # h=0, Q=0: contributes theta_10 * q^{-1/4}
    # h=1/2, Q=1: contributes theta_00 * 1 (with sign from n and (-1)^Q)
    Z_A = 0.5 * th10 * q**(-0.25) * N + 0.5 * th00 * N

    # Mobius strip
    hat_tau_2 = 2 * (1j * t + 0.5)
    hat_q = np.exp(-2 * np.pi * t + 1j * np.pi)

    th10_hat = vartheta_10(hat_tau_2, n_max)
    th00_hat = vartheta_00(hat_tau_2, n_max)

    if np.isnan(np.abs(th10_hat)) or np.isnan(np.abs(th00_hat)):
        return np.nan

    Z_M = np.real(th10_hat * hat_q**(-0.25) - th00_hat)

    # Zero mode subtraction term
    Z_sub = 3 * (np.exp(-2 * np.pi * t) - 1)

    return Z_A + Z_M + Z_sub

def analyze_small_t_behavior():
    """
    Analyze the small-t (UV in open string = IR in closed string) behavior.
    The modular transformation relates:
    Z_open(t) ~ Z_closed(l) where l = 1/(2t) or similar
    """
    print("\nAnalyzing small-t behavior (open string UV / closed string IR):")
    print("-"*60)

    # At small t, use modular transformation:
    # theta_ab(tau) = (-i*tau)^{-1/2} * theta_{ba}(-1/tau)

    for t in [0.01, 0.02, 0.05, 0.1, 0.2, 0.5]:
        z_total = Z_total_openstring(t)

        # Compute asymptotic form using modular transformation
        # At small t: theta_10(2it) ~ (2t)^{-1/2} * theta_01(i/(2t))
        #           theta_00(2it) ~ (2t)^{-1/2} * theta_00(i/(2t))

        # Closed string limit:
        l = 1/(4*t)  # closed string modulus

        tau_cl = 1j * l
        th01_cl = vartheta_01(tau_cl, 100)
        th00_cl = vartheta_00(tau_cl, 100)

        # Asymptotic prediction: Z ~ N * (2t)^{-1/2} * [combination of closed string thetas]
        z_asymp = 24 * (2*t)**(-0.5) * np.abs(th00_cl)

        print(f"t = {t:.3f}: Z_actual = {z_total:10.4f}, Z_asymp ~ N*(2t)^{{-1/2}} = {z_asymp:10.4f}, "
              f"ratio = {z_total/z_asymp if z_asymp != 0 else np.nan:.4f}")

def compute_with_proper_measure():
    """
    Compute the integral with proper measure and regularization.

    The key insight: the divergence at small t corresponds to IR
    in the closed string channel. Tadpole cancellation should
    remove this divergence.

    But we're computing for D(-1) (point-like), where there's no
    volume factor to cancel. The integral may have leftover divergence
    from non-compact directions.
    """
    print("\n" + "="*70)
    print("Computing integral with physical cutoffs")
    print("="*70)

    # Asymptotic constant (zero mode contribution)
    Z_asymp = Z_total_openstring(10.0)
    print(f"\nAsymptotic constant Z(t->infinity) = {Z_asymp:.6f}")

    # The regulated integral
    def integrand_reg(t):
        z = Z_total_openstring(t)
        if np.isnan(z):
            return 0.0
        return (z - Z_asymp) / (2 * t)

    # Compute for various cutoffs
    t_max = 20.0

    print(f"\nIntegral from t_min to t_max = {t_max}:")
    print("-"*50)

    results = []
    for t_min in [0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005]:
        integral, error = integrate.quad(integrand_reg, t_min, t_max, limit=500)
        results.append((t_min, integral))
        print(f"  t_min = {t_min:.4f}: I = {integral:12.4f}")

    # Fit divergence structure: I = A/sqrt(t) + B*log(t) + C
    t_mins = np.array([r[0] for r in results])
    integrals = np.array([r[1] for r in results])

    X = np.column_stack([
        t_mins**(-0.5),
        np.log(t_mins),
        np.ones_like(t_mins)
    ])

    coeffs, _, _, _ = np.linalg.lstsq(X, integrals, rcond=None)
    A, B, C = coeffs

    print(f"\nFit: I = {A:.4f}/sqrt(t) + {B:.4f}*log(t) + {C:.4f}")

    # The physical meaning:
    # A/sqrt(t) term: massless closed string exchange (should cancel with tadpole)
    # B*log(t) term: marginal deformation or threshold correction
    # C: finite part

    # For a well-defined orientifold, A should vanish (tadpole cancellation)
    # But we're getting non-zero A because we're computing only part of the amplitude

    # The full amplitude should include:
    # 1. DN sector (what we computed) - strings from D7s to ED3
    # 2. Contribution from the O-plane directly to ED3

    print(f"\nPhysical interpretation:")
    print(f"  A = {A:.4f}: massless tadpole (should vanish with full amplitude)")
    print(f"  B = {B:.4f}: logarithmic correction")
    print(f"  C = {C:.4f}: finite part (K_0 normalization)")

    # The Mobius strip we computed couples ED3 to the O-plane via D7 strings,
    # but there should also be a direct O-plane contribution.

    # For now, take C as the answer with caveat about missing terms
    K0_estimate = (2 * np.pi)**(-1.5) * np.exp(C)

    print(f"\n" + "="*70)
    print(f"RESULT (with caveats):")
    print(f"  Finite part C = {C:.4f}")
    print(f"  K_0 = (2pi)^{{-3/2}} * exp({C:.4f}) = {K0_estimate:.6e}")
    print(f"  log(K_0) = {np.log(K0_estimate):.4f}")
    print(f"\n  Note: A != 0 indicates incomplete tadpole cancellation.")
    print(f"  The full amplitude may include additional sectors.")
    print("="*70)

    return C

def main():
    print("="*70)
    print("K_0 Computation for C^3/Z_3 Orientifold - Version 5")
    print("Closed string analysis and tadpole structure")
    print("="*70)

    # Test tadpole in closed string channel
    test_tadpole_cancellation()

    # Analyze small-t behavior
    analyze_small_t_behavior()

    # Compute with proper measure
    C_finite = compute_with_proper_measure()

    return C_finite

if __name__ == "__main__":
    main()
