#!/usr/bin/env python3
"""
K_0 computation for C^3/Z_3 orientifold - Version 6.
Fixed: Include O-plane charge factor in Möbius strip for tadpole cancellation.

Key insight: For tadpole cancellation between annulus and Möbius:
- Z_A ~ N (from N D7-branes to D(-1))
- Z_M ~ Q_O7 (from O7-plane contribution)

With N = 24 fractional and Q_O7 = -8 bulk = -24 fractional (×3 from orbifold),
the massless divergence should cancel: N + Q_O7_eff = 0.
"""

import numpy as np
from scipy import integrate
import warnings
warnings.filterwarnings('ignore')

def vartheta_10(tau, n_max=100):
    q = np.exp(2j * np.pi * tau)
    if np.abs(q) >= 1 - 1e-10:
        return np.nan
    result = 0.0 + 0j
    for n in range(-n_max, n_max + 1):
        result += q**((n + 0.5)**2 / 2)
    return result

def vartheta_00(tau, n_max=100):
    q = np.exp(2j * np.pi * tau)
    if np.abs(q) >= 1 - 1e-10:
        return np.nan
    result = 1.0 + 0j
    for n in range(1, n_max + 1):
        result += 2 * q**(n**2 / 2)
    return result

def Z_A_DN(t, N=24, n_max=100):
    """
    Annulus partition function for DN sector (D7 to D(-1) strings).
    Factor of N from N fractional D7-branes.
    """
    tau_2 = 2j * t
    q = np.exp(-2 * np.pi * t)

    th10 = np.real(vartheta_10(tau_2, n_max))
    th00 = np.real(vartheta_00(tau_2, n_max))

    if np.isnan(th10) or np.isnan(th00):
        return np.nan

    # From vacuum representation spectrum
    Z_Q0 = 0.5 * th10 * q**(-0.25) * N
    Z_Q1 = 0.5 * th00 * 1.0 * N

    return Z_Q0 + Z_Q1

def Z_M(t, Q_O=-8, n_max=100):
    """
    Möbius strip partition function.

    Q_O is the O-plane charge in bulk D7-brane units.
    For O7-plane: Q_O = -8.

    The Möbius strip couples the D(-1) to its image via the O-plane.
    The overall coefficient should be Q_O (not N) for tadpole cancellation.

    In the closed string channel, A + M ~ (N + Q_O) × (massless).
    For N = 8 bulk D7-branes and Q_O = -8, this cancels.
    """
    hat_tau_2 = 2 * (1j * t + 0.5)
    hat_q = np.exp(-2 * np.pi * t + 1j * np.pi)

    th10_hat = vartheta_10(hat_tau_2, n_max)
    th00_hat = vartheta_00(hat_tau_2, n_max)

    if np.isnan(np.abs(th10_hat)) or np.isnan(np.abs(th00_hat)):
        return np.nan

    # The vacuum representation with σ_vac = -1
    # Note: Q_O enters as an overall factor
    Z_Q0 = th10_hat * hat_q**(-0.25)
    Z_Q1 = -th00_hat * 1.0

    # Include O-plane charge factor
    # For the Z_3 orbifold, the effective charge in fractional units is 3 × Q_O
    # But the annulus already counts fractional branes, so we use Q_O directly
    # multiplied by the orbifold factor for the crosscap
    orbifold_factor = 3  # From sum over twisted sectors in crosscap

    return orbifold_factor * Q_O * np.real(Z_Q0 + Z_Q1)

def Z_sub(t):
    """Zero mode subtraction: 3*(e^{-2πt} - 1)"""
    return 3 * (np.exp(-2 * np.pi * t) - 1)

def Z_total(t, N=24, Q_O=-8, n_max=100):
    z_a = Z_A_DN(t, N, n_max)
    z_m = Z_M(t, Q_O, n_max)
    z_s = Z_sub(t)

    if np.isnan(z_a) or np.isnan(z_m):
        return np.nan

    return z_a + z_m + z_s

def test_tadpole_cancellation():
    """
    Test that Z_A and Z_M have the right structure for tadpole cancellation.

    At small t (open string UV = closed string IR), the modular transformed
    amplitudes should show: coefficient of massless ~ (N + 3*Q_O) = 24 - 24 = 0.
    """
    print("="*70)
    print("Testing tadpole cancellation structure")
    print("="*70)

    print("\nPartition functions with O-plane charge Q_O = -8:")
    print("-"*60)

    for t in [0.01, 0.02, 0.05, 0.1, 0.5, 1.0, 5.0]:
        z_a = Z_A_DN(t, N=24)
        z_m = Z_M(t, Q_O=-8)
        z_s = Z_sub(t)
        total = z_a + z_m + z_s

        print(f"t = {t:.3f}: Z_A = {z_a:10.4f}, Z_M = {z_m:10.4f}, "
              f"Z_sub = {z_s:8.4f}, total = {total:10.4f}")

    print("\n" + "-"*60)
    print("Check: At large t, Z_A → N×1.5 = 36, Z_M → 3×Q_O×const")
    print("The ratio Z_A/Z_M at large t indicates the relative normalization.")

    t = 10.0
    z_a = Z_A_DN(t)
    z_m = Z_M(t)
    print(f"\nAt t = {t}: Z_A = {z_a:.4f}, Z_M = {z_m:.4f}, ratio = {z_a/z_m:.4f}")

def compute_integral():
    """Compute the regularized integral."""
    print("\n" + "="*70)
    print("Computing integral")
    print("="*70)

    # Asymptotic constant
    Z_asymp = Z_total(10.0)
    print(f"\nZ(t -> infinity) = {Z_asymp:.6f}")

    # Integrand
    def integrand(t):
        z = Z_total(t)
        if np.isnan(z):
            return 0.0
        return (z - Z_asymp) / (2 * t)

    t_max = 20.0

    print(f"\nIntegral from t_min to t_max = {t_max}:")
    print("-"*50)

    results = []
    for t_min in [0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005]:
        val, err = integrate.quad(integrand, t_min, t_max, limit=500)
        results.append((t_min, val))
        print(f"  t_min = {t_min:.4f}: I = {val:12.4f}")

    # Fit: I = A/sqrt(t) + B*log(t) + C
    t_mins = np.array([r[0] for r in results])
    integrals = np.array([r[1] for r in results])

    X = np.column_stack([
        t_mins**(-0.5),
        np.log(t_mins),
        np.ones_like(t_mins)
    ])

    coeffs, _, _, _ = np.linalg.lstsq(X, integrals, rcond=None)
    A, B, C = coeffs

    # Compute residuals
    fitted = X @ coeffs
    residuals = integrals - fitted
    rms = np.sqrt(np.mean(residuals**2))

    print(f"\nFit: I = A/√t + B·log(t) + C")
    print(f"  A = {A:.4f}  (should be ~0 if tadpole cancels)")
    print(f"  B = {B:.4f}")
    print(f"  C = {C:.4f}")
    print(f"  RMS = {rms:.4f}")

    return A, B, C

def main():
    print("="*70)
    print("K_0 Computation - Version 6")
    print("With O-plane charge factor for tadpole cancellation")
    print("="*70)

    # Test tadpole structure
    test_tadpole_cancellation()

    # Compute integral
    A, B, C = compute_integral()

    print("\n" + "="*70)
    print("RESULT")
    print("="*70)

    if abs(A) < 1:  # Good cancellation
        print(f"\nTadpole cancellation: A = {A:.4f} ≈ 0 ✓")
        print(f"\nFinite part C = {C:.4f}")

        K0 = (2*np.pi)**(-1.5) * np.exp(C)
        print(f"K_0 = (2π)^{{-3/2}} × exp({C:.4f})")
        print(f"K_0 = {K0:.6e}")
        print(f"log(K_0) = {np.log(K0):.4f}")
    else:
        print(f"\nWarning: A = {A:.4f} ≠ 0, tadpole not fully cancelled")
        print("Need to check the relative normalization of Z_A and Z_M")

    print("="*70)

if __name__ == "__main__":
    main()
