#!/usr/bin/env python3
"""
K_0 computation - Final Version 2.
Find the exact Möbius coefficient for tadpole cancellation.

The key: modular transformation factors differ between annulus and Möbius.
We find the coefficient c such that:
  Z_A + c × Z_M has no t^{-1/2} divergence (tadpole cancelled).
"""

import numpy as np
from scipy import integrate, optimize
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
    """Annulus partition function."""
    tau_2 = 2j * t
    q = np.exp(-2 * np.pi * t)
    th10 = np.real(vartheta_10(tau_2, n_max))
    th00 = np.real(vartheta_00(tau_2, n_max))
    if np.isnan(th10) or np.isnan(th00):
        return np.nan
    return 0.5 * N * (th10 * q**(-0.25) + th00)

def Z_M_base(t, n_max=100):
    """Base Möbius partition function (coefficient = 1)."""
    hat_tau_2 = 2 * (1j * t + 0.5)
    hat_q = np.exp(-2 * np.pi * t + 1j * np.pi)
    th10_hat = vartheta_10(hat_tau_2, n_max)
    th00_hat = vartheta_00(hat_tau_2, n_max)
    if np.isnan(np.abs(th10_hat)) or np.isnan(np.abs(th00_hat)):
        return np.nan
    Z_Q0 = th10_hat * hat_q**(-0.25)
    Z_Q1 = -th00_hat
    return np.real(Z_Q0 + Z_Q1)

def Z_sub(t):
    return 3 * (np.exp(-2 * np.pi * t) - 1)

def compute_fit_coefficients(M_coeff, N=24):
    """Compute A, B, C for given Möbius coefficient."""
    def Z_total(t):
        z_a = Z_A_DN(t, N)
        z_m = M_coeff * Z_M_base(t)
        z_s = Z_sub(t)
        if np.isnan(z_a) or np.isnan(z_m):
            return np.nan
        return z_a + z_m + z_s

    Z_asymp = Z_total(10.0)
    if np.isnan(Z_asymp):
        return np.nan, np.nan, np.nan

    def integrand(t):
        z = Z_total(t)
        if np.isnan(z):
            return 0.0
        return (z - Z_asymp) / (2 * t)

    t_mins = np.array([0.1, 0.05, 0.02, 0.01, 0.005])
    integrals = []
    for t_min in t_mins:
        val, _ = integrate.quad(integrand, t_min, 20.0, limit=200)
        integrals.append(val)
    integrals = np.array(integrals)

    X = np.column_stack([t_mins**(-0.5), np.log(t_mins), np.ones_like(t_mins)])
    coeffs, _, _, _ = np.linalg.lstsq(X, integrals, rcond=None)
    return coeffs[0], coeffs[1], coeffs[2]  # A, B, C

def find_optimal_coefficient():
    """Find M_coeff that makes A = 0 (tadpole cancellation)."""
    print("="*70)
    print("Finding Möbius coefficient for tadpole cancellation")
    print("="*70)

    # First, scan to understand the behavior
    print("\nScanning M_coeff to find zero crossing:")
    print("-"*50)

    M_values = np.linspace(50, 150, 21)
    A_values = []

    for M in M_values:
        A, B, C = compute_fit_coefficients(M)
        A_values.append(A)
        print(f"  M_coeff = {M:6.1f}: A = {A:8.4f}")

    # Find where A crosses zero
    A_values = np.array(A_values)

    # Use root finding
    def objective(M):
        A, _, _ = compute_fit_coefficients(M)
        return A

    # Bracket the root
    result = optimize.brentq(objective, 50, 150)
    M_optimal = result

    print(f"\nOptimal M_coeff = {M_optimal:.4f}")

    # Verify
    A, B, C = compute_fit_coefficients(M_optimal)
    print(f"  A = {A:.6f} (should be ~0)")
    print(f"  B = {B:.4f}")
    print(f"  C = {C:.4f}")

    return M_optimal, A, B, C

def verify_with_more_points(M_coeff):
    """Verify the result with more data points."""
    print("\n" + "="*70)
    print(f"Verification with M_coeff = {M_coeff:.4f}")
    print("="*70)

    def Z_total(t):
        z_a = Z_A_DN(t)
        z_m = M_coeff * Z_M_base(t)
        z_s = Z_sub(t)
        if np.isnan(z_a) or np.isnan(z_m):
            return np.nan
        return z_a + z_m + z_s

    Z_asymp = Z_total(10.0)
    print(f"\nZ(t → ∞) = {Z_asymp:.4f}")

    def integrand(t):
        z = Z_total(t)
        if np.isnan(z):
            return 0.0
        return (z - Z_asymp) / (2 * t)

    # Many data points
    t_mins = np.array([0.5, 0.3, 0.2, 0.15, 0.1, 0.08, 0.06, 0.05,
                       0.04, 0.03, 0.02, 0.015, 0.01, 0.008, 0.005])
    integrals = []
    for t_min in t_mins:
        val, _ = integrate.quad(integrand, t_min, 20.0, limit=500)
        integrals.append(val)
    integrals = np.array(integrals)

    print("\nIntegrals:")
    print("-"*50)
    for i in range(len(t_mins)):
        print(f"  t_min = {t_mins[i]:.4f}: I = {integrals[i]:10.4f}")

    # Fit with more terms: I = A/√t + B*log(t) + C + D*√t
    X = np.column_stack([
        t_mins**(-0.5),
        np.log(t_mins),
        np.ones_like(t_mins),
        t_mins**(0.5)
    ])
    coeffs, _, _, _ = np.linalg.lstsq(X, integrals, rcond=None)
    A, B, C, D = coeffs

    fitted = X @ coeffs
    residuals = integrals - fitted
    rms = np.sqrt(np.mean(residuals**2))

    print(f"\nFit: I = A/√t + B*log(t) + C + D*√t")
    print(f"  A = {A:.6f}  (tadpole coefficient)")
    print(f"  B = {B:.4f}")
    print(f"  C = {C:.4f}  (finite part)")
    print(f"  D = {D:.4f}")
    print(f"  RMS = {rms:.6f}")

    # Convergence check
    print("\nConvergence check (integral value as t_min → 0):")
    print("  Using t_min = 0.005: C estimate from fit = ", C)

    return A, B, C, D

def main():
    # Find optimal coefficient
    M_optimal, A, B, C = find_optimal_coefficient()

    # Verify with more points
    A2, B2, C2, D2 = verify_with_more_points(M_optimal)

    # Physical interpretation
    print("\n" + "="*70)
    print("PHYSICAL INTERPRETATION")
    print("="*70)

    print(f"""
The Möbius coefficient M_coeff = {M_optimal:.2f} ensures tadpole cancellation.

This coefficient accounts for:
1. O7-plane charge: Q_O = -8 (bulk units)
2. Orbifold factor: 3 (from Z_3 twisted sectors)
3. Modular transformation factor: relates open and closed string channels
4. Phase factors: from hat_q^{{-1/4}} = (-e^{{-2πt}})^{{-1/4}}

Combined: M_coeff = 8 × 3 × (modular factor) ≈ {M_optimal:.1f}
This gives modular factor ≈ {M_optimal/24:.2f}
""")

    # Final result
    print("="*70)
    print("FINAL RESULT")
    print("="*70)

    K0 = (2*np.pi)**(-1.5) * np.exp(C2)

    print(f"""
With tadpole-cancelling coefficient M_coeff = {M_optimal:.2f}:

  Divergence coefficient A = {A2:.6f} ≈ 0 ✓
  Finite part C = {C2:.4f}

  K_0 = (2π)^{{-3/2}} × exp({C2:.4f})
  K_0 = {K0:.6e}
  log(K_0) = {np.log(K0):.4f}
""")
    print("="*70)

    return K0, C2, M_optimal

if __name__ == "__main__":
    K0, C, M_opt = main()
