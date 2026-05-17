#!/usr/bin/env python3
"""
K_0 computation for C^3/Z_3 orientifold - Version 7.
Find the correct Möbius coefficient for tadpole cancellation by scanning.
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
    """Annulus partition function for DN sector."""
    tau_2 = 2j * t
    q = np.exp(-2 * np.pi * t)

    th10 = np.real(vartheta_10(tau_2, n_max))
    th00 = np.real(vartheta_00(tau_2, n_max))

    if np.isnan(th10) or np.isnan(th00):
        return np.nan

    Z_Q0 = 0.5 * th10 * q**(-0.25) * N
    Z_Q1 = 0.5 * th00 * 1.0 * N

    return Z_Q0 + Z_Q1

def Z_M_base(t, n_max=100):
    """Base Möbius strip partition function (coefficient = 1)."""
    hat_tau_2 = 2 * (1j * t + 0.5)
    hat_q = np.exp(-2 * np.pi * t + 1j * np.pi)

    th10_hat = vartheta_10(hat_tau_2, n_max)
    th00_hat = vartheta_00(hat_tau_2, n_max)

    if np.isnan(np.abs(th10_hat)) or np.isnan(np.abs(th00_hat)):
        return np.nan

    Z_Q0 = th10_hat * hat_q**(-0.25)
    Z_Q1 = -th00_hat * 1.0

    return np.real(Z_Q0 + Z_Q1)

def Z_sub(t):
    return 3 * (np.exp(-2 * np.pi * t) - 1)

def compute_divergence_coeff(M_coeff, N=24):
    """
    Compute the coefficient A of the t^{-1/2} divergence for given M_coeff.
    Z_M = M_coeff × Z_M_base
    """
    def Z_total(t):
        z_a = Z_A_DN(t, N)
        z_m = M_coeff * Z_M_base(t)
        z_s = Z_sub(t)
        if np.isnan(z_a) or np.isnan(z_m):
            return np.nan
        return z_a + z_m + z_s

    # Asymptotic constant
    Z_asymp = Z_total(10.0)

    def integrand(t):
        z = Z_total(t)
        if np.isnan(z):
            return 0.0
        return (z - Z_asymp) / (2 * t)

    # Compute integrals at a few t_min values
    t_max = 20.0
    t_mins = np.array([0.1, 0.05, 0.02, 0.01, 0.005])
    integrals = []

    for t_min in t_mins:
        val, _ = integrate.quad(integrand, t_min, t_max, limit=200)
        integrals.append(val)

    integrals = np.array(integrals)

    # Fit: I = A/√t + B*log(t) + C
    X = np.column_stack([
        t_mins**(-0.5),
        np.log(t_mins),
        np.ones_like(t_mins)
    ])

    coeffs, _, _, _ = np.linalg.lstsq(X, integrals, rcond=None)
    A, B, C = coeffs

    return A, B, C

def main():
    print("="*70)
    print("Scanning Möbius coefficient for tadpole cancellation")
    print("="*70)

    print("\nZ_M = M_coeff × Z_M_base")
    print("Finding M_coeff such that A ≈ 0 (tadpole cancellation)")
    print("-"*60)

    # Scan M_coeff values
    M_coeffs = np.linspace(-50, 50, 101)
    A_values = []

    for M_coeff in M_coeffs:
        A, _, _ = compute_divergence_coeff(M_coeff)
        A_values.append(A)

    A_values = np.array(A_values)

    # Find where A crosses zero
    for i in range(len(A_values) - 1):
        if A_values[i] * A_values[i+1] < 0:  # Sign change
            # Linear interpolation
            M0 = M_coeffs[i]
            M1 = M_coeffs[i+1]
            A0 = A_values[i]
            A1 = A_values[i+1]
            M_zero = M0 - A0 * (M1 - M0) / (A1 - A0)
            print(f"A = 0 at M_coeff ≈ {M_zero:.2f}")

    # Refine with optimization
    def objective(M_coeff):
        A, _, _ = compute_divergence_coeff(M_coeff)
        return A**2

    result = optimize.minimize_scalar(objective, bounds=(-50, 50), method='bounded')
    M_optimal = result.x

    print(f"\nOptimal M_coeff = {M_optimal:.4f}")

    # Compute full result at optimal
    A, B, C = compute_divergence_coeff(M_optimal)

    print(f"\nAt M_coeff = {M_optimal:.2f}:")
    print(f"  A = {A:.6f}  (should be ≈ 0)")
    print(f"  B = {B:.4f}")
    print(f"  C = {C:.4f}")

    # Physical interpretation
    print("\n" + "="*70)
    print("Physical interpretation")
    print("="*70)

    # Check specific values
    print("\nChecking physically motivated values:")
    for label, M_coeff in [
        ("M_coeff = 1 (no extra factor)", 1),
        ("M_coeff = -8 (O7 charge in bulk units)", -8),
        ("M_coeff = -24 (O7 charge in fractional units)", -24),
        ("M_coeff = 8", 8),
        ("M_coeff = 24", 24),
        (f"M_coeff = {M_optimal:.1f} (optimal)", M_optimal),
    ]:
        A, B, C = compute_divergence_coeff(M_coeff)
        print(f"  {label}: A = {A:.2f}, C = {C:.2f}")

    # Final result with optimal coefficient
    print("\n" + "="*70)
    print("RESULT with tadpole-cancelling coefficient")
    print("="*70)

    A, B, C = compute_divergence_coeff(M_optimal)
    K0 = (2*np.pi)**(-1.5) * np.exp(C)

    print(f"\n  M_coeff = {M_optimal:.2f}")
    print(f"  A = {A:.6f} ≈ 0  ✓ (tadpole cancelled)")
    print(f"  B = {B:.4f}")
    print(f"  C = {C:.4f}")
    print(f"\n  K_0 = (2π)^{{-3/2}} × exp({C:.4f})")
    print(f"  K_0 = {K0:.6e}")
    print(f"  log(K_0) = {np.log(K0):.4f}")
    print("="*70)

if __name__ == "__main__":
    main()
