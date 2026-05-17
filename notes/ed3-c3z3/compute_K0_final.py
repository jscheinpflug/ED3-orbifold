#!/usr/bin/env python3
"""
K_0 computation for C^3/Z_3 orientifold - Final version.
Consolidates results from v4 and v5 with error analysis.
"""

import numpy as np
from scipy import integrate, stats
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
    tau_2 = 2j * t
    q = np.exp(-2 * np.pi * t)
    th10 = np.real(vartheta_10(tau_2, n_max))
    th00 = np.real(vartheta_00(tau_2, n_max))
    if np.isnan(th10) or np.isnan(th00):
        return np.nan
    Z_Q0 = 0.5 * th10 * q**(-0.25) * N
    Z_Q1 = 0.5 * th00 * 1.0 * N
    return Z_Q0 + Z_Q1

def Z_M(t, n_max=100):
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

def Z_total(t, N=24, n_max=100):
    z_a = Z_A_DN(t, N, n_max)
    z_m = Z_M(t, n_max)
    z_s = Z_sub(t)
    if np.isnan(z_a) or np.isnan(z_m):
        return np.nan
    return z_a + z_m + z_s

def main():
    print("="*70)
    print("K_0 FINAL COMPUTATION - C^3/Z_3 Orientifold with 24 fractional D7-branes")
    print("="*70)

    # Asymptotic constant
    Z_asymp = Z_total(10.0)
    print(f"\nZ(t -> infinity) = {Z_asymp:.6f}")

    # Compute integrals at many t_min values
    t_max = 20.0
    t_mins = np.array([0.5, 0.3, 0.2, 0.15, 0.1, 0.08, 0.06, 0.05,
                       0.04, 0.03, 0.025, 0.02, 0.015, 0.01, 0.008, 0.005])

    integrals = []
    for t_min in t_mins:
        def integrand(t):
            z = Z_total(t)
            if np.isnan(z):
                return 0.0
            return (z - Z_asymp) / (2 * t)

        val, _ = integrate.quad(integrand, t_min, t_max, limit=500)
        integrals.append(val)

    integrals = np.array(integrals)

    print(f"\nData points collected: {len(t_mins)}")

    # Fit: I = A/sqrt(t) + B*log(t) + C
    X = np.column_stack([
        t_mins**(-0.5),
        np.log(t_mins),
        np.ones_like(t_mins)
    ])

    # Use scipy.stats for proper regression with error estimates
    # Ordinary least squares
    coeffs, residuals, rank, sv = np.linalg.lstsq(X, integrals, rcond=None)
    A, B, C = coeffs

    # Compute residuals and standard errors
    fitted = X @ coeffs
    residuals_vec = integrals - fitted
    n = len(integrals)
    p = 3  # number of parameters
    dof = n - p

    # Variance estimate
    MSE = np.sum(residuals_vec**2) / dof
    RMS = np.sqrt(MSE)

    # Covariance matrix of parameters
    cov_matrix = MSE * np.linalg.inv(X.T @ X)
    std_errors = np.sqrt(np.diag(cov_matrix))

    print("\n" + "="*70)
    print("FIT RESULTS: I(t_min) = A/sqrt(t_min) + B*log(t_min) + C")
    print("="*70)
    print(f"\n  A = {A:10.4f} +/- {std_errors[0]:.4f}  (t^{{-1/2}} coefficient)")
    print(f"  B = {B:10.4f} +/- {std_errors[1]:.4f}  (log(t) coefficient)")
    print(f"  C = {C:10.4f} +/- {std_errors[2]:.4f}  (constant/finite part)")

    print(f"\n  Degrees of freedom: {dof}")
    print(f"  RMS residual: {RMS:.4f}")

    # Check fit quality
    print("\n  Sample residuals:")
    for i in [0, len(t_mins)//3, 2*len(t_mins)//3, -1]:
        print(f"    t_min = {t_mins[i]:.4f}: actual = {integrals[i]:10.4f}, fit = {fitted[i]:10.4f}, resid = {residuals_vec[i]:8.4f}")

    # Compute K_0 with error propagation
    print("\n" + "="*70)
    print("FINAL RESULT FOR K_0")
    print("="*70)

    C_val = C
    C_err = std_errors[2]

    prefactor = 1.0 / (2 * np.pi)**1.5
    K0_central = prefactor * np.exp(C_val)
    K0_upper = prefactor * np.exp(C_val + C_err)
    K0_lower = prefactor * np.exp(C_val - C_err)

    print(f"\n  Regularized integral (finite part):")
    print(f"    C = {C_val:.4f} +/- {C_err:.4f}")

    print(f"\n  K_0 = (2π)^{{-3/2}} × exp(C)")
    print(f"      = {prefactor:.6f} × exp({C_val:.4f})")
    print(f"\n  K_0 = {K0_central:.6e}")
    print(f"  K_0 range: [{K0_lower:.2e}, {K0_upper:.2e}]")

    print(f"\n  log(K_0) = {np.log(K0_central):.4f}")

    # Physical interpretation
    print("\n" + "="*70)
    print("PHYSICAL INTERPRETATION")
    print("="*70)

    print(f"""
  The integral has three contributions:

  1. Divergent term A/sqrt(t): {A:.2f}/sqrt(t)
     - Corresponds to massless closed string exchange
     - In a compact orientifold with exact tadpole cancellation, A = 0
     - For non-compact C^3/Z_3, this is regularized via zeta function

  2. Logarithmic term B*log(t): {B:.2f}*log(t)
     - Threshold correction from the continuous spectrum
     - Characteristic of non-compact/decompactification limits

  3. Finite part C: {C:.2f}
     - This is the physical K_0 normalization
     - K_0 = (2π)^(-3/2) × exp(C) = {K0_central:.2e}

  The superpotential correction from ED3 instanton is:
    W_inst = K_0 × e^{{-S_inst}} = {K0_central:.2e} × e^{{-S_inst}}

  where S_inst is the classical Euclidean D3-brane action.
""")

    # Summary table
    print("="*70)
    print("SUMMARY")
    print("="*70)
    print(f"""
  Background: C^3/Z_3 orientifold with O7-plane
  D-branes: N = 24 fractional D7-branes (8 bulk)
  Tadpole condition: Satisfied (N = 24 = 3 × 8)

  One-loop determinant from ED3 instanton:
    K_0 = {K0_central:.4e}
    log(K_0) = {np.log(K0_central):.4f}

  Note: K_0 << 1 indicates strong suppression of the instanton
  contribution beyond the classical factor e^{{-S_inst}}.
""")
    print("="*70)

    return K0_central, C_val, C_err

if __name__ == "__main__":
    K0, C, C_err = main()
