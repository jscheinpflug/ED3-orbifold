#!/usr/bin/env python3
"""
K_0 computation for C^3/Z_3 orientifold - Version 4.
Proper zeta-function regularization following arXiv:2204.02981.

Key insight: The small-t divergence goes as t^{-1/2} (from theta function
modular transformation). Zeta-function regularization:
  ζ(s) = ∫_0^∞ t^{s-1} Z''(t) dt / Γ(s)
The pole at s=1/2 is handled by analytic continuation, giving finite ζ'(0).
"""

import numpy as np
from scipy import integrate, optimize
import warnings
warnings.filterwarnings('ignore')

def vartheta_10(tau, n_max=100):
    """vartheta_{10}(tau) = sum_{n} q^{(n+1/2)^2/2}"""
    q = np.exp(2j * np.pi * tau)
    if np.abs(q) >= 1 - 1e-10:
        return np.nan
    result = 0.0 + 0j
    for n in range(-n_max, n_max + 1):
        result += q**((n + 0.5)**2 / 2)
    return result

def vartheta_00(tau, n_max=100):
    """vartheta_{00}(tau) = sum_{n} q^{n^2/2}"""
    q = np.exp(2j * np.pi * tau)
    if np.abs(q) >= 1 - 1e-10:
        return np.nan
    result = 1.0 + 0j
    for n in range(1, n_max + 1):
        result += 2 * q**(n**2 / 2)
    return result

def Z_A_DN(t, N=24, n_max=100):
    """Annulus partition function Z_{A,DN}."""
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
    """Mobius strip partition function."""
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
    """Zero mode subtraction: 3*(e^{-2πt} - 1)"""
    return 3 * (np.exp(-2 * np.pi * t) - 1)

def Z_total(t, N=24, n_max=100):
    """Total partition function."""
    z_a = Z_A_DN(t, N, n_max)
    z_m = Z_M(t, n_max)
    z_s = Z_sub(t)

    if np.isnan(z_a) or np.isnan(z_m):
        return np.nan

    return z_a + z_m + z_s

def compute_asymptotic_constant(t_test=10.0, N=24, n_max=100):
    """Compute the asymptotic constant (zero mode contribution) at large t."""
    return Z_total(t_test, N, n_max)

def integrand_regularized(t, Z_asymp, N=24, n_max=100):
    """Regularized integrand: [Z(t) - Z_asymp] / (2t)"""
    if t < 1e-12:
        return 0.0

    z = Z_total(t, N, n_max)

    if np.isnan(z):
        return np.nan

    return (z - Z_asymp) / (2 * t)

def compute_raw_integral(t_min, t_max, Z_asymp, N=24, n_max=100):
    """Compute the integral with given cutoffs."""
    def safe_integrand(t):
        result = integrand_regularized(t, Z_asymp, N, n_max)
        if np.isnan(result) or np.isinf(result):
            return 0.0
        return result

    result, error = integrate.quad(safe_integrand, t_min, t_max, limit=500)
    return result, error

def extract_finite_part():
    """
    Extract the finite part of the integral using zeta-function regularization.

    The integral has the form:
    I(t_min) = A * t_min^{-1/2} + B * log(t_min) + C + o(1)

    We fit to multiple data points and extract C (the finite part).
    """
    print("="*70)
    print("Extracting finite part via zeta-function regularization")
    print("="*70)

    # Compute asymptotic constant
    Z_asymp = compute_asymptotic_constant(t_test=10.0)
    print(f"\nAsymptotic constant Z(t→∞) = {Z_asymp:.6f}")

    # Collect data at various t_min values
    t_max = 20.0  # Fixed UV cutoff (converges quickly)
    n_max = 150

    t_mins = [0.5, 0.2, 0.1, 0.08, 0.06, 0.05, 0.04, 0.03, 0.02, 0.015, 0.01]
    integrals = []

    print(f"\nComputing integrals for t_max = {t_max}:")
    print("-"*50)

    for t_min in t_mins:
        integral, _ = compute_raw_integral(t_min, t_max, Z_asymp, n_max=n_max)
        integrals.append(integral)
        print(f"  t_min = {t_min:.4f}: I = {integral:12.4f}, t^{{-1/2}} = {t_min**(-0.5):8.4f}")

    t_mins = np.array(t_mins)
    integrals = np.array(integrals)

    # Fit: I = A * t^{-1/2} + B * log(t) + C
    # Using matrix form: y = X @ [A, B, C]
    X = np.column_stack([
        t_mins**(-0.5),
        np.log(t_mins),
        np.ones_like(t_mins)
    ])

    # Least squares fit
    coeffs, residuals, rank, sv = np.linalg.lstsq(X, integrals, rcond=None)
    A, B, C = coeffs

    print("\n" + "="*70)
    print("Fit: I(t_min) = A * t_min^{-1/2} + B * log(t_min) + C")
    print("="*70)
    print(f"  A (divergent coeff) = {A:.6f}")
    print(f"  B (log coeff)       = {B:.6f}")
    print(f"  C (finite part)     = {C:.6f}")

    # Residuals
    fitted = X @ coeffs
    print(f"\nFit quality:")
    for i, t_min in enumerate(t_mins):
        resid = integrals[i] - fitted[i]
        print(f"  t_min = {t_min:.4f}: actual = {integrals[i]:10.4f}, fit = {fitted[i]:10.4f}, resid = {resid:8.4f}")

    rms_error = np.sqrt(np.mean((integrals - fitted)**2))
    print(f"\n  RMS error = {rms_error:.4f}")

    # Try also with just A*t^{-1/2} + C (simpler model)
    print("\n" + "-"*50)
    print("Simpler fit: I(t_min) = A * t_min^{-1/2} + C")
    print("-"*50)

    X_simple = np.column_stack([t_mins**(-0.5), np.ones_like(t_mins)])
    coeffs_simple, _, _, _ = np.linalg.lstsq(X_simple, integrals, rcond=None)
    A_s, C_s = coeffs_simple

    print(f"  A (divergent coeff) = {A_s:.6f}")
    print(f"  C (finite part)     = {C_s:.6f}")

    fitted_simple = X_simple @ coeffs_simple
    rms_simple = np.sqrt(np.mean((integrals - fitted_simple)**2))
    print(f"  RMS error = {rms_simple:.4f}")

    # Richardson extrapolation approach
    print("\n" + "-"*50)
    print("Richardson extrapolation (assuming I ~ A*t^{-1/2} + C)")
    print("-"*50)

    # Use pairs of points with t_ratio = 2 to eliminate the t^{-1/2} term
    # I(t) = A * t^{-1/2} + C
    # I(t/4) = A * 2*t^{-1/2} + C
    # => C = 2*I(t) - I(t/4)

    richardson_estimates = []
    for i in range(len(t_mins)):
        for j in range(i+1, len(t_mins)):
            t1, t2 = t_mins[i], t_mins[j]
            I1, I2 = integrals[i], integrals[j]

            # I1 = A/sqrt(t1) + C
            # I2 = A/sqrt(t2) + C
            # Solve: C = (I1 * sqrt(t1) - I2 * sqrt(t2)) / (sqrt(t1) - sqrt(t2))
            #       (after some algebra, a better form:)
            # C = I2 - (I1 - I2) * sqrt(t2) / (sqrt(t1) - sqrt(t2))

            r = np.sqrt(t1/t2)  # > 1 since t1 > t2
            C_est = (I2 * r - I1) / (r - 1)
            richardson_estimates.append(C_est)

    # Use most reliable pairs (small t values, close ratios)
    # Take median to be robust
    C_richardson = np.median(richardson_estimates)
    print(f"  Median C estimate from all pairs: {C_richardson:.4f}")

    # Focus on smallest t values
    small_pairs = []
    for i, t1 in enumerate(t_mins[-4:], len(t_mins)-4):
        for j, t2 in enumerate(t_mins[-4:], len(t_mins)-4):
            if i < j:
                I1, I2 = integrals[i], integrals[j]
                r = np.sqrt(t1/t2)
                C_est = (I2 * r - I1) / (r - 1)
                small_pairs.append((t1, t2, C_est))

    print("\n  Estimates from smallest t_min pairs:")
    for t1, t2, C_est in small_pairs:
        print(f"    ({t1:.3f}, {t2:.3f}): C = {C_est:.4f}")

    C_small = np.mean([x[2] for x in small_pairs])
    print(f"\n  Mean C from small-t pairs: {C_small:.4f}")

    return C, C_s, C_richardson, C_small

def compute_K0(integral_value):
    """Compute K_0 from the regularized integral."""
    prefactor = 1.0 / (2 * np.pi)**1.5
    K0 = prefactor * np.exp(integral_value)
    return K0

def main():
    print("="*70)
    print("K_0 Computation for C^3/Z_3 Orientifold - Version 4")
    print("Zeta-function regularization")
    print("="*70)

    # Extract finite part
    C_3param, C_2param, C_richardson, C_small = extract_finite_part()

    print("\n" + "="*70)
    print("FINAL RESULTS")
    print("="*70)

    estimates = {
        '3-param fit (A*t^{-1/2} + B*log(t) + C)': C_3param,
        '2-param fit (A*t^{-1/2} + C)': C_2param,
        'Richardson extrapolation (median)': C_richardson,
        'Richardson extrapolation (small-t pairs)': C_small,
    }

    print("\nFinite part estimates:")
    for name, C in estimates.items():
        K0 = compute_K0(C)
        print(f"\n  {name}:")
        print(f"    Regularized integral = {C:.6f}")
        print(f"    K_0 = {K0:.6e}")
        if K0 > 0:
            print(f"    log(K_0) = {np.log(K0):.6f}")

    # Best estimate: use mean of methods
    C_best = np.mean([C_3param, C_2param, C_small])
    K0_best = compute_K0(C_best)

    print("\n" + "="*70)
    print(f"BEST ESTIMATE (mean of methods):")
    print(f"  Regularized integral = {C_best:.6f}")
    print(f"  K_0 = (2π)^{{-3/2}} exp({C_best:.4f})")
    print(f"  K_0 = {K0_best:.6e}")
    print(f"  log(K_0) = {np.log(K0_best):.6f}")
    print("="*70)

if __name__ == "__main__":
    main()
