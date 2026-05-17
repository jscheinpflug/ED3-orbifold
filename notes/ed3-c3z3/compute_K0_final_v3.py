#!/usr/bin/env python3
"""
K_0 computation - Final Version 3.
Precise determination with proper error analysis.
"""

import numpy as np
from scipy import integrate, optimize
import warnings
warnings.filterwarnings('ignore')

def vartheta_10(tau, n_max=150):
    q = np.exp(2j * np.pi * tau)
    if np.abs(q) >= 1 - 1e-10:
        return np.nan
    result = 0.0 + 0j
    for n in range(-n_max, n_max + 1):
        result += q**((n + 0.5)**2 / 2)
    return result

def vartheta_00(tau, n_max=150):
    q = np.exp(2j * np.pi * tau)
    if np.abs(q) >= 1 - 1e-10:
        return np.nan
    result = 1.0 + 0j
    for n in range(1, n_max + 1):
        result += 2 * q**(n**2 / 2)
    return result

def Z_A_DN(t, N=24, n_max=150):
    tau_2 = 2j * t
    q = np.exp(-2 * np.pi * t)
    th10 = np.real(vartheta_10(tau_2, n_max))
    th00 = np.real(vartheta_00(tau_2, n_max))
    if np.isnan(th10) or np.isnan(th00):
        return np.nan
    return 0.5 * N * (th10 * q**(-0.25) + th00)

def Z_M_base(t, n_max=150):
    hat_tau_2 = 2 * (1j * t + 0.5)
    hat_q = np.exp(-2 * np.pi * t + 1j * np.pi)
    th10_hat = vartheta_10(hat_tau_2, n_max)
    th00_hat = vartheta_00(hat_tau_2, n_max)
    if np.isnan(np.abs(th10_hat)) or np.isnan(np.abs(th00_hat)):
        return np.nan
    return np.real(th10_hat * hat_q**(-0.25) - th00_hat)

def Z_sub(t):
    return 3 * (np.exp(-2 * np.pi * t) - 1)

def compute_integral(M_coeff, t_min, t_max=30.0):
    """Compute the regularized integral for given M_coeff."""
    def Z_total(t):
        z_a = Z_A_DN(t)
        z_m = M_coeff * Z_M_base(t)
        z_s = Z_sub(t)
        if np.isnan(z_a) or np.isnan(z_m):
            return np.nan
        return z_a + z_m + z_s

    Z_asymp = Z_total(15.0)

    def integrand(t):
        z = Z_total(t)
        if np.isnan(z):
            return 0.0
        return (z - Z_asymp) / (2 * t)

    val, err = integrate.quad(integrand, t_min, t_max, limit=1000)
    return val, err

def find_optimal_M():
    """Find M_coeff that gives convergent integral."""
    print("="*70)
    print("Finding optimal Möbius coefficient")
    print("="*70)

    # For each M_coeff, check if integral converges as t_min → 0
    def divergence_rate(M_coeff):
        """Measure how fast integral diverges with t_min."""
        I1, _ = compute_integral(M_coeff, 0.02)
        I2, _ = compute_integral(M_coeff, 0.01)
        I3, _ = compute_integral(M_coeff, 0.005)

        # If converging, differences should decrease
        d1 = I2 - I1  # Change from t=0.02 to t=0.01
        d2 = I3 - I2  # Change from t=0.01 to t=0.005

        # For t^{-1/2} divergence: d ~ t^{-1/2}, so d2/d1 ~ sqrt(2)
        # For log divergence: d ~ log ratio, constant
        # For convergence: d2 < d1

        return d2 - d1  # Should be ~0 for log, negative for convergent

    print("\nSearching for M_coeff where integral converges:")
    print("-"*50)

    for M in [70, 72, 74, 76, 78, 80]:
        I1, _ = compute_integral(M, 0.02)
        I2, _ = compute_integral(M, 0.01)
        I3, _ = compute_integral(M, 0.005)
        d1 = I2 - I1
        d2 = I3 - I2
        ratio = d2/d1 if abs(d1) > 0.01 else np.nan
        print(f"  M = {M}: I(0.02)={I1:.2f}, I(0.01)={I2:.2f}, I(0.005)={I3:.2f}, "
              f"d2/d1 = {ratio:.3f}")

    # For sqrt(t_min) divergence: d2/d1 = sqrt(0.01/0.005) / sqrt(0.02/0.01) = sqrt(2)/sqrt(2) = 1
    # Actually: d = I(t1) - I(t2) ~ -A*(1/sqrt(t2) - 1/sqrt(t1))
    # d1 ~ A*(10 - 7.07) = 2.93*A
    # d2 ~ A*(14.14 - 10) = 4.14*A
    # d2/d1 = 4.14/2.93 = 1.41

    # For log: d1 ~ B*log(2) = 0.69*B, d2 ~ B*log(2) = 0.69*B, ratio = 1

    # Optimal is where the sqrt term vanishes, leaving only log
    print("\n  If d2/d1 ≈ 1.41, divergence is √t")
    print("  If d2/d1 ≈ 1.0, divergence is log(t) only")
    print("  If d2/d1 < 1.0, converging")

    # Fine search
    def objective(M):
        I1, _ = compute_integral(M, 0.02)
        I2, _ = compute_integral(M, 0.01)
        I3, _ = compute_integral(M, 0.005)
        d1 = I2 - I1
        d2 = I3 - I2
        # Target: d2/d1 = 1 (log divergence only)
        return (d2/d1 - 1.0)**2 if abs(d1) > 0.01 else 1e10

    result = optimize.minimize_scalar(objective, bounds=(60, 90), method='bounded')
    M_optimal = result.x

    print(f"\nOptimal M_coeff = {M_optimal:.4f}")

    # Verify
    I1, _ = compute_integral(M_optimal, 0.02)
    I2, _ = compute_integral(M_optimal, 0.01)
    I3, _ = compute_integral(M_optimal, 0.005)
    d1 = I2 - I1
    d2 = I3 - I2
    print(f"  d2/d1 = {d2/d1:.4f} (target: 1.0)")

    return M_optimal

def compute_final_result(M_coeff):
    """Compute K_0 with given M_coeff."""
    print("\n" + "="*70)
    print(f"Computing K_0 with M_coeff = {M_coeff:.4f}")
    print("="*70)

    # Collect data
    t_mins = np.array([0.5, 0.3, 0.2, 0.15, 0.1, 0.08, 0.06, 0.05,
                       0.04, 0.03, 0.025, 0.02, 0.015, 0.01])
    integrals = []
    for t_min in t_mins:
        val, _ = compute_integral(M_coeff, t_min)
        integrals.append(val)
    integrals = np.array(integrals)

    print("\nData:")
    for i in range(len(t_mins)):
        print(f"  t_min = {t_mins[i]:.4f}: I = {integrals[i]:10.4f}")

    # If √t term is cancelled, fit: I = B*log(t) + C
    # But there may still be subleading √t terms from massive modes

    # Fit: I = B*log(t) + C (for t < 0.1 where massive modes are suppressed)
    mask = t_mins <= 0.1
    X = np.column_stack([np.log(t_mins[mask]), np.ones(np.sum(mask))])
    coeffs, residuals, _, _ = np.linalg.lstsq(X, integrals[mask], rcond=None)
    B, C = coeffs

    fitted = X @ coeffs
    rms = np.sqrt(np.mean((integrals[mask] - fitted)**2))

    print(f"\nFit (t_min ≤ 0.1): I = B*log(t) + C")
    print(f"  B = {B:.4f}")
    print(f"  C = {C:.4f}")
    print(f"  RMS = {rms:.6f}")

    # Alternative: extrapolate using Richardson
    # If I(t) = B*log(t) + C + D*√t + ...
    # I(t/4) - 2*I(t/2) + I(t) should give estimate of C corrected for √t

    # Richardson extrapolation for log + √t
    t1, t2, t3 = 0.04, 0.01, 0.0025  # Factor of 4 spacing
    I1, _ = compute_integral(M_coeff, t1)
    I2, _ = compute_integral(M_coeff, t2)

    # For I = B*log(t) + C:
    # I1 = B*log(t1) + C
    # I2 = B*log(t2) + C
    # C = (I1*log(t2) - I2*log(t1)) / (log(t2) - log(t1))
    C_rich = (I1 * np.log(t2) - I2 * np.log(t1)) / (np.log(t2) - np.log(t1))
    B_rich = (I2 - I1) / (np.log(t2) - np.log(t1))

    print(f"\nRichardson extrapolation (assuming I = B*log(t) + C):")
    print(f"  Using t_min = {t1} and {t2}")
    print(f"  B = {B_rich:.4f}")
    print(f"  C = {C_rich:.4f}")

    # Best estimate: average of methods
    C_best = (C + C_rich) / 2
    C_err = abs(C - C_rich) / 2

    print(f"\nBest estimate: C = {C_best:.4f} ± {C_err:.4f}")

    K0 = (2*np.pi)**(-1.5) * np.exp(C_best)
    K0_low = (2*np.pi)**(-1.5) * np.exp(C_best - C_err)
    K0_high = (2*np.pi)**(-1.5) * np.exp(C_best + C_err)

    print("\n" + "="*70)
    print("FINAL RESULT")
    print("="*70)
    print(f"""
  Möbius coefficient for tadpole cancellation: M_coeff = {M_coeff:.2f}

  Regularized integral (finite part):
    C = {C_best:.4f} ± {C_err:.4f}

  K_0 = (2π)^{{-3/2}} × exp(C)
  K_0 = {K0:.6e}
  K_0 range: [{K0_low:.2e}, {K0_high:.2e}]

  log(K_0) = {np.log(K0):.4f} ± {C_err:.4f}
""")
    print("="*70)

    return K0, C_best, C_err

def main():
    M_optimal = find_optimal_M()
    K0, C, C_err = compute_final_result(M_optimal)
    return K0, C, M_optimal

if __name__ == "__main__":
    K0, C, M_opt = main()
