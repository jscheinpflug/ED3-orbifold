#!/usr/bin/env python3
"""
Numerical computation of K_0 for C^3/Z_3 orientifold - Version 2.
More careful implementation following arXiv:2204.02981.
"""

import numpy as np
from scipy import integrate
import warnings
warnings.filterwarnings('ignore')

def theta3(q, n_max=200):
    """theta_3(q) = theta_00 = sum_{n=-inf}^{inf} q^{n^2}
    Note: some conventions use q^{n^2/2}, we use q^{n^2}."""
    if np.abs(q) >= 1:
        return np.nan
    result = 1.0
    for n in range(1, n_max + 1):
        result += 2 * np.real(q**(n**2))
    return result

def theta4(q, n_max=200):
    """theta_4(q) = theta_01 = sum_{n=-inf}^{inf} (-1)^n q^{n^2}"""
    if np.abs(q) >= 1:
        return np.nan
    result = 1.0
    for n in range(1, n_max + 1):
        result += 2 * ((-1)**n) * np.real(q**(n**2))
    return result

def theta2(q, n_max=200):
    """theta_2(q) = theta_10 = sum_{n=-inf}^{inf} q^{(n+1/2)^2}"""
    if np.abs(q) >= 1:
        return np.nan
    result = 0.0
    for n in range(-n_max, n_max + 1):
        result += np.real(q**((n + 0.5)**2))
    return result

def eta_q(q, n_max=200):
    """Dedekind eta in terms of q: eta = q^{1/24} prod_{n=1}^{inf} (1-q^n)
    Here q = e^{2*pi*i*tau}"""
    if np.abs(q) >= 1:
        return np.nan
    result = q**(1/24)
    for n in range(1, n_max + 1):
        result *= (1 - q**n)
    return result

# For the reference conventions: theta_{ab}(tau) uses q = e^{2*pi*i*tau}
# theta_00(tau) = sum_n q^{n^2/2}
# We need to be careful about the power of q

def vartheta_00(tau, n_max=200):
    """vartheta_{00}(tau) = sum_{n} q^{n^2/2} where q = e^{2*pi*i*tau}"""
    q = np.exp(2j * np.pi * tau)
    if np.abs(q) >= 1 - 1e-10:
        return np.nan
    result = 1.0 + 0j
    for n in range(1, n_max + 1):
        result += 2 * q**(n**2 / 2)
    return result

def vartheta_01(tau, n_max=200):
    """vartheta_{01}(tau) = sum_{n} (-1)^n q^{n^2/2}"""
    q = np.exp(2j * np.pi * tau)
    if np.abs(q) >= 1 - 1e-10:
        return np.nan
    result = 1.0 + 0j
    for n in range(1, n_max + 1):
        result += 2 * ((-1)**n) * q**(n**2 / 2)
    return result

def vartheta_10(tau, n_max=200):
    """vartheta_{10}(tau) = sum_{n} q^{(n+1/2)^2/2}"""
    q = np.exp(2j * np.pi * tau)
    if np.abs(q) >= 1 - 1e-10:
        return np.nan
    result = 0.0 + 0j
    for n in range(-n_max, n_max + 1):
        result += q**((n + 0.5)**2 / 2)
    return result

def dedekind_eta(tau, n_max=200):
    """eta(tau) = q^{1/24} prod_{n>=1} (1-q^n), q = e^{2*pi*i*tau}"""
    q = np.exp(2j * np.pi * tau)
    if np.abs(q) >= 1 - 1e-10:
        return np.nan
    result = q**(1/24)
    for n in range(1, n_max + 1):
        result *= (1 - q**n)
    return result

def test_theta():
    """Test theta functions at known values."""
    print("Testing theta functions...")
    tau = 1j  # t = 1
    q = np.exp(-2 * np.pi)

    th00 = vartheta_00(tau)
    th01 = vartheta_01(tau)
    th10 = vartheta_10(tau)
    eta_val = dedekind_eta(tau)

    print(f"tau = i:")
    print(f"  theta_00 = {th00:.10f}")
    print(f"  theta_01 = {th01:.10f}")
    print(f"  theta_10 = {th10:.10f}")
    print(f"  eta = {eta_val:.10f}")

    # Known relation: theta_00^4 = theta_01^4 + theta_10^4 (Jacobi identity)
    jacobi_check = np.abs(th00**4 - th01**4 - th10**4)
    print(f"  Jacobi identity check |th00^4 - th01^4 - th10^4| = {jacobi_check:.2e}")

    # Check at tau = 2i
    tau = 2j
    th00_2 = vartheta_00(tau)
    th10_2 = vartheta_10(tau)
    print(f"\ntau = 2i:")
    print(f"  theta_00(2i) = {th00_2:.10f}")
    print(f"  theta_10(2i) = {th10_2:.10f}")

    return True

def Z_A_DN_simple(t, n_max=100):
    """
    Simplified annulus partition function for testing.

    For the C^3/Z_3 orbifold with 24 fractional D7-branes,
    the annulus contribution from the DN sector.

    Key point from eq. (E.48): The partition function is a sum
    over representations weighted by theta functions.

    For the vacuum rep only (leading contribution):
    Z_{A,DN} = (1/2) * theta_10(2*tau) * q^{-1/4} * N
             - (1/2) * theta_00(2*tau) * (-1) * q^0 * N
             + higher order terms

    where N = 24 (number of fractional D7-branes).
    """
    tau = 1j * t  # purely imaginary
    tau_2 = 2j * t

    q = np.exp(-2 * np.pi * t)

    # Theta functions at 2*tau
    th10 = np.real(vartheta_10(tau_2, n_max))
    th00 = np.real(vartheta_00(tau_2, n_max))

    if np.isnan(th10) or np.isnan(th00):
        return np.nan

    N = 24  # 24 fractional D7-branes for tadpole cancellation

    # From the reference eq. (E.44) and (E.48):
    # The vacuum representation contributes:
    # - (h=0, Q=0) with n=1: contributes to Q=0 sum
    # - (h=1/2, Q=1) with n=-1: contributes to Q=1 sum

    # Z_{A,DN} = sum_Q (-1)^Q/2 * theta_{1-Q,0}(2*tau) * sum_h n_{h,Q} q^{h-(1+Q)/4}

    # Q=0: coeff = +1/2, uses theta_10(2*tau)
    #      h=0 contributes: q^{0-1/4} = q^{-1/4}
    # Q=1: coeff = -1/2, uses theta_00(2*tau)
    #      h=1/2 with n=-1 contributes: -1 * q^{1/2-1/2} = -q^0 = -1

    # Leading terms:
    Z_Q0 = 0.5 * th10 * q**(-0.25) * N
    Z_Q1 = -0.5 * th00 * (-1) * 1.0 * N  # Note: q^0 = 1

    return Z_Q0 + Z_Q1

def Z_M_simple(t, n_max=100):
    """
    Simplified Mobius strip partition function.

    hat_q = -e^{-2*pi*t} = e^{-2*pi*t + i*pi}
    hat_tau = tau + 1/2 = i*t + 1/2

    For the vacuum representation with sigma_vac = -1:
    Z_M involves theta functions at 2*hat_tau = 2*i*t + 1
    """
    hat_tau = 1j * t + 0.5
    hat_tau_2 = 2 * hat_tau  # = 2*i*t + 1

    hat_q = np.exp(-2 * np.pi * t + 1j * np.pi)  # = -e^{-2*pi*t}

    # Theta functions at 2*hat_tau
    th10_hat = vartheta_10(hat_tau_2, n_max)
    th00_hat = vartheta_00(hat_tau_2, n_max)

    if np.isnan(np.abs(th10_hat)) or np.isnan(np.abs(th00_hat)):
        return np.nan

    # From eq. (E.45):
    # Z_M = sum_Q (-1)^{1-Q} theta_{1-Q,0}(2*hat_tau) *
    #       sum_{h,sigma} sigma * n * e^{-i*pi*(h-Q/2)} * hat_q^{h-(1+Q)/4}

    # For vac with sigma=-1:
    # Q=0, h=0: (-1)^1 * theta_10 * (-1) * 1 * e^0 * hat_q^{-1/4}
    #         = theta_10 * hat_q^{-1/4}
    # Q=1, h=1/2, n=-1: (-1)^0 * theta_00 * (-1)*(-1) * e^0 * hat_q^0
    #         = -theta_00

    Z_Q0 = th10_hat * hat_q**(-0.25)
    Z_Q1 = -th00_hat * 1.0  # hat_q^0 = 1

    result = Z_Q0 + Z_Q1

    return np.real(result)

def Z_subtraction(t):
    """Zero mode subtraction: 3 * (e^{-2*pi*t} - 1)"""
    return 3 * (np.exp(-2 * np.pi * t) - 1)

def integrand(t, n_max=100):
    """Full integrand: [Z_A + Z_M + Z_sub] / (2*t)"""
    if t < 1e-12:
        return 0.0

    z_a = Z_A_DN_simple(t, n_max)
    z_m = Z_M_simple(t, n_max)
    z_sub = Z_subtraction(t)

    if np.isnan(z_a) or np.isnan(z_m):
        return np.nan

    total = z_a + z_m + z_sub
    return total / (2 * t)

def debug_partition_functions():
    """Debug the partition functions at various t values."""
    print("\nDebugging partition functions:")
    print("-" * 60)

    for t in [0.1, 0.5, 1.0, 2.0, 5.0]:
        z_a = Z_A_DN_simple(t)
        z_m = Z_M_simple(t)
        z_sub = Z_subtraction(t)
        total = z_a + z_m + z_sub
        intgd = integrand(t)

        print(f"t = {t:.1f}:")
        print(f"  Z_A = {z_a:.6f}")
        print(f"  Z_M = {z_m:.6f}")
        print(f"  Z_sub = {z_sub:.6f}")
        print(f"  Total = {total:.6f}")
        print(f"  Integrand = {intgd:.6f}")

def compute_integral(t_min=0.01, t_max=100.0, n_max=100):
    """Compute the integral from t_min to t_max."""

    def safe_integrand(t):
        result = integrand(t, n_max)
        if np.isnan(result):
            return 0.0
        return result

    result, error = integrate.quad(safe_integrand, t_min, t_max, limit=500)

    return result, error

def main():
    print("="*70)
    print("K_0 Computation for C^3/Z_3 Orientifold")
    print("="*70)

    # Test theta functions first
    test_theta()

    # Debug partition functions
    debug_partition_functions()

    # Convergence tests
    print("\n" + "="*70)
    print("Convergence Tests")
    print("="*70)

    # Test 1: Mode truncation
    print("\n1. Mode truncation convergence (t_min=0.1, t_max=10):")
    for n_max in [50, 100, 150, 200]:
        result, error = compute_integral(0.1, 10.0, n_max)
        print(f"  n_max = {n_max:3d}: integral = {result:12.6f}, error = {error:.2e}")

    # Test 2: t_min convergence (IR cutoff)
    print("\n2. IR cutoff convergence (t_max=100, n_max=150):")
    for t_min in [1.0, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01]:
        result, error = compute_integral(t_min, 100.0, 150)
        print(f"  t_min = {t_min:.3f}: integral = {result:12.6f}")

    # Test 3: t_max convergence (UV cutoff in open string = IR in closed)
    print("\n3. UV cutoff convergence (t_min=0.01, n_max=150):")
    for t_max in [10, 20, 50, 100, 200, 500]:
        result, error = compute_integral(0.01, t_max, 150)
        print(f"  t_max = {t_max:4d}: integral = {result:12.6f}")

    # Final result
    print("\n" + "="*70)
    print("FINAL RESULT")
    print("="*70)

    t_min = 0.01
    t_max = 500
    n_max = 200

    integral, error = compute_integral(t_min, t_max, n_max)

    prefactor = 1.0 / (2 * np.pi)**1.5
    K0 = prefactor * np.exp(integral)

    print(f"\nIntegration parameters:")
    print(f"  t_min = {t_min}")
    print(f"  t_max = {t_max}")
    print(f"  n_max = {n_max}")
    print(f"\nIntegral = {integral:.8f}")
    print(f"Integration error = {error:.2e}")
    print(f"\nPrefactor = 1/(2*pi)^(3/2) = {prefactor:.10f}")
    print(f"\nK_0 = {K0:.10e}")
    print(f"K_0 = {K0:.8f}")

    # Express in useful form
    print(f"\nexp(integral) = {np.exp(integral):.8f}")
    print(f"log(K_0) = {np.log(K0):.6f}")

if __name__ == "__main__":
    main()
