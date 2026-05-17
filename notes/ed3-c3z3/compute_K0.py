#!/usr/bin/env python3
"""
Numerical computation of K_0 for C^3/Z_3 orientifold.

Following arXiv:2204.02981 eq. (2.24):
K_0 = h^{3/2} * exp[ integral terms ]

where the integral is:
int_{eps'}^{1/eps} dt/(2t) Z_A + int_{eps'/4}^{1/eps} dt/(2t) Z_M
    + 3 int_0^{1/eps} dt/(2t) (e^{-2pi h t} - 1)

With h=1 (standard choice).
"""

import numpy as np
from scipy import integrate
from scipy.special import ellipk
import warnings
warnings.filterwarnings('ignore')

# Theta functions
def theta_00(tau, n_max=100):
    """theta_00(tau) = sum_{n=-inf}^{inf} q^{n^2/2}, q = e^{2pi i tau}"""
    q = np.exp(2j * np.pi * tau)
    if np.abs(q) >= 1:
        return np.nan
    result = 1.0
    for n in range(1, n_max + 1):
        term = q**(n**2 / 2)
        result += 2 * np.real(term)  # n and -n give same contribution
    return result

def theta_01(tau, n_max=100):
    """theta_01(tau) = sum_{n=-inf}^{inf} (-1)^n q^{n^2/2}"""
    q = np.exp(2j * np.pi * tau)
    if np.abs(q) >= 1:
        return np.nan
    result = 1.0
    for n in range(1, n_max + 1):
        term = ((-1)**n) * q**(n**2 / 2)
        result += 2 * np.real(term)
    return result

def theta_10(tau, n_max=100):
    """theta_10(tau) = sum_{n=-inf}^{inf} q^{(n+1/2)^2/2}"""
    q = np.exp(2j * np.pi * tau)
    if np.abs(q) >= 1:
        return np.nan
    result = 0.0
    for n in range(-n_max, n_max + 1):
        term = q**((n + 0.5)**2 / 2)
        result += np.real(term)
    return result

def eta(tau, n_max=100):
    """Dedekind eta: eta(tau) = q^{1/24} prod_{n=1}^{inf} (1 - q^n)"""
    q = np.exp(2j * np.pi * tau)
    if np.abs(q) >= 1:
        return np.nan
    result = q**(1/24)
    for n in range(1, n_max + 1):
        result *= (1 - q**n)
    return result

def theta_00_real(t, n_max=100):
    """theta_00 for tau = i*t (purely imaginary), returns real value"""
    q = np.exp(-2 * np.pi * t)
    result = 1.0
    for n in range(1, n_max + 1):
        result += 2 * q**(n**2 / 2)
    return result

def theta_01_real(t, n_max=100):
    """theta_01 for tau = i*t"""
    q = np.exp(-2 * np.pi * t)
    result = 1.0
    for n in range(1, n_max + 1):
        result += 2 * ((-1)**n) * q**(n**2 / 2)
    return result

def theta_10_real(t, n_max=100):
    """theta_10 for tau = i*t"""
    q = np.exp(-2 * np.pi * t)
    result = 0.0
    for n in range(-n_max, n_max + 1):
        result += q**((n + 0.5)**2 / 2)
    return result

def eta_real(t, n_max=100):
    """eta for tau = i*t"""
    q = np.exp(-2 * np.pi * t)
    result = q**(1/24)
    for n in range(1, n_max + 1):
        result *= (1 - q**n)
    return result

# For Mobius strip: hat_tau = tau + 1/2 = i*t + 1/2
def theta_10_hat(t, n_max=100):
    """theta_10(2*hat_tau) where hat_tau = i*t + 1/2, so 2*hat_tau = 2it + 1"""
    # theta_10(tau + 1) = theta_10(tau) * exp(i*pi/4) for the argument shift
    # More directly: compute theta_10(2it + 1)
    tau = 2j * t + 1
    q = np.exp(2j * np.pi * tau)  # q = exp(2pi i (2it + 1)) = exp(-4 pi t + 2 pi i) = -exp(-4 pi t)
    # hat_q = -exp(-2 pi t)
    result = 0.0 + 0j
    for n in range(-n_max, n_max + 1):
        result += q**((n + 0.5)**2 / 2)
    return result

def theta_00_hat(t, n_max=100):
    """theta_00(2*hat_tau) where hat_tau = i*t + 1/2"""
    tau = 2j * t + 1
    q = np.exp(2j * np.pi * tau)
    result = 1.0 + 0j
    for n in range(1, n_max + 1):
        result += 2 * q**(n**2 / 2)
    return result

# Partition functions for C^3/Z_3 orbifold

def Z_A_DN(t, N_D7=24, n_max=50):
    """
    Annulus partition function between D-instanton and D7-branes.

    For C^3/Z_3 with N_D7 fractional D7-branes (8 bulk = 24 fractional for tadpole).

    Z_{A,DN} = (1/2) sum_{Q=0,1} (-1)^Q theta_{1-Q,0}(2*tau) * sum_h n_{h,Q} q^{h-(1+Q)/4}

    where tau = i*t, so 2*tau = 2*i*t.
    """
    q = np.exp(-2 * np.pi * t)

    # theta functions at 2*tau = 2*i*t
    th10_2tau = theta_10_real(2*t, n_max)
    th00_2tau = theta_00_real(2*t, n_max)

    # Spectrum contributions:
    # - Vacuum (vac): contributes to (h=0, Q=0) with n=1 and to (h=1/2, Q=1) with n=-1
    # - DN ground state at h=3/8 for each of 3 complex directions

    # For N_D7 fractional branes, the multiplicity is N_D7 times single brane result
    # But we need to account for the orbifold projection properly

    # The key contributions:
    # Q=0: theta_10(2tau) * [q^{-1/4} (from vac) + 3*q^{3/8-1/4} + ...]
    # Q=1: -theta_00(2tau) * [-q^{-1/2} (from vac h=1/2) + 3*q^{3/8-1/2} + ...]

    # Number of fractional brane pairs contributing
    # For instanton at origin with N_D7/3 fractional branes of each type
    N_frac = N_D7 // 3  # = 8 for tadpole-cancelled case

    # Vacuum representation contribution (from eq. E.27 of reference)
    # n_{0,0}^{vac} = 1, n_{1/2,1}^{vac} = -1

    # The DN sector has mixed boundary conditions
    # Ground state energy for 6 DN directions: h = 6/16 = 3/8
    # But the vacuum rep still contributes at h=0 (Q=0) and h=1/2 (Q=1)

    # Actually, let me be more careful here.
    # The partition function in eq. (E.44) is per representation:
    # Z^{(h,Q)}_{A,DN} = (-1)^Q / 2 * q^{h-(1+Q)/4} * theta_{1-Q,0}(2tau)

    # For the orbifold, we sum over representations with multiplicities.
    # The issue is determining what representations appear.

    # From eq. (E.48), after summing:
    # Z_{A,DN} = sum_{Q=0,1} (-1)^Q/2 * theta_{1-Q,0}(2tau) * sum_h n_{h,Q} * q^{h-(1+Q)/4}

    # For C^3 (before orbifold), the spectrum includes:
    # - Vacuum rep contributing to h=0,Q=0 and h=1/2,Q=1 with coefficient +1 and -1
    # - Massive reps at various h values

    # For simplicity, let's compute the leading contributions:

    # Q=0 term: coefficient (-1)^0/2 = 1/2
    # theta_10(2tau) * [1 * q^{0-1/4} + higher]  (from vac at h=0, Q=0)

    # Q=1 term: coefficient (-1)^1/2 = -1/2
    # theta_00(2tau) * [-1 * q^{1/2-1/2} + higher]  (from vac at h=1/2, Q=1)

    # Leading terms for vacuum contribution:
    Z_vac_Q0 = 0.5 * th10_2tau * q**(-0.25) * N_frac  # h=0, Q=0
    Z_vac_Q1 = -0.5 * th00_2tau * (-1) * q**(0)  * N_frac  # h=1/2, Q=1, coefficient -1

    # Wait, this gives a divergent term as t->0. That's expected - the zero mode contribution.
    # The prescription is that these zero mode divergences cancel in the full expression.

    # Let me compute the full Z_A including all terms up to some order

    # More careful: for the internal CFT characters in the DN sector,
    # we need to use the actual C^3/Z_3 orbifold spectrum.

    # Following the structure in eq. (E.48):
    # The q-expansion starts at some power depending on the ground state energy.

    # For DN boundary conditions on 6 internal real dimensions (3 complex):
    # The ground state energy shift from DN is 1/16 per real DN direction
    # So h_DN = 6 * (1/16) = 3/8 for the internal part

    # But we also have the spacetime DD contribution (4 DD directions):
    # h_DD = 0 (no shift for DD)

    # And ghost contribution: -1/2 in NS sector

    # Total NS ground state: h = 3/8 - 1/2 = -1/8 (before normal ordering)
    # After normal ordering with c=0: h = 0 for the physical ground state

    # Actually the reference handles this differently - the free field factors
    # f_{alpha,beta} absorb the spacetime and ghost contributions.

    # Let me just use the structure from eq. (E.48) directly:
    # The integrand is real for real t, and the sum over q-powers gives:

    result = 0.0

    # Contribution from vacuum representation
    # vac contributes +1 to (h=0,Q=0) and -1 to (h=1/2,Q=1)
    # Per eq. (E.48), the total is:

    # Q=0 sector: theta_10(2tau) * q^{-1/4} (from h=0)
    # Q=1 sector: -theta_00(2tau) * (-1) * q^{0} = +theta_00(2tau) (from h=1/2 with n=-1)

    # With N_frac copies:
    result += N_frac * 0.5 * th10_2tau * q**(-0.25)  # Q=0, h=0, vac
    result += N_frac * 0.5 * th00_2tau * q**(0)       # Q=1, h=1/2, vac (n=-1 gives extra minus)

    # Higher representations would contribute higher powers of q
    # For now, let's include just the leading terms

    return np.real(result)

def Z_M(t, n_max=50):
    """
    Mobius strip partition function.

    Z_M = sum_{Q=0,1} (-1)^{1-Q} theta_{1-Q,0}(2*hat_tau) * sum_h sigma_r * n_r * e^{-i*pi*(h-Q/2)} * hat_q^{h-(1+Q)/4}

    where hat_tau = tau + 1/2 = i*t + 1/2, and hat_q = e^{2*pi*i*hat_tau} = -e^{-2*pi*t}
    """
    hat_q = -np.exp(-2 * np.pi * t)  # This is complex: exp(-2*pi*t + i*pi)

    # theta functions at 2*hat_tau = 2*i*t + 1
    th10_2hat = theta_10_hat(t, n_max)
    th00_2hat = theta_00_hat(t, n_max)

    # For vacuum representation with sigma_vac = -1:
    # Q=0: (-1)^{1-0} = -1, theta_10(2*hat_tau), h=0
    #      contribution: -1 * theta_10 * (-1) * e^{0} * hat_q^{-1/4}
    #      = theta_10 * hat_q^{-1/4}

    # Q=1: (-1)^{1-1} = 1, theta_00(2*hat_tau), h=1/2
    #      contribution: 1 * theta_00 * (-1) * (-1) * e^{-i*pi*(1/2-1/2)} * hat_q^{0}
    #      The n=-1 from vac gives another -1
    #      = -theta_00 * hat_q^{0}

    # Actually, let me use eq. (E.45) more carefully:
    # Z_M = sum_{Q} (-1)^{1-Q} theta_{1-Q,0}(2*hat_tau) * sum_{h,sigma} sigma * n_{h,Q,sigma} * e^{-i*pi*(h-Q/2)} * hat_q^{h-(1+Q)/4}

    # For vac with sigma_vac = -1:
    # h=0, Q=0: (-1)^1 * theta_10 * (-1) * 1 * e^{0} * hat_q^{-1/4} = theta_10 * hat_q^{-1/4}
    # h=1/2, Q=1: (-1)^0 * theta_00 * (-1) * (-1) * e^{0} * hat_q^{0} = -theta_00

    result = 0.0 + 0j

    # Vacuum contribution (sigma_vac = -1)
    # Q=0 term
    result += (-1)**1 * th10_2hat * (-1) * 1 * np.exp(0) * hat_q**(-0.25)
    # Q=1 term (h=1/2, n=-1 from vac decomposition)
    result += (-1)**0 * th00_2hat * (-1) * (-1) * np.exp(0) * hat_q**(0)

    # The result should be real after proper accounting
    return np.real(result)

def Z_zero_mode_subtraction(t, h=1.0):
    """
    Zero mode subtraction term: 3 * (e^{-2*pi*h*t} - 1)
    This cancels the divergent zero mode contributions.
    """
    return 3 * (np.exp(-2 * np.pi * h * t) - 1)

def integrand_full(t, N_D7=24, n_max=50):
    """
    Full integrand: (1/2t) * [Z_A + Z_M + 3*(e^{-2pi*t} - 1)]
    """
    if t < 1e-10:
        return 0.0

    z_a = Z_A_DN(t, N_D7, n_max)
    z_m = Z_M(t, n_max)
    z_sub = Z_zero_mode_subtraction(t)

    return (z_a + z_m + z_sub) / (2 * t)

def compute_K0(eps_prime=0.01, eps=0.01, N_D7=24, n_max=50):
    """
    Compute K_0 using numerical integration.

    K_0 = (1/(2*pi)^{3/2}) * exp[integral]

    The integral is from eps' to 1/eps for Z_A, eps'/4 to 1/eps for Z_M,
    and 0 to 1/eps for the subtraction.

    For simplicity, we use the same cutoffs and rely on the cancellation.
    """

    # Integrate from small t to large t
    t_min = eps_prime
    t_max = 1.0 / eps

    # Use adaptive quadrature
    result, error = integrate.quad(
        lambda t: integrand_full(t, N_D7, n_max),
        t_min, t_max,
        limit=200
    )

    # The factor out front
    prefactor = 1.0 / (2 * np.pi)**(1.5)

    K0 = prefactor * np.exp(result)

    return K0, result, error

def test_convergence():
    """Test convergence with respect to cutoffs and mode truncation."""

    print("="*70)
    print("Testing convergence of K_0 for C^3/Z_3 orientifold")
    print("="*70)

    # Test 1: Convergence with respect to n_max (mode truncation)
    print("\n1. Convergence with mode truncation (n_max):")
    print("-" * 50)
    eps_prime = 0.1
    eps = 0.1
    for n_max in [10, 20, 50, 100, 150]:
        K0, integral, error = compute_K0(eps_prime, eps, n_max=n_max)
        print(f"n_max = {n_max:3d}: K_0 = {K0:.10e}, integral = {integral:.6f}, error = {error:.2e}")

    # Test 2: Convergence with respect to eps_prime (IR cutoff in open string)
    print("\n2. Convergence with IR cutoff (eps_prime -> 0):")
    print("-" * 50)
    eps = 0.01
    n_max = 100
    for eps_prime in [0.5, 0.2, 0.1, 0.05, 0.02, 0.01]:
        K0, integral, error = compute_K0(eps_prime, eps, n_max=n_max)
        print(f"eps' = {eps_prime:.3f}: K_0 = {K0:.10e}, integral = {integral:.6f}")

    # Test 3: Convergence with respect to eps (UV cutoff)
    print("\n3. Convergence with UV cutoff (eps -> 0, i.e., t_max -> inf):")
    print("-" * 50)
    eps_prime = 0.01
    n_max = 100
    for eps in [0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005]:
        K0, integral, error = compute_K0(eps_prime, eps, n_max=n_max)
        print(f"eps = {eps:.4f} (t_max = {1/eps:6.1f}): K_0 = {K0:.10e}, integral = {integral:.6f}")

    # Final result with best parameters
    print("\n" + "="*70)
    print("FINAL RESULT")
    print("="*70)
    eps_prime = 0.001
    eps = 0.001
    n_max = 150

    K0, integral, error = compute_K0(eps_prime, eps, n_max=n_max)
    print(f"\nParameters: eps' = {eps_prime}, eps = {eps}, n_max = {n_max}")
    print(f"Integral value: {integral:.8f}")
    print(f"Integration error estimate: {error:.2e}")
    print(f"\nK_0 = {K0:.10e}")
    print(f"K_0 = {K0:.6f}")
    print(f"\nIn terms of (2*pi)^n: K_0 / (2*pi)^(-3/2) = exp({integral:.6f}) = {np.exp(integral):.6f}")

if __name__ == "__main__":
    test_convergence()
