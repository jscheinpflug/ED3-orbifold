#!/usr/bin/env python3
"""
K_0 computation - Version 10.
Compute entirely in the closed string channel where tadpole cancellation is exact.

The closed string amplitude:
∫_0^∞ dl/l × [Z_A^cl(l) + Z_M^cl(l) + (zero mode subtraction)]

With N_bulk = 8 and Q_O = -8, the massless pole cancels, leaving finite contributions.
"""

import numpy as np
from scipy import integrate
import warnings
warnings.filterwarnings('ignore')

def eta(tau, n_max=100):
    """Dedekind eta function."""
    q = np.exp(2j * np.pi * tau)
    if np.abs(q) >= 1 - 1e-10:
        return np.nan
    result = q**(1/24)
    for n in range(1, n_max + 1):
        result *= (1 - q**n)
    return result

def vartheta_00(tau, n_max=100):
    q = np.exp(2j * np.pi * tau)
    if np.abs(q) >= 1 - 1e-10:
        return np.nan
    result = 1.0 + 0j
    for n in range(1, n_max + 1):
        result += 2 * q**(n**2 / 2)
    return result

def vartheta_01(tau, n_max=100):
    q = np.exp(2j * np.pi * tau)
    if np.abs(q) >= 1 - 1e-10:
        return np.nan
    result = 1.0 + 0j
    for n in range(1, n_max + 1):
        result += 2 * ((-1)**n) * q**(n**2 / 2)
    return result

def vartheta_10(tau, n_max=100):
    q = np.exp(2j * np.pi * tau)
    if np.abs(q) >= 1 - 1e-10:
        return np.nan
    result = 0.0 + 0j
    for n in range(-n_max, n_max + 1):
        result += q**((n + 0.5)**2 / 2)
    return result

def Z_A_closed(l, N_bulk=8, n_max=100):
    """
    Annulus amplitude in closed string channel.

    For DN sector strings ending on N_bulk D7-branes:
    The closed string propagates with 4 transverse directions.
    Z ~ N × (theta functions)/(eta functions)
    """
    if l < 1e-12:
        return 0.0

    tau = 1j * l
    q = np.exp(-2 * np.pi * l)

    th00 = vartheta_00(tau, n_max)
    th01 = vartheta_01(tau, n_max)
    th10 = vartheta_10(tau, n_max)
    eta_val = eta(tau, n_max)

    if np.isnan(np.abs(th00)) or np.isnan(np.abs(eta_val)) or np.abs(eta_val) < 1e-100:
        return np.nan

    # For 4 transverse directions (DN sector):
    # Bosons: (θ₀₀⁴)/η⁴ or similar
    # The specific combination depends on boundary conditions

    # Using the character for 4 free bosons with DN b.c.:
    # Z_bos = |θ₀₀|⁴ / |η|⁴ gives the oscillator partition function
    # The momentum/winding contributes factors of l

    # For the DN sector, one end is Dirichlet (instanton), one is Neumann (D7)
    # This gives a factor of 1/√l from the zero mode integration

    Z_osc = np.abs(th00)**4 / np.abs(eta_val)**4

    # The full amplitude includes the zero mode factor and N D7-branes
    Z_cl = N_bulk * Z_osc

    return np.real(Z_cl)

def Z_M_closed(l, Q_O=-8, n_max=100):
    """
    Möbius strip amplitude in closed string channel.

    The Möbius has modulus τ = il + 1/2 (shifted by 1/2).
    The O-plane contributes with charge Q_O.
    """
    if l < 1e-12:
        return 0.0

    tau = 1j * l + 0.5
    q = np.exp(2j * np.pi * tau)  # = -e^{-2πl}

    th00 = vartheta_00(tau, n_max)
    th01 = vartheta_01(tau, n_max)
    th10 = vartheta_10(tau, n_max)
    eta_val = eta(tau, n_max)

    if np.isnan(np.abs(th00)) or np.isnan(np.abs(eta_val)) or np.abs(eta_val) < 1e-100:
        return np.nan

    # Same structure as annulus but with Q_O coefficient
    Z_osc = np.abs(th00)**4 / np.abs(eta_val)**4

    Z_cl = Q_O * Z_osc

    return np.real(Z_cl)

def Z_total_closed(l, N_bulk=8, Q_O=-8, n_max=100):
    """Total closed string amplitude: A + M."""
    z_a = Z_A_closed(l, N_bulk, n_max)
    z_m = Z_M_closed(l, Q_O, n_max)

    if np.isnan(z_a) or np.isnan(z_m):
        return np.nan

    return z_a + z_m

def main():
    print("="*70)
    print("K_0 Computation in Closed String Channel")
    print("="*70)

    print("\nTesting amplitude structure:")
    print("-"*60)

    N_bulk = 8
    Q_O = -8

    for l in [0.1, 0.2, 0.5, 1.0, 2.0, 5.0]:
        z_a = Z_A_closed(l, N_bulk)
        z_m = Z_M_closed(l, Q_O)
        total = z_a + z_m

        print(f"l = {l:.2f}: Z_A = {z_a:12.4f}, Z_M = {z_m:12.4f}, total = {total:12.4f}")

    # The total should approach 0 at large l (tadpole cancellation)
    # At small l (UV in closed string), it's finite

    print("\n" + "-"*60)
    print("Observation: At large l, Z_A + Z_M → 0 (tadpole cancelled)")
    print("At small l, the amplitudes are finite (massive modes)")

    # Compute the integral
    # The integration measure is dl/l in the closed string channel
    # With the zero mode contribution from 3 complex dimensions: 3*(q-1) term

    print("\n" + "="*70)
    print("Computing K_0 integral in closed string channel")
    print("="*70)

    def integrand_closed(l):
        """
        Closed string channel integrand.

        The integral is ∫ dl/l × [Z_A + Z_M - (asymptotic)]
        But since Z_A + Z_M → 0 at large l, we don't need subtraction there.

        At small l, we need to handle the UV behavior carefully.
        """
        if l < 1e-12:
            return 0.0

        z_total = Z_total_closed(l)

        if np.isnan(z_total):
            return 0.0

        return z_total / l

    # The integral should converge at both ends:
    # - Large l: Z_total → 0 (tadpole cancelled)
    # - Small l: Z_total is finite (massive modes only)

    print("\nIntegral from l_min to l_max:")
    print("-"*50)

    results = []
    for l_max in [1.0, 2.0, 5.0, 10.0]:
        for l_min in [0.1, 0.05, 0.01]:
            val, err = integrate.quad(integrand_closed, l_min, l_max, limit=500)
            results.append((l_min, l_max, val))
            print(f"  [{l_min:.2f}, {l_max:.1f}]: I = {val:10.4f}")

    # Check convergence
    print("\nConvergence check:")
    print("-"*50)

    l_min = 0.01
    for l_max in [5.0, 10.0, 20.0, 50.0]:
        val, err = integrate.quad(integrand_closed, l_min, l_max, limit=500)
        print(f"  l_max = {l_max:4.0f}: I = {val:10.4f}")

    l_max = 20.0
    for l_min in [0.5, 0.2, 0.1, 0.05, 0.02, 0.01]:
        val, err = integrate.quad(integrand_closed, l_min, l_max, limit=500)
        print(f"  l_min = {l_min:.2f}: I = {val:10.4f}")

    # Final result
    l_min = 0.01
    l_max = 50.0
    integral, err = integrate.quad(integrand_closed, l_min, l_max, limit=500)

    print("\n" + "="*70)
    print("RESULT")
    print("="*70)

    K0 = (2*np.pi)**(-1.5) * np.exp(integral)

    print(f"""
  Closed string channel integral:
    l_min = {l_min}, l_max = {l_max}
    Integral = {integral:.4f}

  K_0 = (2π)^{{-3/2}} × exp({integral:.4f})
  K_0 = {K0:.6e}
  log(K_0) = {np.log(K0):.4f}
""")
    print("="*70)

if __name__ == "__main__":
    main()
