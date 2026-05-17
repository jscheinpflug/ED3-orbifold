#!/usr/bin/env python3
"""
K_0 computation - Version 9.
Use bulk normalization (N_bulk = 8, Q_O = -8) for tadpole cancellation.

The closed string channel analysis (v8) showed exact cancellation with these values.
Now apply the same normalization in the open string channel.
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

def Z_A_DN(t, N_bulk=8, n_max=100):
    """
    Annulus partition function using bulk D7-brane count.
    N_bulk = 8 corresponds to 24 fractional branes (×3 from orbifold).
    """
    tau_2 = 2j * t
    q = np.exp(-2 * np.pi * t)

    th10 = np.real(vartheta_10(tau_2, n_max))
    th00 = np.real(vartheta_00(tau_2, n_max))

    if np.isnan(th10) or np.isnan(th00):
        return np.nan

    # Using bulk counting: N_bulk = 8
    # The factor of 3 from orbifold sectors is implicit in the spectrum
    Z_Q0 = 0.5 * th10 * q**(-0.25) * N_bulk
    Z_Q1 = 0.5 * th00 * 1.0 * N_bulk

    return Z_Q0 + Z_Q1

def Z_M(t, Q_O=-8, n_max=100):
    """
    Möbius strip partition function with O-plane charge.
    Q_O = -8 for O7-plane in bulk units.
    """
    hat_tau_2 = 2 * (1j * t + 0.5)
    hat_q = np.exp(-2 * np.pi * t + 1j * np.pi)

    th10_hat = vartheta_10(hat_tau_2, n_max)
    th00_hat = vartheta_00(hat_tau_2, n_max)

    if np.isnan(np.abs(th10_hat)) or np.isnan(np.abs(th00_hat)):
        return np.nan

    Z_Q0 = th10_hat * hat_q**(-0.25)
    Z_Q1 = -th00_hat * 1.0

    # Include O-plane charge
    # Note: The coefficient should be |Q_O| = 8 with the appropriate sign
    # from the crosscap state overlap
    return Q_O * np.real(Z_Q0 + Z_Q1)

def Z_sub(t):
    """Zero mode subtraction term."""
    return 3 * (np.exp(-2 * np.pi * t) - 1)

def Z_total(t, N_bulk=8, Q_O=-8, n_max=100):
    z_a = Z_A_DN(t, N_bulk, n_max)
    z_m = Z_M(t, Q_O, n_max)
    z_s = Z_sub(t)

    if np.isnan(z_a) or np.isnan(z_m):
        return np.nan

    return z_a + z_m + z_s

def compute_and_fit():
    print("="*70)
    print("K_0 Computation with bulk normalization")
    print("N_bulk = 8, Q_O = -8 (tadpole cancelling)")
    print("="*70)

    # Check partition function values
    print("\nPartition functions:")
    print("-"*60)
    for t in [0.01, 0.05, 0.1, 0.5, 1.0, 5.0]:
        z_a = Z_A_DN(t, N_bulk=8)
        z_m = Z_M(t, Q_O=-8)
        z_s = Z_sub(t)
        total = z_a + z_m + z_s
        print(f"t = {t:.3f}: Z_A = {z_a:10.4f}, Z_M = {z_m:10.4f}, "
              f"Z_sub = {z_s:8.4f}, total = {total:10.4f}")

    # Asymptotic constant
    Z_asymp = Z_total(10.0)
    print(f"\nZ(t → ∞) = {Z_asymp:.4f}")

    # Integrand
    def integrand(t):
        z = Z_total(t)
        if np.isnan(z):
            return 0.0
        return (z - Z_asymp) / (2 * t)

    # Compute integrals
    t_max = 20.0
    t_mins = np.array([0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005])
    integrals = []

    print(f"\nIntegrals (t_max = {t_max}):")
    print("-"*50)
    for t_min in t_mins:
        val, _ = integrate.quad(integrand, t_min, t_max, limit=500)
        integrals.append(val)
        print(f"  t_min = {t_min:.4f}: I = {val:12.4f}")

    integrals = np.array(integrals)

    # Fit: I = A/√t + B*log(t) + C
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
    print(f"  A = {A:.4f}  (tadpole coefficient)")
    print(f"  B = {B:.4f}  (log correction)")
    print(f"  C = {C:.4f}  (finite part)")
    print(f"  RMS = {rms:.4f}")

    return A, B, C

def main():
    A, B, C = compute_and_fit()

    print("\n" + "="*70)
    print("RESULT")
    print("="*70)

    K0 = (2*np.pi)**(-1.5) * np.exp(C)

    print(f"""
  Bulk normalization: N_bulk = 8, Q_O = -8

  Divergence coefficient A = {A:.4f}
  {"✓ Tadpole cancelled!" if abs(A) < 1 else "⚠ Tadpole not fully cancelled"}

  Finite part C = {C:.4f}

  K_0 = (2π)^{{-3/2}} × exp({C:.4f})
  K_0 = {K0:.6e}
  log(K_0) = {np.log(K0):.4f}
""")
    print("="*70)

if __name__ == "__main__":
    main()
