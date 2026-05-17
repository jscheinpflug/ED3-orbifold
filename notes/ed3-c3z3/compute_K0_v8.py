#!/usr/bin/env python3
"""
K_0 computation - Version 8.
Compute in the CLOSED string channel where tadpole cancellation is manifest.

The modular transformations:
- Annulus: τ_open = it → τ_closed = i/(2t) = il where l = 1/(2t)
- Möbius: hat_τ_open = it + 1/2 → l = 1/(8t)

In the closed string channel at large l (small t in open), the massless
modes dominate. Tadpole cancellation requires the massless poles to cancel
between annulus and Möbius.
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

def Z_annulus_closed(l, N=24, n_max=100):
    """
    Annulus amplitude in closed string channel.

    For the DN sector (D7 to D(-1) strings):
    - Transform from open string modulus t to closed string modulus l = 1/(2t)
    - The theta functions transform under S: θ(τ) → (-iτ)^{-1/2} θ(-1/τ)

    The closed string channel amplitude at modulus il:
    Z_A^{closed}(l) = 2l × Z_A^{open}(t=1/(2l))

    where the factor 2l comes from the Jacobian dt = -dl/(2l²).
    """
    if l < 1e-12:
        return 0.0

    tau_cl = 1j * l  # Closed string modulus

    th00 = vartheta_00(tau_cl, n_max)
    th01 = vartheta_01(tau_cl, n_max)
    th10 = vartheta_10(tau_cl, n_max)
    eta_val = eta(tau_cl, n_max)

    if np.isnan(np.abs(th00)) or np.isnan(np.abs(eta_val)) or np.abs(eta_val) < 1e-100:
        return np.nan

    # The closed string amplitude for DN sector
    # At large l (massless modes): Z ~ N × l × (constant)
    # The pole in l represents the massless closed string propagator

    # Using the transformation formulas for theta functions:
    # The open string partition function transforms as:
    # Z_open(t) → Z_closed(l) = (1/l) × Z_open(1/(2l))
    # (the factor of 1/l comes from the proper time measure transformation)

    # For the spectrum: 4 transverse bosons, 2 fermion zero modes
    # Z ~ (θ₀₀⁴ - θ₀₁⁴)/η¹² for the bosonic oscillators
    # (this is a simplified expression; the actual formula involves the DN boundary conditions)

    # Massless contribution at large l:
    # The theta functions approach 1 at large l, and η ~ q^{1/24} → 0
    # So Z ~ 1/η¹² → ∞, which represents the massless pole

    # For our purposes, we extract the coefficient of the massless pole:
    # Z = c × l + (massive contributions)

    # The coefficient c must cancel between A and M for tadpole cancellation
    # c_A = N, c_M = Q_O = -8 (in bulk units)

    # For now, return a simplified closed string amplitude
    # that captures the massless pole structure
    Z_cl = N * np.abs(th00)**4 / np.abs(eta_val)**8

    return np.real(Z_cl)

def Z_mobius_closed(l, Q_O=-8, n_max=100):
    """
    Möbius strip amplitude in closed string channel.

    The Möbius strip transforms with modulus l = 1/(8t).
    The crosscap state involves a factor related to the O-plane charge.
    """
    if l < 1e-12:
        return 0.0

    # Möbius strip closed string modulus includes the 1/2 shift
    hat_tau_cl = 1j * l + 0.5

    th00 = vartheta_00(hat_tau_cl, n_max)
    th01 = vartheta_01(hat_tau_cl, n_max)
    th10 = vartheta_10(hat_tau_cl, n_max)
    eta_val = eta(hat_tau_cl, n_max)

    if np.isnan(np.abs(th00)) or np.isnan(np.abs(eta_val)) or np.abs(eta_val) < 1e-100:
        return np.nan

    # The Möbius amplitude with O-plane contribution
    # Crosscap state normalization gives factor Q_O
    Z_cl = Q_O * np.abs(th00)**4 / np.abs(eta_val)**8

    return np.real(Z_cl)

def test_closed_string_tadpole():
    """
    Test tadpole cancellation in the closed string channel.

    At large l, the massless pole should cancel:
    Z_A(l) + Z_M(l) ~ (N + Q_O) × l + O(1)

    For N = 8 bulk (= 24 fractional) and Q_O = -8, this gives 0.
    """
    print("="*70)
    print("Testing tadpole in closed string channel")
    print("="*70)

    print("\nWith N = 24 (fractional), Q_O = -8 (bulk):")
    print("Expect: coefficient of massless pole ~ N_bulk + Q_O = 8 - 8 = 0")
    print("-"*60)

    # Note: N = 24 fractional = 8 bulk
    # The closed string channel should see the bulk charge N_bulk = 8

    N_bulk = 8  # In closed string channel, we use bulk counting

    for l in [0.5, 1.0, 2.0, 5.0, 10.0]:
        z_a = Z_annulus_closed(l, N=N_bulk)  # Use bulk counting
        z_m = Z_mobius_closed(l, Q_O=-8)
        total = z_a + z_m

        # At large l, total/l should approach 0 for tadpole cancellation
        coeff = total / l if l > 0 else 0

        print(f"l = {l:5.1f}: Z_A = {z_a:12.2f}, Z_M = {z_m:12.2f}, "
              f"total = {total:12.2f}, total/l = {coeff:10.2f}")

    print("\n" + "-"*60)
    print("Note: The ratio total/l should approach a constant at large l.")
    print("This constant = (N_bulk + Q_O) × (massless state normalization).")
    print("For tadpole cancellation, we need N_bulk + Q_O = 0, i.e., N_bulk = 8.")

def main():
    test_closed_string_tadpole()

    print("\n" + "="*70)
    print("Conclusion")
    print("="*70)

    print("""
The closed string channel analysis shows that with:
- N = 8 bulk D7-branes (= 24 fractional)
- Q_O = -8 O7-plane charge

The total contribution Z_A + Z_M has:
- A massless pole (proportional to l) that may or may not cancel
- The ratio total/l at large l indicates the tadpole coefficient

For a well-defined K_0 integral, the massless pole must cancel.
If it doesn't cancel exactly with my simplified formulas, more careful
treatment of the orbifold structure and boundary states is needed.
    """)

if __name__ == "__main__":
    main()
