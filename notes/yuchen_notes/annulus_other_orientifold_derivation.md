# Annulus Derivation For The Other Orientifold

This note derives the annulus partition function `Z_A` for the second
orientifold in [superpotential.tex](/Users/yuchenwang/Desktop/D-instanton/Superpotential/superpotential.tex:288),
namely

$$
(z_1,z_2,z_3)\mapsto -(z_1,z_2,z_3).
$$

The O7-plane and the 8 D7-branes wrap the exceptional divisor
$$
E \simeq \mathbb P^2.
$$

The main point is that, at the orbifold point, a space-filling D7 wrapping
`E` should be represented by a fractional `D3`, just as the ED3 wrapping `E`
is represented by a fractional `D(-1)`. So the relevant open-string system is
the fractional `D(-1)`-`D3` annulus, not the old `D(-1)`-flat-`D7` annulus from
the `z_3=0` orientifold.

I keep the fractional Chan-Paton multiplicities
$$
N_0,\qquad N_1,\qquad N_2
$$
general throughout.

## 1. General Form Of The Annulus

Take the ED3 instanton to be fractional type `0`. The D7 stack at the orbifold
point is then represented as
$$
N_0\widetilde D_0+N_1\widetilde D_1+N_2\widetilde D_2.
$$

For this orientifold, a single physical D7-brane wrapping the collapsed
exceptional divisor is represented at the orbifold point by the regular
representation
$$
\widetilde D_0+\widetilde D_1+\widetilde D_2.
$$
Thus eight physical D7-branes correspond to
$$
N_0=N_1=N_2=8,
$$
not to `N_0+N_1+N_2=8`.

The annulus therefore takes the same orbifold-projected form as before:
$$
Z_A[N_0,N_1,N_2]
=
\frac{1}{6}\sum_{k=0}^2
\Bigl(N_0+N_1\omega^{-k}+N_2\omega^k\Bigr)
\sum_{\alpha,\beta=0}^1 Z_{\alpha\beta}^k(\tau),
$$
with
$$
\omega=e^{2\pi i/3},
\qquad
q=e^{2\pi i\tau},
\qquad
\tau=it.
$$

So the task is to compute the `D(-1)`-`D3` building blocks
$$
Z_{00}^k,\quad Z_{01}^k,\quad Z_{10}^k,\quad Z_{11}^k.
$$

## 2. Open-String Spectrum

For the `D(-1)`-`D3` system:

- the four noncompact directions `X^\mu` are `ND`,
- the three internal complex directions `z_1,z_2,z_3` are all `DD`.

So:

- in the `NS` sector, the massless ground states come from the four `ND`
  fermion zero modes and form an `SO(4)` spinor,
- in the `R` sector, the massless ground states come from the three internal
  `DD` fermion zero modes and form an `SO(6)` spinor.

The orbifold acts as
$$
g:\quad (z_1,z_2,z_3)\mapsto (\omega z_1,\omega z_2,\omega z_3).
$$

As in [superpotential.tex](/Users/yuchenwang/Desktop/D-instanton/Superpotential/superpotential.tex:215),
I use the spin-field action
$$
(e^{iH_1/2},e^{iH_2/2},e^{iH_3/2})
\mapsto
(e^{i\pi/3}e^{iH_1/2},e^{i\pi/3}e^{iH_2/2},e^{-2i\pi/3}e^{iH_3/2}),
$$
so an internal spin field with weights `(s_1,s_2,s_3)` picks up phase
$$
\omega^{\,s_1+s_2-2s_3}.
$$

## 3. Definitions Of `Z_{\alpha\beta}^k`

I keep the same definitions as in the first annulus note:

$$
Z_{00}^k =
\operatorname{Tr}_{\mathcal H^{\rm NS}_{-1,3}}
\left[g^k(-1)^{N_{bc}+N_{\beta\gamma}}q^{L_0}\right],
$$

$$
Z_{01}^k =
\operatorname{Tr}_{\mathcal H^{\rm NS}_{-1,3}}
\left[g^k(-1)^{N_{bc}+N_{\beta\gamma}}(-1)^Fq^{L_0}\right],
$$

$$
Z_{10}^k =
-\frac12
\operatorname{Tr}_{\mathcal H^{\rm R}_{-1,3}}
\left[g^k(-1)^{N_{bc}+N_{\beta\gamma}}q^{L_0}\right],
$$

$$
Z_{11}^k =
-\frac12
\operatorname{Tr}_{\mathcal H^{\rm R}_{-1,3}}
\left[g^k(-1)^{N_{bc}+N_{\beta\gamma}}(-1)^Fq^{L_0}\right].
$$

## 4. NS Sector

The `NS` ground states are the `SO(4)` spin fields
$$
ce^{-\phi}S_{\nu_1,\nu_2},
\qquad
\nu_1,\nu_2=\pm \frac12.
$$

They are orbifold-neutral, and their unprojected trace gives
$$
\sum_{\nu_1,\nu_2} 1 = 4.
$$

With an extra `(-1)^F` insertion, the ground-state trace vanishes:
$$
\sum_{\nu_1,\nu_2}(2\nu_1)(2\nu_2)=0.
$$

Hence
$$
Z_{01}^k=0
$$
for all `k`.

The oscillator product in the `NS` sector is
$$
Z_{00}^k
=
4
\prod_{n=1}^\infty
\frac{(1+q^n)^4(1-q^n)^2}{(1-q^{n-\frac12})^4(1+q^{n-\frac12})^2}
\prod_{n=1}^\infty
\frac{(1+\omega^k q^{n-\frac12})^3(1+\omega^{-k} q^{n-\frac12})^3}
{(1-\omega^k q^n)^3(1-\omega^{-k} q^n)^3}.
$$

For `k=0`, this simplifies to
$$
Z_{00}^0
=
4\,
\prod_{n=1}^\infty
\frac{(1+q^n)^4(1+q^{n-\frac12})^4}
{(1-q^n)^4(1-q^{n-\frac12})^4}
=
\frac{4}{\theta_4(q)^4}.
$$

For `k=1,2`, using the product identities from
[annulus_Zabk_derivation.md](/Users/yuchenwang/Desktop/D-instanton/Superpotential/annulus_Zabk_derivation.md:678),
one gets
$$
Z_{00}^{k\neq 0}
=
i(\omega^{k/2}-\omega^{-k/2})^3
\eta(q)^3\,
\frac{\theta_2(q)^2}{\theta_4(q)^2\theta_3(q)}
\frac{\theta_3(\omega^k,q)^3}{\theta_1(\omega^k,q)^3}.
$$

Since the original product is symmetric under `\omega^k\leftrightarrow\omega^{-k}`,
this gives
$$
Z_{00}^1=Z_{00}^2
=
3\sqrt3\,\eta(q)^3\,
\frac{\theta_2(q)^2\theta_3(\omega,q)^3}
{\theta_4(q)^2\theta_3(q)\theta_1(\omega,q)^3}.
$$

## 5. R Sector Without `(-1)^F`

The unprojected `R` ground states are labelled by
$$
(s_1,s_2,s_3)\in\left\{\pm \frac12\right\}^3.
$$

Their orbifold phase is
$$
\omega^{\,s_1+s_2-2s_3}.
$$

So the `R` ground-state trace is
$$
C_k
:=
\sum_{s_1,s_2,s_3}
\omega^{k(s_1+s_2-2s_3)}
=
\left(\omega^{k/2}+\omega^{-k/2}\right)^2
\left(\omega^k+\omega^{-k}\right).
$$

Explicitly,
$$
C_0=8,
\qquad
C_1=C_2=-1.
$$

The full `R`-sector trace without `(-1)^F` is
$$
Z_{10}^k
=
-\frac12 C_k
\prod_{n=1}^\infty
\frac{(1+q^{n-\frac12})^4(1-q^n)^2}{(1-q^{n-\frac12})^4(1+q^n)^2}
\prod_{n=1}^\infty
\frac{(1+\omega^k q^n)^3(1+\omega^{-k} q^n)^3}
{(1-\omega^k q^n)^3(1-\omega^{-k} q^n)^3}.
$$

For `k=0`, this gives
$$
Z_{10}^0
=
-4\,
\prod_{n=1}^\infty
\frac{(1+q^{n-\frac12})^4(1+q^n)^4}
{(1-q^{n-\frac12})^4(1-q^n)^4}
=
-\frac{4}{\theta_4(q)^4}.
$$

So already
$$
Z_{00}^0+Z_{10}^0=0.
$$

For `k=1,2`, the theta-function form is
$$
Z_{10}^{k\neq 0}
=
-C_k\,i
\left(
\frac{\omega^{k/2}-\omega^{-k/2}}
{\omega^{k/2}+\omega^{-k/2}}
\right)^3
\eta(q)^3\,
\frac{\theta_3(q)^2}{\theta_4(q)^2\theta_2(q)}
\frac{\theta_2(\omega^k,q)^3}{\theta_1(\omega^k,q)^3}.
$$

For `k=1,2`, this reduces to
$$
Z_{10}^1=Z_{10}^2
=
3\sqrt3\,\eta(q)^3\,
\frac{\theta_3(q)^2\theta_2(\omega,q)^3}
{\theta_4(q)^2\theta_2(q)\theta_1(\omega,q)^3}.
$$

## 6. R Sector With `(-1)^F`

Now insert
$$
(-1)^F=(2s_1)(2s_2)(2s_3).
$$

The ground-state sum is
$$
D_k
:=
\sum_{s_1,s_2,s_3}
(2s_1)(2s_2)(2s_3)\,
\omega^{k(s_1+s_2-2s_3)}
$$
$$
=
\left(\omega^{k/2}-\omega^{-k/2}\right)^2
\left(\omega^{-k}-\omega^k\right).
$$

Therefore
$$
D_0=0,
\qquad
D_1=3i\sqrt3,
\qquad
D_2=-3i\sqrt3.
$$

With `(-1)^F`, the `R`-sector oscillator factors cancel completely, so
$$
Z_{11}^k=-\frac12 D_k.
$$

Hence
$$
Z_{11}^0=0,
\qquad
Z_{11}^1=-\frac{3i\sqrt3}{2},
\qquad
Z_{11}^2=+\frac{3i\sqrt3}{2}.
$$

## 7. Spin-Structure Sums

Collecting everything:

$$
\sum_{\alpha,\beta} Z_{\alpha\beta}^0=0.
$$

Define
$$
\mathcal F(q)
:=
3\sqrt3\,\eta(q)^3\,
\frac{1}{\theta_4(q)^2\theta_1(\omega,q)^3}
\left[
\frac{\theta_2(q)^2\theta_3(\omega,q)^3}{\theta_3(q)}
+
\frac{\theta_3(q)^2\theta_2(\omega,q)^3}{\theta_2(q)}
\right].
$$

Then
$$
\sum_{\alpha,\beta} Z_{\alpha\beta}^1
=
\mathcal F(q)-\frac{3i\sqrt3}{2},
$$

$$
\sum_{\alpha,\beta} Z_{\alpha\beta}^2
=
\mathcal F(q)+\frac{3i\sqrt3}{2}.
$$

## 8. Final Formula For General `N_0,N_1,N_2`

Substituting into the orbifold-projected annulus gives
$$
Z_A[N_0,N_1,N_2]
=
\frac{1}{6}
\Bigl(N_0+N_1\omega^2+N_2\omega\Bigr)
\left(\mathcal F(q)-\frac{3i\sqrt3}{2}\right)
$$
$$
+\frac{1}{6}
\Bigl(N_0+N_1\omega+N_2\omega^2\Bigr)
\left(\mathcal F(q)+\frac{3i\sqrt3}{2}\right).
$$

Using
$$
\omega+\omega^2=-1,
\qquad
\omega-\omega^2=i\sqrt3,
$$
this simplifies to
$$
\boxed{
Z_A[N_0,N_1,N_2]
=
\frac{2N_0-N_1-N_2}{6}\,\mathcal F(q)
-\frac{3}{4}(N_1-N_2).
}
$$

For the physical eight-D7 stack, `N_0=N_1=N_2=8`, and hence
$$
\boxed{Z_A=0.}
$$

## 9. Immediate Checks

The `k=0` sector cancels identically:
$$
Z_{00}^0+Z_{10}^0+Z_{11}^0=0.
$$

This is the expected cancellation between the `D(-1)`-`D3` bosonic `w` zero
modes and the fermionic `\mu` zero modes before orbifolding.

The full annulus therefore comes entirely from the twisted sectors `k=1,2`,
which is exactly what one should expect for fractional branes on
`\mathbb C^3/\mathbb Z_3`.

The large-`t` limit also matches the zero-mode supertrace. This coefficient is
not the raw number of zero modes: a bosonic zero mode contributes `+1`, while a
fermionic zero mode contributes `-1/2` in the `K_0` normalization.

In the NS sector, the four ground states are orbifold-neutral, and the GSO
projection keeps half of them. Thus only the `N_0` Chan-Paton sector contributes,
with
$$
\frac12(4N_0)=2N_0.
$$
The GSO-even R ground states have phases
$$
1,\qquad \omega,\qquad \omega,\qquad \omega.
$$
Thus the `N_0` sector keeps one fermionic zero mode, the `N_1` sector keeps
three, and the `N_2` sector keeps none. The total open-channel weighted
zero-mode coefficient is therefore
$$
Z_A(t\to\infty)
=2N_0-\frac12(N_0+3N_1)
=\frac32(N_0-N_1).
$$

For the physical regular stack `N_0=N_1=N_2=8`, both the small-`t` and large-`t`
annulus contributions vanish.
