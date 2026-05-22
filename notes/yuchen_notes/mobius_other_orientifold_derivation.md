# Mobius Derivation For The Other Orientifold

This note derives the Möbius-strip partition function `Z_M` for the second
orientifold in [superpotential.tex](/Users/yuchenwang/Desktop/D-instanton/Superpotential/superpotential.tex:288),
namely

$$
(z_1,z_2,z_3)\mapsto -(z_1,z_2,z_3).
$$

I use the same conventions as in
[mobius_ZMabk_derivation.md](/Users/yuchenwang/Desktop/D-instanton/Superpotential/mobius_ZMabk_derivation.md),
but now the geometric `\mathbb Z_2` reflection acts on all three internal
complex coordinates.

For the standard unweighted crosscap, the final answer is

$$
\boxed{
Z_M
=
2\sqrt3\,\eta(\hat q)^3
\left[
\frac{\theta_4(\hat q)^2\theta_3(\omega,\hat q)^3}
{\theta_2(\hat q)^2\theta_3(\hat q)\theta_1(\omega,\hat q)^3}
+
\frac{\theta_3(\hat q)^2\theta_4(\omega,\hat q)^3}
{\theta_2(\hat q)^2\theta_4(\hat q)\theta_1(\omega,\hat q)^3}
\right],
}
$$

with

$$
\hat q=e^{2\pi i(\tau+\frac12)}.
$$

Unlike the first orientifold, there is no untwisted constant term: the `k=0`
sector cancels identically, and the standard Möbius strip is built entirely
from the twisted sectors.

## 1. Definitions

We define

$$
Z_M
=
\frac{1}{12}\sum_{k=0}^2\sum_{\alpha,\beta=0}^1
Z_{M,\alpha\beta}^k(\tau),
$$

with

$$
Z_{M,00}^k
=
\operatorname{Tr}_{\mathcal H_{\rm inst}^{\rm NS}}
\left[g^k\Omega(-1)^{N_{bc}+N_{\beta\gamma}}q^{L_0}\right],
$$

$$
Z_{M,01}^k
=
\operatorname{Tr}_{\mathcal H_{\rm inst}^{\rm NS}}
\left[g^k\Omega(-1)^{N_{bc}+N_{\beta\gamma}}(-1)^Fq^{L_0}\right],
$$

$$
Z_{M,10}^k
=
-\frac12
\operatorname{Tr}_{\mathcal H_{\rm inst}^{\rm R}}
\left[g^k\Omega(-1)^{N_{bc}+N_{\beta\gamma}}q^{L_0}\right],
$$

$$
Z_{M,11}^k
=
-\frac12
\operatorname{Tr}_{\mathcal H_{\rm inst}^{\rm R}}
\left[g^k\Omega(-1)^{N_{bc}+N_{\beta\gamma}}(-1)^Fq^{L_0}\right].
$$

Here

$$
g:\quad (z_1,z_2,z_3)\mapsto (\omega z_1,\omega z_2,\omega z_3),
\qquad
\omega=e^{2\pi i/3}.
$$

As usual, absorb the universal oscillator sign of `\Omega` by shifting

$$
\tau\mapsto \hat\tau:=\tau+\frac12,
\qquad
\hat q:=e^{2\pi i\hat\tau}=-q.
$$

## 2. Residual `\mathbb Z_2` Action After The `\hat\tau` Shift

For `D(-1)` open strings, the universal Dirichlet artifact after the `\hat\tau`
shift contributes a minus sign to all coordinates. The geometric involution of
the present orientifold flips all three internal complex coordinates as well.

Therefore the residual insertion acts on oscillator modes as

$$
(X^\mu,Z^1,Z^2,Z^3)\mapsto (-X^\mu,+Z^1,+Z^2,+Z^3),
$$

while

$$
(b,c,\beta,\gamma)\mapsto (b,c,\beta,\gamma).
$$

This is the key difference from the first orientifold.

## 3. Ground-State Charges

### 3.1 NS Ground State

The NS ground state is

$$
ce^{-\phi}|0\rangle,
\qquad
L_0=-\frac12.
$$

It is odd under the residual `\mathbb Z_2`, neutral under `\mathbb Z_3`, and
has the standard

$$
(-1)^F=-1.
$$

Hence the NS ground-state factors are

$$
G_{\rm NS,00}=-\hat q^{-1/2},
\qquad
G_{\rm NS,01}=+\hat q^{-1/2}.
$$

### 3.2 R Ground States

The R ground states are

$$
ce^{-\phi/2}S_{\mathbf s}|0\rangle,
\qquad
\mathbf s=(\nu_1,\nu_2,s_1,s_2,s_3),
$$

with

$$
\nu_i,s_i=\pm \frac12,
\qquad
L_0=0.
$$

Under the orbifold,

$$
S_{\mathbf s}\mapsto
\omega^{s_1+s_2-2s_3}S_{\mathbf s}.
$$

Since the residual insertion flips only the four noncompact directions, its
action on the spin ground state is proportional to

$$
(2\nu_1)(2\nu_2).
$$

The overall sign of this lift is conventional and will not affect the standard
unweighted Möbius sum, so I leave it implicit.

Then the R ground-state trace without `(-1)^F` vanishes:

$$
\sum_{\nu_1,\nu_2}(2\nu_1)(2\nu_2)=0,
$$

and therefore

$$
Z_{M,10}^k=0.
$$

For `Z_{M,11}^k`, after multiplying by

$$
(-1)^F=(2\nu_1)(2\nu_2)(2s_1)(2s_2)(2s_3),
$$

the `\nu` factors square to `1`, leaving a factor `4` from the sum over
`\nu_1,\nu_2`. So the relevant internal zero-mode sum is

$$
D_k
:=
\sum_{s_1,s_2,s_3}
(2s_1)(2s_2)(2s_3)\,
\omega^{k(s_1+s_2-2s_3)}.
$$

This evaluates to

$$
D_k
=
\left(\omega^{k/2}-\omega^{-k/2}\right)^2
\left(\omega^{-k}-\omega^k\right).
$$

Hence

$$
D_0=0,
\qquad
D_1=3i\sqrt3,
\qquad
D_2=-3i\sqrt3.
$$

Up to the conventional overall sign choice discussed above, this gives

$$
Z_{M,11}^0=0,
\qquad
Z_{M,11}^1=\mp 6i\sqrt3,
\qquad
Z_{M,11}^2=\pm 6i\sqrt3.
$$

For the standard unweighted Möbius strip, these two twisted `R` pieces cancel
between `k=1` and `k=2`, so this sign ambiguity will drop out of `Z_M`.

## 4. Oscillator Products

After the `\hat\tau` shift:

- `X^\mu` has residual charge `-1`,
- each `Z^i` has residual charge `+1`,
- ghosts are neutral.

Therefore the NS-sector products are

$$
Z_{M,00}^k
=
-\hat q^{-1/2}
\prod_{n=1}^\infty
\frac{(1-\hat q^n)^2}
{(1+\hat q^n)^4(1-\omega^k\hat q^n)^3(1-\omega^{-k}\hat q^n)^3}
$$
$$
\times
\prod_{r>0}
\frac{(1-\hat q^r)^4(1+\omega^k\hat q^r)^3(1+\omega^{-k}\hat q^r)^3}
{(1+\hat q^r)^2},
$$

$$
Z_{M,01}^k
=
\hat q^{-1/2}
\prod_{n=1}^\infty
\frac{(1-\hat q^n)^2}
{(1+\hat q^n)^4(1-\omega^k\hat q^n)^3(1-\omega^{-k}\hat q^n)^3}
$$
$$
\times
\prod_{r>0}
\frac{(1+\hat q^r)^4(1-\omega^k\hat q^r)^3(1-\omega^{-k}\hat q^r)^3}
{(1-\hat q^r)^2}.
$$

As already explained,

$$
Z_{M,10}^k=0.
$$

The `Z_{M,11}^k` oscillator product cancels completely, so only the zero-mode
constants above remain.

## 5. The `k=0` Sector

For `k=0`, the two NS products are identical up to the ground-state sign:

$$
Z_{M,01}^0=-Z_{M,00}^0.
$$

Since also

$$
Z_{M,11}^0=0,
$$

the full untwisted contribution vanishes:

$$
\sum_{\alpha,\beta} Z_{M,\alpha\beta}^0=0.
$$

So this orientifold has no `k=0` Möbius contribution.

## 6. Twisted Sectors In Theta-Function Form

From the standard product identities used in
[mobius_ZMabk_derivation.md](/Users/yuchenwang/Desktop/D-instanton/Superpotential/mobius_ZMabk_derivation.md),
for `k=1,2` one finds

$$
Z_{M,00}^1=Z_{M,00}^2
=
12\sqrt3\,
\eta(\hat q)^3
\frac{\theta_4(\hat q)^2\theta_3(\omega,\hat q)^3}
{\theta_2(\hat q)^2\theta_3(\hat q)\theta_1(\omega,\hat q)^3},
$$

and

$$
Z_{M,01}^1=Z_{M,01}^2
=
12\sqrt3\,
\eta(\hat q)^3
\frac{\theta_3(\hat q)^2\theta_4(\omega,\hat q)^3}
{\theta_2(\hat q)^2\theta_4(\hat q)\theta_1(\omega,\hat q)^3}.
$$

The `R` pieces satisfy

$$
Z_{M,11}^1+Z_{M,11}^2=0
$$

for the standard unweighted sum.

Therefore

$$
\sum_{\alpha,\beta} Z_{M,\alpha\beta}^1
=
12\sqrt3\,\eta(\hat q)^3
\left[
\frac{\theta_4(\hat q)^2\theta_3(\omega,\hat q)^3}
{\theta_2(\hat q)^2\theta_3(\hat q)\theta_1(\omega,\hat q)^3}
+
\frac{\theta_3(\hat q)^2\theta_4(\omega,\hat q)^3}
{\theta_2(\hat q)^2\theta_4(\hat q)\theta_1(\omega,\hat q)^3}
\right],
$$

and the same for `k=2`.

## 7. Final Formula For The Standard Crosscap

Since the `k=0` sector vanishes and the `k=1,2` sectors are equal,

$$
Z_M
=
\frac{1}{12}
\left(
\sum_{\alpha,\beta} Z_{M,\alpha\beta}^1
+
\sum_{\alpha,\beta} Z_{M,\alpha\beta}^2
\right)
$$

which gives

$$
\boxed{
Z_M
=
2\sqrt3\,\eta(\hat q)^3
\left[
\frac{\theta_4(\hat q)^2\theta_3(\omega,\hat q)^3}
{\theta_2(\hat q)^2\theta_3(\hat q)\theta_1(\omega,\hat q)^3}
+
\frac{\theta_3(\hat q)^2\theta_4(\omega,\hat q)^3}
{\theta_2(\hat q)^2\theta_4(\hat q)\theta_1(\omega,\hat q)^3}
\right].
}
$$

This is the Möbius-strip partition function for the second orientifold with the
standard unweighted orbifold sum.

## 8. Comparison With The First Orientifold

The structure is qualitatively different from the `z_3=0` orientifold:

- there is no untwisted constant term,
- the standard Möbius strip is built entirely from twisted sectors,
- the `R`-sector constants cancel between `k=1` and `k=2`,
- the final answer is purely the twisted `NS` contribution.

If needed, one can next repeat the crosscap-state analysis and replace the
standard weights `(1,1,1)` by the orbifold crosscap family

$$
(1,\omega^s,\omega^{2s}),
\qquad
s=0,1,2,
$$

exactly as in the first orientifold.

## 9. The Other Two Crosscap Couplings

Now repeat the same orbifold crosscap-state logic as in
[orbifold_crosscap_state_derivation.md](/Users/yuchenwang/Desktop/D-instanton/Superpotential/orbifold_crosscap_state_derivation.md).
The allowed crosscap-family couplings are

$$
(c_0,c_1,c_2)=(1,\omega^s,\omega^{2s}),
\qquad
s=0,1,2.
$$

For the present orientifold, the untwisted sector vanishes:

$$
\sum_{\alpha,\beta} Z_{M,\alpha\beta}^0=0.
$$

So only the twisted sectors matter.

Define

$$
T(t):=Z_{M,00}^1(t)+Z_{M,01}^1(t),
$$

and write the twisted `R` constants as

$$
Z_{M,11}^1=\varepsilon\,6i\sqrt3,
\qquad
Z_{M,11}^2=-\varepsilon\,6i\sqrt3,
\qquad
\varepsilon=\pm 1,
$$

where `\varepsilon` records the conventional sign choice in the lift of the
orientifold action to the R ground states. Then

$$
\sum_{\alpha,\beta} Z_{M,\alpha\beta}^1
=
T(t)+\varepsilon\,6i\sqrt3,
$$

$$
\sum_{\alpha,\beta} Z_{M,\alpha\beta}^2
=
T(t)-\varepsilon\,6i\sqrt3.
$$

Therefore the three allowed crosscaps are

### 9.1 `s=0`

$$
Z_M^{(0)}(t)
=
\frac{1}{12}\Bigl[(T+\varepsilon\,6i\sqrt3)+(T-\varepsilon\,6i\sqrt3)\Bigr]
=
\frac{1}{6}T(t).
$$

This is exactly the standard Möbius strip derived above.

### 9.2 `s=1`

$$
Z_M^{(1)}(t)
=
\frac{1}{12}\Bigl[\omega(T+\varepsilon\,6i\sqrt3)+\omega^2(T-\varepsilon\,6i\sqrt3)\Bigr].
$$

Using

$$
\omega+\omega^2=-1,
\qquad
\omega-\omega^2=i\sqrt3,
$$

we get

$$
\boxed{
Z_M^{(1)}(t)
=
-\frac{1}{12}T(t)-\frac{3\varepsilon}{2}.
}
$$

### 9.3 `s=2`

Similarly,

$$
\boxed{
Z_M^{(2)}(t)
=
-\frac{1}{12}T(t)+\frac{3\varepsilon}{2}.
}
$$

So the two rotated crosscaps differ only by the sign of the constant shift, as
expected.

## 10. Large-`t` Check Against `K_0`

For an `O(1)` ED3 instanton, the universal instanton zero modes contribute the
same weighted count `3` that appears in the original definition of `K_0` in
[superpotential.tex](/Users/yuchenwang/Desktop/D-instanton/Superpotential/superpotential.tex:152).
Since the annulus of this orientifold contains only stretched `D(-1)`-`D3`
strings, that universal `3` must come from the standard Möbius strip itself.
Therefore

$$
Z_M^{(0)}(t\to\infty)\to 3,
$$

hence

$$
T(t\to\infty)\to 18.
$$

It follows that

$$
Z_M^{(1)}(t\to\infty)\to -\frac{3}{2}-\frac{3\varepsilon}{2},
$$

$$
Z_M^{(2)}(t\to\infty)\to -\frac{3}{2}+\frac{3\varepsilon}{2}.
$$

So, up to the conventional sign `\varepsilon`, the two rotated crosscaps have
large-`t` limits

$$
\{0,-3\}.
$$

Now combine this with the annulus large-`t` limit from
[annulus_other_orientifold_derivation.md](/Users/yuchenwang/Desktop/D-instanton/Superpotential/annulus_other_orientifold_derivation.md):

$$
Z_A(t\to\infty)
=
\frac{3}{2}(N_0-N_1).
$$

For the physical eight-D7 stack in this orientifold, the fractional
multiplicities are

$$
N_0=N_1=N_2=8,
$$

and hence

$$
Z_A(t\to\infty)=0.
$$

For `K_0`, the large-`t` numerator is

$$
Z_A(t)+Z_M^{(s)}(t)+3(e^{-2\pi t}-1)
\to
Z_A(\infty)+Z_M^{(s)}(\infty)-3.
$$

For the standard crosscap, `Z_M^{(0)}(\infty)=3`, and the large-`t`
divergence cancels. The rotated crosscaps have large-`t` limits `{0,-3}` and
therefore do not cancel with the physical regular D7 stack.

Therefore:

$$
\boxed{
\text{for the second orientifold, the physical stack }N_0=N_1=N_2=8
\text{ is compatible with the standard crosscap.}
}
$$
