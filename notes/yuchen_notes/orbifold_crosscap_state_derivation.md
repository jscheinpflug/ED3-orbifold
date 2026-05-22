# Orbifold crosscap-state derivation of the Mobius strip

This note derives the allowed Mobius-strip orbifold weights from actual crosscap
states, using the normalization and modular-crossing conventions of
`stringbook-main-2/texsource/string notes.tex`.

The result is:

1. The orbifold does **not** allow arbitrary coefficients `c_k` in the Mobius strip.
2. Just as for fractional boundary states, the orbifold crosscaps come in a
   `Z_3` family related by the quantum symmetry.
3. Therefore the only allowed relative weights are

$$
(c_0,c_1,c_2)=(1,1,1),\quad (1,\omega,\omega^2),\quad (1,\omega^2,\omega).
$$

4. None of these three choices can simultaneously fix, within the original
   `K_0` subtraction:
   - the small-`t` `1/t` term,
   - the small-`t` constant term,
   - and the large-`t` behavior of `K_0`.

So changing `Z_M` within the actual crosscap-state family does not solve the
problem.

## 1. Crosscap normalization from the string book

The string notes define the type-II crosscap state and its modular crossing in
the same style as boundary states.

For the bosonic part, the crosscap state obeys

$$
(\alpha_n^\mu+(-)^n\widetilde\alpha_{-n}^\mu)|\otimes\rangle=0,
$$

and similarly for the ghosts; see
[string notes.tex](/Users/yuchenwang/Desktop/D-instanton/Superpotential/stringbook-main-2/texsource/string%20notes.tex:10278).

For the superstring, the NSNS and RR crosscap states are given in
[string notes.tex](/Users/yuchenwang/Desktop/D-instanton/Superpotential/stringbook-main-2/texsource/string%20notes.tex:10422),
and the Mobius modular-crossing formula is

$$
\operatorname{Tr}_{{\cal H}_{\rm NS}^o}
\frac{1+(-)^F}{2}\,
\Omega\,
(-)^{N_{bc}+N_{\beta\gamma}}
b_0c_0 e^{-\pi tL_0}

=
\frac{-i}{t}
\langle\!\langle B|_{\rm NSNS}
e^{-\frac{\pi}{2t}L_0^+}
b_0c_0
|\otimes\rangle_{\rm NSNS},
$$

from
[string notes.tex](/Users/yuchenwang/Desktop/D-instanton/Superpotential/stringbook-main-2/texsource/string%20notes.tex:10490).

For an O`p` plane in flat space, the crosscap normalization is fixed relative to
the D`p` boundary state by

$$
{\cal N}_{\otimes p,{\rm NSNS/RR}}
=
\pm 2^{p-4}{\cal N}_{Dp,{\rm NSNS/RR}},
$$

where the sign distinguishes the two orientifold projections; see
[string notes.tex](/Users/yuchenwang/Desktop/D-instanton/Superpotential/stringbook-main-2/texsource/string%20notes.tex:10787).

For our purposes, the important point is that once the orientifold type is fixed,
the overall crosscap normalization is not free.

## 2. Orbifold boundary states in the string book

The string book explicitly constructs fractional boundary states in an orbifold
`M/Z_3`.

If `|B,g^\ell\rangle` denotes the boundary state in the `g^\ell`-twisted closed
sector, then the fractional boundary states are

$$
|\widetilde B_m\rangle
=
\frac{1}{\sqrt3}\sum_{\ell=0}^{2}\omega^{m\ell}|B,g^\ell\rangle,
\qquad
\omega=e^{2\pi i/3},
$$

with `m=0,1,2`; see
[string notes.tex](/Users/yuchenwang/Desktop/D-instanton/Superpotential/stringbook-main-2/texsource/string%20notes.tex:11999).

Their cylinder amplitudes project onto open-string states with definite orbifold
eigenvalue. In particular,

$$
{\cal H}_{\widetilde B_m\widetilde B_{m'}}
\simeq
\left({\cal H}_{BB}\right)_{g=\omega^{m-m'}}.
$$

This is the key template we will now reuse for crosscaps.

## 3. Orbifold crosscaps by the same construction

The book does not spell out the orbifold crosscap state separately, but the same
orbifold logic applies.

Start from a parent-theory crosscap that is invariant under the orbifold action.
Let

$$
|\otimes,g^\ell\rangle
$$

denote the crosscap state in the `g^\ell`-twisted closed sector, obtained by
ending the orbifold defect line on the crosscap, exactly analogously to the
boundary construction.

Then the orbifold-invariant crosscap is

$$
|\widetilde\otimes_0\rangle
=
\frac{1}{\sqrt3}\sum_{\ell=0}^{2}|\otimes,g^\ell\rangle.
$$

Because the orbifold theory has the quantum symmetry `\widetilde Z_3`, which acts
on the `g^\ell`-twisted sector by multiplication with `\omega^\ell`, we obtain a
family of orbifold crosscaps

$$
|\widetilde\otimes_m\rangle
=
\frac{1}{\sqrt3}\sum_{\ell=0}^{2}\omega^{m\ell}|\otimes,g^\ell\rangle,
\qquad m=0,1,2.
$$

This is the exact crosscap analogue of the book's fractional boundary-state
construction.

## 4. Mobius strip from fractional brane and orbifold crosscap

Take the fractional ED3 instanton to be type `n`, so

$$
|\widetilde B_n\rangle
=
\frac{1}{\sqrt3}\sum_{r=0}^{2}\omega^{nr}|B,g^r\rangle.
$$

Take the crosscap to be `|\widetilde\otimes_m\rangle`.

Using the same modular-crossing logic as in the book, only equal twisted sectors
pair in the closed channel, so the Mobius overlap becomes

$$
\langle\!\langle \widetilde B_n|
e^{-\pi \ell H_c}
|\widetilde\otimes_m\rangle
=
\frac13\sum_{k=0}^{2}\omega^{(m-n)k}
\langle\!\langle B,g^k|
e^{-\pi \ell H_c}
|\otimes,g^k\rangle.
$$

By the Mobius modular-crossing equation, this is equivalent in the open channel to

$$
Z_M^{(m-n)}
=
\frac1{12}\sum_{k=0}^{2}\omega^{(m-n)k}
\sum_{\alpha,\beta} Z_{M,\alpha\beta}^k.
$$

So the allowed orbifold weights are not arbitrary:

$$
c_k=\omega^{(m-n)k}.
$$

Since shifting both `m` and `n` by the same amount only changes conventions, we
may take the ED3 to be type `n=0`. Then the three allowed crosscap choices are

$$
(c_0,c_1,c_2)=(1,1,1),\quad (1,\omega,\omega^2),\quad (1,\omega^2,\omega).
$$

This is the crosscap-state derivation of the Mobius weights.

## 5. Explicit Mobius-strip formulas for the three crosscaps

From
[mobius_ZMabk_derivation.md](/Users/yuchenwang/Desktop/D-instanton/Superpotential/mobius_ZMabk_derivation.md:864),
define

$$
T(t):=Z_{M,00}^1(t)+Z_{M,01}^1(t),
$$

so that

$$
\sum_{\alpha,\beta} Z_{M,\alpha\beta}^0 = 16,
\qquad
\sum_{\alpha,\beta} Z_{M,\alpha\beta}^1 = T(t)+2i\sqrt3,
\qquad
\sum_{\alpha,\beta} Z_{M,\alpha\beta}^2 = T(t)-2i\sqrt3.
$$

Then:

### 5.1 Untwisted crosscap `m=0`

$$
Z_M^{(0)}(t)
=
\frac1{12}\Bigl(16 + (T+2i\sqrt3)+(T-2i\sqrt3)\Bigr)
=
\frac43+\frac16 T(t).
$$

This is the current formula in the notes.

### 5.2 Quantum-symmetry rotated crosscap `m=1`

$$
Z_M^{(1)}(t)
=
\frac1{12}\Bigl(16 + \omega(T+2i\sqrt3)+\omega^2(T-2i\sqrt3)\Bigr).
$$

Using

$$
\omega+\omega^2=-1,
\qquad
\omega-\omega^2=i\sqrt3,
$$

we get

$$
Z_M^{(1)}(t)=\frac56-\frac1{12}T(t).
$$

### 5.3 Quantum-symmetry rotated crosscap `m=2`

Similarly,

$$
Z_M^{(2)}(t)
=
\frac1{12}\Bigl(16 + \omega^2(T+2i\sqrt3)+\omega(T-2i\sqrt3)\Bigr)
=
\frac{11}{6}-\frac1{12}T(t).
$$

## 6. Their asymptotics

From the current result,

$$
Z_M^{(0)}(t/4)\sim \frac{2}{\sqrt3}\frac1t+\frac43,
\qquad
Z_M^{(0)}(t\to\infty)\to 3,
$$

so necessarily

$$
T(t/4)\sim \frac{4\sqrt3}{t},
\qquad
T(t\to\infty)\to 10.
$$

Therefore:

### 6.1 `m=0`

$$
Z_M^{(0)}(t/4)\sim \frac{2}{\sqrt3}\frac1t+\frac43,
\qquad
Z_M^{(0)}(t\to\infty)\to 3.
$$

### 6.2 `m=1`

$$
Z_M^{(1)}(t/4)\sim -\frac{1}{\sqrt3}\frac1t+\frac56,
\qquad
Z_M^{(1)}(t\to\infty)\to 0.
$$

### 6.3 `m=2`

$$
Z_M^{(2)}(t/4)\sim -\frac{1}{\sqrt3}\frac1t+\frac{11}{6},
\qquad
Z_M^{(2)}(t\to\infty)\to 1.
$$

So the three actual crosscap states give three and only three Mobius asymptotics.

## 7. Combine with the most general annulus boundary state

From
[annulus_boundary_state_redo.md](/Users/yuchenwang/Desktop/D-instanton/Superpotential/annulus_boundary_state_redo.md:199),
the most general D7 fractional stack

$$
N_0\widetilde D_0+N_1\widetilde D_1+N_2\widetilde D_2,
\qquad
N_0+N_1+N_2=8,
$$

has annulus asymptotics

$$
Z_A(t)\sim
\frac{2N_0-N_1-N_2}{4\sqrt3}\frac1t
+
\frac{-2N_0+N_1-5N_2}{12},
$$

and

$$
Z_A(t\to\infty)=-\frac{N_2}{2}.
$$

We now test each actual crosscap.

### 7.1 Crosscap `m=0`

To cancel the small-`t` `1/t` term:

$$
\frac{2N_0-N_1-N_2}{4\sqrt3}+\frac{2}{\sqrt3}=0
\quad\Rightarrow\quad
N_0=0.
$$

Then `N_1+N_2=8`.

To cancel the small-`t` constant:

$$
\frac{-2N_0+N_1-5N_2}{12}+\frac43=0
\quad\Rightarrow\quad
N_1=N_2=4.
$$

Then the large-`t` numerator of the `K_0` integrand is

$$
Z_A+Z_M-3 \to -2+3-3=-2\neq 0.
$$

So `m=0` fails.

### 7.2 Crosscap `m=1`

To cancel the small-`t` `1/t` term:

$$
\frac{2N_0-N_1-N_2}{4\sqrt3}-\frac{1}{\sqrt3}=0
\quad\Rightarrow\quad
N_0=4.
$$

Then `N_1+N_2=4`.

To cancel the small-`t` constant:

$$
\frac{-2N_0+N_1-5N_2}{12}+\frac56=0
\quad\Rightarrow\quad
N_2=1,\quad N_1=3.
$$

Then the large-`t` numerator is

$$
Z_A+Z_M-3 \to -\frac12+0-3=-\frac72\neq 0.
$$

So `m=1` fails.

### 7.3 Crosscap `m=2`

Again the small-`t` `1/t` cancellation gives

$$
N_0=4,
\qquad
N_1+N_2=4.
$$

The small-`t` constant cancellation gives

$$
\frac{-2N_0+N_1-5N_2}{12}+\frac{11}{6}=0
\quad\Rightarrow\quad
N_2=3,\quad N_1=1.
$$

Then the large-`t` numerator is

$$
Z_A+Z_M-3 \to -\frac32+1-3=-\frac72\neq 0.
$$

So `m=2` also fails.

## 8. Conclusion

Using the actual orbifold crosscap states, there are only three allowed Mobius
choices, related by the quantum symmetry of the orbifold:

$$
Z_M^{(0)},\quad Z_M^{(1)},\quad Z_M^{(2)}.
$$

None of them can be paired with any D7 fractional stack

$$
N_0\widetilde D_0+N_1\widetilde D_1+N_2\widetilde D_2
$$

so as to satisfy all three requirements:

1. cancel the small-`t` `1/t` term,
2. cancel the small-`t` constant term,
3. avoid a large-`t` divergence with the original subtraction of only the
   universal ED3 zero modes.

So the conclusion survives the crosscap-state derivation:

Changing `Z_M` within the actual orbifold crosscap family does not fix the
problem if one keeps the original `K_0` formula unchanged.

For the physically natural orientifold-invariant annulus choice

$$
(N_0,N_1,N_2)=(0,4,4),
$$

the small-`t` closed-channel massless exchange is canceled, but the large-`t`
tail signals extra `3-7` Ramond zero modes from the type-2 D7 branes.  In that
situation, the correct conclusion is not that the annulus or Mobius amplitudes
are wrong; it is that the large-`t` subtraction in `K_0` has to be generalized
to include the full zero-mode sector.
