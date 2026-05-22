# Redo of the ED3-D7 annulus boundary-state calculation

This note redoes the annulus boundary-state bookkeeping using the orbifold
boundary-state construction from `stringbook-main-2/texsource/string notes.tex`.
The goal is not to recompute every theta-function identity, but to isolate the
precise assumption behind the current ED3-D7 annulus boundary state and see
where it clashes with the rest of the orientifold setup.

The punchline is:

- The annulus coefficient should be described in terms of a direct sum of
  fractional D7 boundary states, not just an ``SO(8) stack times one character''.
- The choice `8 omega^{-k}` corresponds to putting all 8 D7-branes in a single
  fractional sector.
- That choice cancels the small-`t` `1/t` term against the current Mobius strip,
  but it is not orientifold-invariant and it leaves the small-`t` constant
  uncanceled.
- The orientifold-invariant split `(N_0,N_1,N_2)=(0,4,4)` cancels both small-`t`
  coefficients, but then it necessarily leaves surviving `3-7` Ramond zero
  modes at large `t`.

So the real tension is not that the annulus and Mobius amplitudes are
algebraically inconsistent.  It is that the orientifold-invariant UV
cancellation and the original large-`t` subtraction in `K_0` cannot both hold
without revisiting the `3-7` zero-mode sector.

## 1. Orbifold boundary states from the string notes

Section `D-branes in orbifolds` in the string notes constructs fractional
boundary states for a `Z_3` orbifold as

$$
|\widetilde B_m\rangle
=
\frac{1}{\sqrt 3}\sum_{\ell=0}^{2}\omega^{m\ell}|B,g^\ell\rangle,
\qquad
\omega=e^{2\pi i/3},
$$

with `m=0,1,2`; see
[string notes.tex](/Users/yuchenwang/Desktop/D-instanton/Superpotential/stringbook-main-2/texsource/string%20notes.tex:11999).

The strip Hilbert space between two fractional branes satisfies

$$
{\cal H}_{\widetilde B_m\widetilde B_{m'}}
\simeq
\left({\cal H}_{BB}\right)_{g=\omega^{m-m'}}.
$$

Equivalently, the annulus projection between fractional branes `m` and `m'` is

$$
\frac13\sum_{k=0}^2 \omega^{(m'-m)k} g^k.
$$

So the orbifold phase in the annulus is not an arbitrary decoration: it tells us
which fractional D7 type the ED3 is talking to.

## 2. General D7 stack at the orbifold point

Take the ED3 instanton to be fractional type `m=0`; this is just a convention.
Let the D7 configuration be a direct sum

$$
|D7_{\rm stack}\rangle
=
N_0|\widetilde D_0\rangle
+N_1|\widetilde D_1\rangle
+N_2|\widetilde D_2\rangle,
\qquad
N_0+N_1+N_2=8.
$$

Then the ED3-D7 annulus takes the form

$$
Z_A[N_0,N_1,N_2]
=
\frac16\sum_{k=0}^2
\Bigl(N_0+N_1\omega^{-k}+N_2\omega^k\Bigr)
\sum_{\alpha,\beta} Z_{\alpha\beta}^k.
$$

So the three basic pure choices are:

- `N=(8,0,0)`: trivial character `8`
- `N=(0,8,0)`: character `8 omega^{-k}`
- `N=(0,0,8)`: character `8 omega^{k}`

The current note
[annulus_normalization_from_boundary_states.md](/Users/yuchenwang/Desktop/D-instanton/Superpotential/annulus_normalization_from_boundary_states.md:386)
implicitly picked the second of these.

## 3. What the three pure choices do in the open channel

The open-channel interpretation is immediate from the orbifold projector.

### 3.1 `N=(8,0,0)` projects to `g=1`

This is the trivial character. It keeps the `g=1` eigenspace of the parent `3-7`
Hilbert space.

In the `3-7` system, the NS ground states are massive (`L_0=1/2`), while the R
ground states are massless and carry orbifold eigenvalues

$$
\omega^{-2s_3}=\omega,\omega^2
$$

from
[annulus_Zabk_derivation.md](/Users/yuchenwang/Desktop/D-instanton/Superpotential/annulus_Zabk_derivation.md:153).
So the trivial character keeps no massless R ground state. Therefore

$$
Z_A^{(0)}(t\to\infty)\to 0.
$$

This is the large-`t` behavior of the ``all type-0'' stack.

### 3.2 `N=(0,8,0)` projects to `g=\omega`

This is the current `8 omega^{-k}` choice.

It selects the R ground state with orbifold eigenvalue `omega`, namely the one
with `s_3=-1/2`. But this state is GSO-odd because

$$
(-1)^F = 2s_3 = -1.
$$

So again no `3-7` zero mode survives, and

$$
Z_A^{(1)}(t\to\infty)\to 0.
$$

This is why the current annulus formula still has acceptable large-`t` behavior
despite using a nontrivial orbifold character.

### 3.3 `N=(0,0,8)` projects to `g=\omega^2`

This is the `8 omega^k` choice.

Now the surviving R ground state is the one with `s_3=+1/2`, which is GSO-even:

$$
(-1)^F = 2s_3 = +1.
$$

So this choice does produce `3-7` fermion zero modes. In the annulus trace this
shows up as the weighted zero-mode coefficient

$$
Z_A^{(2)}(t\to\infty)\to -4,
$$

corresponding to eight fermionic zero modes, each weighted by `-1/2`.

## 4. Small-`t` data of the three pure choices

The three pure choices are enough to determine the asymptotics of any direct sum,
because the annulus is linear in the boundary state.

From
[annulus_normalization_from_boundary_states.md](/Users/yuchenwang/Desktop/D-instanton/Superpotential/annulus_normalization_from_boundary_states.md:289)
and
[annulus_normalization_from_boundary_states.md](/Users/yuchenwang/Desktop/D-instanton/Superpotential/annulus_normalization_from_boundary_states.md:345),
we already have

$$
Z_A^{(0)}(t)
\sim
\frac{4}{\sqrt3}\frac1t-\frac43,
$$

$$
Z_A^{(1)}(t)
\sim
-\frac{2}{\sqrt3}\frac1t+\frac23.
$$

Switching `omega^{-k}` to `omega^k` changes only the `Z_{11}^k` piece and shifts
the annulus by `-4`, so

$$
Z_A^{(2)}(t)=Z_A^{(1)}(t)-4
\sim
-\frac{2}{\sqrt3}\frac1t-\frac{10}{3}.
$$

Therefore a general direct sum satisfies

$$
Z_A[N_0,N_1,N_2]
=
\frac{N_0}{8}Z_A^{(0)}
+\frac{N_1}{8}Z_A^{(1)}
+\frac{N_2}{8}Z_A^{(2)}.
$$

So the small-`t` asymptotics are

$$
Z_A[N_0,N_1,N_2](t)
\sim
\frac{2N_0-N_1-N_2}{4\sqrt3}\frac1t
+
\frac{-2N_0+N_1-5N_2}{12}.
$$

And the large-`t` limit is

$$
Z_A[N_0,N_1,N_2](t\to\infty)
=
-\frac{N_2}{2}.
$$

This single formula already explains all three previously discussed special
cases.

## 5. What the current Mobius strip asks for

The current Mobius strip asymptotics are

$$
Z_M(t/4)\sim \frac{2}{\sqrt3}\frac1t+\frac43,
\qquad
Z_M(t\to\infty)\to 3,
$$

from
[annulus_normalization_from_boundary_states.md](/Users/yuchenwang/Desktop/D-instanton/Superpotential/annulus_normalization_from_boundary_states.md:297)
and
[superpotential.tex](/Users/yuchenwang/Desktop/D-instanton/Superpotential/superpotential.tex:270).

So:

- To cancel the small-`t` `1/t` pole of `Z_M(t/4)`, the annulus must satisfy

$$
\frac{2N_0-N_1-N_2}{4\sqrt3}=-\frac{2}{\sqrt3}.
$$

Using `N_0+N_1+N_2=8`, this becomes

$$
N_0=0,\qquad N_1+N_2=8.
$$

So the current crosscap wants a D7 stack made only from type-1 and type-2
fractional branes.

- To also cancel the small-`t` constant term `+4/3` of `Z_M(t/4)`, the annulus
must satisfy

$$
\frac{-2N_0+N_1-5N_2}{12}=-\frac43.
$$

Together with `N_0=0` and `N_1+N_2=8`, this gives

$$
N_1=N_2=4.
$$

So the unique direct sum that cancels both the `1/t` and constant terms is

$$
(N_0,N_1,N_2)=(0,4,4).
$$

But then the large-`t` limit is

$$
Z_A(t\to\infty)=-\frac{N_2}{2}=-2,
$$

which means four `3-7` fermion zero modes are present.  With the original definition of
`K_0`, which subtracts only the universal ED3 zero modes, the large-`t`
integrand is then IR divergent.

## 6. Orientifold invariance

There is an independent consistency condition. Worldsheet parity conjugates the
orbifold representation, so the nontrivial `Z_3` characters are exchanged:

$$
\widetilde D_1 \longleftrightarrow \widetilde D_2.
$$

Therefore an orientifold-invariant D7 stack must satisfy

$$
N_1=N_2.
$$

This is exactly what happened above: the only candidate that cancels both small-`t`
singular pieces is the orientifold-invariant split

$$
(0,4,4).
$$

But that same orientifold-invariant split necessarily contains type-2 branes and
therefore necessarily gives `3-7` zero modes.

By contrast, the current `8 omega^{-k}` choice corresponds to

$$
(N_0,N_1,N_2)=(0,8,0),
$$

which avoids `3-7` zero modes and cancels the small-`t` `1/t` pole, but it is
not orientifold-invariant and leaves the small-`t` constant term uncanceled.

## 7. What is actually going wrong

The boundary-state problem is now visible quite sharply:

1. The current Mobius strip is fixing the crosscap couplings of the O7-plane.
2. The current annulus with `8 omega^{-k}` is solving the twisted tadpole problem
   for a non-orientifold-invariant D7 boundary state `(0,8,0)`.
3. The orientifold-invariant D7 boundary state that would cancel both the `1/t`
   and constant pieces at small `t` is `(0,4,4)`.
4. But `(0,4,4)` inevitably contains the type-2 fractional D7s that generate
   `3-7` zero modes, so it fails at large `t`.

So the annulus and Mobius amplitudes are not in conflict because of a small
algebra mistake in the phase `omega^{-k}` versus `omega^k`.  The real issue is
that the current crosscap data, the orientifold-invariant D7 stack, and the
original assumption of no extra `3-7` zero modes cannot all hold at once.

Stated differently:

- `8` gives the right large-`t` behavior but the wrong small-`t` tadpole.
- `8 omega^{-k}` gives the right small-`t` `1/t` coefficient and the right
  large-`t` behavior, but the wrong small-`t` constant and it is not
  orientifold-invariant.
- `(0,4,4)` gives the right small-`t` behavior, but the wrong large-`t`
  behavior.

This strongly suggests that the part that has to change is not just the annulus
phase, but one of the following:

- the orientifold action on the fractional D7 boundary states,
- the crosscap/Mobius-strip couplings in the twisted sectors,
- the claim that the chosen O7/D7 configuration has no `3-7` zero modes,
- or the formula for `K_0`, which must be generalized once the `3-7` zero-mode
  sector is present.
