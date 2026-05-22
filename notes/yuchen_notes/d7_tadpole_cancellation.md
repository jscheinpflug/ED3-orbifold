# D7 Tadpole Normalization In The `C^3/Z_3` Orientifold

This note fixes the normalization of the D7 Chan-Paton multiplicities
`N_0,N_1,N_2`.

Xi's string notes give the flat-space orientifold normalization
$$
Q_{O_p^-}=-2^{p-5}
$$
in units of a bulk D`p`-brane pair. Equivalently, a D`p`-brane sitting on the
orientifold plane carries half the charge of such a bulk pair. For `p=7`, this
means
$$
Q_{O7^-}=-4\quad\text{bulk-pair units}
        =-8\quad\text{fixed-plane D7 units}.
$$
This is the standard local statement that an `O7^-` plane is cancelled by eight
D7-branes on top of it, with `SO(8)` Chan-Paton factor.

For the first orientifold, the D7-branes are the ordinary D7-branes sitting on
the local O7-plane. Including the `Z_3` orbifold, a D7 stack with equivariant
Chan-Paton action
$$
\gamma_g=\mathrm{diag}
\left(1_{N_0},\omega\,1_{N_1},\omega^2\,1_{N_2}\right),
$$
the orbifold boundary state has the schematic form
$$
|B_{D7}\rangle_{\rm orb}
=\frac{1}{\sqrt3}\sum_{k=0}^2
\mathrm{Tr}\,\gamma_g^k\,|B_{D7};k\rangle .
$$
The crosscap state has the analogous orbifold normalization
$$
|\otimes_{O7}\rangle_{\rm orb}
=\frac{1}{\sqrt3}\sum_{k=0}^2 c_k\,|\otimes_{O7};k\rangle .
$$
Therefore the untwisted D7 tadpole is proportional to
$$
\frac{1}{\sqrt3}\left[
\mathrm{Tr}\,\gamma_g^0-8
\right]
=
\frac{1}{\sqrt3}\left[
N_0+N_1+N_2-8
\right].
$$
The `1/\sqrt3` is common to the boundary and crosscap states, so it does not
turn the condition into `24`. The untwisted D7 tadpole condition is
$$
\boxed{N_0+N_1+N_2=8.}
$$

For the second orientifold, the D7-branes wrap the exceptional divisor. At the
orbifold point, one physical D7-brane of this type is represented by the regular
representation
$$
\widetilde D_0+\widetilde D_1+\widetilde D_2.
$$
Thus the eight physical D7-branes required by the O7-plane correspond to
$$
\boxed{N_0=N_1=N_2=8}
$$
in fractional-brane multiplicities. In this case `N_0+N_1+N_2=24` counts
fractional constituents, not physical D7-branes.

For the first orientifold, the twisted small-`t` cancellation against the
crosscap gives
$$
\frac{2N_0-N_1-N_2}{4\sqrt3}+\frac{2}{\sqrt3}=0,
$$
and
$$
\frac{-2N_0+N_1-5N_2}{12}+\frac43=0.
$$
Together with the untwisted condition `N_0+N_1+N_2=8`, these imply
$$
N_0=0,\qquad N_1=N_2=4.
$$

For the second orientifold, the annulus tadpole is proportional to
$$
2N_0-N_1-N_2
$$
and the constant part is proportional to `N_1-N_2`; both vanish for
`N_0=N_1=N_2=8`.
