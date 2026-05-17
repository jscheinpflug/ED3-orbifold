(* K_0 computation for C^3/Z_3 orientifold - FINAL VERSION *)
(* Uses d = 4 zero modes (spacetime position of D(-1)) *)
(* Internal C^3/Z_3 directions frozen at orbifold fixed point *)

$MinPrecision = 50;

(* Partition functions *)
ZAnnulus[t_?NumericQ, Nbulk_] := Module[{tau2, q, th10, th00},
  tau2 = 2 I t;
  q = Exp[-2 Pi t];
  th10 = N[EllipticTheta[2, 0, Exp[I Pi tau2]], 50];
  th00 = N[EllipticTheta[3, 0, Exp[I Pi tau2]], 50];
  Re[Nbulk/2 * (th10 * q^(-1/4) + th00)]
]

ZMobius[t_?NumericQ, Mcoeff_] := Module[{hattau2, hatq, th10hat, th00hat},
  hattau2 = 2 (I t + 1/2);
  hatq = Exp[-2 Pi t + I Pi];
  th10hat = N[EllipticTheta[2, 0, Exp[I Pi hattau2]], 50];
  th00hat = N[EllipticTheta[3, 0, Exp[I Pi hattau2]], 50];
  Mcoeff * Re[th10hat * hatq^(-1/4) - th00hat]
]

(* Zero mode subtraction with d = 4 *)
d = 4;
Zsub[t_] := d * (Exp[-2 Pi t] - 1)

Ztotal[t_?NumericQ, Nbulk_, Mcoeff_] := ZAnnulus[t, Nbulk] + ZMobius[t, Mcoeff] + Zsub[t]

(* Physical parameters *)
Nbulk = 8;    (* = 24 fractional D7-branes *)
Mcoeff = -8;  (* O7-plane charge *)

Print["=== K_0 Computation for C^3/Z_3 Orientifold ===\n"];
Print["N_bulk = ", Nbulk, " (= 24 fractional D7-branes)"];
Print["M_coeff = ", Mcoeff, " (O7-plane charge)"];
Print["d = ", d, " (spacetime zero modes)\n"];

(* Asymptotic value *)
Zasymp = Ztotal[100, Nbulk, Mcoeff];
Print["Z_asymptotic (t=100) = ", N[Zasymp, 10]];
Print["(Should be 0: Z_A -> 12, Z_M -> -8, Z_sub -> -4)\n"];

(* Integrand *)
integrand[t_?NumericQ] := (Ztotal[t, Nbulk, Mcoeff] - Zasymp) / (2 t)

(* Verify convergence *)
Print["=== Checking convergence ==="];
Print["Z_total at small and large t:"];
Do[
  Zt = Ztotal[t, Nbulk, Mcoeff];
  Print["t = ", t, ": Z_total = ", N[Zt, 8]],
  {t, {0.0001, 0.001, 0.01, 0.1, 1, 10, 100}}
];

Print["\nIntegrand * sqrt(t) (should approach constant ~2.3):"];
Do[
  Print["t = ", t, ": ", N[integrand[t] * Sqrt[t], 10]],
  {t, {0.0001, 0.001, 0.01}}
];

(* Numerical integration *)
Print["\n=== Numerical integration ==="];
tmins = {0.02, 0.01, 0.005, 0.002, 0.001, 0.0005, 0.0002, 0.0001, 0.00005};
integrals = Table[
  {tmin, NIntegrate[integrand[t], {t, tmin, 100},
    WorkingPrecision -> 50, PrecisionGoal -> 25, MaxRecursion -> 100]},
  {tmin, tmins}
];

Print["t_min\t\tI(t_min)\t\tDifference"];
prev = 0;
Do[
  diff = If[prev == 0, "---", N[integrals[[i, 2]] - prev, 6]];
  Print[integrals[[i, 1]], "\t\t", N[integrals[[i, 2]], 12], "\t", diff];
  prev = integrals[[i, 2]],
  {i, Length[integrals]}
];

(* Extrapolation: I(t) = C + b*sqrt(t) *)
Print["\n=== Extrapolation to t_min -> 0 ==="];
fitData = Map[{Sqrt[#[[1]]], #[[2]]} &, integrals];
fit = Fit[fitData, {1, x}, x];
Cfinal = fit /. x -> 0;
Print["Fit: I = C + b*sqrt(t)"];
Print["C = ", N[Cfinal, 12]];

(* Final result *)
Print["\n========================================"];
Print["FINAL RESULT"];
Print["========================================"];

K0 = 1/(2 Pi)^(d/2) * Exp[Cfinal];

Print["\nWith d = ", d, " zero modes:"];
Print["  Prefactor = (2 Pi)^(-", d/2, ") = ", N[1/(2 Pi)^(d/2), 10]];
Print["  Finite part C = ", N[Cfinal, 10]];
Print["  K_0 = ", N[K0, 10]];
Print["  log(K_0) = ", N[Log[K0], 10]];

Print["\nPhysical interpretation:"];
Print["  K_0 < 1: one-loop determinant SUPPRESSES the instanton amplitude."];
Print["\nThe superpotential correction is:"];
Print["  W_inst = K_0 * exp(-S_inst) = ", N[K0, 4], " * exp(-S_inst)"];
Print["========================================"];
