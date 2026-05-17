(* K_0 computation - Version 3 *)
(* Systematic search for M_coeff that gives convergent integral *)
(* And stability analysis *)

SetOptions[NIntegrate, WorkingPrecision -> 30, PrecisionGoal -> 12, MaxRecursion -> 100];

(* Partition functions *)
ZAnnulus[t_?NumericQ, Nbulk_] := Module[{tau2, q, th10, th00},
  tau2 = 2 I t;
  q = Exp[-2 Pi t];
  th10 = EllipticTheta[2, 0, Exp[I Pi tau2]];
  th00 = EllipticTheta[3, 0, Exp[I Pi tau2]];
  Re[Nbulk/2 * (th10 * q^(-1/4) + th00)]
]

ZMobius[t_?NumericQ, Mcoeff_] := Module[{hattau2, hatq, th10hat, th00hat},
  hattau2 = 2 (I t + 1/2);
  hatq = Exp[-2 Pi t + I Pi];
  th10hat = EllipticTheta[2, 0, Exp[I Pi hattau2]];
  th00hat = EllipticTheta[3, 0, Exp[I Pi hattau2]];
  Mcoeff * Re[th10hat * hatq^(-1/4) - th00hat]
]

Zsub[t_?NumericQ] := 3 (Exp[-2 Pi t] - 1)

Ztotal[t_?NumericQ, Nbulk_, Mcoeff_] := ZAnnulus[t, Nbulk] + ZMobius[t, Mcoeff] + Zsub[t]

Nbulk = 8;

(* Find the M_coeff where integral converges *)
Print["=== Finding M_coeff for convergent integral ===\n"];

(* For convergence, the integrand must decay faster than 1/t at small t *)
(* Check: integrand(t) * t -> 0 as t -> 0 *)

checkConvergence[Mcoeff_] := Module[{Zasymp, intgnd},
  Zasymp = Ztotal[50, Nbulk, Mcoeff];
  intgnd[t_] := (Ztotal[t, Nbulk, Mcoeff] - Zasymp) / (2 t);
  (* Check integrand * t at small t *)
  Table[{t, intgnd[t] * t}, {t, {0.001, 0.01, 0.1}}]
]

Print["Checking integrand * t at small t (should -> 0 for convergence):"];
Do[
  data = checkConvergence[M];
  Print["M_coeff = ", M, ": ", Map[{#[[1]], N[#[[2]], 4]} &, data]],
  {M, {-10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10}}
];

(* More refined search *)
Print["\n=== Refined search ==="];

findOptimalM[] := Module[{Zasymp, intgnd, objective},
  objective[Mcoeff_?NumericQ] := Module[{Za},
    Zasymp = Ztotal[50, Nbulk, Mcoeff];
    intgnd[t_] := (Ztotal[t, Nbulk, Mcoeff] - Zasymp) / (2 t);
    (* Minimize |integrand * t| at small t *)
    Abs[intgnd[0.001] * 0.001] + Abs[intgnd[0.01] * 0.01]
  ];
  FindMinimum[objective[M], {M, 0}]
]

result = findOptimalM[];
Moptimal = M /. result[[2]];
Print["Optimal M_coeff: ", N[Moptimal, 10]];
Print["Objective value: ", N[result[[1]], 10]];

(* Verify convergence with optimal M *)
Print["\n=== Verification with optimal M_coeff ==="];
Zasymp = Ztotal[50, Nbulk, Moptimal];
Print["Z_asymptotic = ", N[Zasymp, 10]];

integrand[t_?NumericQ] := (Ztotal[t, Nbulk, Moptimal] - Zasymp) / (2 t)

Print["\nIntegrand * t at small t:"];
Do[Print["t = ", t, ": integrand*t = ", N[integrand[t] * t, 10]], {t, {0.0001, 0.001, 0.01, 0.1}}];

Print["\nIntegrand * sqrt(t) at small t:"];
Do[Print["t = ", t, ": integrand*sqrt(t) = ", N[integrand[t] * Sqrt[t], 10]], {t, {0.0001, 0.001, 0.01, 0.1}}];

(* Numerical integration *)
Print["\n=== Numerical integration ==="];
Print["If integrand ~ c/sqrt(t), integral ~ 2c*sqrt(t_max) - 2c*sqrt(t_min), diverges"];
Print["If integrand ~ constant, integral converges"];

results = Table[
  {tmin, NIntegrate[integrand[t], {t, tmin, 50}]},
  {tmin, {0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001, 0.0005}}
];

Print["\nIntegral values:"];
Do[Print["t_min = ", r[[1]], ": I = ", N[r[[2]], 12]], {r, results}];

(* Check differences *)
Print["\nDifferences (delta_I / delta_sqrt(t)):"];
Do[
  t1 = results[[i, 1]]; t2 = results[[i+1, 1]];
  I1 = results[[i, 2]]; I2 = results[[i+1, 2]];
  deltaI = I2 - I1;
  deltaSqrtT = Sqrt[t1] - Sqrt[t2];
  Print["  ", t1, " -> ", t2, ": deltaI = ", N[deltaI, 6],
        ", deltaI/delta_sqrt_t = ", N[deltaI/deltaSqrtT, 6]],
  {i, Length[results] - 1}
];

(* If deltaI/delta_sqrt_t is constant, we have 1/sqrt(t) divergence *)
(* If deltaI is decreasing, we're converging *)

Print["\n=== Stability against mode truncation ==="];
(* The theta functions are computed with machine precision *)
(* Check that results don't change with higher precision *)

intHighPrec = NIntegrate[integrand[t], {t, 0.001, 50},
  WorkingPrecision -> 50, PrecisionGoal -> 20];
intLowPrec = NIntegrate[integrand[t], {t, 0.001, 50},
  WorkingPrecision -> 20, PrecisionGoal -> 10];

Print["High precision (WP=50): ", N[intHighPrec, 15]];
Print["Low precision (WP=20):  ", N[intLowPrec, 15]];
Print["Difference: ", N[Abs[intHighPrec - intLowPrec], 10]];

(* Final extrapolation *)
Print["\n=== Final result ==="];

(* Use Sequence acceleration on the integral values *)
seqData = Map[#[[2]] &, results];
accelerated = SequenceLimit[seqData];

Print["Raw integral at t_min = 0.0005: ", N[results[[-1, 2]], 12]];
Print["Sequence-accelerated limit: ", N[accelerated, 12]];

Cfinal = accelerated;
K0 = 1/(2 Pi)^(3/2) * Exp[Cfinal];

Print["\n========================================"];
Print["FINAL RESULT"];
Print["========================================"];
Print["M_coeff = ", N[Moptimal, 6]];
Print["C (finite part) = ", N[Cfinal, 10]];
Print["K_0 = ", N[K0, 10]];
Print["log(K_0) = ", N[Log[K0], 10]];
Print["========================================"];
