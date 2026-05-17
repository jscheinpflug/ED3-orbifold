(* K_0 computation - Version 2 *)
(* With M_coeff = -8 (O-plane charge) and proper convergence analysis *)

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

(* Use N_bulk = 8 and M_coeff = -8 for tadpole cancellation *)
Nbulk = 8;
Mcoeff = -8;

Print["=== Using N_bulk = ", Nbulk, ", M_coeff = ", Mcoeff, " ===\n"];

(* Compute asymptotic value *)
Zasymp = Ztotal[50, Nbulk, Mcoeff];
Print["Z_asymptotic (t=50) = ", N[Zasymp, 15]];

(* Integrand *)
integrand[t_?NumericQ] := (Ztotal[t, Nbulk, Mcoeff] - Zasymp) / (2 t)

(* Analyze small-t behavior *)
Print["\n=== Small-t analysis ==="];
Print["Checking integrand structure:"];
smallTdata = Table[{t, integrand[t], integrand[t] * t, integrand[t] * Sqrt[t]},
  {t, {0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1}}];

Print["t\t\tintegrand\t\tintegrand*t\t\tintegrand*sqrt(t)"];
Do[Print[r[[1]], "\t\t", N[r[[2]], 6], "\t\t", N[r[[3]], 6], "\t\t", N[r[[4]], 6]],
   {r, smallTdata}];

(* If integrand ~ A/t + B/sqrt(t) + C + ..., then:
   integrand*t ~ A + B*sqrt(t) + C*t + ...
   integrand*sqrt(t) ~ A/sqrt(t) + B + C*sqrt(t) + ... *)

(* Extract leading behavior *)
Print["\n=== Extracting divergence structure ==="];

(* Fit: integrand = a/t + b/sqrt(t) + c *)
fitData = Table[{t, integrand[t]}, {t, 0.001, 0.1, 0.001}];

(* Use linear regression on transformed data *)
(* Let y = integrand, x = 1/t, z = 1/sqrt(t) *)
(* y = a*x + b*z + c *)

Print["Fitting integrand = a/t + b/sqrt(t) + c"];
fitResult = Fit[Table[{1/t, 1/Sqrt[t], 1, integrand[t]}, {t, 0.005, 0.1, 0.005}],
                {x, y, 1}, {x, y}];
Print["Best fit: ", fitResult /. {x -> "1/t", y -> "1/sqrt(t)"}];

(* Alternative: assume just log divergence (a=0) *)
Print["\nAssuming integrand = b/sqrt(t) + c + d*sqrt(t):"];
fitResult2 = FindFit[Table[{t, integrand[t]}, {t, 0.001, 0.1, 0.002}],
                      b/Sqrt[t] + c + d*Sqrt[t], {b, c, d}, t];
Print["Best fit: b = ", b /. fitResult2, ", c = ", c /. fitResult2, ", d = ", d /. fitResult2];

(* Check if the 1/sqrt(t) term is really cancelled *)
Print["\nChecking 1/sqrt(t) coefficient:"];
bcoeff = (b /. fitResult2);
Print["b (1/sqrt(t) coefficient) = ", N[bcoeff, 10]];

If[Abs[bcoeff] < 0.1,
  Print["Small b -> 1/sqrt(t) divergence is cancelled!"];
  cvalue = (c /. fitResult2);
  Print["Finite part c = ", N[cvalue, 10]],
  Print["b not small -> divergence remains"]
];

(* Direct numerical integration with increasing precision *)
Print["\n=== Direct numerical integration ==="];
Print["Integrating from t_min to t_max = 100:"];

results = {};
Do[
  integral = NIntegrate[integrand[t], {t, tmin, 100},
    WorkingPrecision -> 30, PrecisionGoal -> 15, MaxRecursion -> 100];
  AppendTo[results, {tmin, integral}];
  Print["t_min = ", tmin, ": I = ", N[integral, 12]],
  {tmin, {0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001, 0.0005, 0.0002, 0.0001}}
];

(* Analyze convergence *)
Print["\n=== Convergence analysis ==="];
Print["Successive differences:"];
diffs = Table[results[[i+1, 2]] - results[[i, 2]], {i, Length[results]-1}];
Do[
  Print["Delta(", results[[i, 1]], " -> ", results[[i+1, 1]], ") = ", N[diffs[[i]], 8]],
  {i, Length[diffs]}
];

(* Check if differences are decreasing *)
Print["\nRatio of successive differences (should approach 0 for convergence):"];
ratios = Table[diffs[[i+1]]/diffs[[i]], {i, Length[diffs]-1}];
Do[Print["Ratio ", i, ": ", N[ratios[[i]], 6]], {i, Length[ratios]}];

(* Richardson extrapolation for 1/sqrt(t) divergence *)
Print["\n=== Richardson extrapolation ==="];
(* If I(eps) = A/sqrt(eps) + B*log(eps) + C + O(sqrt(eps)) *)
(* Use three points to extract C *)

(* For I = b/sqrt(t) + c, we have:
   I(t1) = b/sqrt(t1) + c
   I(t2) = b/sqrt(t2) + c
   c = (I1*sqrt(t1) - I2*sqrt(t2)) / (sqrt(t1) - sqrt(t2)) *)

t1 = 0.001; t2 = 0.0001;
I1 = results[[7, 2]];  (* t_min = 0.001 *)
I2 = results[[10, 2]]; (* t_min = 0.0001 *)

Cextrapolated = (I1 * Sqrt[t1] - I2 * Sqrt[t2]) / (Sqrt[t1] - Sqrt[t2]);
Print["Extrapolated C (assuming I = b/sqrt(t) + c): ", N[Cextrapolated, 10]];

(* If I = A*log(t) + C *)
Clog = (I1 * Log[t2] - I2 * Log[t1]) / (Log[t2] - Log[t1]);
Print["Extrapolated C (assuming I = A*log(t) + C): ", N[Clog, 10]];

(* Final result *)
Print["\n========================================"];
Print["FINAL RESULT"];
Print["========================================"];

(* Use the log extrapolation as it's more stable *)
Cfinal = Clog;
K0 = 1/(2 Pi)^(3/2) * Exp[Cfinal];

Print["\nWith N_bulk = 8, M_coeff = -8 (tadpole cancelled):"];
Print["  Finite part C = ", N[Cfinal, 10]];
Print["  K_0 = (2 Pi)^(-3/2) * Exp[", N[Cfinal, 6], "]"];
Print["  K_0 = ", N[K0, 10]];
Print["  log(K_0) = ", N[Log[K0], 10]];
Print["========================================"];
