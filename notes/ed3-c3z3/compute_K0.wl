(* K_0 computation for C^3/Z_3 orientifold *)
(* Clean implementation with proper modular transformation *)

(* Theta functions using Mathematica's built-in *)
theta00[tau_] := EllipticTheta[3, 0, Exp[I Pi tau]]
theta01[tau_] := EllipticTheta[4, 0, Exp[I Pi tau]]
theta10[tau_] := EllipticTheta[2, 0, Exp[I Pi tau]]

(* Verify theta function identities *)
Print["=== Verifying theta functions ==="];
tau = I;
Print["At tau = i:"];
Print["  theta00 = ", N[theta00[tau], 10]];
Print["  theta01 = ", N[theta01[tau], 10]];
Print["  theta10 = ", N[theta10[tau], 10]];
Print["  Jacobi identity |th00^4 - th01^4 - th10^4| = ",
      N[Abs[theta00[tau]^4 - theta01[tau]^4 - theta10[tau]^4], 10]];

(* Partition functions *)
(* Note: Mathematica's EllipticTheta uses nome q = Exp[I Pi tau], not q = Exp[2 I Pi tau] *)
(* So theta_ab(tau) in our convention = EllipticTheta[..., Exp[I Pi tau]] *)

ZAnnulus[t_, Nbulk_] := Module[{tau2, q, th10, th00},
  tau2 = 2 I t;  (* 2*tau for DN sector *)
  q = Exp[-2 Pi t];
  th10 = EllipticTheta[2, 0, Exp[I Pi tau2]];
  th00 = EllipticTheta[3, 0, Exp[I Pi tau2]];
  (* Z_A = N/2 * [th10 * q^(-1/4) + th00] *)
  Re[Nbulk/2 * (th10 * q^(-1/4) + th00)]
]

ZMobius[t_, Mcoeff_] := Module[{hattau2, hatq, th10hat, th00hat},
  hattau2 = 2 (I t + 1/2);  (* 2*hat_tau *)
  hatq = Exp[-2 Pi t + I Pi];  (* = -Exp[-2 Pi t] *)
  th10hat = EllipticTheta[2, 0, Exp[I Pi hattau2]];
  th00hat = EllipticTheta[3, 0, Exp[I Pi hattau2]];
  (* Z_M = Mcoeff * [th10hat * hatq^(-1/4) - th00hat] *)
  Mcoeff * Re[th10hat * hatq^(-1/4) - th00hat]
]

Zsub[t_] := 3 (Exp[-2 Pi t] - 1)

Ztotal[t_, Nbulk_, Mcoeff_] := ZAnnulus[t, Nbulk] + ZMobius[t, Mcoeff] + Zsub[t]

(* Test partition functions *)
Print["\n=== Partition function values ==="];
Print["Using N_bulk = 8 (= 24 fractional), M_coeff = 1:"];
Do[
  Print["t = ", t, ": Z_A = ", N[ZAnnulus[t, 8], 6],
        ", Z_M = ", N[ZMobius[t, 1], 6],
        ", Z_sub = ", N[Zsub[t], 6],
        ", Total = ", N[Ztotal[t, 8, 1], 6]],
  {t, {0.1, 0.5, 1.0, 5.0}}
]

(* Find asymptotic value *)
Print["\n=== Asymptotic behavior ==="];
Zasymptotic[Nbulk_, Mcoeff_] := Ztotal[20, Nbulk, Mcoeff]
Print["Z(t -> infinity) for N=8, M=1: ", N[Zasymptotic[8, 1], 10]];

(* The integrand *)
integrand[t_, Nbulk_, Mcoeff_] := Module[{Zasymp},
  Zasymp = Zasymptotic[Nbulk, Mcoeff];
  (Ztotal[t, Nbulk, Mcoeff] - Zasymp) / (2 t)
]

(* Small-t behavior analysis *)
Print["\n=== Small-t behavior ==="];
Print["Checking Z(t) * sqrt(2t) at small t (should approach constant if Z ~ t^(-1/2)):"];
Do[
  Print["t = ", t, ": Z*sqrt(2t) = ", N[Ztotal[t, 8, 1] * Sqrt[2 t], 6]],
  {t, {0.01, 0.02, 0.05, 0.1}}
]

(* The coefficient of t^(-1/2) should cancel for tadpole cancellation *)
(* At small t: Z_A ~ N * (2t)^(-1/2), Z_M ~ Mcoeff * alpha * (2t)^(-1/2) *)
(* where alpha comes from the phase in hatq^(-1/4) *)

Print["\n=== Finding tadpole-cancelling coefficient ==="];
(* Compute the small-t coefficient *)
smallTcoeff[Nbulk_, Mcoeff_] := Module[{t = 0.001},
  Ztotal[t, Nbulk, Mcoeff] * Sqrt[2 t]
]

Print["Small-t coefficient Z*sqrt(2t) for various M_coeff:"];
Do[
  Print["M_coeff = ", M, ": coeff = ", N[smallTcoeff[8, M], 6]],
  {M, {1, 10, 20, 30, 40, 50, 60, 70, 80}}
]

(* Find M where coefficient = 0 *)
Moptimal = M /. FindRoot[smallTcoeff[8, M] == 0, {M, 50}];
Print["\nOptimal M_coeff for tadpole cancellation: ", N[Moptimal, 10]];

(* Verify *)
Print["Verification at small t with optimal M_coeff:"];
Do[
  Print["t = ", t, ": Z*sqrt(2t) = ", N[Ztotal[t, 8, Moptimal] * Sqrt[2 t], 10]],
  {t, {0.001, 0.002, 0.005, 0.01}}
]

(* Now compute the integral with the optimal coefficient *)
Print["\n=== Computing integral with tadpole-cancelled coefficient ==="];
Zasymp = Zasymptotic[8, Moptimal];
Print["Z_asymptotic = ", N[Zasymp, 10]];

(* The integrand should now be integrable *)
integrandOpt[t_?NumericQ] := (Ztotal[t, 8, Moptimal] - Zasymp) / (2 t)

(* Test integrability *)
Print["\nIntegrand values:"];
Do[
  Print["t = ", t, ": integrand = ", N[integrandOpt[t], 6]],
  {t, {0.01, 0.05, 0.1, 0.5, 1.0, 5.0}}
]

(* Numerical integration *)
Print["\n=== Numerical integration ==="];
Print["Testing convergence as t_min -> 0:"];

results = Table[
  {tmin, NIntegrate[integrandOpt[t], {t, tmin, 30},
    WorkingPrecision -> 20, MaxRecursion -> 50]},
  {tmin, {0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001}}
];

Do[Print["t_min = ", r[[1]], ": integral = ", r[[2]]], {r, results}];

(* Check for convergence *)
Print["\nDifferences (should decrease for convergence):"];
diffs = Table[results[[i+1, 2]] - results[[i, 2]], {i, Length[results]-1}];
Do[Print["Delta I (", results[[i, 1]], " -> ", results[[i+1, 1]], ") = ", diffs[[i]]],
   {i, Length[diffs]}];

(* If converging, extrapolate to t_min = 0 *)
Print["\n=== Final result ==="];
(* Use Richardson extrapolation or limit *)
finalIntegral = Limit[NIntegrate[integrandOpt[t], {t, eps, 30},
  WorkingPrecision -> 15], eps -> 0];
Print["Attempting limit..."];

(* If limit doesn't work, use smallest stable value *)
Print["Using smallest t_min result as estimate:"];
Cvalue = results[[-1, 2]];
Print["C = ", Cvalue];

K0 = 1/(2 Pi)^(3/2) * Exp[Cvalue];
Print["\nK_0 = (2 Pi)^(-3/2) * Exp[", Cvalue, "]"];
Print["K_0 = ", N[K0, 10]];
Print["log(K_0) = ", N[Log[K0], 10]];
