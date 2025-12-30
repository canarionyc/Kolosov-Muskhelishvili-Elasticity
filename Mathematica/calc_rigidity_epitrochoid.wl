(* ::Package:: *)

(* ::Title:: *)
(*Interactive Calculation of Torsional Rigidity for an Epitrochoid*)


(* ::Text:: *)
(*This notebook interactively calculates a key component of the torsional rigidity for a bar with an epitrochoidal cross-section, based on Muskhelishvili's complex variable methods.*)
(*The calculation uses the conformal mapping function \[Omega](z) = b(z + m z^n) to map the unit circle to the epitrochoid shape.*)
(**)
(*The integral term required for the torsional rigidity is calculated using the residue theorem.*)


Manipulate[
  (* Main calculation block *)
  Module[{omega, omegaSigma, omegaBarSigma, dOmegaSigma, integrand, residue, integralVal, z, t, plot},

    (* 1. Define Mapping Function and related expressions symbolically *)
    omega[zeta_] := b (zeta + m zeta^n);
    omegaSigma = omega[sigma];
    omegaBarSigma = b (1/sigma + m (1/sigma)^n); (* Conjugate on unit circle *)
    dOmegaSigma = D[omegaSigma, sigma];
    integrand = omegaSigma^2 * omegaBarSigma * dOmegaSigma;

    (* 2. Calculate Residue and the final Integral Value *)
    (* The residue at sigma=0 is the coefficient of the 1/sigma term in the series expansion of the integrand. *)
    residue = Residue[integrand, {sigma, 0}];
    integralVal = 2 * Pi * I * residue;

    (* 3. Generate the plot of the epitrochoid shape *)
    (* The condition for the map to be one-to-one (non-self-intersecting) is m <= 1/n *)
    plot = ParametricPlot[
      {b (Cos[t] + m Cos[n t]), b (Sin[t] + m Sin[n t])},
      {t, 0, 2 Pi},
      PlotRange -> All,
      ImageSize -> Medium,
      AxesLabel -> {"x", "y"},
      PlotLabel -> Style[StringForm["Epitrochoid Shape (m \[LessEqual] `` is non-intersecting)", NumberForm[1/n, {Infinity, 3}]], 12],
      ColorFunction -> "Rainbow"
    ];

    (* 4. Display the results in a clean grid *)
    Grid[{
      {Style["Inputs", Bold, 14], SpanFromLeft},
      {"Mapping Function, \[Omega](z)", TraditionalForm[omega[z]]},
      {"", ""}, (* Spacer *)
      {Style["Calculation on Unit Circle (\[Sigma])", Bold, 14], SpanFromLeft},
      {Style["Integrand", Italic], TraditionalForm[Defer[Collect[Expand[integrand], sigma]]]},
      {Style["Residue at \[Sigma]=0", Italic], Simplify[residue]},
      {Style["Integral Term Value (2\[Pi]i * Residue)", Bold, Italic], Style[Simplify[integralVal], Bold, Red, 14]},
      {"", ""}, (* Spacer *)
      {Style["Visualization", Bold, 14], SpanFromLeft},
      {plot, SpanFromLeft}
    }, Alignment -> Left, Spacings -> {1, 2}]
  ],

  (* Controls for Manipulate *)
  {{n, 3, "Number of lobes (n)"}, 2, 10, 1, Appearance -> "Labeled", ImageSize -> Small},
  {{b, 1, "Scale (b)"}, 0.1, 2, Appearance -> "Labeled", ImageSize -> Small},
  {{m, 0.2, "Shape factor (m)"}, 0, 0.5, Appearance -> "Labeled", ImageSize -> Small},

  (* Options for Manipulate *)
  ControlPlacement -> Left,
  Alignment -> Center,
  Paneled -> True,
  FrameLabel -> Style["Interactive Torsional Rigidity Calculation", Bold, 16]
]
