(* ::Package:: *)
overallTimeStart = AbsoluteTime[];

SetOptions[$Output, FormatType->OutputForm];

currentPath=DirectoryName @ $InputFileName;
SetDirectory[currentPath];
proj         = $ScriptCommandLine[[2]];
readInput    = ToExpression[$ScriptCommandLine[[3]]];
surfInput    = $ScriptCommandLine[[4]];
outputToVMEC = ToExpression[$ScriptCommandLine[[5]]];
vmecInput    = $ScriptCommandLine[[6]];
nthetaM      = ToExpression[$ScriptCommandLine[[8]]];
nphiM        = ToExpression[$ScriptCommandLine[[9]]];
plotFit      = ToExpression[$ScriptCommandLine[[19]]];
pi = N[Pi,10];

(*Fast Numerical Integrator*)
GaussLegendreQuadrature[f_, {x_, a_, b_}, n_Integer: 10, prec_: MachinePrecision] := 
  Module[{nodes, weights}, {nodes, weights} = Most[NIntegrate`GaussRuleData[n, prec]];
   (b - a) weights.Map[Function[x, f], Rescale[nodes, {0, 1}, {a, b}]]];

(*Obtain Axis and Surface Parameters*)
If[readInput == 1,
    timestart1 = AbsoluteTime[]; WriteString[$Output, "  Read from surf_input..."];
    Get[StringJoin[currentPath,"src/surfinput_read.wls"]];
    WriteString[$Output, " Successful in "<>ToString[AbsoluteTime[]-timestart1]<>"s\n"];
    ,
    timestart1 = AbsoluteTime[]; WriteString[$Output, "  Read from VMEC file...\n"];
    Get[StringJoin[currentPath,"src/vmec_read.wls"]];
    WriteString[$Output, "  Successful in "<>ToString[AbsoluteTime[]-timestart1]<>"s\n"];
]

(*VMEC Axis Frenet-Serret Frame*)
timestart1 = AbsoluteTime[]; WriteString[$Output, "  Computing Frenet-Serret Frame..."]; Off[Simplify::time];
curvVMEC[t_] = Chop[Norm[Cross[closedcurvVMEC'[t], closedcurvVMEC''[t]]]/Norm[closedcurvVMEC'[t]]^3, 10^-6]//Re;
torsVMEC[t_] = Chop[Dot[Cross[closedcurvVMEC'[t], closedcurvVMEC''[t]], closedcurvVMEC'''[t]]/Norm[Cross[closedcurvVMEC'[t], closedcurvVMEC''[t]]]^2, 10^-6];
sprimeVMEC[t_] = Chop[Simplify[Chop[ComplexExpand[Norm[closedcurvVMEC'[t]]], 10^-6], TimeConstraint -> 0.01], 10^-6];
b0VMEC[t_] = Simplify[Chop[closedcurvVMEC'[t]/sprimeVMEC[t], 10^-6], TimeConstraint -> 0.01];
k0VMEC[t_] = Simplify[Chop[b0VMEC'[t]/sprimeVMEC[t], 10^-6], TimeConstraint -> 0.01]/curvVMEC[t]//Re;
t0VMEC[t_] = Chop[Cross[b0VMEC[t], k0VMEC[t]], 10^-6]//Re;
FluxSurfaceVMEC[theta_, phi_] = Simplify[Chop[{RBCVMEC[theta, phi] Cos[phi], RBCVMEC[theta, phi] Sin[phi], ZBSVMEC[theta, phi]}, 10^-6], TimeConstraint -> 0.01]//Re;
WriteString[$Output, " Successful in "<>ToString[AbsoluteTime[]-timestart1]<>"s\n"];

(*Start Kernels for Parallel Computation*)
LaunchKernels[];

(*Obtain Desired Surface in Mercier's Coordinates*)
If[muDefined==0,
  timestart1 = AbsoluteTime[]; WriteString[$Output, "  Fitting VMEC to Mercier's theta and rho...\n"];
  Get[StringJoin[currentPath,"src/fit_to_mercier.wls"]];
  WriteString[$Output, "  Successful in "<>ToString[AbsoluteTime[]-timestart1]<>"s\n"];
  ,
  timestart1 = AbsoluteTime[]; WriteString[$Output, "  Constructing Mercier's theta and rho..."];
  rho[theta_, phi_] = Sqrt[(vmecPSI*Sqrt[1 - mu[phi]^2] )/(2*(B0[phi]) (1 + mu[phi]*Cos[2 (theta + delta[phi])]))];
  FluxSurfacenFitVMEC[theta_, s_] = Chop[closedcurvVMEC[s] + (rho[theta, s]) (Cos[theta] k0VMEC[s] + Sin[theta] t0VMEC[s]), 10^-6]//Quiet;
	RBCFit[phi_, theta_] = Re[Cos[phi]*FluxSurfacenFitVMEC[theta, phi][[1]] + Sin[phi]*FluxSurfacenFitVMEC[theta, phi][[2]]]//Quiet;
	ZBSFit[phi_, theta_] = Re[FluxSurfacenFitVMEC[theta, phi][[3]]]//Quiet;
  WriteString[$Output, " Successful in "<>ToString[AbsoluteTime[]-timestart1]<>"s\n"];
];

(*Compute Iota on Axis from Expression of Near-Axis Expansion*)
If[muDefined==1,
  iota0VMEC1 = GaussLegendreQuadrature[(D[delta[phi], phi] - torsVMEC[phi]) Sqrt[1 - (mu[phi])^2], {phi, 0, 2*pi}, 70, 20]/(2*pi);
  iota0VMEC2 = (delta[2*pi] - delta[0])/(2*pi);
  iotaFit = iota0VMEC1 - iota0VMEC2;
  Print["  Iota from near-axis formula: "<>ToString[iotaFit]];
  ,
  iota0VMEC1 = GaussLegendreQuadrature[(D[deltaVMEC /. fitParams, phi] - torsVMEC[phi]) Sqrt[1 - (muVMEC /. fitParams)^2], {phi, 0, 2*pi}, 70, 20]/(2*pi);
  iota0VMEC2 = (((deltaVMEC /. fitParams) /. phi -> 2*pi) - ((deltaVMEC /. fitParams) /. phi -> 0))/(2*pi);
  iotaFit = iota0VMEC1 - iota0VMEC2;
  If[readInput==1,
    Print["  Iota from near-axis formula: "<>ToString[iotaFit]];
    ,
    Print["  Iota from VMEC: "<>ToString[iotaAxisOut]<>", Iota from near-axis: "<>ToString[iotaFit]];
  ]
];

(*Plot Obtained Surfaces*)
If[plotFit==1,
  timestart1 = AbsoluteTime[]; WriteString[$Output, "  Creating Figures...\n"];
  Get[StringJoin[currentPath,"src/plot_mercier.wls"]];
  WriteString[$Output, "  Successful in "<>ToString[AbsoluteTime[]-timestart1]<>"s\n"];
]

(*Output VMEC file*)
If[outputToVMEC == 1,
   timestart1 = AbsoluteTime[]; WriteString[$Output, "  Mercier to VMEC...\n"]
   Get[StringJoin[currentPath,"src/mercier_to_vmec.wls"]];
   WriteString[$Output, "  Successful in "<>ToString[AbsoluteTime[]-timestart1]<>"s\n"];
]

WriteString[$Output,"  SENAC Mathematica took "<>ToString[AbsoluteTime[]-overallTimeStart]<>"s\n"]

CloseKernels[];

(*Pause[5]*)