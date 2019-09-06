(* ::Package:: *)

(*Start Kernels for Parallel Computation*)
LaunchKernels[];

(*Set directories, form of output to terminal and start timer*)
SetOptions[$Output, FormatType->OutputForm];
currentPath=DirectoryName @ $InputFileName;
SetDirectory[currentPath];
overallTimeStart = AbsoluteTime[];

(*Global Variables*)
proj           = $ScriptCommandLine[[2]];
surfInput      = $ScriptCommandLine[[3]];
readFit        = ToExpression[$ScriptCommandLine[[4]]];
outputToVMEC   = ToExpression[$ScriptCommandLine[[5]]];
vmecInput      = $ScriptCommandLine[[6]];
vmecOutput     = $ScriptCommandLine[[7]];
nsurfaces      = ToExpression[$ScriptCommandLine[[8]]];
nthetaM        = ToExpression[$ScriptCommandLine[[9]]];
nphiM          = ToExpression[$ScriptCommandLine[[10]]];
deltac0        = ToExpression[$ScriptCommandLine[[11]]];
deltal0        = ToExpression[$ScriptCommandLine[[12]]];
deltalmin      = ToExpression[$ScriptCommandLine[[13]]];
deltalmax      = ToExpression[$ScriptCommandLine[[14]]];
muc0           = ToExpression[$ScriptCommandLine[[15]]];
mucMin         = ToExpression[$ScriptCommandLine[[16]]];
mucMax         = ToExpression[$ScriptCommandLine[[17]]];
nModes         = ToExpression[$ScriptCommandLine[[18]]];
maxiterations  = ToExpression[$ScriptCommandLine[[19]]];
plotFit        = ToExpression[$ScriptCommandLine[[20]]];
maxm           = ToExpression[$ScriptCommandLine[[21]]];
maxn           = ToExpression[$ScriptCommandLine[[22]]];
maxRecursTheta = ToExpression[$ScriptCommandLine[[23]]];
maxRecursPhi   = ToExpression[$ScriptCommandLine[[24]]];
pi = N[Pi,10];

(*Fast Numerical Integrator*)
GaussLegendreQuadrature[f_, {x_, a_, b_}, n_Integer: 10, prec_: MachinePrecision] := 
  Module[{nodes, weights}, {nodes, weights} = Most[NIntegrate`GaussRuleData[n, prec]];
   (b - a) weights.Map[Function[x, f], Rescale[nodes, {0, 1}, {a, b}]]];

If[readFit==1,
  timestart1 = AbsoluteTime[]; WriteString[$Output, "  Reading Fit From Text File...\n"];
  Get[StringJoin[currentPath,"data/",ToString[proj],"/",ToString[proj],"_fit_functions.mx"]];
  WriteString[$Output, "  Successful in "<>ToString[AbsoluteTime[]-timestart1]<>"s\n"];
];

(*Obtain Axis and Surface Parameters*)
If[FileExistsQ[vmecOutput],
    timestart1 = AbsoluteTime[]; WriteString[$Output, "  Read from VMEC file...\n"];
    Get[StringJoin[currentPath,"src/vmec_read.wls"]];
    WriteString[$Output, "  Successful in "<>ToString[AbsoluteTime[]-timestart1]<>"s\n"];
    ,
    If[FileExistsQ[surfInput],
      nsurfaces=1;
      timestart1 = AbsoluteTime[]; WriteString[$Output, "  Read from surf_input..."];
      Get[StringJoin[currentPath,"src/surfinput_read.wls"]];
      WriteString[$Output, " Successful in "<>ToString[AbsoluteTime[]-timestart1]<>"s\n"];
      ,
      WriteString[$Output, "No Input File to Read From.\n"];
      CloseKernels[];Exit[];
    ];
];

(*VMEC Axis Frenet-Serret Frame*)
timestart1 = AbsoluteTime[]; WriteString[$Output, "  Computing VMEC's Frenet-Serret Frame and Surface..."]; Off[Simplify::time];
curvVMEC[t_] = Chop[Norm[Cross[closedcurvVMEC'[t], closedcurvVMEC''[t]]]/Norm[closedcurvVMEC'[t]]^3, 10^-6]//Re;
torsVMEC[t_] = Chop[Dot[Cross[closedcurvVMEC'[t], closedcurvVMEC''[t]], closedcurvVMEC'''[t]]/Norm[Cross[closedcurvVMEC'[t], closedcurvVMEC''[t]]]^2, 10^-6];
sprimeVMEC[t_] = Chop[Simplify[Chop[ComplexExpand[Norm[closedcurvVMEC'[t]]], 10^-6], TimeConstraint -> 0.01], 10^-6];
b0VMEC[t_] = Simplify[Chop[closedcurvVMEC'[t]/sprimeVMEC[t], 10^-6], TimeConstraint -> 0.01];
k0VMEC[t_] = Simplify[Chop[b0VMEC'[t]/sprimeVMEC[t], 10^-6], TimeConstraint -> 0.01]/curvVMEC[t]//Re;
t0VMEC[t_] = Chop[Cross[b0VMEC[t], k0VMEC[t]], 10^-6]//Re;
(*VMEC Flux Surface*)
If[muDefined==0,
  Table[
    FluxSurfaceVMEC[theta_, phi_, iradius] = Simplify[Chop[{RBCVMEC[theta, phi, iradius]*Cos[phi], RBCVMEC[theta, phi, iradius]*Sin[phi], ZBSVMEC[theta, phi, iradius]}, 10^-6], TimeConstraint -> 0.01]//Re//Quiet;
  , {iradius,1,nsurfaces}];
]
WriteString[$Output, " Successful in "<>ToString[AbsoluteTime[]-timestart1]<>"s\n"];

(*Obtain Desired Surface in Mercier's Coordinates*)
If[muDefined==0,
  (*Mercier rho function in terms of axis quantities a.k.a. Mercier's theta and phi*)
  muVMEC = muc[0] + Sum[muc[i] Cos[vmecNFP i phi], {i, 1, nModes}];
  deltaVMEC = deltal*phi+deltac[0]+Sum[deltas[i] Sin[vmecNFP i phi],{i,1,nModes}];
  B0VMEC = B0c[0]+ Sum[B0c[i] Cos[vmecNFP i phi],{i,1,nModes}];
  rho[theta_, phi_] = Sqrt[(PSI*Sqrt[1 - muVMEC^2] )/(2*(B0VMEC)*(1 + muVMEC*Cos[2*(theta + deltaVMEC)]))];
  If[readFit==0,
    timestart1 = AbsoluteTime[]; WriteString[$Output, "  Fitting VMEC to Mercier's theta and rho...\n"];
    Get[StringJoin[currentPath,"src/fit_to_mercier.wls"]];
    WriteString[$Output, "  Successful in "<>ToString[AbsoluteTime[]-timestart1]<>"s\n"];
  ]
  ,
  timestart1 = AbsoluteTime[]; WriteString[$Output, "  Constructing Mercier's theta and rho..."];
  rho[theta_, phi_] = Sqrt[(vmecPSI[[nsurfaces]]*Sqrt[1 - mu[phi]^2] )/(2*(B0[phi]) (1 + mu[phi]*Cos[2 (theta + delta[phi])]))];
  FluxSurfacenFitVMEC[theta_, s_, 1] = Chop[closedcurvVMEC[s] + (rho[theta, s]) (Cos[theta] k0VMEC[s] + Sin[theta] t0VMEC[s]), 10^-6]//Quiet;
	RBCFit[phi_, theta_, 1] = Re[Cos[phi]*FluxSurfacenFitVMEC[theta, phi, 1][[1]] + Sin[phi]*FluxSurfacenFitVMEC[theta, phi, 1][[2]]]//Quiet;
	ZBSFit[phi_, theta_, 1] = Re[FluxSurfacenFitVMEC[theta, phi, 1][[3]]]//Quiet;
  WriteString[$Output, " Successful in "<>ToString[AbsoluteTime[]-timestart1]<>"s\n"];
];

(*Compute Iota on Axis from Expression of Near-Axis Expansion*)
timestart1 = AbsoluteTime[]; WriteString[$Output, "  Computing total rotation angle N of the normal vector..."];
quadrant = 0; nNormal = 0; diffQuadrant = 0; qnew=0;
Table[
  normalx = k0VMEC[t][[1]] Cos[t] + k0VMEC[t][[2]] Sin[t];
  normaly = k0VMEC[t][[3]];
  If[Abs[normalx] > Abs[normaly],If[normalx > 0, qnew = 2, qnew = 4],If[normaly > 0, qnew = 1, qnew = 3]];
  If[quadrant == 0,quadrant = qnew, diffQuadrant = qnew - quadrant;];
  If[Abs[diffQuadrant] == 3, diffQuadrant = diffQuadrant/3];nNormal = nNormal + diffQuadrant/2; quadrant = qnew;
, {t, 0, 2*pi, 2*pi/vmecNFP/15}];
WriteString[$Output, " N="<>ToString[nNormal]<>" Successful in "<>ToString[AbsoluteTime[]-timestart1]<>"s\n"];
If[muDefined==1,
  iota1 = GaussLegendreQuadrature[(D[delta[phi], phi] - torsVMEC[phi]) Sqrt[1 - (mu[phi])^2], {phi, 0, 2*pi}, 70, 20]/(2*pi);
  iota2 = (delta[2*pi] - delta[0])/(2*pi);
  iotaFit = iota1 - iota2 - nNormal;
  Print["  Iota from near-axis formula: "<>ToString[iotaFit]];
  ,
  iotaFit=ParallelTable[
    iota1 = GaussLegendreQuadrature[(D[deltaVMEC /. fitParams[iradius], phi] - torsVMEC[phi]) Sqrt[1 - (muVMEC /. fitParams[iradius])^2], {phi, 0, 2*pi}, 70, 20]/(2*pi);
    iota2 = (((deltaVMEC /. fitParams[iradius]) /. phi -> 2*pi) - ((deltaVMEC /. fitParams[iradius]) /. phi -> 0))/(2*pi);
    iota2 - iota1 - nNormal,
  {iradius,1,nsurfaces}];
  If[FileExistsQ[vmecOutput],
    Print["  Iota from VMEC: "<>ToString[iotaAxisOut]<>", Iota from near-axis: "<>ToString[iotaFit]];
    ,
    Print["  Iota from near-axis formula: "<>ToString[iotaFit]];
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