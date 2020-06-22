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
ordern         = ToExpression[$ScriptCommandLine[[8]]];
nsurfaces      = ToExpression[$ScriptCommandLine[[9]]];
nthetaM        = ToExpression[$ScriptCommandLine[[10]]];
nphiM          = ToExpression[$ScriptCommandLine[[11]]];
deltac0        = ToExpression[$ScriptCommandLine[[12]]];
deltal0        = ToExpression[$ScriptCommandLine[[13]]];
deltalmin      = ToExpression[$ScriptCommandLine[[14]]];
deltalmax      = ToExpression[$ScriptCommandLine[[15]]];
muc0           = ToExpression[$ScriptCommandLine[[16]]];
mucMin         = ToExpression[$ScriptCommandLine[[17]]];
mucMax         = ToExpression[$ScriptCommandLine[[18]]];
nModes         = ToExpression[$ScriptCommandLine[[19]]];
maxiterations  = ToExpression[$ScriptCommandLine[[20]]];
plotFit        = ToExpression[$ScriptCommandLine[[21]]];
plotOriginal   = ToExpression[$ScriptCommandLine[[22]]];
maxm           = ToExpression[$ScriptCommandLine[[23]]];
maxn           = ToExpression[$ScriptCommandLine[[24]]];
maxRecursTheta = ToExpression[$ScriptCommandLine[[25]]];
maxRecursPhi   = ToExpression[$ScriptCommandLine[[26]]];
nPlotTheta     = ToExpression[$ScriptCommandLine[[27]]];
nPlotPhi       = ToExpression[$ScriptCommandLine[[28]]];
plotPointsFig  = ToExpression[$ScriptCommandLine[[29]]];
maxRecursPlot  = ToExpression[$ScriptCommandLine[[30]]];
ImageSizePlot  = ToExpression[$ScriptCommandLine[[31]]];
ImageResolutionPlot  = ToExpression[$ScriptCommandLine[[32]]];
nfigsSurf            = ToExpression[$ScriptCommandLine[[33]]];
nPlots               = ToExpression[$ScriptCommandLine[[34]]];
nthetapointsBsurface = ToExpression[$ScriptCommandLine[[35]]];
ntphipointsBsurface  = ToExpression[$ScriptCommandLine[[36]]];
npointsPolPlots      = ToExpression[$ScriptCommandLine[[37]]];
exportBFieldSurface  = ToExpression[$ScriptCommandLine[[38]]];
keepfit              = ToExpression[$ScriptCommandLine[[39]]];
export3DSurface      = ToExpression[$ScriptCommandLine[[40]]];
readlowfit           = ToExpression[$ScriptCommandLine[[41]]];
plotPolFig           = ToExpression[$ScriptCommandLine[[42]]];
quasisymmetry        = ToExpression[$ScriptCommandLine[[43]]];
surfInputQS          = $ScriptCommandLine[[44]];
phiedge              = ToExpression[$ScriptCommandLine[[45]]];
Bzero                = ToExpression[$ScriptCommandLine[[46]]];
If[quasisymmetry==1,surfInput=surfInputQS];
pi = N[Pi,10];
mu0const = 4*pi*10^(-7);
$HistoryLength = 0;

(*Fast Numerical Integrator*)
GaussLegendreQuadrature[f_, {x_, a_, b_}, n_Integer: 10, prec_: MachinePrecision] := 
  Module[{nodes, weights}, {nodes, weights} = Most[NIntegrate`GaussRuleData[n, prec]];
   (b - a) weights.Map[Function[x, f], Rescale[nodes, {0, 1}, {a, b}]]];

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
timestart1 = AbsoluteTime[]; WriteString[$Output, "  Computing VMEC's Frenet-Serret Frame and Surfaces..."]; Off[Simplify::time];
curvVMEC[t_] = Chop[Norm[Cross[closedcurvVMEC'[t], closedcurvVMEC''[t]]]/ComplexExpand[Norm[closedcurvVMEC'[t]]^3, 10^-6]]//Re;
torsVMEC[t_] = Chop[Dot[Cross[closedcurvVMEC'[t], closedcurvVMEC''[t]], closedcurvVMEC'''[t]]/ComplexExpand[Norm[Cross[closedcurvVMEC'[t], closedcurvVMEC''[t]]]]^2, 10^-6];
sprimeVMEC[t_] = Chop[Simplify[Chop[ComplexExpand[Norm[closedcurvVMEC'[t]]], 10^-6], TimeConstraint -> 0.01], 10^-6];
b0VMEC[t_] = Simplify[Chop[closedcurvVMEC'[t]/sprimeVMEC[t], 10^-6], TimeConstraint -> 0.01];
k0VMEC[t_] = Chop[b0VMEC'[t]/sprimeVMEC[t]/curvVMEC[t], 10^-6]//Re;(*Simplify[Chop[b0VMEC'[t]/sprimeVMEC[t], 10^-6], TimeConstraint -> 0.01]/curvVMEC[t]//Re;*)
t0VMEC[t_] = Chop[Cross[b0VMEC[t], k0VMEC[t]], 10^-6]//Re;
(*VMEC Flux Surface*)
If[muDefined==0,
  Table[
    FluxSurfaceVMEC[theta_, phi_, iradius] = Simplify[Chop[{RBCVMEC[theta, phi, iradius]*Cos[phi], RBCVMEC[theta, phi, iradius]*Sin[phi], ZBSVMEC[theta, phi, iradius]}, 10^-6], TimeConstraint -> 0.01]//Re//Quiet;
  , {iradius,1,nsurfaces}];
]
WriteString[$Output, " Successful in "<>ToString[AbsoluteTime[]-timestart1]<>"s\n"];

(*Compute the total rotation angle N of the normal Frenet-Serret vector*)
timestart1 = AbsoluteTime[]; WriteString[$Output, "  Computing total rotation angle N of the normal vector..."];
quadrant = 0; nNormal = 0; diffQuadrant = 0; qnew=0; printv = "";nl = "\n";
Table[
  normalx = k0VMEC[t][[1]]*Cos[t] + k0VMEC[t][[2]]*Sin[t];
  normaly = k0VMEC[t][[3]];
  If[Abs[normalx] > Abs[normaly],If[normalx > 0, qnew = 2, qnew = 4],If[normaly > 0, qnew = 1, qnew = 3]];
  If[quadrant == 0,quadrant = qnew, diffQuadrant = qnew - quadrant;];
  If[Abs[diffQuadrant] == 3, diffQuadrant = diffQuadrant/3];nNormal = nNormal + diffQuadrant/2; quadrant = qnew;
, {t, 0, 2*pi, 2*pi/vmecNFP/15}];
WriteString[$Output, " N="<>ToString[nNormal]<>" Successful in "<>ToString[AbsoluteTime[]-timestart1]<>"s\n"];

(*Obtain Desired Surface in Mercier's Coordinates*)
If[muDefined==0,
  (*Construct Rho of Psi*)
  timestart1 = AbsoluteTime[]; WriteString[$Output, "  Constructing rho of psi..."];
  Get[StringJoin[currentPath,"src/rho_of_psi.wls"]];
  WriteString[$Output, " Successful in "<>ToString[AbsoluteTime[]-timestart1]<>"s\n"];
  (*Perform Fit*)
  timestart1 = AbsoluteTime[]; WriteString[$Output, "  Fitting VMEC to Mercier's theta and rho...\n"];
  Get[StringJoin[currentPath,"src/fit_to_mercier.wls"]];
  WriteString[$Output, "  Successful in "<>ToString[AbsoluteTime[]-timestart1]<>"s\n"];
  ,
  timestart1 = AbsoluteTime[]; WriteString[$Output, "  Constructing Mercier's theta and rho..."];
  psin[theta_, phi_, 2]=pi*B0[phi]*(1 + mu[phi]*Cos[2 (theta + delta[phi])])/Sqrt[1 - mu[phi]^2];
  psin[theta_, phi_, 3]=pi*B0[phi]*(psi31[phi]*Cos[theta]+psi32[phi]*Sin[theta]+psi33[phi]*Cos[3*theta]+psi34[phi]*Sin[3*theta]);
  rho[theta_, phi_]=Sqrt[vmecPSI[[nsurfaces]]/psin[theta, phi, 2]] - vmecPSI[[nsurfaces]]*psin[theta, phi, 3]/(2*psin[theta, phi, 2]^2);
  rhotrue[theta_, phi_]:=x/.Solve[{psin[theta, phi, 2]*x^2+psin[theta, phi, 3]*x^3==vmecPSI[[nsurfaces]], x>0},x,Reals][[1]];
  FluxSurfacenFitVMEC[theta_, s_, 1] := Chop[closedcurvVMEC[s] + (rho[theta, s])*(Cos[theta]*k0VMEC[s] + Sin[theta]*t0VMEC[s]), 10^-6]//Quiet;
  FluxSurfacenFitVMECtrue[theta_, s_, 1] := Chop[closedcurvVMEC[s] + (rhotrue[theta, s])*(Cos[theta]*k0VMEC[s] + Sin[theta]*t0VMEC[s]), 10^-6]//Quiet;
	RBCFit[phi_, theta_, 1] := Re[Cos[phi]*FluxSurfacenFitVMEC[theta, phi, 1][[1]] + Sin[phi]*FluxSurfacenFitVMEC[theta, phi, 1][[2]]]//Quiet;
	ZBSFit[phi_, theta_, 1] := Re[FluxSurfacenFitVMEC[theta, phi, 1][[3]]]//Quiet;
  RBCFittrue[phi_, theta_, 1] := Re[Cos[phi]*FluxSurfacenFitVMECtrue[theta, phi, 1][[1]] + Sin[phi]*FluxSurfacenFitVMECtrue[theta, phi, 1][[2]]]//Quiet;
	ZBSFittrue[phi_, theta_, 1] := Re[FluxSurfacenFitVMECtrue[theta, phi, 1][[3]]]//Quiet;
  WriteString[$Output, " Successful in "<>ToString[AbsoluteTime[]-timestart1]<>"s\n"];
];

(*Compute Iota on Axis from Expression of Near-Axis Expansion*)
If[muDefined==1,
  printv = printv <> "Iota on Axis From Near-Axis Formula" <> nl;
  iotaFit = GaussLegendreQuadrature[(-torCurVMEC*sprimeVMEC[phi]/2-D[delta[phi], phi] + torsVMEC[phi]*sprimeVMEC[phi])*Sqrt[1 - (mu[phi])^2]+D[delta[phi], phi], {phi, 0, 2*pi}, 250, 20]/(2*pi)-nNormal;
  Print["  Iota on Axis from near-axis formula: "<>ToString[iotaFit]];
  printv = printv <> ToString[iotaFit] <> nl;
  If[quasisymmetry==1,
    iotaQS = StringCases[StringDelete[FindList[surfInput, "iota"][[1]], " "], x : NumberString :> ToExpression[x]][[1]];
    Print["  Iota on Axis from QS Garren-Boozer: "<>ToString[-iotaQS]];
  ]
  ,
  iotaFit=ParallelTable[
    GaussLegendreQuadrature[(-torCurVMEC*sprimeVMEC[phi]/2-D[deltaVMEC /. fitParams[iradius], phi] + torsVMEC[phi]*sprimeVMEC[phi])*Sqrt[1 - (muVMEC /. fitParams[iradius])^2]+D[deltaVMEC /. fitParams[iradius], phi], {phi, 0, 2*pi}, 250, 20]/(2*pi)-nNormal,
  {iradius,1,nsurfaces}];
  If[FileExistsQ[vmecOutput],
    Print["  Iota from VMEC:      "<>ToString[iotaAxisOut]];
    Print["  Iota from near-axis: "<>ToString[iotaFit]];
    printv = printv <> "Iota on Axis From Near-Axis Formula" <> nl;
    printv = printv <> ToString[iotaFit] <> nl;
    printv = printv <> "Iota From VMEC" <> nl;
    printv = printv <> ToString[iotaAxisOut] <> nl;
    ,
    Print["  Iota on Axis from near-axis formula: "<>ToString[iotaFit]];
    printv = printv <> "Iota on Axis From Near-Axis Formula" <> nl;
    printv = printv <> ToString[iotaFit] <> nl;
  ]
];

(*Compute Magnetic Well*)
If[FileExistsQ[vmecOutput],
  psia = Last[psiiVMEC];
  shalf = Table[(1. (i - 1/2))/(nns - 1), {i, 1, nns - 1, 1}];
  data = Drop[gmnc[[All, 1]], 1];
  fitRes = Fit[Transpose@{shalf[[1 ;; 15]], data[[1 ;; 15]]}, {1, x, x^2},x][[2]] /. x -> 1;
  Vwell0VMEC = -((2*pi)^2)*fitRes/psia^2;
]

(* Expression with integration by parts*)
(*
VppTemp[s_]=0.5*1/(4 Pi^2 B0R[s]^4 (1 - muR[s]^2)^(5/2)) (3 Pi (-1 + muR[s]^2)^2 (B0R'[s]*(sprimeVMEC[s]^(-1)))^2 + Pi B0R[s]^2 (-(-1 + muR[s]^2)^2 (-lambda^2 + curvVMEC[s]^2 (2 - 2 Cos[2 deltaR[s]] muR[s])
            + muR[s]^2 (lambda + 2 torsVMEC[s] - 2 deltaR'[s]*(sprimeVMEC[s]^(-1))) (lambda - 2 torsVMEC[s] + 2 deltaR'[s]*(sprimeVMEC[s]^(-1)))) + (muR'[s]*(sprimeVMEC[s]^(-1)))^2) +
            4 B0R[s] (-1 + muR[s]^2) (-Sqrt[1 - muR[s]^2] (Pi^2 pprime0 + 2 curvVMEC[s] (Cos[deltaR[s]] Psi31c[s] + 
            Sin[deltaR[s]] Psi31s[s])) + muR[s]^2 Sqrt[1 - muR[s]^2] (Pi^2 pprime0 + curvVMEC[s] (-Cos[deltaR[s]] (Psi31c[s] + 3 Psi33c[s]) - 
            Sin[deltaR[s]] (Psi31s[s] - 3 Psi33s[s]))) + muR[s] (3 curvVMEC[s] Sqrt[1 - muR[s]^2] (Cos[deltaR[s]] (Psi31c[s] + Psi33c[s]) + 
            Sin[deltaR[s]] (-Psi31s[s] + Psi33s[s])) + Pi B0R'[s]*(sprimeVMEC[s]^(-1))*muR'[s]*(sprimeVMEC[s]^(-1)))));
*)

(* Expression without integration by parts*)
VppTemp[s_]=0.5*(1/(4 Pi^2 B0R[s]^3 (1 - muR[s]^2)^(5/2)))((-1 + muR[s]) (1 + muR[s]) (-8 curvVMEC[s] Sqrt[1 - muR[s]^2] (Cos[deltaR[s]] Psi31c[s] + 
        Sin[deltaR[s]] Psi31s[s]) + muR[s] (12 curvVMEC[s] Sqrt[1 - muR[s]^2] (Cos[deltaR[s]] (Psi31c[s] + Psi33c[s]) + Sin[deltaR[s]] (-Psi31s[s] + Psi33s[s])) 
        + Pi B0R'[s]*(sprimeVMEC[s]^(-1)) muR'[s]*(sprimeVMEC[s]^(-1))) - Pi (4 Pi pprime0 Sqrt[1 - muR[s]^2] + D[B0R'[s]*(sprimeVMEC[s]^(-1)),s]*(sprimeVMEC[s]^(-1))) 
        + muR[s]^2 (-4 curvVMEC[s] Sqrt[1 - muR[s]^2] (Cos[deltaR[s]] (Psi31c[s] + 3 Psi33c[s]) + 
           Sin[deltaR[s]] (Psi31s[s] - 3 Psi33s[s])) + Pi (4 Pi pprime0 Sqrt[1 - muR[s]^2] + D[B0R'[s]*(sprimeVMEC[s]^(-1)),s]*(sprimeVMEC[s]^(-1))))) 
           + Pi B0R[s] (lambda^2 + 2 curvVMEC[s]^2 (-1 + Cos[2 deltaR[s]] muR[s]) (-1 + muR[s]^2)^2 + muR[s] (-muR[s] (3 lambda^2 + muR[s]^2 (-3 lambda^2 + 
              8 (torsVMEC[s] - deltaR'[s]*(sprimeVMEC[s]^(-1)))^2) - 4 (torsVMEC[s] - deltaR'[s]*(sprimeVMEC[s]^(-1)))^2 + muR[s]^4 (lambda + 2 torsVMEC[s] - 
              2 deltaR'[s]*(sprimeVMEC[s]^(-1))) (lambda - 2 torsVMEC[s] + 2 deltaR'[s]*(sprimeVMEC[s]^(-1))) + 2 muR'[s]*(sprimeVMEC[s]^(-1))^2) 
              + (-1 + muR[s]^2) D[muR'[s]*(sprimeVMEC[s]^(-1)),s]*(sprimeVMEC[s]^(-1)))));

(*Compute DWell*)
Laxis=GaussLegendreQuadrature[sprimeVMEC[phi], {phi, 0, 2*pi}, 250, 20];
DwellTemp[s_]=(pprime0/((2*pi)^1))*(Laxis/(4*pi*PsiI))*(VppP/Laxis-pprime0/(B0R[s]^3));

(*Integrate them over s to get final answer*)
If[muDefined==1,
  pprime0      = pprime0VMEC;
  lambda       = torCurVMEC;
  B0R[phi_]    = B0[phi];
  deltaR[phi_] = delta[phi];
  muR[phi_]    = mu[phi];
  Psi31c[phi_] = Pi*B0[phi]*(Cos[  deltaR[s]]*psi31[phi] - psi32[phi]*Sin[deltaR[s]]);
  Psi31s[phi_] = Pi*B0[phi]*(Cos[  deltaR[s]]*psi32[phi] + psi31[phi]*Sin[deltaR[s]]);
  Psi33c[phi_] = Pi*B0[phi]*(Cos[3*deltaR[s]]*psi33[phi] - psi34[phi]*Sin[3 deltaR[s]]);
  Psi33s[phi_] = Pi*B0[phi]*(Cos[3*deltaR[s]]*psi34[phi] + psi33[phi]*Sin[3 deltaR[s]]);
  Vpp = GaussLegendreQuadrature[VppTemp[phi]*sprimeVMEC[phi], {phi, 0, 2*pi}, 250, 20];
  Print["  Magnetic well: from near-axis = "<>ToString[Vpp]<>If[FileExistsQ[vmecOutput],", from VMEC = "<>ToString[Vwell0VMEC],""]];
  printv = printv <> "Magnetic well near-axis" <> nl;
  printv = printv <> ToString[Vpp] <> nl;
  If[FileExistsQ[vmecOutput],PsiI = vmecPSI[[1]];
    Dwell = GaussLegendreQuadrature[DwellTemp[phi]*sprimeVMEC[phi], {phi, 0, 2*pi}, 250, 20];
    printv = printv <> "Magnetic well VMEC" <> nl;
    printv = printv <> ToString[Vwell0VMEC] <> nl;
  ]
  ,
  vPPDwell=ParallelTable[
    B0R[phi_]    = B0VMEC /. fitParams[iradius];
    deltaR[phi_] = deltaVMEC /. fitParams[iradius];
    muR[phi_]    = muVMEC /. fitParams[iradius];
    pprime0      = pprime0VMEC;
    lambda       = torCurVMEC;
    PsiI         = vmecPSI[[iradius]];
    If[ordern>2,
      Psi31c[phi_] = (B0R[phi]^(3/2))*Sum[Cos[m*vmecNFP*phi + deltaR[phi]]*psic[3, 1, m] /. fitParams[iradius], {m, 0, nModes}];
      Psi31s[phi_] = (B0R[phi]^(3/2))*Sum[Sin[m*vmecNFP*phi + deltaR[phi]]*psic[3, 1, m] /. fitParams[iradius], {m, 0, nModes}];
      Psi33c[phi_] = (B0R[phi]^(3/2))*Sum[Cos[m*vmecNFP*phi + deltaR[phi]]*psic[3, 3, m] /. fitParams[iradius], {m, 0, nModes}];
      Psi33s[phi_] = (B0R[phi]^(3/2))*Sum[Sin[m*vmecNFP*phi + deltaR[phi]]*psic[3, 3, m] /. fitParams[iradius], {m, 0, nModes}];
      ,Psi31c[phi_] = 0;Psi31s[phi_] = 0;Psi33c[phi_] = 0;Psi33s[phi_] = 0;
    ];
    VpptempI=GaussLegendreQuadrature[VppTemp[phi]*sprimeVMEC[phi], {phi, 0, 2*pi}, 200, 15];
    DwelltempI=GaussLegendreQuadrature[(DwellTemp[phi]/.VppP->VpptempI)*sprimeVMEC[phi], {phi, 0, 2*pi}, 200, 15];
    {VpptempI,DwelltempI}
  ,{iradius,1,nsurfaces}];
  Vpp=vPPDwell[[All,1]];
  Dwell=vPPDwell[[All,2]];
  Print["  Magnetic well: from near-axis = "<>ToString[Vpp]<>", from VMEC = "<>ToString[Vwell0VMEC]];
  printv = printv <> "Magnetic well near-axis" <> nl;
  printv = printv <> ToString[Vpp] <> nl;
  printv = printv <> "Magnetic well VMEC" <> nl;
  printv = printv <> ToString[Vwell0VMEC] <> nl;
  Print["  Dwell: from near-axis = "<>ToString[Dwell]<>", from VMEC = "<>ToString[DWellVMEC]];
  printv = printv <> "Dwell" <> nl;
  printv = printv <> ToString[Dwell] <> nl;
  printv = printv <> "Dwell VMEC" <> nl;
  printv = printv <> ToString[DWellVMEC] <> nl;
];

fitFile = OpenAppend[StringJoin[currentPath,"data/",ToString[proj],"/",ToString[proj],"_fit_parameters_ordern_"<>ToString[ordern]<>".txt"]];
WriteString[fitFile, printv];Close[fitFile];

(*Plot Obtained Surfaces*)
If[plotFit==1 || plotOriginal==1,
  timestart1 = AbsoluteTime[]; WriteString[$Output, "  Creating Figures...\n"];
  Get[StringJoin[currentPath,"src/plot_mercier.wls"]];
  WriteString[$Output, "  Successful in "<>ToString[AbsoluteTime[]-timestart1]<>"s\n"];
];

(*Output VMEC file*)
If[outputToVMEC == 1,
   timestart1 = AbsoluteTime[]; WriteString[$Output, "  Mercier to VMEC...\n"]
   Get[StringJoin[currentPath,"src/mercier_to_vmec.wls"]];
   WriteString[$Output, "  Successful in "<>ToString[AbsoluteTime[]-timestart1]<>"s\n"];
];

(*Finish Evaluation*)
WriteString[$Output,"  SENAC Mathematica took "<>ToString[AbsoluteTime[]-overallTimeStart]<>"s\n"];
CloseKernels[];
Exit[];