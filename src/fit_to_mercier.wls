(* ::Package:: *)

deltac0           = ToExpression[$ScriptCommandLine[[10]]];
deltal0           = ToExpression[$ScriptCommandLine[[11]]];
deltalmin         = ToExpression[$ScriptCommandLine[[12]]];
deltalmax         = ToExpression[$ScriptCommandLine[[13]]];
If[deltalmin==deltalmax, deltalmin=-1.2*vmecNFP;deltalmax=1.2*vmecNFP];
muc0              = ToExpression[$ScriptCommandLine[[14]]];
mucMin            = ToExpression[$ScriptCommandLine[[15]]];
mucMax            = ToExpression[$ScriptCommandLine[[16]]];
nModes            = ToExpression[$ScriptCommandLine[[17]]];
maxiterations     = ToExpression[$ScriptCommandLine[[18]]];
If[readInput!=2,B0Est=1];

(*Compute Mercier Angles Theta and Phi as a function VMEC's theta and phi*)
WriteString[$Output, "    Axis's Phi Evaluation...  "];timestart = AbsoluteTime[];
phiAxis = Interpolation[Flatten[ParallelTable[{{theta, phi}, x /. FindRoot[Dot[b0VMEC[x], FluxSurfaceVMEC[theta, phi] - closedcurvVMEC[x]], {x, phi}]}, {theta, 0, 2*pi, 2*pi/nthetaM}, {phi, 0, 2*pi, 2*pi/nphiM}], 1], Method -> "Spline", InterpolationOrder -> 2];
(**phiAxis[theta_,phi_]=phi;*)
ycomponent[theta_, phi_] = Dot[FluxSurfaceVMEC[theta, phi] - closedcurvVMEC[phiAxis[theta,phi]], t0VMEC[phiAxis[theta,phi]]];
xcomponent[theta_, phi_] = Dot[FluxSurfaceVMEC[theta, phi] - closedcurvVMEC[phiAxis[theta,phi]], k0VMEC[phiAxis[theta,phi]]];
thetaMercierFunc[theta_, phi_] = Chop[ArcTan[xcomponent[theta, phi], ycomponent[theta, phi]], 10^-8];
Print["Successful in "<>ToString[AbsoluteTime[]-timestart]<>"s"];

(*Compute Mercier Angle*)
WriteString[$Output, "    Mercier's Angle Evaluation...  "];timestart = AbsoluteTime[];
thetaMercierTable = ParallelTable[thetaMercierFunc[theta, phi], {theta, 0, 2*pi,2*pi/nthetaM}, {phi, 0, 2*pi, 2*pi/nphiM}]//Quiet;
thetaMercier = ListInterpolation[thetaMercierTable[[All, All]], {{0, 2*pi}, {0, 2*pi}}, Method -> "Spline", InterpolationOrder -> 2];
Print["Successful in "<>ToString[AbsoluteTime[]-timestart]<>"s"];Clear[thetaMercierTable];

rhoSurf[theta_, phi_] = FluxSurfaceVMEC[theta, phi] - closedcurvVMEC[phiAxis[theta,phi]];
rhonVMEC2[theta_, phi_] = Chop[Dot[rhoSurf[theta, phi], k0VMEC[phiAxis[theta,phi]]]^2 + Dot[rhoSurf[theta, phi], t0VMEC[phiAxis[theta,phi]]]^2, 10^-6];
FluxSurfaceVMECn[theta_, phi_] = Chop[closedcurvVMEC[phiAxis[theta,phi]] + Sqrt[rhonVMEC2[theta, phi]]*(Cos[thetaMercier[theta, phi]]*k0VMEC[phiAxis[theta,phi]] + Sin[thetaMercier[theta, phi]]*t0VMEC[phiAxis[theta,phi]])];

(*Flux Surface After Fitting for mu and delta*)
WriteString[$Output, "    Mercier's Rho Evaluation...  "];timestart = AbsoluteTime[];
dataVMEC = Flatten[ParallelTable[{theta, phi, Re[Sqrt[rhonVMEC2[theta, phi]]]}, {theta, 0, 2*pi, 2*pi/nthetaM}, {phi, 0, 2*pi, 2*pi/nphiM}], 1];
Print["Successful in "<>ToString[AbsoluteTime[]-timestart]<>"s"];

(*Mercier rho function in terms of axis quantities a.k.a. Mercier's theta and phi*)
muVMEC = muc[0] + Sum[muc[i] Cos[vmecNFP i phi], {i, 1, nModes}];
deltaVMEC = deltal*phi+deltac[0]+Sum[deltas[i] Sin[vmecNFP i phi],{i,1,nModes}];
B0VMEC = B0c[0]+ Sum[B0c[i] Cos[vmecNFP i phi],{i,1,nModes}];
rho[theta_, phi_] = Sqrt[(vmecPSI*Sqrt[1 - muVMEC^2] )/(2*(B0VMEC)*(1 + muVMEC*Cos[2*(theta + deltaVMEC)]))];

WriteString[$Output, "    Computing Fit to Surface...  "];timestart = AbsoluteTime[];
Off[InterpolatingFunction::dmval]; (*Hack not to specify that the arguments should vary between 0 and 2 pi *)
modelVMEC = rho[thetaMercier[theta, phi], phiAxis[theta,phi]];
nlm = NonlinearModelFit[dataVMEC, {modelVMEC, mucMin < Sum[muc[i], {i, 0, nModes}] < mucMax, B0c[0] > 0., -pi <= deltac[0] <= pi,  deltalmin <= deltal <= deltalmax}, 
   Flatten[{{{deltal, deltal0}, {B0c[0], B0Est}, {muc[0], muc0}, {deltac[0], deltac0}}, 
     Flatten[Table[{{B0c[i],0.01},{muc[i], 0.01},{deltas[i],0.01}}, {i, 1, nModes}], 1]}, 1],
   {theta,phi}, Method -> {"NMinimize"}, MaxIterations -> maxiterations
];
Clear[dataVMEC];

fitParams = nlm["BestFitParameters"];
WriteString[$Output, "Successful in "<>ToString[AbsoluteTime[]-timestart]<>"s\n"];
WriteString[$Output, "    Best fit parameters "<>ToString[fitParams]<>"\n"];
Export[StringJoin[currentPath,"data/",ToString[proj],"/",ToString[proj],"_fit_parameters.txt"],fitParams];

rho[theta_, phi_] = rho[theta, phi] /. fitParams;
FluxSurfacenFitVMEC[theta_, phi_] = Chop[closedcurvVMEC[phi] + (rho[theta, phi])*(Cos[theta]*k0VMEC[phi] + Sin[theta]*t0VMEC[phi]), 10^-6]//Quiet;

(*The inverse transform phi_axis -> phi_boundary if needed*)
(*phiBoundary = Interpolation[Flatten[ParallelTable[{{theta, phia}, phi /. FindRoot[phiAxis[theta, phi] - phia, {phi, phia, 0, 2*pi}]}, {theta, 0, 2*pi, 2*pi/nthetaM}, {phia, 0, 2*pi, 2*pi/nphiM}], 1]];*)
RBCFit[phi_, theta_] = Re[Cos[phi]*FluxSurfacenFitVMEC[theta, phiAxis[theta,phi]][[1]] + Sin[phi]*FluxSurfacenFitVMEC[theta, phiAxis[theta,phi]][[2]]]//Quiet;
ZBSFit[phi_, theta_] = Re[FluxSurfacenFitVMEC[theta, phiAxis[theta,phi]][[3]]]//Quiet;