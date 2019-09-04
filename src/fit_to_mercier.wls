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

ycomponent[theta_, phi_] = Dot[FluxSurfaceVMEC[theta, phi] - closedcurvVMEC[phi], t0VMEC[phi]];
xcomponent[theta_, phi_] = Dot[FluxSurfaceVMEC[theta, phi] - closedcurvVMEC[phi], k0VMEC[phi]];
thetaMercierFunc[theta_, phi_] = Chop[ArcTan[xcomponent[theta, phi], ycomponent[theta, phi]], 10^-6];

(*Compute Mercier Angle*)
WriteString[$Output, "    Mercier's Angle Evaluation...  "];timestart = AbsoluteTime[];
thetaMercierTable = ParallelTable[thetaMercierFunc[theta, s], {theta, 0, 2*pi,2*pi/nthetaM}, {s, 0, 2*pi, 2*pi/nphiM}]//Quiet;
thetaMercier = ListInterpolation[thetaMercierTable[[All, All]], {{0, 2*pi}, {0, 2*pi}}];
Print["Successful in "<>ToString[AbsoluteTime[]-timestart]<>"s"];
Clear[thetaMercierTable];

rhoSurf[theta_, phi_] = FluxSurfaceVMEC[theta, phi] - closedcurvVMEC[phi];
rhonVMEC2[theta_, phi_] = Chop[Dot[rhoSurf[theta, phi], k0VMEC[phi]]^2 + Dot[rhoSurf[theta, phi], t0VMEC[phi]]^2, 10^-6];
FluxSurfaceVMECn[theta_, s_] = Chop[closedcurvVMEC[s] + Sqrt[rhonVMEC2[theta, s]] (Cos[thetaMercier[ theta, s]] k0VMEC[s] + Sin[thetaMercier[ theta, s]] t0VMEC[s])];

(*Flux Surface After Fitting for mu and delta*)
WriteString[$Output, "    Mercier's Rho Evaluation...  "];timestart = AbsoluteTime[];
dataVMEC = Flatten[ParallelTable[{theta, phi, Re[Sqrt[rhonVMEC2[theta, phi]]]}, {theta, 0, 2*pi, 2*pi/nthetaM}, {phi, 0, 2*pi, 2*pi/nphiM}], 1];
Print["Successful in "<>ToString[AbsoluteTime[]-timestart]<>"s"];

muVMEC = muc[0] + Sum[muc[i] Cos[vmecNFP i phi]+ mus[i] Sin[vmecNFP i phi], {i, 1, nModes}];
deltaVMEC = deltal*phi + deltac[0]+Sum[deltac[i] Cos[vmecNFP i phi]+deltas[i] Sin[vmecNFP i phi],{i,1,nModes}];
B0VMEC = B0c[0]+ Sum[B0c[i] Cos[vmecNFP i phi]+B0s[i] Sin[vmecNFP i phi],{i,1,nModes}];
rho[theta_, phi_] = Sqrt[(vmecPSI*Sqrt[1 - muVMEC^2] )/(2*(B0VMEC)*(1 + muVMEC*Cos[2 (theta + deltaVMEC)]))];
modelVMEC = rho[thetaMercier[theta, phi], phi];

WriteString[$Output, "    Computing Fit to Surface...  "];timestart = AbsoluteTime[];
nlm = NonlinearModelFit[dataVMEC, {modelVMEC, mucMin < Sum[muc[i], {i, 0, nModes}] < mucMax, B0c[0] > 0., -pi <= deltac[0] <= pi,  deltalmin <= deltal <= deltalmax}, 
   Flatten[{{{deltal, deltal0}, {B0c[0], B0Est}, {muc[0], muc0}, {deltac[0], deltac0}}, 
     Flatten[Table[{{B0c[i],0.01},{B0s[i],0.01},{muc[i], 0.01}, {mus[i], 0.0001},{deltac[i],0.01},{deltas[i],0.01}}, {i, 1, nModes}], 1]}, 1],
   {theta,phi}, Method -> {"NMinimize"}, MaxIterations -> maxiterations
];
Clear[dataVMEC];

fitParams = nlm["BestFitParameters"];
WriteString[$Output, "Successful in "<>ToString[AbsoluteTime[]-timestart]<>"s\n"];
WriteString[$Output, "    Best fit parameters "<>ToString[fitParams]<>"\n"];
Export[StringJoin[currentPath,"data/",ToString[proj],"/",ToString[proj],"_fit_parameters.txt"],fitParams];

FluxSurfacenFitVMEC[theta_, s_] = Chop[closedcurvVMEC[s] + (rho[theta, s] /. fitParams) (Cos[theta] k0VMEC[s] + Sin[theta] t0VMEC[s]), 10^-6]//Quiet;
RBCFit[phi_, theta_] = Re[Cos[phi] FluxSurfacenFitVMEC[theta, phi][[1]] + Sin[phi] FluxSurfacenFitVMEC[theta, phi][[2]]]//Quiet;
ZBSFit[phi_, theta_] = Re[FluxSurfacenFitVMEC[theta, phi][[3]]]//Quiet;

If[plotFit==1 && readInput==2,
    timestart = AbsoluteTime[];
    Off[ParametricPlot::precw,NMinimize::precw,NMaximize::precw];
    rangePlot = {
        {0.9*NMinimize[{RBCVMEC[phi, theta], 0 < phi <= 2*pi, 0 < theta <= 2*pi}, {phi, theta}][[1]],
        1.1*NMaximize[{RBCVMEC[phi, theta], 0 < phi <= 2*pi, 0 < theta <= 2*pi}, {phi, theta}][[1]]},
        {0.9*NMinimize[{ZBSVMEC[phi, theta], 0 < phi <= 2*pi,0 < theta <= 2*pi}, {phi, theta}][[1]],
        1.1*NMaximize[{ZBSVMEC[phi, theta], 0 < phi <= 2*pi, 0 < theta <= 2*pi}, {phi, theta}][[1]]}
    }//Quiet;
    nPlots = 8;
    fig=Grid[{ParallelTable[ParametricPlot[{RBCFit[phi, theta], ZBSFit[phi, theta]}, {theta, 0, 2*pi}, AspectRatio -> 1, PlotRange -> rangePlot, PlotPoints -> 15, MaxRecursion -> 1, ImageSize->Large], {phi, 0, 2*pi, 2*pi/nPlots}], 
        ParallelTable[ParametricPlot[{RBCVMEC[theta, phi], ZBSVMEC[theta, phi]}, {theta, 0, 2*pi}, AspectRatio -> 1, PlotRange -> rangePlot, PlotPoints -> 15, MaxRecursion -> 1, ImageSize->Large], {phi, 0, 2*pi, 2*pi/nPlots}]
    }]//Quiet;
    Export[StringJoin[currentPath,"data/",ToString[proj],"/Figures/",ToString[proj],"_fit_plotcompare.pdf"],fig, "AllowRasterization" -> True, ImageSize -> 1000, ImageResolution -> 600]
    Print["    Fit plot took "<>ToString[AbsoluteTime[]-timestart]<>"s"];
    Clear[fig,posvecpsiTable,figVMEC,figFit,Nposvecpsitemp];
]