(* ::Package:: *)

factorPhi=1.0;
(*Initialize Variables*)
timestartPhiAxis=0;timestartMercierAngle=0;timestartrhoVMEC=0;timestartFit=0;
Table[
  WriteString[$Output, "    Fiting Surface "<>ToString[iradius]<>" of "<>ToString[nsurfaces]<>"...  "];timestartFitIteration=AbsoluteTime[];
  (*Compute Mercier Angles Theta and Phi as a function VMEC's theta and phi*)
  timestart=AbsoluteTime[];
  phiAxis = Interpolation[Flatten[ParallelTable[{{theta, phi}, x /. FindRoot[Dot[b0VMEC[x], FluxSurfaceVMEC[theta, phi, iradius] - closedcurvVMEC[x]], {x, phi}]}, {theta, 0, 2*pi, 2*pi/nthetaM/1.2}, {phi, 0, factorPhi*2*pi/vmecNFP, factorPhi*2*pi/nphiM/vmecNFP/1.4}], 1], Method -> "Spline", InterpolationOrder -> 3];
  timestartPhiAxis=timestartPhiAxis+AbsoluteTime[]-timestart;

  (*Compute Mercier Angle*)
  timestart = AbsoluteTime[];
  ycomponent[theta_, phi_] = Dot[FluxSurfaceVMEC[theta, phi, iradius] - closedcurvVMEC[phiAxis[theta,phi]], t0VMEC[phiAxis[theta,phi]]];
  xcomponent[theta_, phi_] = Dot[FluxSurfaceVMEC[theta, phi, iradius] - closedcurvVMEC[phiAxis[theta,phi]], k0VMEC[phiAxis[theta,phi]]];
  thetaMercierFunc[theta_, phi_] = Chop[ArcTan[xcomponent[theta, phi], ycomponent[theta, phi]], 10^-5];
  timestartMercierAngle=timestartMercierAngle+AbsoluteTime[]-timestart;Clear[thetaMercierTable];
  If[iradius>1 && keepfit==1,
    fitParams[iradius] = fitParams[iradius-1];
  ,
    (*Compute Rho from VMEC*)
    rhoSurf[theta_, phi_] = FluxSurfaceVMEC[theta, phi, iradius] - closedcurvVMEC[phiAxis[theta,phi]];
    rhonVMEC2[theta_, phi_] = Chop[Dot[rhoSurf[theta, phi], k0VMEC[phiAxis[theta,phi]]]^2 + Dot[rhoSurf[theta, phi], t0VMEC[phiAxis[theta,phi]]]^2, 10^-4];

    (*Flux Surface After Fitting for mu and delta*)
    timestart = AbsoluteTime[];
    dataVMEC = Flatten[ParallelTable[{theta, phi, Re[Sqrt[rhonVMEC2[theta, phi]]]}, {theta, 0, 2*pi, 2*pi/nthetaM}, {phi, 0, factorPhi*2*pi/vmecNFP, factorPhi*2*pi/nphiM/vmecNFP}], 1];
    timestartrhoVMEC=timestartrhoVMEC+AbsoluteTime[]-timestart;

    (*Performing the fit*)
    timestart = AbsoluteTime[];
    Off[InterpolatingFunction::dmval]; (*Hack not to specify that the arguments should vary between 0 and 2 pi *)
    modelVMEC = Re[rho[thetaMercierFunc[theta, phi], phiAxis[theta,phi]]/.PSI->vmecPSI[[iradius]]];
    If[!FileExistsQ[vmecOutput],B0Est=1];If[deltalmin==deltalmax, deltalmin=-1.2*vmecNFP;deltalmax=1.2*vmecNFP];
    If[iradius>1,
      deltalmin = Round[deltal/.fitParams[iradius-1]]-0.01*Abs[Round[deltal/.fitParams[iradius-1]]];
      deltalmax = Round[deltal/.fitParams[iradius-1]]+0.01*Abs[Round[deltal/.fitParams[iradius-1]]];
      deltal0   = Round[deltal/.fitParams[iradius-1]];
      muc0      = muc[0]/.fitParams[iradius-1];
      B0Est     = B0c[0]/.fitParams[iradius-1];
      deltac0   = deltac[0]/.fitParams[iradius-1];
    ];

    nlm = NonlinearModelFit[dataVMEC,
      DeleteCases[Flatten[
        {modelVMEC, mucMin < (muVMEC/.phi->0) < mucMax, mucMin < (muVMEC/.phi->pi) < mucMax, 0.01<B0c[0], 0.01 < (B0VMEC/.phi->0), 0.01 < (B0VMEC/.phi->pi), -pi <= deltac[0] <= pi,  deltalmin <= deltal <= deltalmax, 
        Flatten[
          Table[{(*{mucMin < muc[i] < mucMax}, {-B0Est < B0c[i] < B0Est},*) {-pi < deltas[i] < pi}},{i,1,nModes}]
        ](*,(*Limit the values of psic and psis*)
        If[ordern > 2, Flatten[
          Table[Table[If[EvenQ[n], {{-10<psic[n, 2 i]<10}, {-10<psis[n, 2 i]<10}}, {{-10<psic[n, 2 i + 1]<10}, {-10<psis[n, 2 i + 1]<10}}], {i,0, Floor[n/2]}], {n, 3, ordern}]
        ],xa]*)
        }
      ,1],xa],
      DeleteCases[Flatten[{
        {{deltal, deltal0}, {B0c[0], B0Est}, {muc[0], muc0}, {deltac[0], deltac0}}, 
        Flatten[Table[{{B0c[i],0.01*B0Est},{muc[i],0.01*Min[Abs[mucMin],Abs[mucMax]]},{deltas[i],0.01}}, {i, 1, nModes}], 1],
        If[ordern > 2, Flatten[Table[Table[Table[If[EvenQ[n], {{psic[n, 2 i, j], 0.1}, {psis[n, 2 i, j], 0.1}}, {{psic[n, 2 i + 1, j], 0.1}, {psis[n, 2 i + 1, j], 0.1}}], {j,0,2*nModes}], {i,0, Floor[n/2]}], {n, 3, ordern}], 3],xa]
      }, 1],xa],
    {theta,phi}, Method -> {"NMinimize"}, MaxIterations -> maxiterations, AccuracyGoal -> 5, PrecisionGoal -> 5
    ];
    fitParams[iradius] = nlm["BestFitParameters"];
    WriteString[$Output, "    Best fit parameters of surface "<>ToString[iradius]<>" - "<>ToString[fitParams[iradius]]<>"\n"];
    timestartFit=timestartFit+AbsoluteTime[]-timestart;
  ];
  (*The inverse transform phi_axis -> phi_boundary if needed*)
  (*phiBoundary = Interpolation[Flatten[ParallelTable[{{theta, phia}, phi /. FindRoot[phiAxis[theta, phi] - phia, {phi, phia, 0, 2*pi}]}, {theta, 0, 2*pi, 2*pi/nthetaM}, {phia, 0, 2*pi, 2*pi/nphiM}], 1]];*)

  (*Construct Obtained Surface*)
  FluxSurfacenFitVMEC[theta_, phi_, iradius] = Chop[closedcurvVMEC[phi] + (rho[theta, phi] /. fitParams[iradius] /.PSI->vmecPSI[[iradius]])*(Cos[theta]*k0VMEC[phi] + Sin[theta]*t0VMEC[phi]), 10^-5]//Quiet; (*Mercier's Coordinates*)
  RBCFit[phi_, theta_, iradius] = Re[Cos[phi]*FluxSurfacenFitVMEC[thetaMercierFunc[theta,phi], phiAxis[theta,phi], iradius][[1]] + Sin[phi]*FluxSurfacenFitVMEC[thetaMercierFunc[theta,phi], phiAxis[theta,phi], iradius][[2]]]//Quiet; (*VMEC's Coordinates*)
  ZBSFit[phi_, theta_, iradius] = Re[FluxSurfacenFitVMEC[thetaMercierFunc[theta,phi], phiAxis[theta,phi], iradius][[3]]]//Quiet; (*VMEC's Coordinates*)
  WriteString[$Output, "Successful in "<>ToString[AbsoluteTime[]-timestartFitIteration]<>"s\n"];

  iotaFitSurf = GaussLegendreQuadrature[(-D[deltaVMEC /. fitParams[iradius], phi] + torsVMEC[phi]*sprimeVMEC[phi])*Sqrt[1 - (muVMEC /. fitParams[iradius])^2]+D[deltaVMEC /. fitParams[iradius], phi], {phi, 0, 2*pi}, 250, 5]/(2*pi)-nNormal;
  If[FileExistsQ[vmecOutput],WriteString[$Output, "Iota From VMEC "<>ToString[iotaAxisOut[[iradius]]]<>", "]];
  WriteString[$Output, "Iota from Near-Axis:  "<>ToString[iotaFitSurf]<>"\n"];

,{iradius,1,nsurfaces}];Clear[dataVMEC];

WriteString[$Output, "    Total time computing phiAxis "<>ToString[timestartPhiAxis]<>"s\n"];
WriteString[$Output, "    Total time computing rhoVMEC data "<>ToString[timestartrhoVMEC]<>"s\n"];
WriteString[$Output, "    Total time computing fits "<>ToString[timestartFit]<>"s\n"];

printv = "";nl = "\n";Do[printv = printv <> ToString[fitParams[i]] <> nl, {i, 1, nsurfaces}];
Export[StringJoin[currentPath,"data/",ToString[proj],"/",ToString[proj],"_fit_parameters_ordern_"<>ToString[ordern]<>".txt"],printv];
(*Save[StringJoin[currentPath,"data/",ToString[proj],"/",ToString[proj],"_fit_functions.mx"], {RBCFit,ZBSFit,FluxSurfacenFitVMEC,nsurfaces,fitParams}];*)