(* ::Package:: *)

(*Initialize Variables*)
timestartPhiAxis=0;timestartMercierAngle=0;timestartrhoVMEC=0;timestartFit=0;

Table[
  WriteString[$Output, "    Fiting Surface "<>ToString[iradius]<>" of "<>ToString[nsurfaces]<>"...  "];timestartFitIteration=AbsoluteTime[];
  
  (*Compute Mercier Angles Theta and Phi as a function VMEC's theta and phi*)
  timestart=AbsoluteTime[];
  phiAxis = Interpolation[Flatten[ParallelTable[{{theta, phi}, x /. FindRoot[Dot[b0VMEC[x], FluxSurfaceVMEC[theta, phi, iradius] - closedcurvVMEC[x]], {x, phi}]}, {theta, 0, 2*pi, 2*pi/nthetaM}, {phi, 0, 2*pi, 2*pi/nphiM}], 1], Method -> "Spline", InterpolationOrder -> 2];
  timestartPhiAxis=timestartPhiAxis+AbsoluteTime[]-timestart;

  (*Compute Mercier Angle*)
  timestart = AbsoluteTime[];
  ycomponent[theta_, phi_] = Dot[FluxSurfaceVMEC[theta, phi, iradius] - closedcurvVMEC[phiAxis[theta,phi]], t0VMEC[phiAxis[theta,phi]]];
  xcomponent[theta_, phi_] = Dot[FluxSurfaceVMEC[theta, phi, iradius] - closedcurvVMEC[phiAxis[theta,phi]], k0VMEC[phiAxis[theta,phi]]];
  thetaMercierFunc[theta_, phi_] = Chop[ArcTan[xcomponent[theta, phi], ycomponent[theta, phi]], 10^-8];
  thetaMercierTable = ParallelTable[thetaMercierFunc[theta, phi], {theta, 0, 2*pi,2*pi/nthetaM}, {phi, 0, 2*pi, 2*pi/nphiM}]//Quiet;
  thetaMercier = ListInterpolation[thetaMercierTable[[All, All]], {{0, 2*pi}, {0, 2*pi}}, Method -> "Spline", InterpolationOrder -> 2];
  timestartMercierAngle=timestartMercierAngle+AbsoluteTime[]-timestart;Clear[thetaMercierTable];

  (*Compute Rho from VMEC*)
  rhoSurf[theta_, phi_] = FluxSurfaceVMEC[theta, phi, iradius] - closedcurvVMEC[phiAxis[theta,phi]];
  rhonVMEC2[theta_, phi_] = Chop[Dot[rhoSurf[theta, phi], k0VMEC[phiAxis[theta,phi]]]^2 + Dot[rhoSurf[theta, phi], t0VMEC[phiAxis[theta,phi]]]^2, 10^-6];
  FluxSurfaceVMECn[theta_, phi_] = Chop[closedcurvVMEC[phiAxis[theta,phi]] + Sqrt[rhonVMEC2[theta, phi]]*(Cos[thetaMercier[theta, phi]]*k0VMEC[phiAxis[theta,phi]] + Sin[thetaMercier[theta, phi]]*t0VMEC[phiAxis[theta,phi]])];

  (*Flux Surface After Fitting for mu and delta*)
  timestart = AbsoluteTime[];
  dataVMEC = Flatten[ParallelTable[{theta, phi, Re[Sqrt[rhonVMEC2[theta, phi]]]}, {theta, 0, 2*pi, 2*pi/nthetaM}, {phi, 0, 1.4*2*pi/vmecNFP, 2*pi/nphiM}], 1];
  timestartrhoVMEC=timestartrhoVMEC+AbsoluteTime[]-timestart;

  (*Performing the fit*)
  timestart = AbsoluteTime[];
  Off[InterpolatingFunction::dmval]; (*Hack not to specify that the arguments should vary between 0 and 2 pi *)
  modelVMEC = rho[thetaMercier[theta, phi], phiAxis[theta,phi]]/.PSI->vmecPSI[[iradius]];
  If[!FileExistsQ[vmecOutput],B0Est=1];If[deltalmin==deltalmax, deltalmin=-1.2*vmecNFP;deltalmax=1.2*vmecNFP];
  If[iradius>1,
    If[iradius>2,
      deltalmin = Round[deltal/.fitParams[iradius-1]]-0.01*Abs[Round[deltal/.fitParams[iradius-1]]];
      deltalmax = Round[deltal/.fitParams[iradius-1]]+0.01*Abs[Round[deltal/.fitParams[iradius-1]]];
    ];
    deltal0   = Round[deltal/.fitParams[iradius-1]];
    muc0      = muc[0]/.fitParams[iradius-1];
    B0Est     = B0c[0]/.fitParams[iradius-1];
    deltac0   = deltac[0]/.fitParams[iradius-1];
  ];

  nlm = NonlinearModelFit[dataVMEC,
    Flatten[
      {modelVMEC, mucMin < muc[0] < mucMax, 0.05*B0Est<B0c[0], -pi <= deltac[0] <= pi,  deltalmin <= deltal <= deltalmax, 
      Flatten[
        Table[{{-0.1*Min[Abs[mucMin],Abs[mucMax]] < muc[i] < 0.1*Min[Abs[mucMin],Abs[mucMax]]}, {-0.01*B0Est < B0c[i] < 0.01*B0Est}, {-0.7 < deltas[i] < 0.7}},{i,1,nModes}]
      ],
      If[ordern > 2, Flatten[
        Table[Table[If[EvenQ[n], {{-1<psic[n, 2 i]<1}, {-1<psis[n, 2 i]<1}}, {{-1<psic[n, 2 i + 1]<1}, {-1<psis[n, 2 i + 1]<1}}], {i,0, Floor[n/2]}], {n, 3, ordern}]
      ],xa]
      }
    ,1],
    DeleteCases[Flatten[{
      {{deltal, deltal0}, {B0c[0], B0Est}, {muc[0], muc0}, {deltac[0], deltac0}}, 
      Flatten[Table[{{B0c[i],0.001*B0Est},{muc[i],0.001*Min[Abs[mucMin],Abs[mucMax]]},{deltas[i],0.01}}, {i, 1, nModes}], 1],
      If[ordern > 2, Flatten[Table[Table[If[EvenQ[n], {{psic[n, 2 i], 0.1}, {psis[n, 2 i], 0.1}}, {{psic[n, 2 i + 1], 0.1}, {psis[n, 2 i + 1], 0.1}}], {i,0, Floor[n/2]}], {n, 3, ordern}], 2],xa]
    }, 1],xa],
   {theta,phi}, Method -> {"NMinimize"}, MaxIterations -> maxiterations
  ];
  fitParams[iradius] = nlm["BestFitParameters"];
  timestartFit=timestartFit+AbsoluteTime[]-timestart;

  (*The inverse transform phi_axis -> phi_boundary if needed*)
  (*phiBoundary = Interpolation[Flatten[ParallelTable[{{theta, phia}, phi /. FindRoot[phiAxis[theta, phi] - phia, {phi, phia, 0, 2*pi}]}, {theta, 0, 2*pi, 2*pi/nthetaM}, {phia, 0, 2*pi, 2*pi/nphiM}], 1]];*)

  (*Construct Obtained Surface*)
  FluxSurfacenFitVMEC[theta_, phi_, iradius] = Chop[closedcurvVMEC[phi] + (rho[theta, phi] /. fitParams[iradius] /.PSI->vmecPSI[[iradius]])*(Cos[theta]*k0VMEC[phi] + Sin[theta]*t0VMEC[phi]), 10^-7]//Quiet; (*Mercier's Coordinates*)
  RBCFit[phi_, theta_, iradius] = Re[Cos[phi]*FluxSurfacenFitVMEC[thetaMercierFunc[theta,phi], phiAxis[theta,phi], iradius][[1]] + Sin[phi]*FluxSurfacenFitVMEC[thetaMercierFunc[theta,phi], phiAxis[theta,phi], iradius][[2]]]//Quiet; (*VMEC's Coordinates*)
  ZBSFit[phi_, theta_, iradius] = Re[FluxSurfacenFitVMEC[thetaMercierFunc[theta,phi], phiAxis[theta,phi], iradius][[3]]]//Quiet; (*VMEC's Coordinates*)
  WriteString[$Output, "Successful in "<>ToString[AbsoluteTime[]-timestartFitIteration]<>"s\n"];
,{iradius,1,nsurfaces}];Clear[dataVMEC];

WriteString[$Output, "    Best fit parameters of first surface "<>ToString[fitParams[1]]<>"\n"];
WriteString[$Output, "    Best fit parameters of last  surface "<>ToString[fitParams[nsurfaces]]<>"\n"];
WriteString[$Output, "    Total time computing phiAxis "<>ToString[timestartPhiAxis]<>"s\n"];
WriteString[$Output, "    Total time computing thetaMercier "<>ToString[timestartMercierAngle]<>"s\n"];
WriteString[$Output, "    Total time computing rhoVMEC data "<>ToString[timestartrhoVMEC]<>"s\n"];
WriteString[$Output, "    Total time computing fits "<>ToString[timestartFit]<>"s\n"];

printv = "";nl = "\n";Do[printv = printv <> ToString[fitParams[i]] <> nl, {i, 1, nsurfaces}];
Export[StringJoin[currentPath,"data/",ToString[proj],"/",ToString[proj],"_fit_parameters.txt"],printv];
(*Save[StringJoin[currentPath,"data/",ToString[proj],"/",ToString[proj],"_fit_functions.mx"], {RBCFit,ZBSFit,FluxSurfacenFitVMEC,nsurfaces,fitParams}];*)