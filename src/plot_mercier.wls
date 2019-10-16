(* ::Package:: *)

(*PLOT MU AND DELTA AS A FUNCTION OF PHI*)

(*Plotting Options*)
If[muDefined==1,nfigsSurf=1];
nfigsSurf=Min[nfigsSurf,nsurfaces];
PlotOptions = Flatten[{LabelStyle -> Directive[FontFamily -> "Latin Modern Roman", FontSize -> 30], ImageSize -> ImageSizePlot, Frame -> True}];
phimax = Table[phi + 2*pi (1 - (nfigsSurf - 1)/nfigsSurf/vmecNFP), {phi, 2*pi/vmecNFP (1 - 1/nfigsSurf), 0, -2*pi/nfigsSurf/vmecNFP/2.5}]-0.0;
phimin = Table[phi, {phi, 0, 2*pi/vmecNFP (1 - 1/nfigsSurf), 2*pi/nfigsSurf/vmecNFP/2.5}]+0.0;

(*Plot Original Surface*)
If[muDefined==0 && plotOriginal==1,
  figVMEC=Array[Null,nfigsSurf];bvmecFig=Array[Null,nfigsSurf];
  timestart = AbsoluteTime[]; WriteString[$Output, "    Plotting Original Surface...\n"];
  j=0;Table[j=j+1;
    (*timestart2 = AbsoluteTime[]; WriteString[$Output, "      Creating Grid for Figure "<>ToString[j]<>" of "<>ToString[nfigsSurf],"..."];
    VMECposvecpsiTable = Table[ParallelTable[FluxSurfaceVMEC[theta, phi, Floor[iradius]][[i]], {theta, 0, 2*pi, 2*pi/(nPlotTheta)}, {phi, phimin[[j]], phimax[[j]], phimax[[j]]/(nPlotPhi)}], {i, 1, 3}];
    VMECNposvecpsitemp = ParallelTable[ListInterpolation[VMECposvecpsiTable[[i]], {{0, 2*pi}, {phimin[[j]], phimax[[j]]}}], {i, 1, 3}];
    VMECNposvecpsi[theta_, phi_] = Table[VMECNposvecpsitemp[[i]][theta, phi], {i, 1, 3}];
    WriteString[$Output, " Successful in "<>ToString[AbsoluteTime[]-timestart2]<>"s\n"];*)
    If[FileExistsQ[vmecOutput],
      If[exportBFieldSurface==1,
        timestart2 = AbsoluteTime[]; WriteString[$Output, "      Exporting magnetic field at surface..."];
        bvmecFig[[j]]=DensityPlot[BVMEC[theta, phi, Floor[iradius]], {phi, 0, 2*pi/vmecNFP}, {theta, 0, 2*pi}, FrameLabel -> {"Toroidal Angle \[CapitalPhi]", HoldForm["Poloidal Angle " Subscript[\[Theta], V]]}, Evaluate@PlotOptions, MaxRecursion -> 1, PlotPoints -> {nthetapointsBsurface,ntphipointsBsurface}, ColorFunction -> "Rainbow", PlotLegends -> BarLegend[Automatic, LegendMarkerSize -> 750, LegendLabel -> HoldForm["B(" Subscript[\[Theta], V] ",\[CapitalPhi])"]]];
        Export[StringJoin[currentPath,"data/",ToString[proj],"/Figures/",ToString[proj],"_VMEC_magnetic_field_surface_"<>ToString[Floor[iradius]]<>"_.pdf"], bvmecFig[[j]], ImageSize->ImageSizePlot, ImageResolution->ImageResolutionPlot];
        WriteString[$Output, " Successful in "<>ToString[AbsoluteTime[]-timestart2]<>"s\n"];
      ];
      timestart2 = AbsoluteTime[]; WriteString[$Output, "      Computing Maximum and Minimum B..."];
      minB0VMEC[j] = FindMinValue[{BVMEC[theta,phi,Floor[iradius]], 0 <= theta <= 2*pi && phimin[[j]] <= phi <= phimax[[j]]}, {theta, phi}] // Quiet;
      maxB0VMEC[j] = FindMaxValue[{BVMEC[theta,phi,Floor[iradius]], 0 <= theta <= 2*pi && phimin[[j]] <= phi <= phimax[[j]]}, {theta, phi}] // Quiet;
      WriteString[$Output, " Successful in "<>ToString[AbsoluteTime[]-timestart2]<>"s\n"];
      timestart2 = AbsoluteTime[]; WriteString[$Output, "      Creating Plot..."];
      figVMEC[[j]]  = ParametricPlot3D[(*VMECNposvecpsi[theta, phi]*)FluxSurfaceVMEC[theta, phi, Floor[iradius]], {theta, 0, 2*pi}, {phi, phimin[[j]], phimax[[j]]}, Boxed -> False, Axes -> False, ImageSize -> Large, Mesh -> None, PlotStyle -> Directive[(*Opacity[0.3+0.7*j/nfigsSurf],*) Specularity[White, 20]], PlotPoints -> {30, plotPointsFig}, MaxRecursion -> maxRecursPlot
        , ColorFunctionScaling -> False, ColorFunction -> (ColorData[{"Rainbow", {minB0VMEC[j], maxB0VMEC[j]}}][(BVMEC[#4,#5,Floor[iradius]])] &)];
      WriteString[$Output, " Successful in "<>ToString[AbsoluteTime[]-timestart2]<>"s\n"];
      ,
      timestart2 = AbsoluteTime[]; WriteString[$Output, "      Computing Maximum and Minimum B..."];
      BFunc[theta_, phi_]=1 + (rho[theta, phi]/.fitParams[Floor[iradius]]/.PSI->vmecPSI[[Floor[iradius]]])*Cos[theta]*curvVMEC[phi];
      minB0 = FindMinValue[{BFunc[theta,phi], 0 <= theta <= 2*pi && phimin[[j]] <= phi <= phimax[[j]]}, {theta, phi}] // Quiet;
      maxB0 = FindMaxValue[{BFunc[theta,phi], 0 <= theta <= 2*pi && phimin[[j]] <= phi <= phimax[[j]]}, {theta, phi}] // Quiet;
      WriteString[$Output, " Successful in "<>ToString[AbsoluteTime[]-timestart2]<>"s\n"];
      timestart2 = AbsoluteTime[]; WriteString[$Output, "      Creating Plot..."];
      figVMEC[[j]]  = ParametricPlot3D[(*VMECNposvecpsi[theta, phi]*)FluxSurfaceVMEC[theta, phi, Floor[iradius]], {theta, 0, 2*pi}, {phi, phimin[[j]], phimax[[j]]}, Boxed -> False, Axes -> False, ImageSize -> Large, Mesh -> None, PlotStyle -> Directive[(*Opacity[0.3+0.7*j/nfigsSurf],*) Specularity[White, 20]], PlotPoints -> {30, plotPointsFig}, MaxRecursion -> maxRecursPlot
        , ColorFunctionScaling -> False, ColorFunction -> (ColorData[{"Rainbow", { minB0, maxB0}}][(BFunc[#4,#5])] &)];
      WriteString[$Output, " Successful in "<>ToString[AbsoluteTime[]-timestart2]<>"s\n"];
    ];
  ,{iradius,1,nsurfaces,If[nfigsSurf==1,1,(nsurfaces-1)/(nfigsSurf-1)]}];
  timestart2 = AbsoluteTime[]; WriteString[$Output, "      Exporting Figure..."];
  figVMECplot=Show[figVMEC(*,ParametricPlot3D[closedcurvVMEC[phi],{phi,0,2*pi},PlotStyle->Black]*)];
  Export[StringJoin[currentPath,"data/",ToString[proj],"/Figures/",ToString[proj],"_VMEC_3Dsurface_all.pdf"], figVMECplot, ImageSize->ImageSizePlot, ImageResolution->ImageResolutionPlot];
  WriteString[$Output, " Successful in "<>ToString[AbsoluteTime[]-timestart2]<>"s\n"];
  timestart2 = AbsoluteTime[]; WriteString[$Output, "      Exporting Figure..."];
  figVMECplot=Show[figVMEC[[1;;Floor[0.75*nfigsSurf]]](*,ParametricPlot3D[closedcurvVMEC[phi],{phi,0,2*pi},PlotStyle->Black]*)];
  Export[StringJoin[currentPath,"data/",ToString[proj],"/Figures/",ToString[proj],"_VMEC_3Dsurface_half.pdf"], figVMECplot, ImageSize->ImageSizePlot, ImageResolution->ImageResolutionPlot];
  WriteString[$Output, " Successful in "<>ToString[AbsoluteTime[]-timestart2]<>"s\n"];
  WriteString[$Output, "    Successful in "<>ToString[AbsoluteTime[]-timestart]<>"s\n"];
];

If[plotFit==1,
(*Plot Iota Profiles*)
If[muDefined==0 && nsurfaces>1,
  WriteString[$Output, "    Plotting Iota Profiles...  "];timestart = AbsoluteTime[];
  If[FileExistsQ[vmecOutput],
    iotaFitPlot=Table[{psiTable[[j]],iotaFit[[j]]},{j,1,Length[iotaFit]}];
    iotaAxisOutPlot=Table[{psiTable[[j]],iotaAxisOut[[j]]},{j,1,Length[iotaAxisOut]}];
    fig1=ListPlot[{iotaAxisOutPlot, iotaFitPlot}, PlotOptions, AspectRatio -> 1, FrameLabel -> {HoldForm["\[Psi]/" Subscript[\[Psi], b]], "Rotational Transform \[Iota]"}, PlotLegends -> {"VMEC \[Iota]", "Near-Axis \[Iota] on Axis"}, PlotStyle -> PointSize[0.018], PlotRange->{{0.,Automatic},{All,All}}];
    ,
    fig1=ListPlot[iotaFit, AspectRatio -> 1, PlotOptions, FrameLabel -> {"Surface", "Rotational Transform \[Iota]"}, PlotStyle -> PointSize[0.018]];
  ];
  Export[StringJoin[currentPath,"data/",ToString[proj],"/Figures/",ToString[proj],"_iota_profile_ordern_"<>ToString[ordern]<>".pdf"],fig1, "AllowRasterization" -> True, ImageResolution -> 200];
  WriteString[$Output, "Successful in "<>ToString[AbsoluteTime[]-timestart]<>"s\n"];
];

(*Compare VMEC and Fitted Surface*)
If[muDefined==0 && plotPolFig==1,
  If[nfigsSurf>1,
    WriteString[$Output, "    Plotting Comparison at Different Poloidal Planes 1...\n"];timestart = AbsoluteTime[];
    j=0;Table[j=j+1;
      timestart2 = AbsoluteTime[]; WriteString[$Output, "      Creating Plot for figure "<>ToString[j]<>" of "<>ToString[nPlots]<>"..."];
      fig2=Show[{
        ParametricPlot[Evaluate[{RBCFit[phi,  theta, #], ZBSFit[phi,  theta, #]} & /@ Table[iradius, {iradius, 1, nsurfaces}][[1 ;; nsurfaces]]], {theta, 0, 2*pi}, PlotPoints->npointsPolPlots, MaxRecursion->1, PlotStyle -> Dashed]//Quiet, 
        ParametricPlot[Evaluate[{RBCVMEC[theta, phi, #], ZBSVMEC[theta, phi, #]} & /@ Table[iradius, {iradius, 1, nsurfaces}][[1 ;; nsurfaces]]], {theta, 0, 2*pi}, PlotPoints->npointsPolPlots, MaxRecursion->1]//Quiet},
        PlotRange -> All, AspectRatio -> Automatic, PlotOptions, FrameLabel -> {"R [meters]", "Z [meters]"}];
      Export[StringJoin[currentPath,"data/",ToString[proj],"/Figures/",ToString[proj],"_poloidalPlane_"<>ToString[j]<>"_ordern_"<>ToString[ordern]<>".pdf"],fig2, "AllowRasterization" -> True];
      WriteString[$Output, " Successful in "<>ToString[AbsoluteTime[]-timestart2]<>"s\n"];
    ,{phi, 0, 2*pi/vmecNFP*(1 - 1/nPlots), 2*pi/nPlots/vmecNFP}];
    WriteString[$Output, "    Successful in "<>ToString[AbsoluteTime[]-timestart]<>"s\n"];
  ];
  WriteString[$Output, "    Plotting Comparison at Different Poloidal Planes 2...  \n"];timestart = AbsoluteTime[];
  j=0;Table[j=j+1;
    timestart2 = AbsoluteTime[]; WriteString[$Output, "      Creating Plot for figure "<>ToString[j]<>" of "<>ToString[nsurfaces]<>"..."];
    fig3[j]=Show[{
      ParametricPlot[Evaluate[{RBCFit[#,  theta, Floor[iradius]], ZBSFit[#,  theta, Floor[iradius]]} & /@ Table[phi, {phi, 0, 2*pi/vmecNFP, 2*pi/nPlots/vmecNFP}][[1 ;; nPlots]]], {theta, 0, 2*pi}, PlotPoints->npointsPolPlots, MaxRecursion->1, PlotStyle -> Dashed]//Quiet, 
      ParametricPlot[Evaluate[{RBCVMEC[theta, #, Floor[iradius]], ZBSVMEC[theta, #, Floor[iradius]]} & /@ Table[phi, {phi, 0, 2*pi/vmecNFP, 2*pi/nPlots/vmecNFP}][[1 ;; nPlots]]], {theta, 0, 2*pi}, PlotPoints->npointsPolPlots, MaxRecursion->1]//Quiet},
    PlotRange -> All, AspectRatio -> Automatic, PlotOptions, FrameLabel -> {"R [meters]", "Z [meters]"}];
    WriteString[$Output, " Successful in "<>ToString[AbsoluteTime[]-timestart2]<>"s\n"];
    timestart2 = AbsoluteTime[]; WriteString[$Output, "      Exporting Figure..."];
    Export[StringJoin[currentPath,"data/",ToString[proj],"/Figures/",ToString[proj],"_poloidalPlane_Surface_"<>ToString[Floor[iradius]]<>"_ordern_"<>ToString[ordern]<>".pdf"],fig3[j], "AllowRasterization" -> True];
    WriteString[$Output, " Successful in "<>ToString[AbsoluteTime[]-timestart2]<>"s\n"];
  ,{iradius,1,nsurfaces,If[nfigsSurf==1,1,(nsurfaces-1)/(nfigsSurf-1)]}];
  WriteString[$Output, "    Successful in "<>ToString[AbsoluteTime[]-timestart]<>"s\n"];
];

(*Plot Obtained Surface*)
If[export3DSurface==1 || exportBFieldSurface==1,
  timestartplotfit = AbsoluteTime[]; WriteString[$Output, "    Plotting Obtained Surface...\n"];
  figFit=Array[Null,nfigsSurf];
  j=0;Table[j=j+1;WriteString[$Output, "     Surface "<>ToString[j]<>" of "<>ToString[nfigsSurf],"...\n"];
  If[muDefined==0,
    phiMmin=phiAxis[j][0,phimin[[j]]];
    phiMmax=phiAxis[j][0,phimax[[j]]];
    If[export3DSurface==1,
      (**)timestart2 = AbsoluteTime[]; WriteString[$Output, "      Creating Grid for Figure "<>ToString[j]<>" of "<>ToString[nfigsSurf],"..."];
      posvecpsiTable = Table[ParallelTable[FluxSurfacenFitVMEC[-theta, phi, Floor[iradius]][[i]], {theta, 0, 2*pi, 2*pi/nPlotTheta}, {phi, phiMmin, phiMmax, (phiMmax-phiMmin)/nPlotPhi}], {i, 1, 3}];
      Nposvecpsitemp = ParallelTable[ListInterpolation[posvecpsiTable[[i]], {{0, 2*pi}, {phiMmin, phiMmax}}], {i, 1, 3}];
      Nposvecpsi[theta_, phi_] = Table[Nposvecpsitemp[[i]][theta, phi], {i, 1, 3}];
      WriteString[$Output, " Successful in "<>ToString[AbsoluteTime[]-timestart2]<>"s\n"];(**)
      (*Nposvecpsi[theta_, phi_] = FluxSurfacenFitVMEC[-theta, phi, Floor[iradius]];*)
    ];
    muN = muVMEC /. fitParams[Floor[iradius]];deltaN = deltaVMEC /. fitParams[Floor[iradius]];
    B0N = B0VMEC /. fitParams[Floor[iradius]];dmuN = D[muN,phi]/sprimeVMEC[phi];ddmuN = D[dmuN,phi]/sprimeVMEC[phi];
    ddeltaN = D[deltaN,phi]/sprimeVMEC[phi];dddeltaN = D[ddeltaN,phi]/sprimeVMEC[phi];
    dB0N = D[B0N,phi]/sprimeVMEC[phi];dtorsVMEC = D[torsVMEC[phi],phi]/sprimeVMEC[phi];dcurvVMEC = D[curvVMEC[phi],phi]/sprimeVMEC[phi];
    rhoN[theta_,phi_] = rho[theta, phi] /. fitParams[Floor[iradius]] /.PSI->vmecPSI[[Floor[iradius]]];
    BFunc[theta_, phi_]=Chop[(1 + rhoN[theta, phi]*Cos[theta]*curvVMEC[phi])*(B0N),10^-6];
    If[ordern==3, BFunc[theta_, phi_] = Chop[BFunc[theta,phi] + 1/2 B0N rhoN[theta, phi]^2 (-((dmuN^2 muN Cos[2 (deltaN + theta)])/(-1 + muN^2)^2) + (ddmuN Cos[2 (deltaN + theta)])/(2 (-1 + muN^2)) + (dB0N dmuN Cos[2 (deltaN + theta)])/(2 B0N (-1 + muN^2)) + 2 Cos[theta]^2 curvVMEC[phi]^2 + (dddeltaN - dtorsVMEC) muN Sin[2 (deltaN + theta)] + (-((B0N dmuN Sin[2 (deltaN + theta)])/(2 (-1 + muN^2))) + B0N muN Cos[2 (deltaN + theta)] (ddeltaN - torsVMEC[phi]))^2/B0N^2 + ((B0N dmuN Cos[2 (deltaN + theta)])/(2 (-1 + muN^2)) + B0N muN Sin[2 (deltaN + theta)] (ddeltaN - torsVMEC[phi]))^2/B0N^2 + dmuN Sin[2 (deltaN + theta)] (ddeltaN - torsVMEC[phi]) + (dB0N muN Sin[2 (deltaN + theta)] (ddeltaN - torsVMEC[phi]))/B0N + 2 muN Cos[2 (deltaN + theta)] (ddeltaN - torsVMEC[phi])^2 + (dmuN Sin[2 (deltaN + theta)] (-ddeltaN + torsVMEC[phi]))/(-1 + muN^2)),10^-6]];
  ,
    phiMmin=0;phiMmax=2*pi;
    If[export3DSurface==1,
      Nposvecpsi[theta_, phi_] = FluxSurfacenFitVMEC[theta, phi, Floor[iradius]];
    ];
    BFunc[theta_, phi_]=(1 + rho[theta, phi]*Cos[theta]*curvVMEC[phi])*B0[phi];
    BFuncExportB[theta_, phi_]=BFunc[theta, phi];
  ];
  If[exportBFieldSurface==1,
    timestart2 = AbsoluteTime[]; WriteString[$Output, "      Exporting magnetic field at surface..."];
    (*timestart=AbsoluteTime[];WriteString[$Output, "\n        Computing thetaMercierInterp..."];
    thetaMercierInterp=ListInterpolation[Table[thetaMercierFunc[theta, phi, Floor[iradius]], {theta, 0, 2*pi, 2*pi/nthetapointsBsurface}, {phi, 0, 2*pi/vmecNFP, 2*pi/vmecNFP/ntphipointsBsurface}], {{0,2*pi}, {0,2*pi/vmecNFP}}];
    WriteString[$Output, " Successful in "<>ToString[AbsoluteTime[]-timestart]<>"s\n"];*)
    timestart=AbsoluteTime[];WriteString[$Output, "\n        Computing BFuncExportBp..."];
    BFuncExportB = ListInterpolation[Table[BFunc[thetaMercierFunc[theta, phi, Floor[iradius]], phiAxis[Floor[iradius]][theta,phi]], {theta, 0, 2*pi, 2*pi/nthetapointsBsurface}, {phi, 0, 2*pi/vmecNFP, 2*pi/vmecNFP/ntphipointsBsurface}], {{0,2*pi}, {0,2*pi/vmecNFP}}];
    WriteString[$Output, " Successful in "<>ToString[AbsoluteTime[]-timestart]<>"s\n"];
    timestart=AbsoluteTime[];WriteString[$Output, "        Creating figure..."];
    bfitFig=DensityPlot[BFuncExportB[theta, phi], {phi, 0, 2*pi/vmecNFP}, {theta, 0, 2*pi}, FrameLabel -> {"Toroidal Angle \[CapitalPhi]", HoldForm["Poloidal Angle " Subscript[\[Theta], V]]}, Evaluate@PlotOptions, MaxRecursion -> 1, PlotRange->All, PlotPoints -> {nthetapointsBsurface,ntphipointsBsurface}, ColorFunction -> "Rainbow", PlotLegends -> BarLegend[Automatic, LegendMarkerSize -> ImageSizePlot, LegendLabel -> HoldForm["B(" Subscript[\[Theta], V] ",\[CapitalPhi])"]]];
    WriteString[$Output, " Successful in "<>ToString[AbsoluteTime[]-timestart]<>"s\n"];
    timestart=AbsoluteTime[];WriteString[$Output, "        Exporting BField figure..."];
    Export[StringJoin[currentPath,"data/",ToString[proj],"/Figures/",ToString[proj],"_fit_magnetic_field_surface_"<>ToString[Floor[iradius]]<>"_ordern_"<>ToString[ordern]<>"_.pdf"], bfitFig, ImageSize->ImageSizePlot, ImageResolution->ImageResolutionPlot];
    WriteString[$Output, " Successful in "<>ToString[AbsoluteTime[]-timestart]<>"s\n"];
    WriteString[$Output, "      Successful in "<>ToString[AbsoluteTime[]-timestart2]<>"s\n"];
  ];
  If[export3DSurface==1,
    If[muDefined==1,signB=1,signB=-1];
    timestart2 = AbsoluteTime[]; WriteString[$Output, "      Computing Maximum and Minimum B..."];
    minB0 = FindMinValue[{BFunc[theta,phi], 0 <= theta <= 2*pi && phiMmin <= phi <= phiMmax}, {theta, phi}] // Quiet;
    maxB0 = FindMaxValue[{BFunc[theta,phi], 0 <= theta <= 2*pi && phiMmin <= phi <= phiMmax}, {theta, phi}] // Quiet;
    WriteString[$Output, " Successful in "<>ToString[AbsoluteTime[]-timestart2]<>"s\n"];
    timestart2 = AbsoluteTime[]; WriteString[$Output, "      Creating Plot..."];
    figFit[[j]]  = ParametricPlot3D[Nposvecpsi[theta, phi], {theta, 0, 2*pi}, {phi, phiMmin, phiMmax}, Boxed -> False, Axes -> False, ImageSize -> Large, Mesh -> None, PlotStyle -> Directive[(*Opacity[0.3+0.7*j/nfigsSurf],*) Specularity[White, 20]], PlotPoints -> {30, plotPointsFig}, MaxRecursion -> maxRecursPlot
    , ColorFunctionScaling -> False, ColorFunction -> (ColorData[{"Rainbow", { minB0, maxB0}}][(BFunc[#4*signB,#5])] &)];
    WriteString[$Output, " Successful in "<>ToString[AbsoluteTime[]-timestart2]<>"s\n"];
  ];
  ,{iradius,1,nsurfaces,If[nfigsSurf==1,1,(nsurfaces-1)/(nfigsSurf-1)]}];
  If[export3DSurface==1,
    timestart2 = AbsoluteTime[]; WriteString[$Output, "      Exporting Figure..."];
    figFitPlot=Show[figFit(*,ParametricPlot3D[closedcurvVMEC[phi],{phi,0,2*pi},PlotStyle->Black]*)];
    Export[StringJoin[currentPath,"data/",ToString[proj],"/Figures/",ToString[proj],"_senac_3Dsurface_all_ordern_"<>ToString[ordern]<>".pdf"], figFitPlot, ImageSize->ImageSizePlot, ImageResolution->ImageResolutionPlot];
    WriteString[$Output, " Successful in "<>ToString[AbsoluteTime[]-timestart2]<>"s\n"];
    If[nfigsSurf>1,
      timestart2 = AbsoluteTime[]; WriteString[$Output, "      Exporting Figure..."];
      figFitPlot=Show[figFit[[1 ;; Floor[0.75*nfigsSurf]]](*,ParametricPlot3D[closedcurvVMEC[phi],{phi,0,2*pi},PlotStyle->Black]*)];
      Export[StringJoin[currentPath,"data/",ToString[proj],"/Figures/",ToString[proj],"_senac_3Dsurface_half_ordern_"<>ToString[ordern]<>".pdf"], figFitPlot, ImageSize->ImageSizePlot, ImageResolution->ImageResolutionPlot];
      WriteString[$Output, " Successful in "<>ToString[AbsoluteTime[]-timestart2]<>"s\n"];
    ];
  ];
  WriteString[$Output, "    Successful in "<>ToString[AbsoluteTime[]-timestartplotfit]<>"s\n"];
];
];
Clear[fig1,fig2,posvecpsiTable,figVMEC,figFit,figVMEC,Nposvecpsitemp];