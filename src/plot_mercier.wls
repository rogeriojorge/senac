(* ::Package:: *)

(*Plotting Options*)
nPlotTheta=30;nPlotPhi=50;plotPointsFig=50;maxRecursPlot=2;ImageSizePlot=1000;ImageResolutionPlot=500;nfigsSurf=Min[4,nsurfaces];phimax=0.9*2*pi;
PlotOptions = Flatten[{LabelStyle -> Directive[FontFamily -> "Latin Modern Roman", FontSize -> 30], ImageSize -> ImageSizePlot, Frame -> True, AspectRatio -> 1}];

(*Plot Original Surface and Axis*)
If[muDefined==0,
  figVMEC=Array[Null,nfigsSurf];
  timestart = AbsoluteTime[]; WriteString[$Output, "    Plotting Original Surface..."];
  j=0;Table[j=j+1;
    VMECposvecpsiTable = Table[ParallelTable[FluxSurfaceVMEC[theta, phi, Floor[iradius]][[i]], {theta, 0, 2*pi, 2*pi/nPlotTheta}, {phi, 0, 2*pi, 2*pi/nPlotPhi}], {i, 1, 3}];
    VMECNposvecpsitemp = ParallelTable[ListInterpolation[VMECposvecpsiTable[[i]], {{0, 2*pi}, {0, 2*pi}}], {i, 1, 3}];
    VMECNposvecpsi[theta_, phi_] = Table[VMECNposvecpsitemp[[i]][theta, phi], {i, 1, 3}];
    If[FileExistsQ[vmecOutput],
      minB0VMEC = FindMinValue[{BVMEC[theta,phi,Floor[iradius]], 0 <= theta <= 2*pi && 0 <= phi <= 2*pi}, {theta, phi}] // Quiet;
      maxB0VMEC = FindMaxValue[{BVMEC[theta,phi,Floor[iradius]], 0 <= theta <= 2*pi && 0 <= phi <= 2*pi}, {theta, phi}] // Quiet;
      figVMEC[[j]]  = ParametricPlot3D[VMECNposvecpsi[theta, phi], {theta, 0, 2*pi}, {phi, 0, phimax}, Boxed -> False, Axes -> False, ImageSize -> Large, Mesh -> None, PlotStyle -> Directive[Opacity[0.8], Specularity[White, 20]], PlotPoints -> plotPointsFig, MaxRecursion -> maxRecursPlot,
          ColorFunctionScaling -> False, ColorFunction -> (ColorData[{"Rainbow", {minB0VMEC, maxB0VMEC}}][(BVMEC[#4,#5,Floor[iradius]])] &)];
      ,
      BFunc[theta_, phi_]=1 + (rho[theta, phi]/.fitParams[Floor[iradius]])*Cos[theta]*curvVMEC[phi];
      minB0 = FindMinValue[{BFunc[theta,phi], 0 <= theta <= 2*pi && 0 <= phi <= 2*pi}, {theta, phi}] // Quiet;
      maxB0 = FindMaxValue[{BFunc[theta,phi], 0 <= theta <= 2*pi && 0 <= phi <= 2*pi}, {theta, phi}] // Quiet;
      figVMEC[[j]]  = ParametricPlot3D[VMECNposvecpsi[theta, phi], {theta, 0, 2*pi}, {phi, 0, 2*pi}, Boxed -> False, Axes -> False, ImageSize -> Large, Mesh -> None, PlotStyle -> Directive[Opacity[0.8], Specularity[White, 20]], PlotPoints -> plotPointsFig, MaxRecursion -> maxRecursPlot,
          ColorFunctionScaling -> False, ColorFunction -> (ColorData[{"Rainbow", { minB0, maxB0}}][(BFunc[#4,#5])] &)];
    ];
  ,{iradius,1,nsurfaces,If[nfigsSurf==1,1,(nsurfaces-1)/(nfigsSurf-1)]}];
  figVMECplot=Show[figVMEC,ParametricPlot3D[closedcurvVMEC[phi],{phi,0,2*pi},PlotStyle->Black]];
  Export[StringJoin[currentPath,"data/",ToString[proj],"/Figures/",ToString[proj],"_VMEC_3Dsurface.pdf"], figVMECplot, ImageSize->ImageSizePlot, ImageResolution->ImageResolutionPlot];
  WriteString[$Output, " Successful in "<>ToString[AbsoluteTime[]-timestart]<>"s\n"];
];

(*Plot Obtained Surface and Axis*)
timestart = AbsoluteTime[]; WriteString[$Output, "    Plotting Obtained Surface..."];
figFit=Array[Null,nfigsSurf];
j=0;Table[j=j+1;
  posvecpsiTable = Table[ParallelTable[FluxSurfacenFitVMEC[theta, phi, Floor[iradius]][[i]], {theta, 0, 2*pi, 2*pi/nPlotTheta}, {phi, 0, 2*pi, 2*pi/nPlotPhi}], {i, 1, 3}];
  Nposvecpsitemp = ParallelTable[ListInterpolation[posvecpsiTable[[i]], {{0, 2*pi}, {0, 2*pi}}], {i, 1, 3}];
  Nposvecpsi[theta_, phi_] = Table[Nposvecpsitemp[[i]][theta, phi], {i, 1, 3}];
  If[muDefined==0,rho[theta_, phi_]=rho[theta, phi]/.fitParams[Floor[iradius]]/.PSI->vmecPSI[[iradius]]];
  BFunc[theta_, phi_]=1 + rho[theta, phi]*Cos[theta]*curvVMEC[phi];
  minB0 = FindMinValue[{BFunc[theta,phi], 0 <= theta <= 2*pi && 0 <= phi <= 2*pi}, {theta, phi}] // Quiet;
  maxB0 = FindMaxValue[{BFunc[theta,phi], 0 <= theta <= 2*pi && 0 <= phi <= 2*pi}, {theta, phi}] // Quiet;
  figFit[[j]]  = ParametricPlot3D[Nposvecpsi[theta, phi], {theta, 0, 2*pi}, {phi, 0, phimax}, Boxed -> False, Axes -> False, ImageSize -> Large, Mesh -> None, PlotStyle -> Directive[Opacity[0.8], Specularity[White, 20]], PlotPoints -> plotPointsFig, MaxRecursion -> maxRecursPlot,
      ColorFunctionScaling -> False, ColorFunction -> (ColorData[{"Rainbow", { minB0, maxB0}}][(BFunc[#4,#5])] &)];
,{iradius,1,nsurfaces,If[nfigsSurf==1,1,(nsurfaces-1)/(nfigsSurf-1)]}];
figFitPlot=Show[figFit,ParametricPlot3D[closedcurvVMEC[phi],{phi,0,2*pi},PlotStyle->Black]];
Export[StringJoin[currentPath,"data/",ToString[proj],"/Figures/",ToString[proj],"_senac_3Dsurface.pdf"], figFitPlot, ImageSize->ImageSizePlot, ImageResolution->ImageResolutionPlot];
WriteString[$Output, " Successful in "<>ToString[AbsoluteTime[]-timestart]<>"s\n"];

(*Plot Iota Profiles*)
If[muDefined==0,
  WriteString[$Output, "    Plotting Iota Profiles...  "];timestart = AbsoluteTime[];
  If[FileExistsQ[vmecOutput],
    fig=ListPlot[{iotaAxisOut, iotaFit}, PlotOptions, FrameLabel -> {"Surface", "Rotational Transform \[Iota]"}, PlotLegends -> {"Original \[Iota]", "Near-Axis \[Iota]"}, PlotStyle -> PointSize[0.018], PlotRange->{{0.5,Automatic},{Automatic,Automatic}}];
    ,
    fig=ListPlot[iotaFit, PlotOptions, FrameLabel -> {"Surface", "Rotational Transform \[Iota]"}, PlotStyle -> PointSize[0.018]];
  ];
  Export[StringJoin[currentPath,"data/",ToString[proj],"/Figures/",ToString[proj],"_iota_profile.pdf"],fig, "AllowRasterization" -> True, ImageSize -> 500, ImageResolution -> 200];
  WriteString[$Output, "Successful in "<>ToString[AbsoluteTime[]-timestart]<>"s\n"];
];

(*Compare VMEC and Fitted Surface*)
If[muDefined==0,
  nPlots = 4;
  WriteString[$Output, "    Plotting Comparison at Different Poloidal Planes...  "];timestart = AbsoluteTime[];
  ParallelTable[
    Off[ParametricPlot::precw,NMinimize::precw,NMaximize::precw];
    rangePlot = {
      {NMinimize[{RBCVMEC[phi, theta, Floor[iradius]], 0 < phi <= 2*pi, 0 < theta <= 2*pi}, {phi, theta}][[1]],
      NMaximize[{ RBCVMEC[phi, theta, Floor[iradius]], 0 < phi <= 2*pi, 0 < theta <= 2*pi}, {phi, theta}][[1]]},
      {NMinimize[{ZBSVMEC[phi, theta, Floor[iradius]], 0 < phi <= 2*pi, 0 < theta <= 2*pi}, {phi, theta}][[1]],
      NMaximize[{ ZBSVMEC[phi, theta, Floor[iradius]], 0 < phi <= 2*pi, 0 < theta <= 2*pi}, {phi, theta}][[1]]}
    }//Quiet;
    fig=Show[{
      ParametricPlot[Evaluate[{RBCFit[#, theta, Floor[iradius]],  ZBSFit[#, theta, Floor[iradius]]}  & /@ Table[phi, {phi, 0, 2*pi/vmecNFP, 2*pi/nPlots/vmecNFP}]], {theta, 0, 2*pi}, PlotPoints->20, MaxRecursion->1, PlotStyle -> Dashed], 
      ParametricPlot[Evaluate[{RBCVMEC[theta, #, Floor[iradius]], ZBSVMEC[theta, #, Floor[iradius]]} & /@ Table[phi, {phi, 0, 2*pi/vmecNFP, 2*pi/nPlots/vmecNFP}]], {theta, 0, 2*pi}, PlotPoints->20, MaxRecursion->1]},
    PlotRange -> rangePlot, PlotPoints -> 60, MaxRecursion -> 1, PlotOptions, FrameLabel -> {"R [meters]", "Z [meters]"}];
    Export[StringJoin[currentPath,"data/",ToString[proj],"/Figures/",ToString[proj],"_poloidalPlane_"<>ToString[Floor[iradius]]<>".pdf"],fig, "AllowRasterization" -> True, ImageSize -> 500, ImageResolution -> 200];
  ,{iradius,1,nsurfaces,If[nfigsSurf==1,1,(nsurfaces-1)/(nfigsSurf-1)]}];
  Clear[fig,posvecpsiTable,figVMEC,figFit,Nposvecpsitemp];
  WriteString[$Output, "Successful in "<>ToString[AbsoluteTime[]-timestart]<>"s\n"];
]
