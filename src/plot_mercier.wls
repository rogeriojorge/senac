(* ::Package:: *)


(*Plotting Options*)
nPlotTheta=30;nPlotPhi=60;plotPointsFig=40;maxRecursPlot=2;ImageSizePlot=700;ImageResolutionPlot=300;
PlotOptions = Flatten[{LabelStyle -> Directive[FontFamily -> (*"Latin Modern Roman"*)"Helvetica", FontSize -> 30], ImageSize -> ImageSizePlot, PlotStyle -> Thickness[0.007], Frame -> True, AspectRatio -> 1}];
BFunc[theta_, phi_]=1 + rho[theta, phi]*Cos[theta]*curvVMEC[phi];
minB0 = FindMinValue[{BFunc[theta,phi], 0 <= theta <= 2*pi && 0 <= phi <= 2*pi}, {theta, phi}] // Quiet;
maxB0 = FindMaxValue[{BFunc[theta,phi], 0 <= theta <= 2*pi && 0 <= phi <= 2*pi}, {theta, phi}] // Quiet;

(*Plot Original Surface and Axis*)
timestart = AbsoluteTime[]; WriteString[$Output, "    Plotting Original Surface..."];
VMECposvecpsiTable = Table[ParallelTable[FluxSurfaceVMEC[theta, phi][[i]], {theta, 0, 2*pi, 2*pi/nPlotTheta}, {phi, 0, 2*pi, 2*pi/nPlotPhi}], {i, 1, 3}];
VMECNposvecpsitemp = ParallelTable[ListInterpolation[VMECposvecpsiTable[[i]], {{0, 2*pi}, {0, 2*pi}}], {i, 1, 3}];
VMECNposvecpsi[theta_, phi_] = Table[VMECNposvecpsitemp[[i]][theta, phi], {i, 1, 3}];
figVMEC  = Show[{
  ParametricPlot3D[VMECNposvecpsi[theta, phi], {theta, 0, 2*pi}, {phi, 0, 2*pi}, Boxed -> False, Axes -> False, ImageSize -> Large, Mesh -> None, PlotStyle -> Directive[Opacity[0.8], Specularity[White, 20]], PlotPoints -> plotPointsFig, MaxRecursion -> maxRecursPlot,
    ColorFunctionScaling -> False, ColorFunction -> (ColorData[{"Rainbow", { minB0, maxB0}}][(BFunc[#4,#5])] &)],
  ParametricPlot3D[closedcurvVMEC[phi],{phi,0,2*pi}]
}];
Export[StringJoin[currentPath,"data/",ToString[proj],"/Figures/",ToString[proj],"_VMEC_3Dsurface.pdf"], figVMEC, ImageSize->ImageSizePlot, ImageResolution->ImageResolutionPlot];
WriteString[$Output, " Successful in "<>ToString[AbsoluteTime[]-timestart]<>"s\n"];

(*Plot Obtained Surface and Axis*)
timestart = AbsoluteTime[]; WriteString[$Output, "    Plotting Obtained Surface..."];
posvecpsiTable = Table[ParallelTable[FluxSurfacenFitVMEC[theta, phi][[i]], {theta, 0, 2*pi, 2*pi/nPlotTheta}, {phi, 0, 2*pi, 2*pi/nPlotPhi}], {i, 1, 3}];
Nposvecpsitemp = ParallelTable[ListInterpolation[posvecpsiTable[[i]], {{0, 2*pi}, {0, 2*pi}}], {i, 1, 3}];
Nposvecpsi[theta_, phi_] = Table[Nposvecpsitemp[[i]][theta, phi], {i, 1, 3}];
figFit  = Show[{
  ParametricPlot3D[Nposvecpsi[theta, phi], {theta, 0, 2*pi}, {phi, 0, 2*pi}, Boxed -> False, Axes -> False, ImageSize -> Large, Mesh -> None, PlotStyle -> Directive[Opacity[0.8], Specularity[White, 20]], PlotPoints -> plotPointsFig, MaxRecursion -> maxRecursPlot,
    ColorFunctionScaling -> False, ColorFunction -> (ColorData[{"Rainbow", { minB0, maxB0}}][(BFunc[#4,#5])] &)],
  ParametricPlot3D[closedcurvVMEC[phi],{phi,0,2*pi}]
}];
Export[StringJoin[currentPath,"data/",ToString[proj],"/Figures/",ToString[proj],"_senac_3Dsurface.pdf"], figFit, ImageSize->ImageSizePlot, ImageResolution->ImageResolutionPlot];
WriteString[$Output, " Successful in "<>ToString[AbsoluteTime[]-timestart]<>"s\n"];

(*Compare VMEC and Fitted Surface*)
If[muDefined==0,
  nPlots = 4;
  WriteString[$Output, "    Plotting comparison with Surface...  "];timestart = AbsoluteTime[];
  Off[ParametricPlot::precw,NMinimize::precw,NMaximize::precw];
  rangePlot = {
    {NMinimize[{RBCVMEC[phi, theta], 0 < phi <= 2*pi, 0 < theta <= 2*pi}, {phi, theta}][[1]],
    NMaximize[{ RBCVMEC[phi, theta], 0 < phi <= 2*pi, 0 < theta <= 2*pi}, {phi, theta}][[1]]},
    {NMinimize[{ZBSVMEC[phi, theta], 0 < phi <= 2*pi, 0 < theta <= 2*pi}, {phi, theta}][[1]],
    NMaximize[{ ZBSVMEC[phi, theta], 0 < phi <= 2*pi, 0 < theta <= 2*pi}, {phi, theta}][[1]]}
  }//Quiet;
  fig=Show[{
    ParametricPlot[Evaluate[{RBCFit[#, theta],  ZBSFit[#, theta]} & /@ Table[phi, {phi, 0, 2*pi/vmecNFP, 2*pi/nPlots/vmecNFP}]], {theta, 0, 2*pi}, PlotStyle -> Dashed], 
    ParametricPlot[Evaluate[{RBCVMEC[theta, #], ZBSVMEC[theta, #]} & /@ Table[phi, {phi, 0, 2*pi/vmecNFP, 2*pi/nPlots/vmecNFP}]], {theta, 0, 2*pi}]},
  PlotRange -> rangePlot, PlotPoints -> 60, MaxRecursion -> 1, PlotOptions, FrameLabel -> {"R [meters]", "Z [meters]"}];
  Export[StringJoin[currentPath,"data/",ToString[proj],"/Figures/",ToString[proj],"_fit_plotcompare.pdf"],fig, "AllowRasterization" -> True, ImageSize -> 500, ImageResolution -> 200]
  WriteString[$Output, "Successful in "<>ToString[AbsoluteTime[]-timestart]<>"s\n"];
  Clear[fig,posvecpsiTable,figVMEC,figFit,Nposvecpsitemp];
]
