(* ::Package:: *)

(*Plotting Options*)
nPlotTheta=80;nPlotPhi=120;plotPointsFig=60;maxRecursPlot=2;ImageSizePlot=700;ImageResolutionPlot=500;nfigsSurf=Min[3,nsurfaces];
PlotOptions = Flatten[{LabelStyle -> Directive[FontFamily -> "Latin Modern Roman", FontSize -> 30], ImageSize -> ImageSizePlot, Frame -> True, AspectRatio -> 1}];
phimax = Table[phi + 2*pi (1 - (nfigsSurf - 1)/nfigsSurf/vmecNFP), {phi, 2*pi/vmecNFP (1 - 1/nfigsSurf), 0, -2*pi/nfigsSurf/vmecNFP/2.6}]-0.03;
phimin = Table[phi, {phi, 0, 2*pi/vmecNFP (1 - 1/nfigsSurf), 2*pi/nfigsSurf/vmecNFP/2.6}]+0.03;

(*Plot Original Surface and Axis*)
If[muDefined==0 && plotOriginal==1,
  figVMEC=Array[Null,nfigsSurf];
  timestart = AbsoluteTime[]; WriteString[$Output, "    Plotting Original Surface...\n"];
  j=0;Table[j=j+1;
    timestart2 = AbsoluteTime[]; WriteString[$Output, "      Creating Grid for Figure "<>ToString[j]<>" of "<>ToString[nfigsSurf],"..."];
    VMECposvecpsiTable = Table[ParallelTable[FluxSurfaceVMEC[theta, phi, Floor[iradius]][[i]], {theta, 0, 2*pi, 2*pi/(nPlotTheta/2)}, {phi, phimin[[j]], phimax[[j]], phimax[[j]]/(nPlotPhi/2)}], {i, 1, 3}];
    VMECNposvecpsitemp = ParallelTable[ListInterpolation[VMECposvecpsiTable[[i]], {{0, 2*pi}, {phimin[[j]], phimax[[j]]}}], {i, 1, 3}];
    VMECNposvecpsi[theta_, phi_] = Table[VMECNposvecpsitemp[[i]][theta, phi], {i, 1, 3}];
    WriteString[$Output, " Successful in "<>ToString[AbsoluteTime[]-timestart2]<>"s\n"];
    If[FileExistsQ[vmecOutput],
      timestart2 = AbsoluteTime[]; WriteString[$Output, "      Computing Maximum and Minimum B..."];
      minB0VMEC = FindMinValue[{BVMEC[theta,phi,Floor[iradius]], 0 <= theta <= 2*pi && phimin[[j]] <= phi <= phimax[[j]]}, {theta, phi}] // Quiet;
      maxB0VMEC = FindMaxValue[{BVMEC[theta,phi,Floor[iradius]], 0 <= theta <= 2*pi && phimin[[j]] <= phi <= phimax[[j]]}, {theta, phi}] // Quiet;
      WriteString[$Output, " Successful in "<>ToString[AbsoluteTime[]-timestart2]<>"s\n"];
      timestart2 = AbsoluteTime[]; WriteString[$Output, "      Creating Plot..."];
      figVMEC[[j]]  = ParametricPlot3D[VMECNposvecpsi[theta, phi], {theta, 0, 2*pi}, {phi, phimin[[j]], phimax[[j]]}, Boxed -> False, Axes -> False, ImageSize -> Large, Mesh -> None, PlotStyle -> Directive[Opacity[0.3+0.7*j/nsurfaces], Specularity[White, 20]], PlotPoints -> plotPointsFig, MaxRecursion -> maxRecursPlot,
          ColorFunctionScaling -> False, ColorFunction -> (ColorData[{"Rainbow", {minB0VMEC, maxB0VMEC}}][(BVMEC[#4,#5,Floor[iradius]])] &)];
      WriteString[$Output, " Successful in "<>ToString[AbsoluteTime[]-timestart2]<>"s\n"];
      ,
      timestart2 = AbsoluteTime[]; WriteString[$Output, "      Computing Maximum and Minimum B..."];
      BFunc[theta_, phi_]=1 + (rho[theta, phi]/.fitParams[Floor[iradius]])*Cos[theta]*curvVMEC[phi];
      minB0 = FindMinValue[{BFunc[theta,phi], 0 <= theta <= 2*pi && phimin[[j]] <= phi <= phimax[[j]]}, {theta, phi}] // Quiet;
      maxB0 = FindMaxValue[{BFunc[theta,phi], 0 <= theta <= 2*pi && phimin[[j]] <= phi <= phimax[[j]]}, {theta, phi}] // Quiet;
      WriteString[$Output, " Successful in "<>ToString[AbsoluteTime[]-timestart2]<>"s\n"];
      timestart2 = AbsoluteTime[]; WriteString[$Output, "      Creating Plot..."];
      figVMEC[[j]]  = ParametricPlot3D[VMECNposvecpsi[theta, phi], {theta, 0, 2*pi}, {phi, phimin[[j]], phimax[[j]]}, Boxed -> False, Axes -> False, ImageSize -> Large, Mesh -> None, PlotStyle -> Directive[Opacity[0.3+0.7*j/nsurfaces], Specularity[White, 20]], PlotPoints -> plotPointsFig, MaxRecursion -> maxRecursPlot,
          ColorFunctionScaling -> False, ColorFunction -> (ColorData[{"Rainbow", { minB0, maxB0}}][(BFunc[#4,#5])] &)];
      WriteString[$Output, " Successful in "<>ToString[AbsoluteTime[]-timestart2]<>"s\n"];
    ];
  ,{iradius,1,nsurfaces,If[nfigsSurf==1,1,(nsurfaces-1)/(nfigsSurf-1)]}];
  timestart2 = AbsoluteTime[]; WriteString[$Output, "      Exporting Figure..."];
  figVMECplot=Show[figVMEC,ParametricPlot3D[closedcurvVMEC[phi],{phi,0,2*pi},PlotStyle->Black]];
  Export[StringJoin[currentPath,"data/",ToString[proj],"/Figures/",ToString[proj],"_VMEC_3Dsurface.pdf"], figVMECplot, ImageSize->ImageSizePlot, ImageResolution->ImageResolutionPlot];
  WriteString[$Output, " Successful in "<>ToString[AbsoluteTime[]-timestart2]<>"s\n"];
  WriteString[$Output, "    Successful in "<>ToString[AbsoluteTime[]-timestart]<>"s\n"];
];

(*Plot Obtained Surface and Axis*)
timestart = AbsoluteTime[]; WriteString[$Output, "    Plotting Obtained Surface...\n"];
figFit=Array[Null,nfigsSurf];
j=0;Table[j=j+1;
  timestart2 = AbsoluteTime[]; WriteString[$Output, "      Creating Grid for Figure "<>ToString[j]<>" of "<>ToString[nfigsSurf],"..."];
  posvecpsiTable = Table[ParallelTable[FluxSurfacenFitVMEC[theta, phi, Floor[iradius]][[i]], {theta, 0, 2*pi, 2*pi/nPlotTheta}, {phi, phimin[[j]], phimax[[j]], phimax[[j]]/nPlotPhi}], {i, 1, 3}];
  Nposvecpsitemp = ParallelTable[ListInterpolation[posvecpsiTable[[i]], {{0, 2*pi}, {0, phimax[[j]]}}], {i, 1, 3}];
  Nposvecpsi[theta_, phi_] = Table[Nposvecpsitemp[[i]][theta, phi], {i, 1, 3}];
  WriteString[$Output, " Successful in "<>ToString[AbsoluteTime[]-timestart2]<>"s\n"];
  timestart2 = AbsoluteTime[]; WriteString[$Output, "      Computing Maximum and Minimum B..."];
  If[muDefined==0,rhoPlot[theta_, phi_]=rho[theta, phi]/.fitParams[Floor[iradius]]/.PSI->vmecPSI[[Floor[iradius]]]];
  BFunc[theta_, phi_]=1 + rhoPlot[theta, phi]*Cos[theta]*curvVMEC[phi];
  minB0 = FindMinValue[{BFunc[theta,phi], 0 <= theta <= 2*pi && phimin[[j]] <= phi <= phimax[[j]]}, {theta, phi}] // Quiet;
  maxB0 = FindMaxValue[{BFunc[theta,phi], 0 <= theta <= 2*pi && phimin[[j]] <= phi <= phimax[[j]]}, {theta, phi}] // Quiet;
  WriteString[$Output, " Successful in "<>ToString[AbsoluteTime[]-timestart2]<>"s\n"];
  timestart2 = AbsoluteTime[]; WriteString[$Output, "      Creating Plot..."];
  figFit[[j]]  = ParametricPlot3D[Nposvecpsi[theta, phi], {theta, 0, 2*pi}, {phi, phimin[[j]], phimax[[j]]}, Boxed -> False, Axes -> False, ImageSize -> Large, Mesh -> None, PlotStyle -> Directive[Opacity[0.3+0.7*j/nsurfaces], Specularity[White, 20]], PlotPoints -> plotPointsFig, MaxRecursion -> maxRecursPlot,
      ColorFunctionScaling -> False, ColorFunction -> (ColorData[{"Rainbow", { minB0, maxB0}}][(BFunc[#4,#5])] &)];
  WriteString[$Output, " Successful in "<>ToString[AbsoluteTime[]-timestart2]<>"s\n"];
,{iradius,1,nsurfaces,If[nfigsSurf==1,1,(nsurfaces-1)/(nfigsSurf-1)]}];
figFitPlot=Show[figFit,ParametricPlot3D[closedcurvVMEC[phi],{phi,0,2*pi},PlotStyle->Black]];
timestart2 = AbsoluteTime[]; WriteString[$Output, "      Exporting Figure..."];
Export[StringJoin[currentPath,"data/",ToString[proj],"/Figures/",ToString[proj],"_senac_3Dsurface.pdf"], figFitPlot, ImageSize->ImageSizePlot, ImageResolution->ImageResolutionPlot];
WriteString[$Output, " Successful in "<>ToString[AbsoluteTime[]-timestart2]<>"s\n"];
WriteString[$Output, "    Successful in "<>ToString[AbsoluteTime[]-timestart]<>"s\n"];

(*Plot Iota Profiles*)
If[muDefined==0,
  WriteString[$Output, "    Plotting Iota Profiles...  "];timestart = AbsoluteTime[];
  If[FileExistsQ[vmecOutput],
    fig1=ListPlot[{iotaAxisOut, iotaFit}, PlotOptions, FrameLabel -> {"Surface", "Rotational Transform \[Iota]"}, PlotLegends -> {"Original \[Iota]", "Near-Axis \[Iota]"}, PlotStyle -> PointSize[0.018], PlotRange->{{0.9,Automatic},{Automatic,Automatic}}];
    ,
    fig1=ListPlot[iotaFit, PlotOptions, FrameLabel -> {"Surface", "Rotational Transform \[Iota]"}, PlotStyle -> PointSize[0.018]];
  ];
  Export[StringJoin[currentPath,"data/",ToString[proj],"/Figures/",ToString[proj],"_iota_profile.pdf"],fig1, "AllowRasterization" -> True, ImageResolution -> 200];
  WriteString[$Output, "Successful in "<>ToString[AbsoluteTime[]-timestart]<>"s\n"];
];

(*Compare VMEC and Fitted Surface*)
If[muDefined==0,
  nPlots = 4;
  If[nfigsSurf==1,
    WriteString[$Output, "    Plotting Comparison at Different Poloidal Planes...  \n"];timestart = AbsoluteTime[];
    j=0;Table[j=j+1;
      timestart2 = AbsoluteTime[]; WriteString[$Output, "      Computing Maximum and Minimum B for Figure "<>ToString[j]<>" of "<>ToString[nfigsSurf],"..."];
      WriteString[$Output, " Successful in "<>ToString[AbsoluteTime[]-timestart2]<>"s\n"];
      timestart2 = AbsoluteTime[]; WriteString[$Output, "      Creating Plot..."];
      fig2[Floor[iradius]]=Show[{
        ParametricPlot[Evaluate[{RBCFit[#,  theta, Floor[iradius]], ZBSFit[#,  theta, Floor[iradius]]} & /@ Table[phi, {phi, 0, 2*pi/vmecNFP, 2*pi/nPlots/vmecNFP}][[1 ;; nPlots]]], {theta, 0, 2*pi}, PlotPoints->30, MaxRecursion->1, PlotStyle -> Dashed]//Quiet, 
        ParametricPlot[Evaluate[{RBCVMEC[theta, #, Floor[iradius]], ZBSVMEC[theta, #, Floor[iradius]]} & /@ Table[phi, {phi, 0, 2*pi/vmecNFP, 2*pi/nPlots/vmecNFP}][[1 ;; nPlots]]], {theta, 0, 2*pi}, PlotPoints->30, MaxRecursion->1]//Quiet},
      PlotRange -> All, PlotOptions, FrameLabel -> {"R [meters]", "Z [meters]"}];
      WriteString[$Output, " Successful in "<>ToString[AbsoluteTime[]-timestart2]<>"s\n"];
      timestart2 = AbsoluteTime[]; WriteString[$Output, "      Exporting Figure..."];
      Export[StringJoin[currentPath,"data/",ToString[proj],"/Figures/",ToString[proj],"_poloidalPlane_Surface_"<>ToString[Floor[iradius]]<>".pdf"],fig2[Floor[iradius]], "AllowRasterization" -> True];
      WriteString[$Output, " Successful in "<>ToString[AbsoluteTime[]-timestart2]<>"s\n"];
    ,{iradius,1,nsurfaces,If[nfigsSurf==1,1,(nsurfaces-1)/(nfigsSurf-1)]}];
    WriteString[$Output, "    Successful in "<>ToString[AbsoluteTime[]-timestart]<>"s\n"];
    ,
    WriteString[$Output, "    Plotting Comparison at Different Poloidal Planes..."];timestart = AbsoluteTime[];
    j=0;Table[j=j+1;
      fig2=Show[{
        ParametricPlot[Evaluate[{RBCFit[phi,  theta, #], ZBSFit[phi,  theta, #]} & /@ Table[iradius, {iradius, 1, nsurfaces}][[1 ;; nsurfaces]]], {theta, 0, 2*pi}, PlotPoints->30, MaxRecursion->1, PlotStyle -> Dashed]//Quiet, 
        ParametricPlot[Evaluate[{RBCVMEC[theta, phi, #], ZBSVMEC[theta, phi, #]} & /@ Table[iradius, {iradius, 1, nsurfaces}][[1 ;; nsurfaces]]], {theta, 0, 2*pi}, PlotPoints->30, MaxRecursion->1]//Quiet},
        PlotRange -> All, PlotOptions, FrameLabel -> {"R [meters]", "Z [meters]"}];
      Export[StringJoin[currentPath,"data/",ToString[proj],"/Figures/",ToString[proj],"_poloidalPlane_"<>ToString[j]<>".pdf"],fig2, "AllowRasterization" -> True];
    ,{phi, 0, 2*pi/vmecNFP*(1 - 1/nPlots), 2*pi/nPlots/vmecNFP}];
    WriteString[$Output, " Successful in "<>ToString[AbsoluteTime[]-timestart]<>"s\n"];
  ];
];
Clear[fig1,fig2,posvecpsiTable,figVMEC,figFit,figVMEC,Nposvecpsitemp];