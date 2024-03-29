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
ordern         = ToExpression[$ScriptCommandLine[[3]]];
nModes         = ToExpression[$ScriptCommandLine[[4]]];
ImageSizePlot  = ToExpression[$ScriptCommandLine[[5]]];
ImageResolutionPlot  = ToExpression[$ScriptCommandLine[[6]]];
coilsPerHalfPeriod   = ToExpression[$ScriptCommandLine[[7]]];
numHalfPeriodsToPlot = ToExpression[$ScriptCommandLine[[8]]];
plotRegcoilSENAC     = ToExpression[$ScriptCommandLine[[9]]];
npointsContourPlot = ToExpression[$ScriptCommandLine[[10]]];
npointsInterpolationCoilPos = ToExpression[$ScriptCommandLine[[11]]];
interpolationOrderCoilPos = ToExpression[$ScriptCommandLine[[12]]];
plotPointsFig = ToExpression[$ScriptCommandLine[[13]]];
maxRecursPlot = ToExpression[$ScriptCommandLine[[14]]];
coilthickness = ToExpression[$ScriptCommandLine[[15]]];
plotpointsCoil = ToExpression[$ScriptCommandLine[[16]]];

pi = N[Pi,10];

(*Plot Regcoil Results*)
regcoilOriginalFile = StringJoin[currentPath,"../vmec/",ToString[proj],"/regcoil_out_",ToString[proj],".nc"];
regcoilFitFile = StringJoin[currentPath,"../data/",ToString[proj],"/regcoil_out_",ToString[proj],"_senac_ordern",ToString[ordern],"_nmodes",ToString[nModes],".nc"];

If[plotRegcoilSENAC==1,
  If[FileExistsQ[regcoilOriginalFile],
    regcoilFile = regcoilOriginalFile;,
    WriteString[$Output, "    No REGCOIL VMEC file\n"];
    CloseKernels[];Exit[];
  ];
];

If[plotRegcoilSENAC==2,
  If[FileExistsQ[regcoilFitFile],
    regcoilFile = regcoilFitFile;,
    WriteString[$Output, "    No REGCOIL SENAC file\n"];
    CloseKernels[];Exit[];
  ];
];

(*Plotting Options*)
regcoilDatasets = Import[regcoilFile, {"Datasets", {"xm_plasma", "xn_plasma", "rmnc_plasma", "zmns_plasma", "xm_coil", "xn_coil", "rmnc_coil", "zmns_coil", "xm_potential", "xn_potential", "single_valued_current_potential_mn", "nfp", "net_poloidal_current_Amperes", "net_toroidal_current_Amperes"}}];
xmplasma = regcoilDatasets[[1]];
xnplasma = regcoilDatasets[[2]];
rmncplasma = regcoilDatasets[[3]];
zmnsplasma = regcoilDatasets[[4]];
xmcoil = regcoilDatasets[[5]];
xncoil = regcoilDatasets[[6]];
rmnccoil = regcoilDatasets[[7]];
zmnscoil = regcoilDatasets[[8]];
xmpotential = regcoilDatasets[[9]];
xnpotential = regcoilDatasets[[10]];
svcurrentpotentialmn = regcoilDatasets[[11]];
nfpcoil = regcoilDatasets[[12]];
netpcurrentA = regcoilDatasets[[13]];
nettcurrentA = regcoilDatasets[[14]];

If[numHalfPeriodsToPlot == 0, numHalfPeriodsToPlot = nfpcoil*2];

currentpotentialmn = svcurrentpotentialmn[[Length[svcurrentpotentialmn[[All, 1]]], All]];
currentpotentialTemp[theta_, phi_] = 0;
Table[angle[theta_, phi_] = xmpotential[[im]]*theta - xnpotential[[im]]*phi;
currentpotentialTemp[theta_, phi_] = currentpotentialTemp[theta, phi] + currentpotentialmn[[im]]*Sin[angle[theta, phi]];, {im, 1, Length[xmpotential]}];
currentpotential = Compile[{{theta, _Real}, {phi, _Real}}, currentpotentialTemp[theta, phi] + netpcurrentA*phi/(2*pi) + nettcurrentA*theta/(2*pi)];
Off[CompiledFunction::cfsa];

minCurr = FindMinValue[{currentpotential[theta, phi], 0 <= theta <= 2*pi && 0 <= phi <= 2*pi}, {theta, phi}];
maxCurr = FindMaxValue[{currentpotential[theta, phi], 0 <= theta <= 2*pi && 0 <= phi <= 2*pi}, {theta, phi}];
coilContours = (ParallelTable[(
  Interpolation[
    ContourPlot[currentpotential[theta, phi] == i, {theta, 0, 2*pi}, {phi, 0, 2*pi}, PlotPoints -> npointsContourPlot, MaxRecursion -> 1][[1, 1, 1]],
   Method -> "Spline", InterpolationOrder -> interpolationOrderCoilPos])
  // Quiet, {i, minCurr, maxCurr, (maxCurr - minCurr)/(coilsPerHalfPeriod*2*nfpcoil + 1)}] // Quiet)[[2 ;; coilsPerHalfPeriod*2*nfpcoil + 1]] // Quiet;

If[!NumberQ[coilContours[[1]][1]],
  WriteString[$Output, "    Disconnected contours, closing.\n"];
  CloseKernels[];Exit[];
]

If[Dimensions[Position[rmncplasma, _?(10^-13 < # < 10^-10 &)]][[1]] > 7 && NumberQ[Position[rmncplasma, _?(10^-13 < # < 10^-10 &)][[1, 1]]], 
  imMaxP = Position[rmncplasma, _?(10^-13 < # < 10^-10 &)][[8, 1]], 
  imMaxP = Length[rmncplasma]
];

plasmaposTemp[theta_, phi_]=ParallelSum[
  angle[theta_, phi_] = xmplasma[[im]]*theta - xnplasma[[im]]*phi;
  {rmncplasma[[im]]*Cos[angle[theta, phi]]*Cos[phi], rmncplasma[[im]]*Cos[angle[theta, phi]]*Sin[phi], zmnsplasma[[im]]*Sin[angle[theta, phi]]}
, {im, 1, imMaxP}];

If[Dimensions[Position[rmnccoil, _?(10^(-13) < # < 10^(-12) &)]][[1]] > 7 && NumberQ[Position[rmnccoil, _?(10^(-13) < # < 10^(-12) &)][[1, 1]]], 
  imMax = Position[rmnccoil, _?(10^(-13) < # < 10^(-12) &)][[8, 1]], 
  imMax = Length[rmnccoil]
]//Quiet;
If[imMax < 300 && Length[rmnccoil] > 300, imMax = 300];

WSposTemp[theta_, phi_] = ParallelSum[
   angle[theta_, phi_] = xmcoil[[im]]*theta - xncoil[[im]]*phi;
   {rmnccoil[[im]]*Cos[angle[theta, phi]]*Cos[phi], rmnccoil[[im]]*Cos[angle[theta, phi]]*Sin[phi], zmnscoil[[im]]*Sin[angle[theta, phi]]}
, {im, 1, imMax}];
Table[
  coilposTemp[theta_, j] = ParallelSum[
    angle2[theta_] = coilContours[[j]][theta];
    anglec[theta_] = xmcoil[[im]]*theta - xncoil[[im]]*angle2[theta];
    {rmnccoil[[im]]*Cos[anglec[theta]]*Cos[angle2[theta]], rmnccoil[[im]]*Cos[anglec[theta]]*Sin[angle2[theta]], zmnscoil[[im]]*Sin[anglec[theta]]}
    , {im, 1, imMax}]
, {j, 1, coilsPerHalfPeriod*numHalfPeriodsToPlot}];

plasmapos[theta_, phi_] = Chop[plasmaposTemp[theta, phi], 10^(-chopResolution)]; 
WSpos[theta_, phi_] = Chop[WSposTemp[theta, phi], 10^(-chopResolution)];
Table[coilpos[theta_, j] = Chop[coilposTemp[theta, j], 10^(-chopResolution)], {j, 1, coilsPerHalfPeriod*numHalfPeriodsToPlot}];

Table[
  coilposTable = Table[{{theta, coilposTemp[theta, j][[1]]}, {theta, coilpos[theta, j][[2]]}, {theta, coilpos[theta, j][[3]]}}, {theta, 0, 2*pi, 2*pi/npointsInterpolationCoilPos}];
  coilposInterp = {Interpolation[coilposTable[[All, 1]], Method -> "Spline", InterpolationOrder -> interpolationOrderCoilPos], Interpolation[coilposTable[[All, 2]], Method -> "Spline", InterpolationOrder -> interpolationOrderCoilPos], Interpolation[coilposTable[[All, 3]], Method -> "Spline", InterpolationOrder -> interpolationOrderCoilPos]};
  coilposFunc[theta_, j] = {coilposInterp[[1]][theta], coilposInterp[[2]][theta], coilposInterp[[3]][theta]};
, {j, 1, coilsPerHalfPeriod*numHalfPeriodsToPlot}];

(*Surface of Constant Toroidal Flux*)
plasmaPlot = ParametricPlot3D[plasmapos[theta, phi], {theta, 0, 2*pi}, {phi, 0, 2*pi}, PlotPoints -> {25, plotPointsFig}, MaxRecursion -> maxRecursPlot, Mesh -> None, PlotStyle -> Directive[Opacity[0.7], Specularity[White, 20]]];
(*Plot of the Coils*)
(*coilsPlot = Show[Table[ParametricPlot3D[coilposFunc[theta, j], {theta, 0, 2*pi}, PlotPoints -> 70, MaxRecursion -> 2, PlotStyle -> {Thickness[0.004], RGBColor[56/255, 122/255, 213/255], Directive[Opacity[1.0], Specularity[White, 20]]}], {j, 1, coilsPerHalfPeriod*numHalfPeriodsToPlot}]];*)
(*Winding Surface*)
(*windingPlot=ParametricPlot3D[WSpos[theta,phi],{theta,0,2*pi},{phi,0,2*pi},PlotPoints->{25,plotPointsFig},MaxRecursion->1,Mesh->None,PlotStyle->Directive[LightBlue,Opacity[0.3]]];*)

WSposTable =  ParallelTable[{{theta, phi}, WSpos[theta, phi]}, {theta, 0, 2*pi, 2*pi/npointsInterpolationCoilPos}, {phi, 0, 2*pi, 2*pi/npointsInterpolationCoilPos/1.2}];
WSposInterp = Interpolation[Flatten[WSposTable, 1], Method -> "Spline", InterpolationOrder -> interpolationOrderCoilPos];

regcoilPlotFig = Show[plasmaPlot, Table[
  tangentTheta[theta_]  = D[WSposInterp[theta, phi], theta] /. phi -> coilContours[[j]][theta];
  utangentTheta[theta_] = tangentTheta[theta]/Norm[tangentTheta[theta]];
  tangentPhi[theta_]    = D[WSposInterp[theta, phi], phi] /. phi -> coilContours[[j]][theta];
  utangentPhi[theta_]   = tangentPhi[theta]/Norm[tangentPhi[theta]];
  unormalWSpos[theta_]  = Cross[utangentPhi[theta], utangentTheta[theta]]/Norm[Cross[utangentPhi[theta], utangentTheta[theta]]];
  ParametricPlot3D[{
    coilposFunc[theta, j] + x*utangentPhi[theta], 
    coilposFunc[theta, j] + coilthickness utangentPhi[theta] + x*unormalWSpos[theta], 
    coilposFunc[theta, j] + x*unormalWSpos[theta], 
    coilposFunc[theta, j] + x*utangentPhi[theta] + coilthickness*unormalWSpos[theta]
    }, {theta, 0, 2 Pi}, {x, 0, coilthickness}, PlotPoints -> {plotpointsCoil, 6}, MaxRecursion -> maxRecursPlot, 
    Mesh -> None, PlotStyle -> Directive[Gray, Specularity[White, 20]]]
, {j, 1, coilsPerHalfPeriod*numHalfPeriodsToPlot}], Boxed -> False, ImageSize -> Large, Axes -> False, PlotRange->All];

If[plotRegcoilSENAC==1,
  Export[StringJoin[currentPath,"../data/",ToString[proj],"/Figures/",ToString[proj],"_Regcoil_Plot_Original.pdf"], regcoilPlotFig, ImageSize->ImageSizePlot, ImageResolution->ImageResolutionPlot];
  ,
  Export[StringJoin[currentPath,"../data/",ToString[proj],"/Figures/",ToString[proj],"_Regcoil_Plot_Fit_ordern",ToString[ordern],"_nmodes",ToString[nModes],".pdf"], regcoilPlotFig, ImageSize->ImageSizePlot, ImageResolution->ImageResolutionPlot];
];

(*Finish Evaluation*)
WriteString[$Output,"  SENAC plot REGCOIL took "<>ToString[AbsoluteTime[]-overallTimeStart]<>"s\n"]
CloseKernels[];