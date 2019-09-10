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
nfigsSurf      = ToExpression[$ScriptCommandLine[[33]]];
nPlots         = ToExpression[$ScriptCommandLine[[34]]];
plotRegcoilFit = ToExpression[$ScriptCommandLine[[35]]];
plotRegcoilOriginal  = ToExpression[$ScriptCommandLine[[36]]];
coilsPerHalfPeriod   = ToExpression[$ScriptCommandLine[[37]]];
numHalfPeriodsToPlot = ToExpression[$ScriptCommandLine[[38]]];
plotRegcoilSENAC = ToExpression[$ScriptCommandLine[[39]]];

pi = N[Pi,10];

(*Plot Regcoil Results*)
regcoilOriginalFile = StringJoin[currentPath,"../vmec/",ToString[proj],"/regcoil_out_",ToString[proj],".nc"];
regcoilFitFile = StringJoin[currentPath,"../data/",ToString[proj],"/regcoil_out_",ToString[proj],"_senac.nc"];

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
    ContourPlot[currentpotential[theta, phi] == i, {theta, 0, 2*pi}, {phi, 0, 2*pi}, PlotPoints -> 40, MaxRecursion -> 1][[1, 1, 1]],
   Method -> "Spline", InterpolationOrder -> 2])
  // Quiet, {i, minCurr, maxCurr, (maxCurr - minCurr)/(coilsPerHalfPeriod*2*nfpcoil + 1)}] // Quiet)[[2 ;; coilsPerHalfPeriod*2*nfpcoil + 1]] // Quiet;

plasmaposTemp[theta_, phi_] = {0, 0, 0};
Table[
  angle[theta_, phi_] = xmplasma[[im]]*theta - xnplasma[[im]]*phi;
  plasmaposTemp[theta_, phi_] = plasmaposTemp[theta, phi] + {rmncplasma[[im]]*Cos[angle[theta, phi]]*Cos[phi], rmncplasma[[im]]*Cos[angle[theta, phi]]*Sin[phi], zmnsplasma[[im]]*Sin[angle[theta, phi]]};
  , {im, 1, Length[xmplasma]}];

imMax = Position[rmnccoil, _?(10^-12 < # < 10^-8 &)][[3, 1]];
WSposTemp[theta_, phi_] = {0, 0, 0};
Table[coilposTemp[theta_, j] = {0, 0, 0}, {j, 1, coilsPerHalfPeriod*numHalfPeriodsToPlot}];
Table[
  angle[theta_, phi_] = xmcoil[[im]]*theta - xncoil[[im]]*phi;
  WSposTemp[theta_, phi_] = WSposTemp[theta, phi] + {rmnccoil[[im]]*Cos[angle[theta, phi]]*Cos[phi], rmnccoil[[im]]*Cos[angle[theta, phi]]*Sin[phi], zmnscoil[[im]]*Sin[angle[theta, phi]]};
  Table[
   angle2[theta_] = coilContours[[j]][theta];
   anglec[theta_] = xmcoil[[im]]*theta - xncoil[[im]]*angle2[theta];
   coilposTemp[theta_, j] = coilposTemp[ theta, j] + {rmnccoil[[im]]*Cos[anglec[theta]]*Cos[angle2[theta]], rmnccoil[[im]]*Cos[anglec[theta]]*Sin[angle2[theta]], zmnscoil[[im]]*Sin[anglec[theta]]};
   , {j, 1, coilsPerHalfPeriod*numHalfPeriodsToPlot}];
  , {im, 1, imMax}];

plasmapos[theta_, phi_] = Chop[plasmaposTemp[theta, phi], 10^-5]; 
WSpos[theta_, phi_] = Chop[WSposTemp[theta, phi], 10^-5];
Table[coilpos[theta_, j] = Chop[coilposTemp[theta, j], 10^-5], {j, 1, coilsPerHalfPeriod*numHalfPeriodsToPlot}];

npointsCoilPos = 50;
Table[
  coilposTable = Table[{{theta, coilposTemp[theta, j][[1]]}, {theta, coilpos[theta, j][[2]]}, {theta, coilpos[theta, j][[3]]}}, {theta, 0, 2*pi, 2*pi/npointsCoilPos}];
  coilposInterp = {Interpolation[coilposTable[[All, 1]], Method -> "Spline", InterpolationOrder -> 2], Interpolation[coilposTable[[All, 2]], Method -> "Spline", InterpolationOrder -> 2], Interpolation[coilposTable[[All, 3]], Method -> "Spline", InterpolationOrder -> 2]};
  coilposFunc[theta_, j] = {coilposInterp[[1]][theta], coilposInterp[[2]][theta], coilposInterp[[3]][theta]};
  , {j, 1, coilsPerHalfPeriod*numHalfPeriodsToPlot}];

(*Surface of Constant Toroidal Flux*)
plasmaPlot = ParametricPlot3D[plasmapos[theta, phi], {theta, 0, 2*pi}, {phi, 0, 2*pi}, PlotPoints -> {30, 50}, MaxRecursion -> 2, Mesh -> None, PlotStyle -> Directive[Opacity[0.8], Specularity[White, 20]]];
(*Plot of the Coils*)
coilsPlot = Show[Table[ParametricPlot3D[coilposFunc[theta, j], {theta, 0, 2*pi}, PlotPoints -> 50, MaxRecursion -> 2, PlotStyle -> {Thickness[0.004], RGBColor[56/255, 122/255, 213/255], Directive[Opacity[1.0], Specularity[White, 20]]}], {j, 1, coilsPerHalfPeriod*numHalfPeriodsToPlot}]];
(*Winding Surface*)
(*windingPlot=ParametricPlot3D[WSpos[theta,phi],{theta,0,2*pi},{phi,0,2*pi},PlotPoints->{15,50},MaxRecursion->1,Mesh->None,PlotStyle->Directive[LightBlue,Opacity[0.3]]];*)
fitRegcoil=Show[plasmaPlot,(*windingPlot,*)coilsPlot, Boxed -> False, Axes -> False];

If[plotRegcoilSENAC==1,
  Export[StringJoin[currentPath,"../data/",ToString[proj],"/Figures/",ToString[proj],"_Regcoil_Plot_Original.pdf"], fitRegcoil, ImageSize->ImageSizePlot, ImageResolution->ImageResolutionPlot];
  ,
  Export[StringJoin[currentPath,"../data/",ToString[proj],"/Figures/",ToString[proj],"_Regcoil_Plot_Fit.pdf"], fitRegcoil, ImageSize->ImageSizePlot, ImageResolution->ImageResolutionPlot];
];

(*Finish Evaluation*)
WriteString[$Output,"  SENAC Mathematica took "<>ToString[AbsoluteTime[]-overallTimeStart]<>"s\n"]
CloseKernels[];

(*Along Normal and Binormal Directions*)
(*For Future Work*)
(*
Table[
  coilposj[t_] = coilposFunc[t, j];
  tangentCoil[t_] = FrenetSerretSystem[coilposFunc[t, j], t][[2, 1]];
  normalCoil[t_] = FrenetSerretSystem[coilposFunc[t, j], t][[2, 2]];
  binormalCoil[t_] = Cross[tangentCoil[t], normalCoil[t]];
  thickCoils1[j] = 
   Interpolation[
    Table[{t, 
      coilposj[t] + 
       coilthickness*(normalCoil[t] + binormalCoilTable[t])}, {t, 0, 
      2*pi, 2*pi/30}], Method -> "Spline", InterpolationOrder -> 2];
  thickCoils2[j] = 
   Interpolation[
    Table[{t, 
      coilposj[t] + 
       coilthickness*(normalCoil[t] - binormalCoil[t])}, {t, 0, 2*pi, 
      2*pi/30}], Method -> "Spline", InterpolationOrder -> 2];
  thickCoils3[j] = 
   Interpolation[
    Table[{t, 
      coilposj[t] + 
       coilthickness*(-normalCoil[t] - binormalCoil[t])}, {t, 0, 2*pi,
       2*pi/30}], Method -> "Spline", InterpolationOrder -> 2];
  thickCoils4[j] = 
   Interpolation[
    Table[{t, 
      coilposj[t] + 
       coilthickness*(-normalCoil[t] + binormalCoil[t])}, {t, 0, 2*pi,
       2*pi/30}], Method -> "Spline", InterpolationOrder -> 2];
  , {j, 1, coilsPerHalfPeriod*numHalfPeriodsToPlot}];
  *)