(* ::Package:: *)

timestart = AbsoluteTime[]; WriteString[$Output, "    Creating RBC and ZBS interpolation functions... "];
RBCtoVMECTable=ParallelTable[{{theta, phi}, Cos[phi]*FluxSurfacenFitVMEC[theta, phi, nsurfaces][[1]] + Sin[phi]*FluxSurfacenFitVMEC[theta, phi, nsurfaces][[2]]}, {theta, 0, 2*pi, 2*pi/maxRecursTheta/1.5}, {phi, 0, 2*pi, 2*pi/maxRecursPhi/1.5}];
ZBStoVMECTable=ParallelTable[{{theta, phi}, FluxSurfacenFitVMEC[theta, phi, nsurfaces][[3]]}, {theta, 0, 2*pi, 2*pi/maxRecursTheta/1.2}, {phi, 0, 2*pi, 2*pi/maxRecursPhi/1.2}];
RBCtoVMECInterp=Interpolation[Flatten[RBCtoVMECTable, 1], Method -> "Spline", InterpolationOrder -> 2];
ZBStoVMECInterp=Interpolation[Flatten[ZBStoVMECTable, 1], Method -> "Spline", InterpolationOrder -> 2];
WriteString[$Output, "Successful in "<>ToString[AbsoluteTime[]-timestart]<>"s."];

timestart = AbsoluteTime[]; WriteString[$Output, "\n    RBCmn Evaluation... "];
RBCmn = Quiet[Table[
  rbstemp[theta_,phi_] = RBCtoVMECInterp[theta,phi]*Cos[mm*theta - vmecNFP*nn*phi];(*Interpolation[Flatten[ParallelTable[{{theta, phi}, RBCtoVMEC[theta, phi, mm, nn]}, {theta, 0, 2*pi, 2*pi/maxRecursTheta/1.1}, {phi, 0, 2*pi, 2*pi/maxRecursPhi/1.1}], 1], Method -> "Spline", InterpolationOrder -> 3];*)
  GaussLegendreQuadrature[GaussLegendreQuadrature[
  rbstemp[theta, phi], {phi, 0, 2*pi}, maxRecursPhi, 10], {theta, 0, 2*pi}, maxRecursTheta, 10]/If[mm == 0 && nn == 0, 4*pi^2, 2*pi^2],
{mm, 0, maxm}, {nn, -maxn, maxn}]];
WriteString[$Output, "Successful in "<>ToString[AbsoluteTime[]-timestart]<>"s."];

timestart = AbsoluteTime[]; WriteString[$Output, "\n    ZBSmn Evaluation... "];
ZBSmn = Quiet[Table[
  zbstemp[theta_,phi_] = ZBStoVMECInterp[theta,phi]*Sin[mm*theta - vmecNFP*nn*phi];(*Interpolation[Flatten[ParallelTable[{{theta, phi}, ZBStoVMEC[theta, phi, mm, nn]}, {theta, 0, 2*pi, 2*pi/maxRecursTheta/1.1}, {phi, 0, 2*pi, 2*pi/maxRecursPhi/1.1}], 1], Method -> "Spline", InterpolationOrder -> 3];*)
  GaussLegendreQuadrature[GaussLegendreQuadrature[
  zbstemp[theta, phi], {phi, 0, 2*pi}, maxRecursPhi, 10], {theta, 0, 2*pi}, maxRecursTheta, 10]/If[mm == 0 && nn == 0, 1, 2*pi^2],
{mm, 0, maxm}, {nn, -maxn, maxn}]];
WriteString[$Output, "Successful in "<>ToString[AbsoluteTime[]-timestart]<>"s.\n"];

(*Print to VMEC input file*)
printv = "";nl = "\n";
printv = printv <> " NFP = "<>ToString[vmecNFP];
printv = printv <> nl;
printv = printv <> " PHIEDGE = "<>ToString[DecimalForm[vmecPSI[[nsurfaces]]]];
printv = printv <> nl;
rbcs[n_, m_] := " RBC(" <> ToString[n] <> "," <> ToString[m] <> ")=";
zbcs[n_, m_] := "  ZBS(" <> ToString[n] <> "," <> ToString[m] <> ")=";
printv = printv <> " RAXIS = ";
If[FileExistsQ[vmecOutput],
  Do[printv = printv <> " " <> ToString[DecimalForm[vmecOutputRaxis[[i]]]], {i, 1, Length[vmecOutputRaxis]}];
  printv = printv <> nl;
  printv = printv <> " ZAXIS = ";
  Do[printv = printv <> " " <> ToString[DecimalForm[vmecOutputZaxis[[i]]]], {i, 1, Length[vmecOutputZaxis]}];
  printv = printv <> nl;
,
  Do[printv = printv <> " " <> ToString[DecimalForm[Internal`StringToMReal@vmecRAXIS[[i + 3]]]], {i, 0, Dimensions[vmecRAXIS][[1]] - 3}];
  printv = printv <> nl;
  printv = printv <> " ZAXIS = ";
  Do[printv = printv <> " " <> ToString[DecimalForm[Internal`StringToMReal@vmecZAXIS[[i + 3]]]], {i, 0, Dimensions[vmecZAXIS][[1]] - 3}];
  printv = printv <> nl;
]

Do[If[m == 0 && n < 0 || Chop[RBCmn[[m + 1, n + maxn + 1]], 10^(-chopResolution)] == 0 && Chop[ZBSmn[[m + 1, n + maxn + 1]], 10^(-chopResolution)] == 0,
  Null,
  printv = printv <>
    rbcs[n, m] <> 
    ToString[DecimalForm[RBCmn[[m + 1, n + maxn + 1]]]] <>
    zbcs[n, m] <> 
    ToString[DecimalForm[ZBSmn[[m + 1, n + maxn + 1]]]] <>
    nl
  ],{m, 0, maxm}, {n, -maxn, maxn}
];
printv = printv <> "/" <> nl;

CopyFile[vmecInput,StringJoin[currentPath,"data/",ToString[proj],"/input.",ToString[proj],"_SENACtoVMEC_ordern",ToString[ordern],"_nmodes",ToString[nModes],".txt"], OverwriteTarget -> True];
newVMECfile = OpenAppend[StringJoin[currentPath,"data/",ToString[proj],"/input.",ToString[proj],"_SENACtoVMEC_ordern",ToString[ordern],"_nmodes",ToString[nModes],".txt"]];
WriteString[newVMECfile, printv];
Close[newVMECfile];