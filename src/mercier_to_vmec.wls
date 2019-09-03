(* ::Package:: *)

maxm           = ToExpression[$ScriptCommandLine[[19]]];
maxn           = ToExpression[$ScriptCommandLine[[20]]];
maxRecursTheta = ToExpression[$ScriptCommandLine[[21]]];
maxRecursPhi   = ToExpression[$ScriptCommandLine[[22]]];

timestart = AbsoluteTime[]; 
Print["Starting RBC Evaluation"];
RBCmn = Quiet[ParallelTable[
     GaussLegendreQuadrature[
       GaussLegendreQuadrature[
        RBCFit[phi, theta] Cos[mm theta - vmecNFP nn phi], {phi, 0, 2*pi}, maxRecursPhi, 20], {theta, 0, 2*pi}, maxRecursTheta, 10]/If[mm == 0 && nn == 0, 4*pi^2, 2*pi^2], {mm, 0, maxm}, {nn, -maxn, maxn}]];
Print["RBCmn from fit took "<>ToString[AbsoluteTime[]-timestart]<>"s"];

timestart = AbsoluteTime[];
Print["Starting ZBS Evaluation"];
ZBSmn =
  Quiet[ParallelTable[
     GaussLegendreQuadrature[
       GaussLegendreQuadrature[
        ZBSFit[phi, theta] Sin[mm theta - vmecNFP nn phi], {phi, 0, 2*pi}, maxRecursPhi, 20], {theta, 0, 2*pi}, maxRecursTheta, 10]/If[mm == 0 && nn == 0, 1, 2*pi^2], {mm, 0, maxm}, {nn, -maxn, maxn}]];
Print["ZBSmn from fit took "<>ToString[AbsoluteTime[]-timestart]<>"s"];

(*Print to VMEC input file*)
printv = "";nl = "\n";
printv = printv <> "NFP = "<>ToString[vmecNFP];
printv = printv <> nl;
printv = printv <> "PHIEDGE = "<>ToString[vmecPSI];
printv = printv <> nl;
rbcs[n_, m_] := "RBC(" <> ToString[n] <> "," <> ToString[m] <> ")=";
zbcs[n_, m_] := "  ZBS(" <> ToString[n] <> "," <> ToString[m] <> ")=";
printv = printv <> "RAXIS = ";
Do[printv = printv <> " " <> ToString[DecimalForm[vmecOutputRaxis[[i]]]], {i, 1, Length[vmecOutputRaxis]}];
printv = printv <> nl;
printv = printv <> "ZAXIS = ";
Do[printv = printv <> " " <> ToString[DecimalForm[vmecOutputZaxis[[i]]]], {i, 1, Length[vmecOutputZaxis]}];
printv = printv <> nl;
Do[If[m == 0 && n < 0 || 
    Chop[RBCmn[[m + 1, n + maxn + 1]], 10^-8] == 0 && 
     Chop[ZBSmn[[m + 1, n + maxn + 1]], 10^-8] == 0, Null,
   printv = printv <> rbcs[n, m] <> 
     ToString[DecimalForm[RBCmn[[m + 1, n + maxn + 1]]]] <> 
     zbcs[n, m] <> 
     ToString[DecimalForm[ZBSmn[[m + 1, n + maxn + 1]]]] <> nl
   ], {m, 0, maxm}, {n, -maxn, maxn}];
printv = printv <> "/" <> nl;

CopyFile[vmecInput,StringJoin[currentPath,"data/",ToString[proj],"/",ToString[proj],"_SENACtoVMEC_input.txt"], OverwriteTarget -> True];
newVMECfile = OpenAppend[StringJoin[currentPath,"data/",ToString[proj],"/",ToString[proj],"_SENACtoVMEC_input.txt"]];
WriteString[newVMECfile, printv];
Close[newVMECfile];