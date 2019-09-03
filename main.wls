(* ::Package:: *)

SetOptions[$Output, FormatType->OutputForm];

currentPath=DirectoryName @ $InputFileName;
SetDirectory[currentPath];
readVMEC = ToExpression[$ScriptCommandLine[[2]]];
outputToVMEC = ToExpression[$ScriptCommandLine[[17]]];
pi = N[Pi,10];

(*Fast Numerical Integration*)
GaussLegendreQuadrature[f_, {x_, a_, b_}, n_Integer: 10, prec_: MachinePrecision] := 
  Module[{nodes, weights}, {nodes, weights} = Most[NIntegrate`GaussRuleData[n, prec]];
   (b - a) weights.Map[Function[x, f], Rescale[nodes, {0, 1}, {a, b}]]];

LaunchKernels[];

If[readVMEC == 1,
    WriteString[$Output, "  VMEC read\n"]
    Get[StringJoin[currentPath,"src/vmec_read.wls"]];
    WriteString[$Output, "  Fitting VMEC to Mercier's theta and rho\n"]
    Get[StringJoin[currentPath,"src/fit_to_mercier.wls"]];
    If[outputToVMEC == 1,
      WriteString[$Output, "  Mercier to VMEC\n"]
      Get[StringJoin[currentPath,"src/mercier_to_vmec.wls"]];
    ]
]
WriteString[$Output,"  SENAC took "<>ToString[TimeUsed[]]<>"s\n"]

CloseKernels[];

(*Pause[5]*)