(* ::Package:: *)

SetOptions[$Output, FormatType->OutputForm];

currentPath=DirectoryName @ $InputFileName;
SetDirectory[currentPath];
readVMEC = ToExpression[$ScriptCommandLine[[2]]];
pi = N[Pi,10];

LaunchKernels[];

If[readVMEC == 1,
    WriteString[$Output, "  VMEC read\n"]
    Get[StringJoin[currentPath,"src/vmec_read.wls"]];
    Get[StringJoin[currentPath,"src/fit_to_mercier.wls"]];
    Get[StringJoin[currentPath,"src/mercier_to_vmec.wls"]];
]
WriteString[$Output,"  SENAC took "<>ToString[TimeUsed[]]<>"s\n"]

CloseKernels[];

(*Pause[5]*)