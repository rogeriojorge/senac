(* ::Package:: *)

vmecInput      = $ScriptCommandLine[[4]];
vmecOutput     = $ScriptCommandLine[[5]];

If[FileExistsQ[vmecOutput],
    VMECdatasets = Import[vmecOutput, {"Datasets", {"iotaf", "rmnc", "zmns","ns","nfp","phi","bmnc","xm","xn","raxis_cc","zaxis_cs"}}];
    iotaAxisOut = VMECdatasets[[1]][[1]];
    rmnc = VMECdatasets[[2]];
    zmns = VMECdatasets[[3]];
    nns  = VMECdatasets[[4]];
    vmecNFP = VMECdatasets[[5]];
 	vmecPSI = Abs[VMECdatasets[[6]]//Last];
 	B0Est = Max[VMECdatasets[[7]]];
 	xm = VMECdatasets[[8]];
 	xn = VMECdatasets[[9]];
 	vmecOutputRaxis = VMECdatasets[[10]];
 	vmecOutputZaxis = VMECdatasets[[11]];
    iradius = nns;
 	vmecOutputnRAXIS[phi_] = Chop[Sum[vmecOutputRaxis[[nn]] Cos[-vmecNFP (nn - 1) phi], {nn, 1, Length[vmecOutputRaxis]}], 10^-4];
 	vmecOutputnZAXIS[phi_] = Chop[Sum[vmecOutputZaxis[[nn]] Sin[-vmecNFP (nn - 1) phi], {nn, 1, Length[vmecOutputZaxis]}], 10^-4];
 	closedcurvVMEC[phi_] = {vmecOutputnRAXIS[phi] Cos[phi], vmecOutputnRAXIS[phi] Sin[phi], vmecOutputnZAXIS[phi]};
 	RBCVMEC[theta_, phi_] = Sum[rmnc[[iradius, i]] Cos[xm[[i]] theta - xn[[i]] phi], {i, 1, Length[xn]}];
 	ZBSVMEC[theta_, phi_] = Sum[zmns[[iradius, i]] Sin[xm[[i]] theta - xn[[i]] phi], {i, 1, Length[xn]}];
,
 	WriteString[$Output, "  No VMEC output file\n"]
	CloseKernels[];
	Abort[];
	Run["exit 1"];
]

(*VMEC Axis Frenet-Serret Frame*)
Off[Simplify::time]
curvVMEC[t_] = Chop[Norm[Cross[closedcurvVMEC'[t], closedcurvVMEC''[t]]]/Norm[closedcurvVMEC'[t]]^3, 10^-6]//Re;
torsVMEC[t_] = Chop[Dot[Cross[closedcurvVMEC'[t], closedcurvVMEC''[t]], closedcurvVMEC'''[t]]/Norm[Cross[closedcurvVMEC'[t], closedcurvVMEC''[t]]]^2, 10^-6];
sprimeVMEC[t_] = Chop[Simplify[Chop[ComplexExpand[Norm[closedcurvVMEC'[t]]], 10^-6], TimeConstraint -> 0.01], 10^-6];
b0VMEC[t_] = Simplify[Chop[closedcurvVMEC'[t]/sprimeVMEC[t], 10^-6], TimeConstraint -> 0.01];
k0VMEC[t_] = Simplify[Chop[b0VMEC'[t]/sprimeVMEC[t], 10^-6], TimeConstraint -> 0.01]/curvVMEC[t]//Re;
t0VMEC[t_] = Chop[Cross[b0VMEC[t], k0VMEC[t]], 10^-6]//Re;

(*Flux Surface of VMEC*)
FluxSurfaceVMEC[theta_, phi_] = Simplify[Chop[{RBCVMEC[theta, phi] Cos[phi], RBCVMEC[theta, phi] Sin[phi], ZBSVMEC[theta, phi]}, 10^-6], TimeConstraint -> 0.01]//Re;