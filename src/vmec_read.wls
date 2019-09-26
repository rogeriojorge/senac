(* ::Package:: *)

muDefined       = 0;
VMECdatasets    = Import[vmecOutput, {"Datasets", {"iotaf", "rmnc", "zmns","ns","nfp","phi","bmnc","xm","xn","raxis_cc","zaxis_cs","xn_nyq","xm_nyq"}}];
rmnc            = VMECdatasets[[2]];
zmns            = VMECdatasets[[3]];
nns             = VMECdatasets[[4]]-1;
vmecNFP         = VMECdatasets[[5]];
bmnc            = VMECdatasets[[7]];
B0Est           = Max[bmnc];
xm              = VMECdatasets[[8]];
xn              = VMECdatasets[[9]];
vmecOutputRaxis = VMECdatasets[[10]];
vmecOutputZaxis = VMECdatasets[[11]];
xnnyq          = VMECdatasets[[12]];
xmnyq          = VMECdatasets[[13]];

vmecOutputnRAXIS[phi_] = Chop[Sum[vmecOutputRaxis[[nn]] Cos[-vmecNFP (nn - 1) phi], {nn, 1, Length[vmecOutputRaxis]}], 10^-7];
vmecOutputnZAXIS[phi_] = Chop[Sum[vmecOutputZaxis[[nn]] Sin[-vmecNFP (nn - 1) phi], {nn, 1, Length[vmecOutputZaxis]}], 10^-7];
closedcurvVMEC[phi_] = Chop[{vmecOutputnRAXIS[phi] Cos[phi], vmecOutputnRAXIS[phi] Sin[phi], vmecOutputnZAXIS[phi]}, 10^-7];
If[nsurfaces>nns,WriteString[$Output, "    The number of surfaces asked is larger than the number of surfaces in VMEC file. Proceeding with "<>ToString[nns]<>" surfaces.\n"];nsurfaces=nns];
iotaAxisOut=Array[Null,nsurfaces];vmecPSI=Array[Null,nsurfaces];
j=0;Table[j=j+1;
	RBCVMEC[theta_, phi_, j] = Chop[Sum[rmnc[[Floor[iradius], i]]*Cos[xm[[i]] theta - xn[[i]] phi], {i, 1, Length[xn]}], 10^-7];
	ZBSVMEC[theta_, phi_, j] = Chop[Sum[zmns[[Floor[iradius], i]]*Sin[xm[[i]] theta - xn[[i]] phi], {i, 1, Length[xn]}], 10^-7];
	BVMEC  [theta_, phi_, j] = Chop[Sum[bmnc[[Floor[iradius], i]]*Cos[xmnyq[[i]] theta - xnnyq[[i]] phi], {i, 1, Length[xnnyq]}], 10^-7];
	iotaAxisOut[[j]]         = VMECdatasets[[1]][[Floor[iradius]]];
	vmecPSI[[j]]             = Abs[VMECdatasets[[6]][[Floor[iradius]]]];
, {iradius,4,If[nsurfaces==1,4,nns],If[nsurfaces==1,1,(nns-4)/(nsurfaces-1)]}];
WriteString[$Output, "Psi on the surfaces = "<>ToString[vmecPSI]<>".\n"];
WriteString[$Output, "Sqrt(psi/psi_b)     = "<>ToString[Sqrt[vmecPSI/(vmecPSI//Last)]]<>".\n"];