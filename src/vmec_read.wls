(* ::Package:: *)

muDefined       = 0;
VMECdatasets    = Import[vmecOutput, {"Datasets", {"iotaf", "rmnc", "zmns","ns","nfp","phi","bmnc","xm","xn","raxis_cc","zaxis_cs","xn_nyq","xm_nyq"}}];
rmnc            = VMECdatasets[[2]];
zmns            = VMECdatasets[[3]];
nns             = VMECdatasets[[4]];
vmecNFP         = VMECdatasets[[5]];
vmecPSI         = Abs[VMECdatasets[[6]]//Last];
bmnc            = VMECdatasets[[7]];
B0Est           = Max[bmnc];
xm              = VMECdatasets[[8]];
xn              = VMECdatasets[[9]];
vmecOutputRaxis = VMECdatasets[[10]];
vmecOutputZaxis = VMECdatasets[[11]];
xnnyq          = VMECdatasets[[12]];
xmnyq          = VMECdatasets[[13]];

vmecOutputnRAXIS[phi_] = Chop[Sum[vmecOutputRaxis[[nn]] Cos[-vmecNFP (nn - 1) phi], {nn, 1, Length[vmecOutputRaxis]}], 10^-4];
vmecOutputnZAXIS[phi_] = Chop[Sum[vmecOutputZaxis[[nn]] Sin[-vmecNFP (nn - 1) phi], {nn, 1, Length[vmecOutputZaxis]}], 10^-4];
closedcurvVMEC[phi_] = {vmecOutputnRAXIS[phi] Cos[phi], vmecOutputnRAXIS[phi] Sin[phi], vmecOutputnZAXIS[phi]};
If[nsurfaces>nns,WriteString[$Output, "    The number of surfaces asked is larger than the number of surfaces in VMEC file. Proceeding with "<>ToString[nns]<>" surfaces.\n"];nsurfaces=nns];
iotaAxisOut=Array[Null,nsurfaces];
j=0;Table[j=j+1;
	RBCVMEC[theta_, phi_, j] = Chop[Sum[rmnc[[Floor[iradius], i]]*Cos[xm[[i]] theta - xn[[i]] phi], {i, 1, Length[xn]}], 10^-8];
	ZBSVMEC[theta_, phi_, j] = Chop[Sum[zmns[[Floor[iradius], i]]*Sin[xm[[i]] theta - xn[[i]] phi], {i, 1, Length[xn]}], 10^-8];
	BVMEC  [theta_, phi_, j] = Chop[Sum[bmnc[[Floor[iradius], i]]*Cos[xmnyq[[i]] theta - xnnyq[[i]] phi], {i, 1, Length[xnnyq]}], 10^-8];
	iotaAxisOut[[j]] = VMECdatasets[[1]][[Floor[iradius]]];
, {iradius,2,nns,(nns-2)/(nsurfaces-1)}];