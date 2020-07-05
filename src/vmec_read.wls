(* ::Package:: *)

muDefined       = 0;
VMECdatasets    = Import[vmecOutput, {"Datasets", {"iotaf", "rmnc", "zmns","ns","nfp","phi","bmnc","xm","xn","raxis_cc","zaxis_cs","xn_nyq","xm_nyq","ctor","am","ac","ai","DWell","gmnc"}}];
rmnc            = VMECdatasets[[2]];
zmns            = VMECdatasets[[3]];
nns             = VMECdatasets[[4]]-1;
vmecNFP         = VMECdatasets[[5]];
psiiVMEC        = VMECdatasets[[6]];
bmnc            = VMECdatasets[[7]];
B0Est           = Max[bmnc];
xm              = VMECdatasets[[8]];
xn              = VMECdatasets[[9]];
vmecOutputRaxis = VMECdatasets[[10]];
vmecOutputZaxis = VMECdatasets[[11]];
xnnyq           = VMECdatasets[[12]];
xmnyq           = VMECdatasets[[13]];
ctor            = VMECdatasets[[14]];
pprime0VMEC     = mu0const*VMECdatasets[[15]][[2]];
torCurVMEC      = ctor*mu0const*VMECdatasets[[16]][[1]]/Last[psiiVMEC];
(*ai              = VMECdatasets[[17]];*)
DWellTemp       = VMECdatasets[[18]];
gmnc            = VMECdatasets[[19]];

surface1=15;

vmecOutputnRAXIS[phi_] = Chop[Sum[vmecOutputRaxis[[nn]] Cos[-vmecNFP (nn - 1) phi], {nn, 1, Length[vmecOutputRaxis]}], 10^-6];
vmecOutputnZAXIS[phi_] = Chop[Sum[vmecOutputZaxis[[nn]] Sin[-vmecNFP (nn - 1) phi], {nn, 1, Length[vmecOutputZaxis]}], 10^-6];
closedcurvVMEC[phi_] = Chop[{vmecOutputnRAXIS[phi] Cos[phi], vmecOutputnRAXIS[phi] Sin[phi], vmecOutputnZAXIS[phi]}, 10^-6];
If[nsurfaces>nns,WriteString[$Output, "    The number of surfaces asked is larger than the number of surfaces in VMEC file. Proceeding with "<>ToString[nns]<>" surfaces.\n"];nsurfaces=nns];
iotaAxisOut=Array[Null,nsurfaces];vmecPSI=Array[Null,nsurfaces];DWellVMEC=Array[Null,nsurfaces];
j=0;Table[j=j+1;
	RBCVMEC[theta_, phi_, j] = Chop[Sum[rmnc[[Floor[iradius], i]]*Cos[xm[[i]] theta - xn[[i]] phi], {i, 1, Length[xn]}], 10^-6];
	ZBSVMEC[theta_, phi_, j] = Chop[Sum[zmns[[Floor[iradius], i]]*Sin[xm[[i]] theta - xn[[i]] phi], {i, 1, Length[xn]}], 10^-6];
	BVMEC  [theta_, phi_, j] = Chop[Sum[bmnc[[Floor[iradius], i]]*Cos[xmnyq[[i]] theta - xnnyq[[i]] phi], {i, 1, Length[xnnyq]}], 10^-6];
	iotaAxisOut[[j]]         = VMECdatasets[[1]][[Floor[iradius]]];
	vmecPSI[[j]]             = Abs[psiiVMEC[[Floor[iradius]]]];
	DWellVMEC[[j]]           = DWellTemp[[Floor[iradius]]];
, {iradius,surface1,If[nsurfaces==1,surface1,nns],If[nsurfaces==1,surface1,(nns-1-surface1)/(nsurfaces-1)]}];
psiTable=vmecPSI/(vmecPSI//Last);
WriteString[$Output, "  Psi on the surfaces = "<>ToString[vmecPSI]<>".\n"];
WriteString[$Output, "  psi/psi_b           = "<>ToString[psiTable]<>".\n"];