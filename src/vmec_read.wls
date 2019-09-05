(* ::Package:: *)

muDefined=0;
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