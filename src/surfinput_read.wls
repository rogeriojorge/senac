(* ::Package:: *)

If[FileExistsQ[surfInput],
	(*Calculate Axis Parameters*)
	vmecNFP = StringCases[StringDelete[FindList[surfInput, "NFP"][[1]], " "], x : NumberString :> ToExpression[x]][[1]];
	If[quasisymmetry==1,vmecPSI={phiedge},vmecPSI = {Abs[StringCases[StringDelete[FindList[surfInput, "PHIEDGE"][[1]], " "], x : NumberString :> ToExpression[x]][[1]]]}];
	vmecRAXIS = StringSplit[FindList[surfInput, "RAXIS"][[1]]];
	vmecZAXIS = StringSplit[FindList[surfInput, "ZAXIS"][[1]]];
	torCurVMEC = StringCases[StringDelete[FindList[surfInput, "torCurr"][[1]], " "], x : NumberString :> ToExpression[x]][[1]];
	pprime0VMEC = StringCases[StringDelete[FindList[surfInput, "pressGr"][[1]], " "], x : NumberString :> ToExpression[x]][[1]];
	nRAXIS[phi_] = Sum[Internal`StringToDouble@vmecRAXIS[[nn + 3]]*Cos[vmecNFP nn phi], {nn, 0, Dimensions[vmecRAXIS][[1]] - 3}];
	nZAXIS[phi_] = Sum[Internal`StringToDouble@vmecZAXIS[[nn + 3]]*Sin[vmecNFP nn phi], {nn, 0, Dimensions[vmecZAXIS][[1]] - 3}];
	closedcurvVMEC[phi_] = {nRAXIS[phi] Cos[phi], nRAXIS[phi] Sin[phi], nZAXIS[phi]};
	(*Check if surface parameters or Fourier components are given*)
	If[SameQ[Find[OpenRead[surfInput],"mu"], EndOfFile], 
		muDefined=0; WriteString[$Output, " Surface Fourier Components Defined in file..."];
		tmaxN = 40;
		vmecRBC = ConstantArray[0, {tmaxN, tmaxN}];
		vmecZBS = ConstantArray[0, {tmaxN, tmaxN}];
		vmecSurf = DeleteCases[Flatten[Table[StringSplit[FindList[surfInput, "RBC(" <> ToString[i] <> "," <> ToString[j] <> ")"][[1]]]//Quiet, {i, -tmaxN, tmaxN}, {j, 0, tmaxN}]], StringSplit[{}[[1]]] // Quiet];
		Table[
			ToExpression[StringJoin[StringReplace[vmecSurf[[i]], {"RBC" -> "vmecRBC", "ZBS" -> "vmecZBS", "(" -> "[[" <> ToString[Floor[tmaxN/2]] <> "+", ")" -> "+" <> ToString[Floor[tmaxN/2]] <> "]]", "e" -> "\[Times]10^("}], ")"]],
		{i, 1, Length[vmecSurf]}];
		nsurfaces=1;nns=1;
		RBCVMEC[theta_, phi_, 1] = Chop[Sum[vmecRBC[[mm, nn]] Cos[(mm - Floor[tmaxN/2]) theta - vmecNFP (nn - Floor[tmaxN/2]) phi], {mm, 1, Length[vmecRBC]}, {nn, 1, Length[vmecRBC]}], 10^(-chopResolution)];
		ZBSVMEC[theta_, phi_, 1] = Chop[Sum[vmecZBS[[mm, nn]] Sin[(mm - Floor[tmaxN/2]) theta - vmecNFP (nn - Floor[tmaxN/2]) phi], {mm, 1, Length[vmecRBC]}, {nn, 1, Length[vmecRBC]}], 10^(-chopResolution)];
		,
		muDefined=1; WriteString[$Output, " Mercier's Parameters Defined in file..."];
		muTable    = StringSplit[FindList[surfInput, "mu"]][[1]];
		deltaTable = StringSplit[FindList[surfInput, "delta"][[1]]];
		mu[phi_] = Sum[Internal`StringToDouble@muTable[[nn + 3]] Cos[vmecNFP nn phi], {nn, 0, Dimensions[muTable][[1]] - 3}];
		delta[phi_] = Internal`StringToDouble@deltaTable[[3]]*phi/2+Sum[Internal`StringToDouble@deltaTable[[nn + 3]] Sin[vmecNFP nn phi], {nn, 1, Dimensions[deltaTable][[1]] - 5}];
		If[quasisymmetry==1,B0[phi_]=Bzero,B0Table    = StringSplit[FindList[surfInput, "B0"][[1]]];B0[phi_] = Sum[Internal`StringToDouble@B0Table[[nn + 3]] Cos[vmecNFP nn phi], {nn, 0, Dimensions[B0Table][[1]] - 3}]];
		psi31Table    = StringSplit[FindList[surfInput, "psi31"]][[1]];
		psi32Table    = StringSplit[FindList[surfInput, "psi32"]][[1]];
		psi33Table    = StringSplit[FindList[surfInput, "psi33"]][[1]];
		psi34Table    = StringSplit[FindList[surfInput, "psi34"]][[1]];
		psi31[phi_] = Sum[Internal`StringToDouble@psi31Table[[nn + 3]] Cos[vmecNFP nn phi], {nn, 0, Dimensions[psi31Table][[1]] - 3}];
		psi32[phi_] = Sum[Internal`StringToDouble@psi32Table[[nn + 3]] Sin[vmecNFP nn phi], {nn, 0, Dimensions[psi32Table][[1]] - 3}];
		psi33[phi_] = Sum[Internal`StringToDouble@psi33Table[[nn + 3]] Cos[vmecNFP nn phi], {nn, 0, Dimensions[psi33Table][[1]] - 3}];
		psi34[phi_] = Sum[Internal`StringToDouble@psi34Table[[nn + 3]] Sin[vmecNFP nn phi], {nn, 0, Dimensions[psi34Table][[1]] - 3}];
	];
	,
 	WriteString[$Output, "  No surf_input file\n"]
	CloseKernels[];
	Abort[];
	Run["exit 1"];
];