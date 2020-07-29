(* ::Package:: *)

Clear[psin, rhon, rho, rhoSeries, positiveRhoPos]

(*Order 0 and 1 are zero*)
psin[theta_, phi_, 0] = 0; 
psin[theta_, phi_, 1] = 0;
rhon[theta_, phi_, 0] = 0;

(*Expression for rho of psi*)
psiSeries[rhoo_, theta_, phi_, orderN_] = Sum[psin[theta, phi, n]*rhoo^n, {n, 0, orderN}];
rhoSeries[psii_, theta_, phi_, orderN_] = Sum[rhon[theta, phi, n]*psii^(n/2), {n, 0, orderN - 1}];
rho[theta_, phi_] = rhoSeries[psii, theta, phi, ordern] /. 
  Solve[CoefficientList[Assuming[{rhoo > 0, psii > 0}, 
    Normal[Series[
      psiSeries[rhoSeries[psii^2, theta, phi, ordern], theta, phi, ordern]
    , {psii, 0, ordern + 1}]]] - psii^2, psii][[3 ;; ordern + 1]] == 0, 
  Flatten[Table[{rhon[theta, phi, n]}, {n, 1, ordern - 1}]]];
positiveRhoPos = Position[Coefficient[rho[1, 1]/.psii->1, 1/Sqrt[psin[1, 1, 2]]], 1][[1, 1]];
rho[theta_, phi_] = rho[theta, phi][[positiveRhoPos]]/.psii->PSI;

(*Order n psi*)
Table[psin[theta_, phi_, n] = ((B0VMEC)^(n/2)) * Sum[If[EvenQ[n], 
  Sum[psic[n, 2 i  , j]*Cos[2*i*theta       - vmecNFP*j*phi]+0*psis[n, 2 i,   j]*Sin[2*i*theta       - vmecNFP*j*phi], {j, 0, nModes}], 
  Sum[psic[n, 2 i+1, j]*Cos[(2*i + 1)*theta - vmecNFP*j*phi]+0*psis[n, 2 i+1, j]*Sin[(2*i + 1)*theta - vmecNFP*j*phi], {j, 0, nModes}]]
,{i, 0, Floor[n/2]}], {n, 0, ordern}];

(*Order 0 and 1 are zero*)
psin[theta_, phi_, 0] = 0; 
psin[theta_, phi_, 1] = 0;

(**)
deltal0=vmecNFP;
deltac0=0;

Switch[proj,
  "TJII", deltal0=0,
  "ARIES", deltal0=3,
  "NCSX", deltal0=3,
  "ESTELL", deltac0=pi/2; deltal0=2
];
(***)

(*Order 2 is an ellipse*)
muVMEC = muc[0] + Sum[muc[i]*Cos[vmecNFP i phi], {i, 1, nModes}];
deltaVMEC = -deltal0*phi/2 + deltac0 + Sum[deltas[i]*Sin[vmecNFP i phi], {i, 1, nModes}];
B0VMEC = B0c[0] + Sum[B0c[i]*Cos[vmecNFP i phi], {i, 1, nModes}];
psin[theta_, phi_, 2] = pi*(B0VMEC)*(1 + muVMEC*Cos[2*(theta + deltaVMEC)])/Sqrt[1 - (muVMEC)^2];