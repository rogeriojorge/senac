eta[phi_]    = ArcTanh[mu[phi]];
phi33c[phi_] = 0;
phi33s[phi_] = 0;

psiMatA[phi_]=
{
{-2*mu[phi]*Tan[2*delta[phi]]*(delta'[phi]/sprimeVMEC[phi]-torsVMEC[phi])+eta'[phi]/sprimeVMEC[phi]                                  ,torsVMEC[phi]*Sec[2*delta[phi]]-2*mu[phi]*(delta'[phi]/sprimeVMEC[phi]-torsVMEC[phi])-Tan[2*delta[phi]]*eta'[phi]/sprimeVMEC[phi] ,3/2 (-2  mu[phi] Tan[2 delta[phi]] (delta'[phi]/sprimeVMEC[phi]-torsVMEC[phi])+eta'[phi]/sprimeVMEC[phi]) ,-(3/2) (2  mu[phi] (delta'[phi]/sprimeVMEC[phi]-torsVMEC[phi])+Tan[2 delta[phi]] eta'[phi]/sprimeVMEC[phi])},
{-torsVMEC[phi] Sec[2 delta[phi]]-2  mu[phi] (delta'[phi]/sprimeVMEC[phi]-torsVMEC[phi])-Tan[2 delta[phi]] eta'[phi]/sprimeVMEC[phi] ,2  mu[phi] Tan[2 delta[phi]] (delta'[phi]/sprimeVMEC[phi]-torsVMEC[phi])-eta'[phi]/sprimeVMEC[phi]                                ,3/2 (2  mu[phi] (delta'[phi]/sprimeVMEC[phi]-torsVMEC[phi])+Tan[2 delta[phi]] eta'[phi]/sprimeVMEC[phi])  ,3/2 (-2  mu[phi] Tan[2 delta[phi]] (delta'[phi]/sprimeVMEC[phi]-torsVMEC[phi])+eta'[phi]/sprimeVMEC[phi])},
{1/2 (-2  mu[phi] Tan[2 delta[phi]] (delta'[phi]/sprimeVMEC[phi]-torsVMEC[phi])+eta'[phi]/sprimeVMEC[phi])                           ,1/2 (2  mu[phi] (delta'[phi]/sprimeVMEC[phi]-torsVMEC[phi])+Tan[2 delta[phi]] eta'[phi]/sprimeVMEC[phi])                          ,0                                                                                                         ,3 torsVMEC[phi] Sec[2 delta[phi]]},
{1/2 (-2  mu[phi] (delta'[phi]/sprimeVMEC[phi]-torsVMEC[phi])-Tan[2 delta[phi]] eta'[phi]/sprimeVMEC[phi])                           ,1/2 (-2  mu[phi] Tan[2 delta[phi]] (delta'[phi]/sprimeVMEC[phi]-torsVMEC[phi])+eta'[phi]/sprimeVMEC[phi])                         ,-3 torsVMEC[phi] Sec[2 delta[phi]]                                                                        ,0}
};

psiMatB[phi_]=
{
1/8 pi B0[phi] (6 Cosh[eta[phi]] curvVMEC'[phi]/sprimeVMEC[phi]+2 Sin[2 delta[phi]] Sinh[eta[phi]] (24 phi33s[phi]+curvVMEC[phi] (torsVMEC[phi]-3 (torsVMEC[phi]+(delta'[phi]/sprimeVMEC[phi]-torsVMEC[phi]))))+2 curvVMEC[phi] Sinh[eta[phi]] eta'[phi]/sprimeVMEC[phi]+Cos[2 delta[phi]] (4 Sinh[eta[phi]] (-12 phi33c[phi]+curvVMEC'[phi]/sprimeVMEC[phi])+3 Cosh[eta[phi]] curvVMEC[phi] eta'[phi]/sprimeVMEC[phi]))-pi curvVMEC[phi] (B0[phi] (-2 Sin[2 delta[phi]] Sinh[eta[phi]] (delta'[phi]/sprimeVMEC[phi]-torsVMEC[phi])+(Cos[2 delta[phi]] Cosh[eta[phi]]+2 Sinh[eta[phi]]) eta'[phi]/sprimeVMEC[phi])+(2 Cosh[eta[phi]]+Cos[2 delta[phi]] Sinh[eta[phi]]) (B0'[phi]/sprimeVMEC[phi])),
1/8 pi B0[phi] (2 curvVMEC[phi] (3 torsVMEC[phi] Cosh[eta[phi]]+2 Sinh[eta[phi]] Tanh[eta[phi]] (delta'[phi]/sprimeVMEC[phi]-torsVMEC[phi]))-2 Cos[2 delta[phi]] Sinh[eta[phi]] (24 phi33s[phi]+curvVMEC[phi] (2 torsVMEC[phi]+3 (delta'[phi]/sprimeVMEC[phi]-torsVMEC[phi])))-Sin[2 delta[phi]] (4 Sinh[eta[phi]] (12 phi33c[phi]+curvVMEC'[phi]/sprimeVMEC[phi])+3 Cosh[eta[phi]] curvVMEC[phi] eta'[phi]/sprimeVMEC[phi]))+pi curvVMEC[phi] (2 B0[phi] Cos[2 delta[phi]] Sinh[eta[phi]] (delta'[phi]/sprimeVMEC[phi]-torsVMEC[phi])+Sin[2 delta[phi]] (B0[phi] Cosh[eta[phi]] eta'[phi]/sprimeVMEC[phi]+Sinh[eta[phi]] (B0'[phi]/sprimeVMEC[phi]))),
1/8 pi (8 B0[phi] curvVMEC[phi] (2 Sin[2 delta[phi]] Sinh[eta[phi]] (delta'[phi]/sprimeVMEC[phi]-torsVMEC[phi])-Cos[2 delta[phi]] Cosh[eta[phi]] eta'[phi]/sprimeVMEC[phi])+B0[phi] (-48 Cosh[eta[phi]] phi33c[phi]+2 Cos[2 delta[phi]] Sinh[eta[phi]] curvVMEC'[phi]/sprimeVMEC[phi]+curvVMEC[phi] Sinh[eta[phi]] (2 torsVMEC[phi] Sin[2 delta[phi]]-2  mu[phi] Sin[4 delta[phi]] (delta'[phi]/sprimeVMEC[phi]-torsVMEC[phi])+Cos[4 delta[phi]] eta'[phi]/sprimeVMEC[phi]))-8 Cos[2 delta[phi]] curvVMEC[phi] Sinh[eta[phi]] (B0'[phi]/sprimeVMEC[phi])),
pi B0[phi] curvVMEC[phi] (2 Cos[2 delta[phi]] Sinh[eta[phi]] (delta'[phi]/sprimeVMEC[phi]-torsVMEC[phi])+Cosh[eta[phi]] Sin[2 delta[phi]] eta'[phi]/sprimeVMEC[phi])-1/8 pi B0[phi] (48 Cosh[eta[phi]] phi33s[phi]+2 Sin[2 delta[phi]] Sinh[eta[phi]] curvVMEC'[phi]/sprimeVMEC[phi]+curvVMEC[phi] Sinh[eta[phi]] (-2 torsVMEC[phi] Cos[2 delta[phi]]+2  mu[phi] Cos[4 delta[phi]] (delta'[phi]/sprimeVMEC[phi]-torsVMEC[phi])+Sin[4 delta[phi]] eta'[phi]/sprimeVMEC[phi]))+pi curvVMEC[phi] Sin[2 delta[phi]] Sinh[eta[phi]] (B0'[phi]/sprimeVMEC[phi])
};

optsCT3 = {Table[psi[i] -> psii[[i]], {i, 1, 4}]} // Quiet // Flatten;

fpsi3 = Compile[{{phi, _Real}, {psii, _Real, 1}},  Evaluate[(-sprimeVMEC[phi] (psiMatA[phi].{psi[1], psi[2], psi[3], psi[4]} + psiMatB[phi])) /. optsCT3]];

dfpsi3 = Compile[{{phi, _Real}, {psii, _Real, 1}}, Evaluate[-sprimeVMEC[phi] (psiMatA[phi]) /. optsCT3]];

qsresidualpsi3[Dmat_, psii_, ns_] := 
  Module[{delta = (2 pi)/ns },
   fpsi3Mat  = Table[fpsi3[ delta i, psii[[i]]], {i, 1, ns}];
   dfpsi3Mat = Table[dfpsi3[delta i, psii[[i]]], {i, 1, ns}];
   f1Mat = {
      Dmat.psii[[All, 1]] + fpsi3Mat[[All, 1]],
      Dmat.psii[[All, 2]] + fpsi3Mat[[All, 2]],
      Dmat.psii[[All, 3]] + fpsi3Mat[[All, 3]],
      Dmat.psii[[All, 4]] + fpsi3Mat[[All, 4]]
      } // Flatten;
   zeroarray = DiagonalMatrix[ConstantArray[0, ns]];
   myArrayFlatten = Flatten /@ Flatten[#, {{1, 3}}] &;
   Jacf1Mat = myArrayFlatten@{
      {Dmat + DiagonalMatrix[dfpsi3Mat[[All, 1, 1]]], 
       DiagonalMatrix[dfpsi3Mat[[All, 1, 2]]], 
       DiagonalMatrix[dfpsi3Mat[[All, 1, 3]]], 
       DiagonalMatrix[dfpsi3Mat[[All, 1, 4]]]},
      {DiagonalMatrix[dfpsi3Mat[[All, 2, 1]]], 
       Dmat + DiagonalMatrix[dfpsi3Mat[[All, 2, 2]]], 
       DiagonalMatrix[dfpsi3Mat[[All, 2, 3]]], 
       DiagonalMatrix[dfpsi3Mat[[All, 2, 4]]]},
      {DiagonalMatrix[dfpsi3Mat[[All, 3, 1]]], 
       DiagonalMatrix[dfpsi3Mat[[All, 3, 2]]], 
       Dmat + DiagonalMatrix[dfpsi3Mat[[All, 3, 3]]], 
       DiagonalMatrix[dfpsi3Mat[[All, 3, 4]]]},
      {DiagonalMatrix[dfpsi3Mat[[All, 4, 1]]], 
       DiagonalMatrix[dfpsi3Mat[[All, 4, 2]]], 
       DiagonalMatrix[dfpsi3Mat[[All, 4, 3]]], 
       Dmat + DiagonalMatrix[dfpsi3Mat[[All, 4, 4]]]}
      };
   {f1Mat, Jacf1Mat}
   ];

solvepsi3[Dmat_, nIts_, tol_, nLineS_, psiA_, ns_] := Module[
  {deltaold = 0, delta = (2 pi)/ns },
  psii = psiA;
  For[i = 1, i <= nIts, i++,
   {f1, Jacf1} = qsresidualpsi3[Dmat, psii, ns];
   deltay = -LinearSolve[Jacf1, f1] // Quiet;
   res[i] = Norm[deltay];
   If[(Abs[Norm[deltay] - Norm[deltaold]])/Sqrt[ns] < tol, Break[]];
   For[j = 1, j < nLineS, j++,
    psiAtemp = psii; 
    psiAtemp[[All, 1]] = psiAtemp[[All, 1]] + deltay[[1 ;; ns]]; 
    psiAtemp[[All, 2]] = psiAtemp[[All, 2]] + deltay[[ns + 1 ;; 2 ns]];
    psiAtemp[[All, 3]] = psiAtemp[[All, 3]] + deltay[[2 ns + 1 ;; 3 ns]]; 
    psiAtemp[[All, 4]] = psiAtemp[[All, 4]] + deltay[[3 ns + 1 ;; 4 ns]];
    {f2, Jacf2} = qsresidualpsi3[Dmat, psiAtemp, ns];
    If[Abs[Norm[f2]] < Abs[Norm[f1]], Break[], deltay = deltay/2];
   ];
   psii = psiAtemp;
   deltaold = deltay;
  ];
psii
]

finalPsi3[nIts_, tol_, nLineS_, psiinit_, ns_] :=
  Module[{delta = (2 Pi)/ns},
   psii = Table[psiinit[delta i], {i, 1, ns}];
   (*Derivative Matrix*)
   column1 = Prepend[ Table[ 0.5 (-1)^n Cot[n/2 delta], {n, 1, ns - 1}], 0];
   column2 = RotateRight[Reverse[column1]];
   Dmat = ToeplitzMatrix[column1, column2];
   (*Solve Psi3*)
   psisol2 = solvepsi3[Dmat, nIts, tol, nLineS, psii, ns]
];

psiinit[t_] = {0.01 Cos[vmecNFP t], 0.01 Sin[vmecNFP t], 0.01 Cos[vmecNFP t], 0.01 Sin[vmecNFP t]};
ns = 220; nIts = 18; nLineS = 5; tol = 10^-6;

sol=finalPsi3[nIts,tol,nLineS,psiinit,ns];

psi31 = ListInterpolation[sol[[All, 1]], {{0, 2 pi}}];
psi32 = ListInterpolation[sol[[All, 2]], {{0, 2 pi}}];
psi33 = ListInterpolation[sol[[All, 3]], {{0, 2 pi}}];
psi34 = ListInterpolation[sol[[All, 4]], {{0, 2 pi}}];