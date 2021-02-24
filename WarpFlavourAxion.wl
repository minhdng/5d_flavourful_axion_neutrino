(* ::Package:: *)

(* :Title: WarpFlavourAxion *)
(* :Context: WarpFlavourAxion` *)
(* :Author: Minh D. Nguyen *)
(* :Summary: Supplementary Mathematica package *)
(* :Copyright: Copyright 2020, Minh D. Nguyen *)
(* :Mathematica Version: 12.0 *)
(* :History:
  V1.00 
  V1.01 Include leptons
  V1.02 Update constraint 3 for PMNS
  V1.03 Include running to Axion scale (10^10 GeV)
  V1.04 Update constraint 1 for PMNS
  V2.00 Finalize for first paper
 *)
(* :Notes:

 *)

(****************************************************************
 * Begin package
 ****************************************************************)


BeginPackage["WarpFlavourAxion`"];


(****************************************************************
 * Messages
 ****************************************************************)


(* Usage messages *)


(* constants *)
$UpMass::usage = "$UpMass return a vector of {up, charm, top} quark masses of in GeV";
$UpMassError::usage = "$UpMassError returns uncertainty of UpMass in GeV";
$DownMass::usage = "$DownMass[1], $DownMass[2], $DownMass[3] return masses of down, strange and bottom quark in GeV";
$DownMassError::usage = "$DownMassError returns uncertainty of DownMass in GeV";
$LeptonMass::usage = "$LeptonMass[1], $LeptonMass[2], $LeptonMass[3] return masses of electron, muon and tau in GeV";
$LeptonMassError::usage = "$LeptonMassError returns uncertainty of LeptonMass in GeV";
$Yu::usage = "$Yu returns a vector of {up, charm, top} Yukawa couplings at 10^10 GeV for tan beta = 3";
$YuError::usage = "$YuError returns a vector of {up, charm, top} Yukawa coupling uncertainty at 10^10 GeV for tan beta = 3";
$Yd::usage = "$Yd returns a vector of {down, strange, bottom} Yukawa couplings at 10^10 GeV for tan beta = 3";
$YdError::usage = "$YdError returns a vector of {down, strange, bottom} Yukawa coupling uncertainty at 10^10 GeV for tan beta = 3";
$Ye::usage = "$Ye returns a vector of {electron, mu, tau} Yukawa couplings at 10^10 GeV for tan beta = 3";
$YeError::usage = "$YeError returns a vector of {electron, mu, tau} Yukawa coupling uncertainty at 10^10 GeV for tan beta = 3";
$CkmTheta12::usage = "0.22848";
$CkmTheta12Sigma::usage = "0.00050"; 
$CkmTheta12Error::usage = "{0.00050, -0.00049}";
$CkmTheta13::usage = "0.00361"; 
$CkmTheta13Sigma::usage = "0.00010"; 
$CkmTheta13Error::usage ="{0.00011, -0.00009}";
$CkmTheta23::usage = "0.04054"; 
$CkmTheta23Sigma::usage = "0.00072"; 
$CkmTheta23Error::usage ="{0.0083; -0.00061}";
$CkmDelta::usage = "1.196"; 
$CkmDeltaSigma::usage = "0.044"; 
$CkmDeltaError::usage = "{0.045, -0.043}";
$CkmLambda::usage = "0.22453"; 
$CkmLambdaError::usage = "0.00044";
$CkmA::usage = "0.836"; 
$CkmAError::usage = "0.015";
$CkmRhoBar::usage = "0.122"; 
$CkmRhoBarError::usage = "0.018";
$CkmEtaBar::usage = "0.355"; 
$CkmEtaBarError::usage = "0.012";
$PmnsTheta12::usage = "33.62 * Pi/180"; 
$PmnsTheta12Error::usage = "{0.78, -0.76} * Pi/180";
$PmnsTheta23::usage = "47.2 * Pi/180"; 
$PmnsTheta23Error::usage = "{1.9, -3.9} * Pi/180";
$PmnsTheta13::usage = "8.52 * Pi/180"; 
$PmnsTheta13Error::usage = "{0.15, -0.15} * Pi/180";
$PmnsDelta::usage = "234 * Pi/180"; 
$PmnsDeltaError::usage = "{43,-31} * Pi/180";

(* Random Functions *)

RandomMatrix::usage = "RandomMatrix[dim_, min_:0.001, max_:3.] generates a random square matrix of dimension dim with complex entries with uniformly distributed absolute value between [min, max]";
RandomUnitaryMatrix::usage = "RandomUnitaryMatrix[dim_, min_:0.001, max_:3.] generates a unitary random square matrix of dimension dim with complex entries with uniformly distributed absolute value between [min, max]";

(* Fermion profiles and overlaps *)

FermionProfile::usage = "FermionProfile[c_?NumericQ, z_, zir_:10^8]";
FermionProfileUVOverlap::usage = "FermionProfileUVOverlap[cL_?NumericQ, cR_?NumericQ, zir_:10^8]";
FermionProfileBulkOverlap::usage = "FermionProfileBulkOverlap[cL_?NumericQ, cR_?NumericQ, zir_:10^8]";
FermionProfileTilde::usage = "FermionProfileTilde[cL_?NumericQ, zir_:10^8, zuv_:1.01]"
TurnRight::usage = "TurnRight[listLR_]";
TurnLeft::usage = "TurnLeft[listPM_]";

(* Constraint 1: CKM and PMNS constraints *)

CkmRhoEtaBarPred::usage = "CkmRhoEtaBarPred[yU_,yD_,yUMinor_,yDMinor_]";
CkmQ::usage = "CkmQ[yU, yD, yUMinor, yDMinor]";

PmnsTanTheta12Pred::usage = "PmnsTanTheta12Pred[y\[Nu]_,ye_,M\[Nu]_,Me_,r13_,r23_]";
PmnsTanTheta13Pred::usage = "PmnsTanTheta13Pred[y\[Nu]_,ye_,M\[Nu]_,Me_,r13_,r23_]";
PmnsTanTheta23Pred::usage = "PmnsTanTheta23Pred[y\[Nu]_,ye_,M\[Nu]_,Me_,r13_,r23_]";
PmnsDeltaPred::usage = "PmnsDeltaPred[y\[Nu]_,ye_,M\[Nu]_,Me_,r13_,r23_]";
PmnsQ::usage = "PmnsQ[yU_,yD_,yUMinor_,yDMinor_]";

(* Constraint 2 *)

QuarkEffYukawa::usage = "QuarkEffYukawa[yU_, yD_, yUMinor_, yDMinor_]";
LeptonEffYukawa::usage = "LeptonEffYukawa[yN_, yE_, yNMinor_, yEMinor_, r13_, r23_]";

QuarkProfileBoundedQ::usage = "QuarkProfileBoundedQ[yU_,yD_,yUMinor_,yDMinor_,v_,\[Beta]_,threshold_]";
LeptonProfileBoundedQ::usage = "LeptonProfileBoundedQ[yU_,yD_,yUMinor_,yDMinor_, r13_, r23_ ,threshold_]";

(* Constraint 3 *)
solveM::usage = "solveM[cM_,quark_,  yU_, yD_, yUMinor_, yDMinor_,seed_:0.1,zir_:10^10,zuv_:1.01]";
solveQuark::usage = "solveQuark[cM_, yU_, yD_, yUMinor_, yDMinor_, seed_:0.1,zir_:10^10,zuv_:1.01]";
QuarkAMatrices::usage = "QuarkAMatrices[yU_, yD_, yUMinor_, yDMinor";

solveMLepton::usage = "solveMLepton[cM_, lepton_, yU_, yD_, yUMinor_, yDMinor_, seed_:0.1, zir_:10^10, zuv_:1.01]";
LeptonAMatrices::usage = "LeptonAMatrices[y\[Nu]_,ye_,M\[Nu]_,Me_,r13_,r23_]:";
UnitaryQ::usage = "UnitaryQ[mat_,thres_]";
QuarkAMatricesUnitaryQ::usage = "QuarkAMatricesUnitaryQ[yU_,yD_,yUMinor_,yDMinor_,thres_, v_, \[Beta]_]";
LeptonAMatricesUnitaryQ::usage = "LeptonAMatricesUnitaryQ[yN_,yE_,yNMinor_,yEMinor_,thres_, v_, \[Beta]_]";
(* QuarkConstraintSummary::usage = "QuarkConstraintSummary[yU_, yD_, yUMinor_, yDMinor_]"; *)

(* Analytical profile and profile overlap functions *)

FermionProfileUVOverlapCM::usage = "FermionProfileUVOverlapCM[cP_?NumericQ, cM_?NumericQ, zir_:10^8]";
FermionProfileBulkOverlapCM::usage = "FermionProfileBulkOverlapCM[cP_?NumericQ, cM_?NumericQ, zir_:10^8]";
FermionAxionOverlap::usage = "FermionAxionOverlap[c_,\[CapitalDelta]_,zir_:10^8]";

AxionProfileUnnormalizedExact::usage = "AxionProfileUnnormalizedExact[z_,\[Sigma]0_,\[CapitalDelta]_,zIR_:10^8,g_:1]";
AxionNormalization::usage = "AxionNormalization[\[Sigma]0_,\[CapitalDelta]_,zIR_:10^8,g_:1]";
AxionProfileExact::usage = "AxionProfileExact[z_,\[Sigma]0_,\[CapitalDelta]_,zIR_:10^8,g_:1]";
FermionAxionOverlap2Unnormalized::usage = "FermionAxionOverlap2Unnormalized[c_,\[CapitalDelta]_,\[Sigma]0_,zIR_:10^8,g_:1]";
FermionAxionOverlap2::usage = "FermionAxionOverlap2[c_,\[CapitalDelta]_,\[Sigma]0_,zIR_:10^8,g_:1]";

(* Plot Utilities *)

FixedRange::usage = "FixedRange[min_, max_, length_]";
ListMirror::usage = "ListMirror[list_]";


(****************************************************************
 * Begin private context
 ****************************************************************)


Begin["`Private`"];


(****************************************************************
 * Constants and experimental values
 ****************************************************************)


(* Quark and lepton masses in GeV *)


$UpMass = {0.0023, 1.275, 173.21};
$UpMassError = {0.0012, 0.025, 1.22};
$DownMass = {0.0048 , 0.095, 4.18}; 
$DownMassError = {0.0008, 0.005, 0.03};
$LeptonMass = {0.00051099895000, 0.1056583755, 1.77686};
$LeptonMassError = {0.00000000000015, 0.0000000023, 0.00012};


(* Quark and lepton Yukawa at 10^10 GeV, tan \[Beta] = 3 
	See running.nb for detail calculation
*)


$Yu = {3.29456*10^-6, 0.00165737, 0.497757};
$YuError = {1.41*10^-9, 7.4*10^-7, 2.*10^-6};
$Yd = {0.0000244146, 0.000486184, 0.0237974};
$YdError = {5.*10^-10, 2.21*10^-7, 4.*10^-7};
$Ye = {2.96535*10^-6, 0.000624451, 0.0106073};
$YeError = {3.*10^-11, 3.*10^-9, 1.*10^-7};


$CkmTheta12 = 0.22848; $CkmTheta12Sigma = 0.00050; $CkmTheta12Error = {0.00050, -0.00049};
$CkmTheta13 = 0.00361; $CkmTheta13Sigma = 0.00010; $CkmTheta13Error ={0.00011, -0.00009};
$CkmTheta23 = 0.04054; $CkmTheta23Sigma = 0.00072; $CkmTheta23Error ={0.0083; -0.00061};
$CkmDelta = 1.196; $CkmDeltaSigma = 0.044; $CkmDeltaError = {0.045, -0.043};


(* Wolfenstefan parameters for CKM *)


$CkmLambda = 0.22453; $CkmLambdaError = 0.00044;
$CkmA = 0.836; $CkmAError = 0.015;
$CkmRhoBar = 0.122; $CkmRhoBarError = 0.018;
$CkmEtaBar = 0.355; $CkmEtaBarError = 0.012;


(* PMNS angle parameters, normal hierarchy *)


$PmnsTheta12 = 33.62 * Pi/180; $PmnsTheta12Error = {0.78, -0.76} * Pi/180;
$PmnsTheta23 = 47.2 * Pi/180; $PmnsTheta23Error = {1.9, -3.9} * Pi/180;
$PmnsTheta13 = 8.52 * Pi/180; $PmnsTheta13Error = {0.15, -0.15} * Pi/180;
$PmnsDelta = 234 * Pi/180; $PmnsDeltaError = {43,-31} * Pi/180;


(****************************************************************
 * Fermion profiles and their overlaps
 ****************************************************************)


FermionProfile[c_?NumericQ, z_, zir_:1. 10^8, zuv_:1.01]:=
	Piecewise[{
		{Sqrt[(1-2c)/(2(zir^(1-2c) - zuv^(1-2c)))] z^(2-c), c != 1/2},
		{1/Sqrt[2(Log[zir]-Log[zuv])], c == 1/2}
	}];


FermionProfileUVOverlap[cL_?NumericQ, cR_?NumericQ, zir_:10^8, zuv_:1.01]:= FermionProfile[cL, zuv, zir] FermionProfile[cR, zuv, zir]


FermionProfileBulkOverlap[cL_?NumericQ, cR_?NumericQ, zir_:10^8, zuv_:1.01]:=
	Piecewise[{{Sqrt[(1-2cL)/(zir^(1-2cL) - zuv^(1-2cL))],cL != 1/2}, {Sqrt[1/Log[zir]], cL == 1/2}}]*
	Piecewise[{{Sqrt[(1-2cR)/(zir^(1-2cR) - zuv^(1-2cL))],cR != 1/2}, {Sqrt[1/Log[zir]], cR == 1/2}}]* 
	Piecewise[{{(1-zir^(-(cL+cR)))/(cL+cR), cL+cR != 0},{Log[zir], cL+cR == 0}}
];


FermionProfileTilde[cL_?NumericQ, zir_:10^8, zuv_:1.01]:= Piecewise[{{Sqrt[(1-2cL)/(zir^(1-2cL) - zuv^(1-2cL))], cL != 1/2}, {Sqrt[1/Log[zir]], cL == 1/2}}];


TurnRight[listLR_]:=Module[{cL=Transpose[listLR][[1]], cR=Transpose[listLR][[2]]}, 
	Transpose[{cL-cR, cL+cR}]
];


TurnLeft[listPM_]:=Module[{cM=Transpose[listPM][[1]], cP=Transpose[listPM][[2]]},
	Transpose[{cP+cM, cP-cM}]/2
];


(****************************************************************
 * Constraint functions
 ****************************************************************)


RandomMatrix[dim_, min_:0.001, max_:3.] := RandomReal[{min, max},{dim, dim}] Exp[I RandomReal[{0., 2\[Pi]},{dim, dim}]]


RandomUnitaryMatrix[dim_, min_:0.001, max_:3.] := QRDecomposition[RandomMatrix[dim, min, max]][[1]];


(* Quark sector, Constraints 1. *)


CkmRhoEtaBarPred[yU_,yD_,yUMinor_,yDMinor_]:= 
Through[
	{Re,Im}
	[ 
		(yD[[3,3]] yUMinor[[3,1]]-yD[[2,3]] yUMinor[[2,1]]+yD[[1,3]] yUMinor[[1,1]])/
		(
			yD[[3,3]] yUMinor[[1,1]] 
			(yD[[2,3]]/yD[[3,3]]-yU[[2,3]]/yU[[3,3]]) 
			(yDMinor[[2,1]]/yDMinor[[1,1]]-yUMinor[[2,1]]/yUMinor[[1,1]])
		) 
	]
];


CkmQ[yU_,yD_,yUMinor_,yDMinor_]:=(
(Norm[ CkmRhoEtaBarPred[yU, yD, yUMinor, yDMinor][[1]] - $CkmRhoBar] < $CkmRhoBarError) 
&& 
(Norm[ CkmRhoEtaBarPred[yU, yD, yUMinor, yDMinor][[2]] + $CkmEtaBar] < $CkmEtaBarError)
);


(* Lepton Sector, Constraint 1 *)


PmnsTanTheta12Pred[ye_, Me_, r13_, r23_]:=With[
{
NeL1 = Sqrt[r23^2 Abs[Me[[1,1]]]^2+r13^2 Abs[Me[[2,1]]]^2+r13^2 r23^2 Abs[Me[[3,1]]]^2],
NeL3 = Sqrt[r13^2 Abs[ye[[1,3]]]^2+r23^2 Abs[ye[[2,3]]]^2+Abs[ye[[3,3]]]^2]
},
Abs[Me[[2, 1]] / Me[[1, 1]]] (r13/r23)
];


PmnsTanTheta23Pred[ye_, Me_, r13_, r23_]:=With[
{
NeL1 = Sqrt[r23^2 Abs[Me[[1,1]]]^2+r13^2 Abs[Me[[2,1]]]^2+r13^2 r23^2 Abs[Me[[3,1]]]^2],
NeL3 = Sqrt[r13^2 Abs[ye[[1,3]]]^2+r23^2 Abs[ye[[2,3]]]^2+Abs[ye[[3,3]]]^2]
},
Abs[ Conjugate[Me[[1, 1]]] ye[[2, 3]] r23^2 + Conjugate[Me[[2, 1]]] ye[[1, 3]] r13^2 ] / (NeL1 Abs[ye[[3, 3]]])
];


PmnsTanTheta13Pred[ye_, Me_, r13_, r23_]:=With[
{
NeL1=Sqrt[r23^2 Abs[Me[[1,1]]]^2+r13^2 Abs[Me[[2,1]]]^2+r13^2 r23^2 Abs[Me[[3,1]]]^2],
NeL3=Sqrt[r13^2 Abs[ye[[1,3]]]^2+r23^2 Abs[ye[[2,3]]]^2+Abs[ye[[3,3]]]^2]
},
Abs[Me[[3, 1]] r13 r23 / NeL1]
];


PmnsDeltaPred[ye_, Me_, r13_, r23_]:=With[
{
NeL1=Sqrt[r23^2 Abs[Me[[1,1]]]^2+r13^2 Abs[Me[[2,1]]]^2+r13^2 r23^2 Abs[Me[[3,1]]]^2],
NeL3=Sqrt[r13^2 Abs[ye[[1,3]]]^2+r23^2 Abs[ye[[2,3]]]^2+Abs[ye[[3,3]]]^2]
},
Mod[2Pi - Arg[Me[[3, 1]] r13 r23 / NeL1], 2Pi]//N
];


(* Different from the CkmQ function, the charged-lepton sector require r13 and r23 for the left-handed profiles *)


PmnsQ[yE_, yEMinor_,r13_,r23_,rho_:1]:=(
(Tan[$PmnsTheta12 +rho $PmnsTheta12Error[[2]]] < PmnsTanTheta12Pred[yE, yEMinor, r13, r23] < Tan[$PmnsTheta12 +rho $PmnsTheta12Error[[1]]])
&&
(Tan[$PmnsTheta13 +rho $PmnsTheta13Error[[2]]] < PmnsTanTheta13Pred[yE, yEMinor, r13, r23] < Tan[$PmnsTheta13 +rho $PmnsTheta13Error[[1]]])
&&
(Tan[$PmnsTheta23 +rho $PmnsTheta23Error[[2]]] < PmnsTanTheta23Pred[yE, yEMinor, r13, r23] < Tan[$PmnsTheta23 +rho $PmnsTheta23Error[[1]]])
&&
($PmnsDelta +rho $PmnsDeltaError[[2]] < PmnsDeltaPred[yE, yEMinor, r13, r23] < $PmnsDelta +rho $PmnsDeltaError[[1]])
);


(* Constraint 2.  *)


QuarkEffYukawa[yU_, yD_, yUMinor_, yDMinor_]:= <|
"u"-> $Yu[[1]] Norm[yUMinor[[1,1]]]/Norm[Det[yU]],
"c"-> $Yu[[2]] Norm[yU[[3,3]]]/Norm[yUMinor[[1,1]]],
"t"-> $Yu[[3]]/Norm[yU[[3,3]]],
"d"-> $Yd[[1]] Norm[yDMinor[[1,1]]]/Norm[Det[yD]],
"s"-> $Yd[[2]] Norm[yD[[3,3]]]/Norm[yDMinor[[1,1]]],
"b"-> $Yd[[3]]/Norm[yD[[3,3]]]
|>;


LeptonEffYukawa[ye_, Me_, r13_, r23_]:= Module[{NeL1, NeL3},
NeL1=Sqrt[r23^2 Abs[Me[[1,1]]]^2 + r13^2 Abs[Me[[2,1]]]^2 + r13^2 r23^2 Abs[Me[[3,1]]]^2];
NeL3=Sqrt[r13^2 Abs[ye[[1,3]]]^2 + r23^2 Abs[ye[[2,3]]]^2 + Abs[ye[[3,3]]]^2];
 <|
"e"-> $Ye[[1]] NeL1 / (Norm[Det[ye]] r23),
"mu"-> $Ye[[2]] r23 NeL3 / NeL1,
"tau"-> $Ye[[3]] / NeL3
|>];


QuarkProfileBoundedQ[yU_, yD_, yUMinor_, yDMinor_, threshold_:0.95]:= 
	AllTrue[ Values[QuarkEffYukawa[yU, yD, yUMinor, yDMinor]], (#<threshold)& ];
LeptonProfileBoundedQ[yE_, yEMinor_,  r13_, r23_,threshold_:0.95]:= 
	AllTrue[ Values[LeptonEffYukawa[yE, yEMinor, r13, r23]], (#<threshold)& ];


(* Constraint 3 *)


solveM[cM_,quark_,  yU_, yD_, yUMinor_, yDMinor_,seed_:0.1,zir_:10^10,zuv_:1.01]:=Module[{cR,cL1, cR1},
If[cM>=0,
cR1=cR/.FindRoot[Log[FermionProfileTilde[cM + cR,zir,zuv]FermionProfileTilde[cR,zir,zuv]]==Log[QuarkEffYukawa[yU, yD, yUMinor, yDMinor][quark]],{cR,seed}];
cL1 = cM + cR1;
{cL1, cR1},
{cL1, cR1} = solveM[-cM, quark, yU, yD, yUMinor, yDMinor, seed];
{cR1, cL1}
]
]


solveQuark[cM_, yU_, yD_, yUMinor_, yDMinor_, seed_:0.1,zir_:10^10,zuv_:1.01]:=Module[{cL,cR,cL1, cL2, cL3, cu1, cu2, cu3, cd1, cd2, cd3, fQ1, fQ2, fQ3, fu1, fu2, fu3, fd1, fd2,fd3},
(* First generation *)
{cL1, cu1} = solveM[cM, "u",yU, yD, yUMinor, yDMinor, seed];
cd1 =cR/.FindRoot[Log[FermionProfileTilde[cR,zir,zuv]]==Log[QuarkEffYukawa[yU, yD, yUMinor, yDMinor]["d"]/FermionProfileTilde[cL1,zir,zuv]],{cR, seed}];
fQ1=FermionProfileTilde[cL1];

(* Second generation *)
fQ2=(fQ1 / $CkmLambda) Norm[yDMinor[[2, 1]]/yDMinor[[1, 1]]-yUMinor[[2, 1]]/yUMinor[[1, 1]]];
cL2 =cL/.FindRoot[ Log[FermionProfileTilde[cL,zir,zuv]] == Log[fQ2], {cL,seed}];
fu2 =QuarkEffYukawa[yU, yD, yUMinor, yDMinor]["c"]/fQ2;
cu2=cR/.FindRoot[ Log[FermionProfileTilde[cR,zir,zuv]] == Log[fu2], {cR,seed}];
fd2 =QuarkEffYukawa[yU, yD, yUMinor, yDMinor]["s"]/fQ2;
cd2=cR/.FindRoot[ Log[FermionProfileTilde[cR,zir,zuv]] == Log[fd2], {cR,seed}];

(* Third generation *)
fQ3=(fQ2 / ($CkmLambda^2 $CkmA)) Norm[yD[[2,3]]/yD[[3,3]]-yU[[2,3]]/yU[[3,3]]];
cL3 =cL/.FindRoot[ Log[FermionProfileTilde[cL,zir,zuv]] == Log[fQ3], {cL,seed}];
fu3 =QuarkEffYukawa[yU, yD, yUMinor, yDMinor]["t"]/fQ3;
cu3=cR/.FindRoot[ Log[FermionProfileTilde[cR,zir,zuv]] == Log[fu3], {cR,seed}];
fd3 =QuarkEffYukawa[yU, yD, yUMinor, yDMinor]["b"]/fQ3;
cd3=cR/.FindRoot[ Log[FermionProfileTilde[cR,zir,zuv]] == Log[fd3], {cR,seed}];

(* Output *)
{cL1,cL2,cL3, cu1,cu2,cu3, cd1,cd2,cd3}
];


QuarkAMatrices[yU_, yD_, yUMinor_, yDMinor_]:= Module[{fQ, fu, fd, AuL, AuR, AdL, AdR, phaseU, phaseD},
fQ = {$CkmLambda/Norm[yDMinor[[2,1]]/yDMinor[[1,1]]-yUMinor[[2,1]]/yUMinor[[1,1]]], 1, Norm[yD[[2,3]]/yD[[3,3]]-yU[[2,3]]/yU[[3,3]]]/($CkmA $CkmLambda^2)};
fu = (QuarkEffYukawa[yU, yD, yUMinor, yDMinor][#]&/@{"u","c","t"})/fQ;
fd = (QuarkEffYukawa[yU, yD, yUMinor, yDMinor][#]&/@{"d","s","b"})/fQ;
phaseU = ({
 {E^(-I(Arg[Det[yU]]-Arg[yUMinor[[1,1]]])), 0, 0},
 {0, E^(-I(Arg[yUMinor[[1,1]]]- Arg[yU[[3,3]]])), 0},
 {0, 0, E^(-I Arg[yU[[3,3]]])}
});
phaseD = ({
 {E^(-I(Arg[Det[yD]]-Arg[yDMinor[[1,1]]])), 0, 0},
 {0, E^(-I(Arg[yDMinor[[1,1]]]- Arg[yD[[3,3]]])), 0},
 {0, 0, E^(-I Arg[yD[[3,3]]])}
});
AuL = ({
 {1, yUMinor[[2,1]]/yUMinor[[1,1]] fQ[[1]]/fQ[[2]], yU[[1,3]]/yU[[3,3]] fQ[[1]]/fQ[[3]]},
 {-((yUMinor[[2,1]]//Conjugate)/(yUMinor[[1,1]]//Conjugate)) fQ[[1]]/fQ[[2]], 1, yU[[2,3]]/yU[[3,3]] fQ[[2]]/fQ[[3]]},
 {(yUMinor[[3,1]]//Conjugate)/(yUMinor[[1,1]]//Conjugate) fQ[[1]]/fQ[[3]], -((yU[[2,3]]//Conjugate)/(yU[[3,3]]//Conjugate)) fQ[[2]]/fQ[[3]], 1}
});
AuR = phaseU . ({
 {1, (yUMinor[[1,2]]//Conjugate)/(yUMinor[[1,1]]//Conjugate) fu[[1]]/fu[[2]], (yU[[3,1]]//Conjugate)/(yU[[3,3]]//Conjugate) fu[[1]]/fu[[3]]},
 {-(yUMinor[[1,2]]/yUMinor[[1,1]]) fu[[1]]/fu[[2]], 1, (yU[[3,2]]//Conjugate)/(yU[[3,3]]//Conjugate) fu[[2]]/fu[[3]]},
 {yUMinor[[1,3]]/yUMinor[[1,1]] fu[[1]]/fu[[3]], -(yU[[3,2]]/yU[[3,3]]) fu[[2]]/fu[[3]], 1}
});
AdL = ({
 {1, yDMinor[[2,1]]/yDMinor[[1,1]] fQ[[1]]/fQ[[2]], yD[[1,3]]/yD[[3,3]] fQ[[1]]/fQ[[3]]},
 {-((yDMinor[[2,1]]//Conjugate)/(yDMinor[[1,1]]//Conjugate)) fQ[[1]]/fQ[[2]], 1, yD[[2,3]]/yD[[3,3]] fQ[[2]]/fQ[[3]]},
 {(yDMinor[[3,1]]//Conjugate)/(yDMinor[[1,1]]//Conjugate) fQ[[1]]/fQ[[3]], -((yD[[2,3]]//Conjugate)/(yD[[3,3]]//Conjugate)) fQ[[2]]/fQ[[3]], 1}
});
AdR = phaseD . ({
 {1, (yDMinor[[1,2]]//Conjugate)/(yDMinor[[1,1]]//Conjugate) fd[[1]]/fd[[2]], (yD[[3,1]]//Conjugate)/(yD[[3,3]]//Conjugate) fd[[1]]/fd[[3]]},
 {-(yDMinor[[1,2]]/yDMinor[[1,1]]) fd[[1]]/fd[[2]], 1, (yD[[3,2]]//Conjugate)/(yD[[3,3]]//Conjugate) fd[[2]]/fd[[3]]},
 {yDMinor[[1,3]]/yDMinor[[1,1]] fd[[1]]/fd[[3]], -(yD[[3,2]]/yD[[3,3]]) fd[[2]]/fd[[3]], 1}
});
{AuL, AuR, AdL, AdR}
];


solveMLepton[cM_, lepton_, yE_, yEMinor_, r13_, r23_, seed_:0.1, zir_:10^10, zuv_:1.01] := Module[{cR, cL1, cR1},
If[cM >= 0,
cR1 = cR/.FindRoot[Log[FermionProfileTilde[cM + cR, zir, zuv]FermionProfileTilde[cR, zir, zuv]]==Log[LeptonEffYukawa[yE, yEMinor, r13, r23][lepton]],{cR, seed}];
cL1 = cM + cR1;
{cL1, cR1},
{cL1, cR1} = solveMLepton[-cM, lepton, yE, yEMinor, r13, r23, seed, zir, zuv];
{cR1, cL1}
]
]


LeptonAMatrices[ye_, Me_, r13_, r23_]:=Module[{fL, fe, N\[Nu]L1, N\[Nu]L3, NeL1, NeL3, phaseE, AeL, AeR},
fL = {r13, r23, 1};
fe = (LeptonEffYukawa[ye, Me, r13, r23][#]&/@{"e","mu","tau"}) / fL;
NeL1=Sqrt[r23^2 Abs[Me[[1,1]]]^2 + r13^2 Abs[Me[[2,1]]]^2 + r13^2 r23^2 Abs[Me[[3,1]]]^2];
NeL3=Sqrt[r13^2 Abs[ye[[1,3]]]^2 + r23^2 Abs[ye[[2,3]]]^2 + Abs[ye[[3,3]]]^2];
phaseE = DiagonalMatrix[{Exp[- I Arg[Det[ye]]], 1, 1}];
AeL={
{(r23 Conjugate[Me[[1,1]]])/NeL1,(r13 (Me[[3,1]] r23^2 Conjugate[ye[[2,3]]]+Me[[2,1]] Conjugate[ye[[3,3]]]))/(NeL1 NeL3),(r13 ye[[1,3]])/NeL3},
{-((r13 Conjugate[Me[[2,1]]])/NeL1),r23 (-Me[[3,1]] r13^2 Conjugate[ye[[1,3]]]+Me[[1,1]] Conjugate[ye[[3,3]]])/(NeL1 NeL3),(r23 ye[[2,3]])/NeL3},
{(r13 r23 Conjugate[Me[[3,1]]])/NeL1,-((Me[[2,1]] r13^2 Conjugate[ye[[1,3]]]+Me[[1,1]] r23^2 Conjugate[ye[[2,3]]])/(NeL1 NeL3)),ye[[3,3]]/NeL3}
};
AeR={
{1, (Me[[1,2]]//Conjugate)/(Me[[1,1]]//Conjugate) fe[[1]]/fe[[2]], (ye[[3,1]]//Conjugate)/(ye[[3,3]]//Conjugate) fe[[1]]/fe[[3]]},
{-(Me[[1,2]]/Me[[1,1]]) fe[[1]]/fe[[2]], 1, (ye[[3,2]]//Conjugate)/(ye[[3,3]]//Conjugate) fe[[2]]/fe[[3]]},
{Me[[1,3]]/Me[[1,1]] fe[[1]]/fe[[3]], -(ye[[3,2]]/ye[[3,3]]) fe[[2]]/fe[[3]], 1}
} . phaseE;
{AeL,AeR}
]


(* Different from UnitaryMatrixQ of Mathematica *)
UnitaryQ[mat_,thres_]:= AllTrue[ Abs[Flatten[mat . ConjugateTranspose[mat]-IdentityMatrix[3]]], (#<thres)& ];


QuarkAMatricesUnitaryQ[yU_, yD_, yUMinor_, yDMinor_, thres_:0.2]:= 
	AllTrue[ {1,2,3,4}, UnitaryQ[QuarkAMatrices[yU,yD,yUMinor,yDMinor][[#]], thres]& ];


LeptonAMatricesUnitaryQ[yE_, yEMinor_, r13_, r23_, thres_:0.2]:= 
	AllTrue[ {1,2}, UnitaryQ[LeptonAMatrices[yE, yEMinor, r13, r23][[#]], thres]& ];


(* Script running through all three constraint *)


(* QuarkConstraintSummary[yU_, yD_, yUMinor_, yDMinor_]:= Module[{AdR},
	Print["Constraint 1 result for quark sector"];
	Print["\!\(\*SubscriptBox[\(Y\), \(u\)]\) = ",MatrixForm[yU]];
	Print["\!\(\*SubscriptBox[\(Y\), \(d\)]\) = ",MatrixForm[yD]];
	Print["The resulted \[Rho]bar = ",RhoEtaBar[yU,yD,yUMinor,yDMinor][[1]], ", \[Eta]bar = ", - RhoEtaBar[yU,yD,yUMinor,yDMinor][[2]]];
	Print["Experimental \[Rho]bar = ",$CkmRhoBar, " \[PlusMinus] ", $CkmRhoBarError,", \[Eta]bar = ",$CkmEtaBar, " \[PlusMinus] ", $CkmEtaBarError];
	Print["CkmQ[] returns ",CkmQ[yU,yD,yUMinor,yDMinor]];
	Print["Constraint 2 result for quark sector"];
	Print["\!\(\*SubscriptBox[\(Y\), \(u\)]\) = ",MatrixForm[yU]];
	Print["\!\(\*SubscriptBox[\(Y\), \(d\)]\) = ",MatrixForm[yD]];
	Print["The resulted effective quark masses ", QuarkEffYukawa[yU,yD,yUMinor,yDMinor]];
	Print["QuarkProfileBoundedQ[] returns ",QuarkProfileBoundedQ[yU,yD,yUMinor,yDMinor]];
	Print["Constraint 3 result for quark sector"];
	Print["\!\(\*SubscriptBox[\(Y\), \(u\)]\) = ",MatrixForm[yU]];
	Print["\!\(\*SubscriptBox[\(Y\), \(d\)]\) = ",MatrixForm[yD]];
	AdR = QuarkAMatrices[yU, yD, yUMinor, yDMinor][[4]];
	Print["The resulted \!\(\*SubscriptBox[SuperscriptBox[\(A\), \(d\)], \(R\)]\) = ", MatrixForm[AdR]];
	Print[ " with norm ", MatrixForm[AdR.ConjugateTranspose[AdR]]];
	Print["QuarkAMatricesUnitaryQ[] returns ", QuarkAMatricesUnitaryQ[yU, yD, yUMinor, yDMinor]];
]
*)


(****************************************************************
 * Axion exact and approximate profiles (Massless case)
 ****************************************************************)


AxionProfileUnnormalizedExact[z_,\[Sigma]0_,\[CapitalDelta]_,zIR_:10^8,g_:1]:=(z/zIR) ((g \[Sigma]0)/(2\[CapitalDelta]))^(1/\[CapitalDelta])(-(BesselI[(\[CapitalDelta]-1)/\[CapitalDelta],(g \[Sigma]0)/\[CapitalDelta]]/BesselI[-((\[CapitalDelta]-1)/\[CapitalDelta]),(g \[Sigma]0)/\[CapitalDelta]]) BesselI[1/\[CapitalDelta],(g \[Sigma]0)/\[CapitalDelta] (z/zIR)^\[CapitalDelta]]+BesselI[-(1/\[CapitalDelta]),(g \[Sigma]0)/\[CapitalDelta] (z/zIR)^\[CapitalDelta]] );


AxionNormalization[\[Sigma]0_,\[CapitalDelta]_,zIR_:10^8,g_:1]:=(zIR/\[Sigma]0)/Sqrt[
\[CapitalDelta]^2/((-1+\[CapitalDelta]) Gamma[-(1/\[CapitalDelta])]^2) HypergeometricPFQ[{1/2-1/\[CapitalDelta]},{1-2/\[CapitalDelta],2-1/\[CapitalDelta]},((g \[Sigma]0)/\[CapitalDelta])^2]+
\[CapitalDelta]^4/( (-1+2 \[CapitalDelta])(\[CapitalDelta]-1)^2 Gamma[-(1/\[CapitalDelta])]^2) ((g \[Sigma]0)/(2\[CapitalDelta]))^2 HypergeometricPFQ[{3/2-1/\[CapitalDelta]},{3-2/\[CapitalDelta],3-1/\[CapitalDelta]},((g \[Sigma]0)/\[CapitalDelta])^2]+
(\[CapitalDelta]-1)/(\[CapitalDelta] Gamma[1/\[CapitalDelta]] Gamma[-(1/\[CapitalDelta])]) ((g \[Sigma]0)/(2\[CapitalDelta]))^(-2+2/\[CapitalDelta]) BesselI[(\[CapitalDelta]-1)/\[CapitalDelta],(g \[Sigma]0)/\[CapitalDelta]]/BesselI[-((\[CapitalDelta]-1)/\[CapitalDelta]),(g \[Sigma]0)/\[CapitalDelta]]  (HypergeometricPFQ[{-(1/2)},{1-1/\[CapitalDelta],-1+1/\[CapitalDelta]},((g \[Sigma]0)/\[CapitalDelta])^2]-1)-
Sin[\[Pi]/\[CapitalDelta]]/(\[CapitalDelta]^2 \[Pi]) ((g \[Sigma]0)/(2\[CapitalDelta]))^(-2+2/\[CapitalDelta]) BesselI[(\[CapitalDelta]-1)/\[CapitalDelta],(g \[Sigma]0)/\[CapitalDelta]]/BesselI[-((\[CapitalDelta]-1)/\[CapitalDelta]),(g \[Sigma]0)/\[CapitalDelta]]  (HypergeometricPFQ[{-(1/2)},{-(1/\[CapitalDelta]),1/\[CapitalDelta]},((g \[Sigma]0)/\[CapitalDelta])^2]-1) +
1/(4(1+\[CapitalDelta])  Gamma[1/\[CapitalDelta]]^2) ((g \[Sigma]0)/(2\[CapitalDelta]))^(-2+4/\[CapitalDelta]) (BesselI[(\[CapitalDelta]-1)/\[CapitalDelta],(g \[Sigma]0)/\[CapitalDelta]]/BesselI[-((\[CapitalDelta]-1)/\[CapitalDelta]),(g \[Sigma]0)/\[CapitalDelta]] )^2 (4  (1+\[CapitalDelta]) HypergeometricPFQ[{-(1/2)+1/\[CapitalDelta]},{1+1/\[CapitalDelta],-1+2/\[CapitalDelta]},((g \[Sigma]0)/\[CapitalDelta])^2]+g^2 \[Sigma]0^2  HypergeometricPFQ[{1/2+1/\[CapitalDelta]},{2+1/\[CapitalDelta],1+2/\[CapitalDelta]},((g \[Sigma]0)/\[CapitalDelta])^2])
];


AxionProfileExact[z_,\[Sigma]0_,\[CapitalDelta]_,zIR_:10^8,g_:1]:=AxionNormalization[\[Sigma]0, \[CapitalDelta], zIR, g]AxionProfileUnnormalizedExact[z, \[Sigma]0, \[CapitalDelta], zIR, g]


(* Axion-fermion-fermion overlap integral (non-flat part) 
	c either cL or -cR
*)


FermionAxionOverlap[c_,\[CapitalDelta]_,zir_:10^8]:=(1-2c)/(zir^(1-2c) - 1) 1/(4\[CapitalDelta](\[CapitalDelta]-1)) ((zir^(1-2 c)-zir^(-2 \[CapitalDelta]))/(1-2 c+2 \[CapitalDelta])+\[CapitalDelta] (1/zir^2-zir^(1-2 c))/(3-2 c));


FermionAxionOverlap2Unnormalized[c_,\[CapitalDelta]_,\[Sigma]0_,zIR_:10^8,g_:1]:=( ((g \[Sigma]0)/(2\[CapitalDelta]))^(2/\[CapitalDelta]) BesselI[(-1+\[CapitalDelta])/\[CapitalDelta],(g \[Sigma]0)/\[CapitalDelta]])/(2 zIR^2 BesselI[-1+1/\[CapitalDelta],(g \[Sigma]0)/\[CapitalDelta]] Gamma[1+1/\[CapitalDelta]]) ((1-2c)/(3-2c) (zIR^(3-2c)-1.)/(zIR^(1-2c) - 1.)-1)


FermionAxionOverlap2[c_,\[CapitalDelta]_,\[Sigma]0_,zIR_:10^8,g_:1]:= FermionAxionOverlap2Unnormalized[c, \[CapitalDelta], \[Sigma]0, zIR, g] AxionNormalization[\[Sigma]0, \[CapitalDelta], zIR, g]


(* fermion overlap functions in Subscript[c, P], Subscript[c, M] *)


FermionProfileUVOverlapCM[cP_?NumericQ, cM_?NumericQ, zir_:10^8]:= Module[{zuv=1.1}, 
	FermionProfile[(cP-cM)/2, zuv, zir] FermionProfile[(cP+cM)/2, zuv, zir]
];


FermionProfileBulkOverlapCM[cP_?NumericQ, cM_?NumericQ, zir_:10^8]:=
	Piecewise[{{Sqrt[(1-(cP+cM))/(zir^(1-(cP+cM)) - 1)],cP+cM!=1},{Sqrt[1/Log[zir]],cP+cM==1}}]*
	Piecewise[{{Sqrt[(1-(cP-cM))/(zir^(1-(cP-cM)) - 1)],cP-cM!=1},{Sqrt[1/Log[zir]],cP-cM==1}}]* 
	Piecewise[{{(1-zir^-cP)/cP,cP!=0},{Log[zir],cP==0}}
];


FixedRange[min_, max_, length_]:=Module[{step=(max-min)/(length-1)},Range[min, max+1/2 step, step]];


ListMirror[list_]:= Reverse[Reverse/@list]~Join~list;


(****************************************************************
 * Graphing functions
 ****************************************************************)


(****************************************************************
 * End private context
 ****************************************************************)


End[];


(****************************************************************
 * End package
 ****************************************************************)


EndPackage[];
