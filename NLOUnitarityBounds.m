(* ::Package:: *)

(* ::Section:: *)
(*Readme*)


(* ::Text:: *)
(*NLOUnitarityBounds v0.1*)
(*Feb 2017*)
(**)
(*This package implements the algorithm for finding the (approximate) NLO unitarity bounds for the quartic couplings in a general renormalizable theory. It*)
(*was made using Mathematica x11.0.1.0.*)
(**)
(*Examples and some documentation of the algorithm are shown *)
(*in the supplementary file NLOUnitarityBounds-example.nb.*)
(**)
(*This is the 1st attempt to implement this algorithm in a Mathematica package; no doubt it could be done better. Please send comments *)
(*and suggestions to the email address below.*)
(**)
(*Copyright 2017*)
(*Chris Murphy		cmurphy at quark dot phys dot bnl dot gov*)
(**)


(* ::Section:: *)
(*History*)


(* ::Text:: *)
(*This is the first attempt at making a Mathematica package to implement the algorithm for finding the (approximate) NLO unitarity bounds for the quartic couplings in a general renormalizable theory*)


(* ::Section:: *)
(*License*)


(* ::Text:: *)
(*This package consists of the files NLOUnitarityBounds.m and NLOUnitarityBounds-example.nb. *)
(*It may be freely distributed and modified under the terms & conditions of the MIT License.*)


(* ::Section:: *)
(*Package*)


(* ::Subsection:: *)
(*Intro*)


BeginPackage["NLOUnitarityBounds`"];


NLOUnitarityBounds::usage =
"NLOUnitarityBounds[ <partialwavematrix> , <quarticcouplings>, <betafunctions> ]:
Returns a list whose entries are {\!\(\*TemplateBox[{\"a\",\"0\",RowBox[{\"(\", \"0\", \")\"}]},\n\"Subsuperscript\"]\), \!\(\*TemplateBox[{\"a\",\"0\",RowBox[{\"(\", \"1\", \")\"}]},\n\"Subsuperscript\"]\)}, the LO and (approximate) NLO contributions to each eigenvalue;\[IndentingNewLine]<partialwavematrix> the tree level partial wave matrix, 1x1 'matrices' should still be entered as {{x}} rather than x;\[IndentingNewLine]<quarticcouplings> a list of each quartic coupling appear in <partialwavematrix>, complex conjugates must also be listed (real and imaginary parts would work as well);\[IndentingNewLine]<betafunctions> a list that gives the beta function for each quartic coupling in the previous list;
";


Begin["`Private`"]


(* ::Subsection:: *)
(*NLOUnitarityBounds function*)


NLOUnitarityBounds[partialwavematrix_,quarticcouplings_,betafunctions_]:=








(* ::Subsection:: *)
(*Module:*)


Module[{evals0,evecL,evecs0,betareps,i1,evals\[Sigma],i2,pwmbeta,evals\[Beta],i3,output,i4},








(* ::Text:: *)
(*LO eigenvectors and eigenvalues:*)


evals0=Eigenvalues[partialwavematrix];
evecL=Length[evals0];
evecs0=Eigenvectors[partialwavematrix];


(* ::Text:: *)
(*"\[Sigma]-terms" of NLO eigenvalues:*)


evals\[Sigma]=Table[(I -1/\[Pi])evals0[[i1]]^2,{i1,evecL}];


(* ::Text:: *)
(*"\[Beta]-terms" of NLO eigenvalues:*)


betareps=Table[quarticcouplings[[i2]]->betafunctions[[i2]],{i2,Length[quarticcouplings]}];
pwmbeta=-(3/2)partialwavematrix/.betareps;
evals\[Beta]=Table[Transpose[evecs0][[i3]].pwmbeta.evecs0[[i3]],{i3,evecL}];


(* ::Text:: *)
(*Output and finish:*)


output=Simplify[Table[{evals0[[i4]],evals\[Sigma][[i4]]+evals\[Beta][[i4]]},{i4,evecL}]]
];


(* ::Subsection:: *)
(*End*)


End[];
EndPackage[];
