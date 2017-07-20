(* ::Package:: *)

(* ::Section:: *)
(*Readme*)


(* ::Text:: *)
(*NLOUnitarityBounds v0.2*)
(*July 2017*)
(**)
(*This package implements the algorithm for finding the (approximate) NLO unitarity bounds for the quartic couplings in a general renormalizable theory. It*)
(*was made using Mathematica x11.0.1.0. Please cite arXiv:1702.08511 if you use this package.*)
(**)
(*Examples and some documentation of the algorithm are shown *)
(*in the supplementary file NLOUnitarityBounds-example.nb.*)
(**)
(*Please send comments *)
(*and suggestions to the email address below.*)
(**)
(*Copyright 2017*)
(*Chris Murphy		cmurphy at quark dot phys dot bnl dot gov*)
(**)


(* ::Section:: *)
(*History*)


(* ::Text:: *)
(*v0.1*)
(*This is the first attempt at making a Mathematica package to implement the algorithm for finding the (approximate) NLO unitarity bounds for the quartic couplings in a general renormalizable theory*)
(**)
(*v0.2*)
(*various minor improvements made*)


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
"NLOUnitarityBounds[ <partialwavematrix> ,<betapartialwavematrix> , <quarticcouplings> <betafunctions> ]:
Returns a list whose entries are {\!\(\*TemplateBox[{\"a\",\"0\",RowBox[{\"(\", \"0\", \")\"}]},\n\"Subsuperscript\"]\), \!\(\*TemplateBox[{\"a\",\"0\",RowBox[{\"(\", \"1\", \")\"}]},\n\"Subsuperscript\"]\)}, the LO and (approximate) NLO contributions to each eigenvalue;\[IndentingNewLine]<partialwavematrix> the tree level partial wave matrix, 1x1 'matrices' should still be entered as {{x}} rather than x;
<betapartialwavematrix> the beta function contribution to the partial wave matrix;\[IndentingNewLine]<quarticcouplings> a list of each quartic coupling appear in <partialwavematrix>, complex conjugates must also be listed (real and imaginary parts would work as well);\[IndentingNewLine]<betafunctions> a list that gives the beta function for each quartic coupling in the previous list. S;
";


Begin["`Private`"]


(* ::Subsection:: *)
(*NLOUnitarityBounds function*)


NLOUnitarityBounds[partialwavematrix_,betapartialwavematrix_,quarticcouplings_,betafunctions_]:=




(* ::Subsection:: *)
(*Module:*)


Module[{evals0,evecs0,evals\[Sigma],evals\[Beta],it},



(* ::Text:: *)
(*LO eigenvectors and eigenvalues:*)


{evals0,evecs0}=Eigensystem[partialwavematrix];


(* ::Text:: *)
(*"\[Sigma]-terms" of NLO eigenvalues:*)


evals\[Sigma]=(I -1/\[Pi])evals0^2;


(* ::Text:: *)
(*"\[Beta]-terms" of NLO eigenvalues:*)


evals\[Beta]=Table[evecs0[[it]].betapartialwavematrix.evecs0[[it]],{it,Length[evecs0]}];


(* ::Text:: *)
(*Output and finish:*)


Transpose[{evals0,evals\[Sigma]+evals\[Beta]}]
];


(* ::Subsection:: *)
(*End*)


End[];
EndPackage[];
