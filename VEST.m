(* ::Package:: *)

BeginPackage["VEST`"]


ClearAll[IObjectQ,IndexQ,DIndexQ,UnitVecQ,ConstQ,nFIndexQ,delta,leviC,x,v,PD,Xd,Vd]

IObjectQ::usage = "IObjectQ[obj] 
Returns True if a obj is an Einstein notated object";
IndexQ::usage = "IndexQ[ind]
Returns True if ind is an index";
DIndexQ::usage = "DIndexQ[ind]
Returns true if ind is a dummy index in the current expression";
UnitVecQ::usage = "UnitVecQ[obj]
Returns True is obj should be identified as a unit vector";
ConstQ::usage = "ConstQ[obj]
Returns True if object is a constant";
nFIndexQ::usage = "nIndexQ[expr,n]
Checks if expr has n free indices";

delta::usage = "delta[i,j] is the identity operator";
leviC::usage = "leviC[i,j,k] is the totally antisymmetric Levi-Civita tensor";
v::usage = "v[j] is the co-ordinate v";
x::usage = "x[j] is the co-ordinate x";

PD::usage = "PD[obj, Xd[i,j,..], Vd[k,l,..]]
Partial derivative of obj with respect to i,j,.. in x and k,l,.. in v";
Xd::usage = "Xd[inds] is a wrapper for inds inside PD[]";
Vd::usage = "Vd[inds] is a wrapper for inds inside PD[]";


ClearAll[DefObject,unDefObject,DefIndices,DefConst,setUnitVec];

DefObject::usage = "DefObject[obj, n, dep, printas]
Define new Einstein object, obj, with n indices, and dependencies 
dep = {True/False, True/False} specifying {x,v} dependence of obj.
printas (string) defines how the object will be printed in output";

unDefObject::usage = "unDefObject[obj]
Deletes obj and all associated definitions";

DefIndices::usage = "DefIndices[ind]
Sets symbols in the list ind to be used as indices.
Options:
useForPrinting->True  - If False, defined indices will not be used for printing (sometimes
			necessary to avoid dummy index name conflicts)";
useForPrinting::usage = "Option for DefIndices, default is True.";
	
DefConst::usage = "DefConst[obj]
Defines obj as a constant to be ignored by PD";

setUnitVec::usage = "setUnitVec[obj]
Sets obj to have the properties of a unit vector";


ClearAll[FindDummies];

FindDummies::usage = "FindDummies[expr]
Checks that expr is a valid expression and finds dummy indice, replacing them with
unique variables. Equivalent terms will not be summed together after application of
FindDummies since unique indices are used for every dummy pair.
FindDummies[expr,\"indexlist\"] exports list of indices to indexlist, default is $IndexList.

Options:
writeDummyDefs->True  - If False, $IndexList will only contain the free index. In addition
			dummies will not be assigned DIndexQ = True UpValue. Option is mainly
			for internal use.";
writeDummyDefs::usage = "Option for FindDummies, default is True";


ClearAll[SymmetryCanonicalize];

SymmetryCanonicalize::usage = "SymmetryCanonicalize[expr]
Canonicalizes every term in a sum based on tensor symmetries and anti-symmetries. The function
works by converting expression into Mathematica's TensorContract notation and applying the
TensorReduce routine. 

All dummies must be identified by DIndexQ meaning FindDummies must be run first.
SymmetryCanonicalize[expr,\"exportlist\"] stores dummies in exportlist, default is $IndexList.";


ClearAll[freeIndexFind];

$IndexList::usage = "Default storage for the list of indices in an expression";
$IndexedObjectsList::usage = "List of all the defined indexed objects";

freeIndexFind::usage = "freeIndexFind[expr]
Returns the free index in expr";

d$m::usage = "Placeholder for indices in canonicalization. Made public since it may occasionally 
be printed";

$clearInternalDummiesAt::usage = "Once the internal dummy count gets to $clearInternalDummiesAt,
it will be restarted (at the end of the next ToCanonical run) to save storing too many definitions.

Default is 10000, very occasionally this may need to be increased if working with very large 
expressions with large substitutions and applying parts of ToCanonical seperately.";


ClearAll[printDummies];

$PrintNiceIndices::usage = "If set to False, prints indexed objects in input notation";

printDummies::usage = "printDummies[expr,\"exportlist\"]
Internal function to replace dummy indices with \"nice looking\" indices (defined 
in DefIndices) and apply pair-wise coloring so printed expressions are more easily read";

$PrintNiceDummies::usage = "If False, dummies will be printed with their true variable
names (d$1,d$2,...). printDummies[expr] must be run first";

$DummyColorList::usage = "List of colors to cycle through in printing dummy indices. Default is
{Blue, Red, Black, Purple, Brown, Gray, Green, Orange, Magenta, Pink}";


ClearAll[vectorForm];

vectorForm::usage = "vectorForm[expr]
Attempts to convert indexed expr into standard vector notation.
Currenly only works with expr derivatives up to fist order (vector 
notation can become cumbersome with high order derivative tensors)";



ClearAll[$vecFreeI,dot,cross,grad,div,curl];

$vecFreeI::usage = "Default free index to use in vector expressions. On loading the package,
this is set to n";

dot::usage = "dot[expr1,expr2]
Dot product for vector input. Works with one rank 2 tensor, but will give a 
warning if this is not a gradient tensor.
For shorthand input use \[EscapeKey]v.\[EscapeKey].
dot[expr1,expr2,ind] sets ind as the free index (if expr1.expr2 has a free index).";

cross::usage = "cross[expr1,expr2]
Cross product for vector input.
For shorthand input use \[EscapeKey]v*\[EscapeKey].
cross[expr1,expr2,ind] sets the free index to ind";

grad::usage = "grad[expr]
Gradient for vector input. Works with both scalar and vector objects.
For shorthand input use \[EscapeKey]vg\[EscapeKey].
grad[expr,ind] sets free index to ind (this should be a list of length 2 if expr is a vector).";

gradv::usage = "gradv[expr]
Gradient in v co-ordinate for vector input. Works with both scalar and vector objects.
For shorthand input use \[EscapeKey]vgv\[EscapeKey].
gradv[expr,ind] sets free index to ind (this should be a list of length 2 if expr is a vector).";


div::usage = "div[expr]
Divergence for vector input."

curl::usage = "curl[expr]
Curl for vector input. 
curl[expr,ind] sets free index to ind.";



ClearAll[PDExpand];

PDExpand::usage = "PDExpand[expr]
Expands out products in derivatives and concatenates nested PD expressions.";


ClearAll[deltaSimplify,leviCSimplify,unitVecReduce,applyUserRules,defineUserRules];

deltaSimplify::usage = "deltaSimplify[expr]
Canonicalizes \!\(\*SubscriptBox[\(\[Delta]\), \(ij\)]\) in expr by applying the summation convention rules associated with \!\(\*SubscriptBox[\(\[Delta]\), \(ij\)]\)";

leviCSimplify::usage = "leviCSimplify[expr]
Canonicalizes \!\(\*SubscriptBox[\(\[CurlyEpsilon]\), \(ijk\)]\) in expr by applying the summation convention rules associated with \!\(\*SubscriptBox[\(\[CurlyEpsilon]\), \(ijk\)]\)";

unitVecReduce::usage = "unitVecReduce[expr]
Applies rules to reduce expressions involving unit vectors (identified by UnitVecQ).
Only works up to second order derivatives at the moment.";

applyUserRules::usage = "applyUserRules[expr]
Applies rules as defined with defineUserRules. Mainly just a convenience so rules can be
included in ToCanonical.";

defineUserRules::usage = "defineUserRules[rules]
Define rules to be applied in applyUserRules (as part of ToCanonical). Rules should be a 
single rule or a list of rules. If extra dummies are needed on the right hand side of the 
rule, any defined indices may be used.";

$UserRuleList::usage = "List of the currently defined user rules";


ClearAll[ToCanonical];

ToCanonical::usage = "ToCanonical[expr]
Canonicalizes indexed expression expr. All objects must have been previously defined (with 
DefObject) and all indices with DefIndices.

ToCanonical is an application of
	 printDummies[
		SymmetryCanonicalize[
			unitVecReduce[
				applyUserRules[
					deltaSimplify[
						leviCSimplify[
							FindDummies[
								PDExpand[expr]]]]]]]

Options:
printProgress->False - Set to True for screen output showing compuation progress.
exportlist->\"$IndexList\" - Varaiable to store list of indices (passed as string).";	
printProgress::usage = "Option for ToCanonical, default is False";
exportlist::usage = "Option for ToCanonical, default is \"$IndexList\"";										



ClearAll[setIndEqual];

setIndEqual::usage = "setIndEqual[lhs,expr]
Sets lhs equal to expr, where expr may contain indices. lhs[ind] will then evaluate to
whatever specific ind is given at the time of execution with unique indices for all 
dummies.
\[EscapeKey]i=\[EscapeKey] produces \!\(\*OverscriptBox[\(=\), \(ind\)]\)  which is used for shorthand input as
lhs \!\(\*OverscriptBox[\(=\), \(ind\)]\) expr";


ClearAll[FullSimplifyVectorForm];

FullSimplifyVectorForm::usage = "FullSimplifyVectorForm[expr]
Attempts to find a simpler form for expr, judged by total number of terms in the sum.

The general technique is to apply \!\(\*SubscriptBox[\(\[CurlyEpsilon]\), \(irs\)]\)\!\(\*SubscriptBox[\(\[CurlyEpsilon]\), \(jrs\)]\)=2\!\(\*SubscriptBox[\(\[Delta]\), \(ij\)]\) identity to each index of each term in
the sum, expanding \[CurlyEpsilon] tensors in different orders from that in which they were initially 
substituted. This process generates a series of equivalent forms for each term. These
forms are then substituted into expression in all combinations that have a possibility 
of generating cancellations (up to a user specified maximum). The shortest expression
obtained in this substitution process is then passed back in and the process whole repeated.
This iteration stops when the expression no longer changes.

FullSimplifyVectorForm does not seem to change any expression that does not contain \!\(\*SubscriptBox[\(\[CurlyEpsilon]\), \(ijk\)]\) tensors.
Work is ongoing to understand how these expressions can be simplified in a systematic way. 

Options:
unitVec->False - Include a unit vector in simplification process, e.g. unitVec->b. At the present
		time this option multiplies expr by unitvec^2=1 and performs the expansion. If option
		numLeviPairs is more than one, unitvec^2=1 is also inserted multiple times.
		Currently only works for one unit vector but could easily be extended.
printSubs->False - If True, prints out all possible equivalent forms found for each term in expr
numLeviPairs->1  - Number of \!\(\*SubscriptBox[\(\[CurlyEpsilon]\), \(ijk\)]\) pairs to be inserted and re-expanded. Increasing this will often
		   cause many more possible forms to be found at the cost of compuation time. All
		   possible expansions of the \!\(\*SubscriptBox[\(\[CurlyEpsilon]\), \(ijk\)]\) tensors are tried on every possible combination
		   of indices.
maxPatternsToSub->100000 - Maximum number of substitutions to try for any given group of terms.
randomFinalForm->0 - There may be multiple possible forms of the final expression, each with the 
			 same number of terms. Set this option to >0 to try some of the other forms 
			 generated. Idea is to possibly remove optimatization procedure from a local
			 minimum, but so far haven't found a concrete example where this  is helpful.";

unitVec::usage="Option for FullSimplifyVectorForm, default is False";
printSubs::usage="Option for FullSimplifyVectorForm, default is False";
numLeviPairs::usage="Option for FullSimplifyVectorForm, default is 1";
maxPatternsToSub::usage="Option for FullSimplifyVectorForm, default is 100000";
randomFinalForm::usage="Option for FullSimplifyVectorForm, default is 0";


ClearAll[CheckTensorZero];

CheckTensorZero::usage = "CheckTensorZero[expr]
Checks if expr evaluates to zero by explicitly substituting in arrays and summing over indicies.
Unit vectors identities are automatically applied for applicable objects.
Returns True if expression is zero and expanded expression otherwise. 
Options:
simplifyFunction->Simplify - specifies function used to simplify the final expression. 
				Together, Expand or PossibleZeroQ may be substatially faster than
				default. Use N for numerical expressions. 
userFormSpec\[Rule]{} - allows expressions to be substituted in by the user to account for rules or 
		  particular vector forms. 
		  Specify as {{obj1, (replacement function)1},{obj2, (replacement function)2},...}, 
		  e.g. {{b[i],#/.#\[LeftDoubleBracket]3\[RightDoubleBracket]->Sqrt[1-#\[LeftDoubleBracket]2\[RightDoubleBracket]^2-#\[LeftDoubleBracket]1\[RightDoubleBracket]^2]&},...} explicitly specifies b as a unit 
		  vector, {{PD[b[i],Xd[j],Vd[]],#/.#->RandomReal[{-5,5},{3,3}]&},...} specifies numerical
		  \!\(\*SubscriptBox[\(b\), \(i, j\)]\). Note that any indices may be used in obj (i.e., PD[b[k],Xd[k],Vd[]] is 
		  equivalent in the previous example.
subsetSearch->False - If numerical forms are given (with userFormSpec) searches the resulting numerical
			  list for subsets that are zero. Specification is given as for the 
			  Mathematica Subset function e.g., 8 (all subsets up to length 8), or 
			  {2,4} (subsets of lengths 2-4). Will throw an error for non-numeric 
			  lists.";
simplifyFunction::usage="Option for CheckTensorZero, default is Simplify";
userFormSpec::usage="Option for CheckTensorZero, default is {}";
subsetSearch::usage="Option for CheckTensorZero, default is False";


ClearAll[myTiming];
myTiming::usage = "myTiming[func]
Prints out func name and time for execution and returns output of func";


(* ::Section:: *)
(*Main package*)


Begin["`Private`"]


(* ::Subsection::Closed:: *)
(*Loading packages, version checking and welcome message*)


Block[{Notation`AutoLoadNotationPalette},
	Notation`AutoLoadNotationPalette=False;
	Needs["Notation`"];
]

If[TrueQ[$VersionNumber<9.0],
	Print["VEST requires Mathematica 9.0 or later!"];Abort[]];

Print[
"-------------------------------------------------------
VEST: a Mathematica vector calculus package
--------------------------------------------------------

jsquire@princeton.edu
"]


(* ::Subsection::Closed:: *)
(*General private functions for minor tasks*)


(*Asigns func[var]=True to each of the variables in vars, function returns False for everything else*)
ClearAll[addTFUpValues];
addTFUpValues[vars_List, func_] :=(#/:func[#]=True)&/@vars;
  	
(* Remove UpValues, check that it exists first *)
(* Used to have a warning here if an UpValue did not exist. Removed this as it didn't
seem to be a big deal*)
ClearAll[removeTFUpValues];
removeTFUpValues[vars_List, func_] := Module[{nvars,noexist},
  	noexist = Cases[func /@ vars, Except[True]];
  	If[Length[noexist]!=0 ,nvars=Select[vars,func[#]===True &],nvars=vars];
  	(#/:func[#]=.)&/@nvars;
  ]

(* Given any sequence of variables symN[i,j,..] prints them out concatenated into
a single string*)
ClearAll[symN];
symN[expr__] := StringJoin[SymbolName /@ {expr}];
symN[] := Sequence[];

(*Create RowBox with indices*)
ClearAll[indBoxCreate];
indBoxCreate[inds__]:=If[And@@NumericQ/@{inds},RowBox[{StringJoin[Flatten[{ToString[#]&/@{inds}}]]}],
		RowBox[MapThread[StyleBox[#1,FontColor->#2]&,{SymbolName/@{inds},indColor/@{inds}}]]]
indBoxCreate[]=RowBox[{}];

(*Create string of indices, used if $PrintNiceIndices is set to off*)
indsToString=("["<>If[StringLength@#!=0,StringDrop[#,-1],#]&@\
	StringJoin[If[Head[#]===Symbol,SymbolName[#],ToString[#]]<>","&/@#]<>"]")&

ClearAll[indColor]
indColor[_]=Black;

(* In outer expression, puts sum into a list structure*)
(* Expand can be very slow, test if expression is expanded first*)
ClearAll[plusToList];
plusToList[expr_] := Block[{pluslist},
	pluslist=Position[expr,Plus,Heads->True];
	If[pluslist==={{0}},List@@expr,\
	If[pluslist==={},{expr},List@@Expand[expr]]]
	];
ClearAll[timesToList];
timesToList[expr_]:= If[Head[expr] === Times,List@@expr, {expr}];


(*Internal variable to keep track of maximum dummy for setIndEqual*)
If[Head[$maxDummyVal]===Symbol,$maxDummyVal=1];

(*Generates a given number of dummy indices in order*)
ClearAll[dumUnique];
dumUnique[n_,start_:1] := Table[ToExpression["d$"<>ToString[i]],{i,start,start+n-1}]

(* Annoying to write out ConstantArray every time*)
ClearAll[conA];
conA[var_,num_]:=ConstantArray[var,num];


General::generalTensor="No general tensor expression support. Some aspects 
of the canonicalizer will not function correctly.";


(* ::Subsection:: *)
(*Initialization and definition functions*)


(*Functions to define objects and indices *) 
DefObject::vecsym="Symmetry cannot be defined for a scalar or vector object";
DefObject::allinsym="All cannot be used as a keyword in VEST symmetry specification, use a list of slots";

DefObject[obj_,n_Integer,dep_List,printas_,syms_:False]:=Module[{},
unDefObject[obj];
obj/:IObjectQ[obj]=True;
obj/:IObjectQ[obj@@conA[_?IndexQ,n]]=True;
(*Could be nice to add in a check for number of arguments*)
(* Partial derivative dependencies *)
If[Not[dep[[1]]],obj/:PD[obj[conA[_?IndexQ,n]/.List->Sequence],\
	Xd[i__/;(And@@NumericQ/@{i})=!=True],Vd[___]]=0];
If[Not[dep[[2]]],obj/:PD[obj[conA[_?IndexQ,n]/.List->Sequence],\
	Xd[___],Vd[i__/;(And@@NumericQ/@{i})=!=True]]=0];
If[n==0 && Not[dep[[1]]],obj/:PD[obj,Xd[__],Vd[___]]=0];
If[n==0 && Not[dep[[2]]],obj/:PD[obj,Xd[___],Vd[__]]=0];
(* Rank checking functions *)
obj/:IObjectRank[obj]:=n;
obj/:IObjectRank[obj@@conA[_?IndexQ,n]]:=n;
(* Printing *)
obj/:MakeBoxes[obj[ind___],StandardForm]:=If[((And@@IndexQ/@{ind}===True||And@@NumericQ/@{ind}===True||\
		$InitObjects)&&IObjectRank[obj]===Length[{ind}])&&$PrintNiceIndices,\
	RowBox[{SubscriptBox[printas,indBoxCreate[ind]]}],\
	RowBox[{ToString[obj]<>indsToString[{ind}]}]];
If[n!=1,obj/:objPrintAs[obj]=printas,obj/:objPrintAs[obj]=Style[printas,Bold]];
(*Symmetries for tensor objects*)
If[syms=!=False,If[n<=1,Message[DefObject::vecsym],\
	If[Length[Position[syms,All]]!=0,Message[DefObject::allinsym],obj/:objectSymmetry[obj]=syms]]];
Print["Defined ",  Subscript[printas, StringJoin[Take[{"i","j","k","l","p","q","r","s"},n]]]];
If[!$InitObjects,
	AppendTo[$IndexedObjectsList,{obj,n}];
(*Add to $Assumptions for Mathematica tensor operations, go to second order, code will add
more when necessary*)
	addObjDomainsToAssumptions[{obj,n},0];
	addObjDomainsToAssumptions[{obj,n},1];
	addObjDomainsToAssumptions[{obj,n},2];]
];
$MaxStoredTensorDeriv=2;

objectSymmetry[_]=False;

unDefObject[obj_]:=Module[{},
	$IndexedObjectsList=DeleteCases[$IndexedObjectsList,{obj,_}];
	$Assumptions=DeleteCases[$Assumptions,obj[___]\[Element]_];
	$Assumptions=DeleteCases[$Assumptions,PD[obj[___],Xd[___],Vd[___]]\[Element]_];
ClearAll[obj];
]

(*Remove old object definitions for situations in which package is just reloaded*)
(*$IndexedObjectsList is a list of the indexed objects*)
If[Head[$IndexedObjectsList]=!=Symbol,
	unDefObject[#[[1]]]&/@$IndexedObjectsList,$IndexedObjectsList={}];



$PrintNiceIndices=True;
$NumDummyColors=10; (*Max is 10*)

(*$InitObjects changes printing and disables various parts of DefObject*)
$InitObjects=True;

If[Head[$FullIndexList]=!=Symbol,
		ClearAll/@$FullIndexList];
$FullIndexList={};

DefIndices::disall="Disallowed index `1`";
Options[DefIndices]={useForPrinting->True};
(*For temporary indices that wish to be added (that are not nice symbols) set to False*)
SetAttributes[DefIndices,Listable];
DefIndices[inds_,OptionsPattern[]]:=Module[{},
	If[Head[inds]=!=Symbol,Message[DefIndices::disall]];
	addTFUpValues[{inds},IndexQ];
(*For use in printing with dummy indices*)
	If[OptionValue[useForPrinting],$FullIndexList=Join[$FullIndexList,{inds}]];
];

(*Defining constant objects, all numbers are automatically constant*)
SetAttributes[DefConst,Listable];
DefConst[obj_]:=Module[{},
If[!TrueQ[ConstQ[obj]],obj/:ConstQ[obj]=True];
$Assumptions=Join[$Assumptions,{obj\[Element]Arrays[{},Reals,{}]}];
Print[obj, " defined as constant"];
];
ConstQ[_?NumericQ]=True;


Protect[d$m];(*Used for marking indices*)


(* ::Subsubsection::Closed:: *)
(*Mathematica tensor specifications*)


$Assumptions={};
(*Adds given derivatives of object at given order to $Assumptions list*)
ClearAll[addObjDomainsToAssumptions];
addObjDomainsToAssumptions[obj_List,derivorder_]:=Module[{ind,xvdep,derivlist},
(*Decide x, v dependence*)
DefIndices[ind,useForPrinting->False];xvdep={0,0};
If[PD[obj[[1]]@@conA[ind,obj[[2]]],Xd[],Vd[ind]]===0,xvdep[[2]]=False,xvdep[[2]]=True];
If[PD[obj[[1]]@@conA[ind,obj[[2]]],Xd[ind],Vd[]]===0,xvdep[[1]]=False,xvdep[[1]]=True];
Which[xvdep=={True,True},
		derivlist=Select[Flatten[Outer[List,Range[0,derivorder],\
			Range[0,derivorder]],1],#[[1]]+#[[2]]==derivorder &],
	xvdep=={True,False}||(xvdep=={False,False}&&derivorder==0),\
		derivlist={{derivorder,0}},
	xvdep=={False,True},
		derivlist={{0,derivorder}},
	xvdep=={False,False},
		derivlist={}];
$Assumptions=Join[$Assumptions,objTensorDomains[obj,#]&/@derivlist];
]


ClearAll[objTensorDomains];
objTensorDomains[obj_List,derivs_List]:=Module[{objS,n,symprop,xindend,vindend,outelem},
objS=obj[[1]];
n=obj[[2]];
symprop=If[objectSymmetry[objS]===False,Hold[Sequence[]],objectSymmetry[objS]];
If[derivs=!={},
	xindend=n+derivs[[1]];vindend=xindend+derivs[[2]];
	outelem=PD[objS@@conA[d$m,n],Xd@@conA[d$m,derivs[[1]]],Vd@@conA[d$m,derivs[[2]]]]\[Element]\
		Arrays[conA[3,n+Total[derivs]],Reals,{ReleaseHold[symprop],Symmetric[Range[n+1,xindend]],\
		Symmetric[Range[xindend+1,vindend]]}],
(*For objects with no dependence on anything*)
outelem={}];

outelem]


(*Adding higher order derivatives if an expression is encountered,
Just add for all currently defined objects*)
ClearAll[IncreaseStoredTensorDeriv];
IncreaseStoredTensorDeriv[maxderiv_]:=Module[{},
Do[addObjDomainsToAssumptions[#,orderit]&/@$IndexedObjectsList,\
	{orderit,$MaxStoredTensorDeriv+1,maxderiv}];
$MaxStoredTensorDeriv=maxderiv;
Print["Higher derivatives added to tensor assumptions list"];
]




(* ::Subsubsection:: *)
(*Partial derivative properties*)


SetAttributes[Xd,Orderless];
SetAttributes[Vd,Orderless];
PD[obj_,Xd[],Vd[]]:=obj;
(*PD is linear*)
PD[obj1_+obj2_,Xd[ind1___],Vd[ind2___]]:=\
	PD[obj1,Xd[ind1],Vd[ind2]]+PD[obj2,Xd[ind1],Vd[ind2]];
(*Behaviour with constants*)
PD[obj_?ConstQ^pow_.,Xd[___],Vd[___]]:=0;
PD[objc_?ConstQ^pow_. obj_,Xd[ind1___],Vd[ind2___]]:=objc^pow PD[obj,Xd[ind1],Vd[ind2]];


(* Partial derivative output form *)
PD/:MakeBoxes[PD[obj1_[oind___],Xd[indx__],Vd[indv__]],StandardForm]:=
	If[$PrintNiceIndices&&((And@@IndexQ/@{indx,indv}===True||\
			And@@NumericQ/@{indx,indv}===True)&&IObjectRank[obj1]===Length[{oind}]),\
		RowBox[{SubscriptBox[MakeBoxes[obj1[oind],StandardForm],RowBox[{",",indBoxCreate[indx],";",indBoxCreate[indv]}]]}],\
		If[obj1===Times||obj1===Plus,RowBox[{ToString[PD[obj1[oind],Xd[indx],Vd[indv]]]}],
			RowBox[{"PD["<>ToString[obj1]<>indsToString[{oind}]<>",Xd"<>indsToString[{indx}]<>",Vd"<>indsToString[{indv}]<>"]"}]]
];
PD/:MakeBoxes[PD[obj1_[oind___],Xd[indx__],Vd[]],StandardForm]:=
	If[$PrintNiceIndices&&((And@@IndexQ/@{indx}===True||\
			And@@NumericQ/@{indx}===True)&&IObjectRank[obj1]===Length[{oind}]),\
		RowBox[{SubscriptBox[MakeBoxes[obj1[oind],StandardForm],RowBox[{",",indBoxCreate[indx]}]]}],
		If[obj1===Times||obj1===Plus,RowBox[{ToString[PD[obj1[oind],Xd[indx],Vd[]]]}],
			RowBox[{"PD["<>ToString[obj1]<>indsToString[{oind}]<>",Xd"<>indsToString[{indx}]<>",Vd[]]"}]]
];
PD/:MakeBoxes[PD[obj1_[oind___],Xd[],Vd[indv__]],StandardForm]:=
	If[$PrintNiceIndices&&((And@@IndexQ/@{indv}===True||\
			And@@NumericQ/@{indv}===True)&&IObjectRank[obj1]===Length[{oind}]),\
		RowBox[{SubscriptBox[MakeBoxes[obj1[oind],StandardForm],RowBox[{";",indBoxCreate[indv]}]]}],\
		If[obj1===Times||obj1===Plus,RowBox[{ToString[PD[obj1[oind],Xd[],Vd[indv]]]}],
			RowBox[{"PD["<>ToString[obj1]<>indsToString[{oind}]<>",Xd[],Vd"<>indsToString[{indv}]<>"]"}]]
];



(* ::Subsubsection::Closed:: *)
(*Pre-loaded object definitions*)


(* delta properties *)
DefObject[delta,2,{False,False},"\[Delta]"]
SetAttributes[delta,Orderless]
delta[ind_?IndexQ,ind_?IndexQ]=3;
$Assumptions=Join[$Assumptions,{delta[d$m,d$m]\[Element]Arrays[{3,3},Reals,Symmetric[{1,2}]]}];
(*Have to re-assign usage messages as they are cleared by DefObject*)
delta::usage = "delta[i,j] is the identity operator";

(* leviC properties *)
DefObject[leviC,3,{False,False},"\[CurlyEpsilon]"]
leviC[ind1_,ind2_?IndexQ,ind2_?IndexQ]=0;
leviC[ind2_?IndexQ,ind1_,ind2_?IndexQ]=0;
leviC[ind2_?IndexQ,ind2_?IndexQ,ind1_]=0;
$Assumptions=Join[$Assumptions,{leviC[d$m,d$m,d$m]\[Element]Arrays[{3,3,3},Reals,Antisymmetric[{1,2,3}]]}];
leviC::usage = "leviC[i,j,k] is the totally antisymmetric Levi-Civita tensor";

(* v properties *)
DefObject[v,1,{False,True},"v"];
v/:PD[v[j_],Xd[],Vd[k_]]:=delta[j,k];
v/:PD[v[j_],Xd[],Vd[k_,l__]]:=0;
$Assumptions=Join[$Assumptions,{v[d$m]\[Element]Arrays[{3},Reals,{}]}];
v::usage = "v[j] is the co-ordinate v";

(* x properties *)
DefObject[x,1,{True,False},"x"];
x/:PD[x[j_],Xd[k_],Vd[]]:=delta[j,k];
x/:PD[x[j_],Xd[k_,l__],Vd[]]:=0;
$Assumptions=Join[$Assumptions,{x[d$m]\[Element]Arrays[{3},Reals,{}]}];
x::usage = "x[j] is the co-ordinate x";



$InitObjects=False;


(* ::Subsection::Closed:: *)
(*Debug/Timing functions*)


(*Written by Leonid Shiffrin, published on Mathematica stack exchange*)

ClearAll[debug];
SetAttributes[debug, HoldAll];
debug[code_] :=
 Internal`InheritedBlock[{Message},
   Module[{inMessage},
     Unprotect[Message];        
     Message[args___] /; ! MatchQ[First[Hold[args]], _$Off] :=
       Block[{inMessage = True},
         Print[{
            Shallow /@ Replace[#, HoldForm[f_[___]] :> HoldForm[f], 1],
            Style[Map[Short, Last[#], {2}], Red]
           } &@Drop[Drop[Stack[_], -7], 4]
         ];
         Message[args];
         Throw[$Failed, Message];
       ] /; ! TrueQ[inMessage];
    Protect[Message];
   ];
   Catch[StackComplete[code], Message]]


(* Useful as a timing call as it can be used anywhere without disrupting evaluation*)
(* Prints out function name and timing *)
SetAttributes[myTiming,HoldAll];
myTiming[f_]:=Module[{outlist},
	$myTimingIter=0;
	outlist=Timing[Reap[f,timingtag]];
	Print[Head[Unevaluated[f]],": ",outlist[[1]] "s"];
	If[Length[outlist[[2,2]]]!=0,
		Print[Flatten[outlist[[2,2]]][[1]]," total: ",Total[Drop[Flatten[outlist[[2,2]]],1]]]];
	outlist[[2,1]]
];

(*Incude a myTimingSum call within a myTiming call and it will print out the sum
of timings for all calls of a particular function*)
ClearAll[myTimingSum];
SetAttributes[myTimingSum,HoldAll]
myTimingSum[f_]:=Module[{outlist},
	If[$myTimingIter==0,Sow[Head[Unevaluated[f]],timingtag]];
	$myTimingIter++;
	outlist=Timing[f];
	Sow[outlist[[1]],timingtag];
	outlist[[2]]
];


(* ::Subsection::Closed:: *)
(*General clean-up functions*)


General::IndexList="IndexList variable should be passed as a string";
(* Cleans up indexlist, removing all dummies and deleting indlist *)
ClearAll[CleanIndexList];
CleanIndexList[indliststr_]:=Module[{dummyunion},
If[!AtomQ[ToExpression[indliststr]],
	If[!StringQ[indliststr],Message[CleanIndexList::IndexList]];
	dummyunion=DeleteDuplicates[Flatten[ToExpression[indliststr<>"[[1]]"]]];
	removeTFUpValues[dummyunion,DIndexQ];
	(*Print["All dummy definitions reset"];*)
	ToExpression["Clear["<>indliststr<>"]"];
]
];
ClearAll["$IndexList"];

 (*Value at which to apply clearHighDummyDefs*)
$clearInternalDummiesAt=10000;

ClearAll[clearHighDummyDefs];
(*To save memory delete d$ definitions above some value. Called from within ToCanonical*)
clearHighDummyDefs[nmax_]:=Remove/@Flatten[StringCases[Names["d$*"],"d$"~~n__/;ToExpression[n]>nmax]];


(* ::Subsection::Closed:: *)
(*Dummy index identification and checking of expression validity*)


(* Find the dummy indices and check that x is a valid expression *)

(*If FindDummies is used as part of setIndEqual, don't write dummies to DIndexQ or $IndexList*)
(*If False, exportlist is just a list of the free indices*)
ClearAll[FindDummiesSingle];
FindDummies::indmatch="Index terms do not match across sum!";
FindDummies::dummies="Too many instances of a dummy index `1`";
Options[FindDummies]={writeDummyDefs->True};

FindDummies[expr_,exportlist_:"$IndexList",OptionsPattern[]]:= 
Module[{xplist,dNstart,indexlist,indexlistunique,xwithdummies,\
	dummylistunique,alldummies,freelistunique},
(*Cleaning previous definitions, exportlist variable is passed around as a string*)
If[OptionValue[writeDummyDefs],
	If[Length[ToExpression[exportlist]]!=0,CleanIndexList[exportlist]]];
(*Want to generate unique indices fo FindDummies to avoid problems when substituting
expressions*)
If[TrueQ[OptionValue[writeDummyDefs]],dNstart=1,dNstart=$maxDummyVal+1];
(*Expand and form list*)
xplist=expandForFindDummies[expr];
(*These two lines allow general vector input by switching obj1\[Rule]obj1[n] and obj0\[Rule]obj0[]*)
xplist=(xplist/.obj_/;IObjectRank[obj]===0:>obj[])/.obj_[][]:>obj[];
xplist=((xplist/.obj_/;(IObjectRank[obj]===1&&!MatchQ[obj,_[___]]):>\
		obj[$vecFreeI])/.obj_[i1_][i2_]:>obj[i2]);

xplist=plusToList[xplist];
(* Find dummies for each term *)
indexlist=Map[(FindDummiesSingle[#,exportlist]&),xplist];

If[Not[SameQ@@(Transpose[indexlist][[2]])],\
	Message[FindDummies::indmatch];CleanIndexList[exportlist];Abort[]];
(*Generate unique dummy indices*)
indexlistunique=Map[(#/.{l1_List,l2_List}:>{dumUnique[\
	Length[l1],(dNstart+=Length[l1])-Length[l1]],l2}&),indexlist];
(*Storing the maximum value of any dummy index in the entire session*)
If[$maxDummyVal<dNstart,$maxDummyVal=dNstart-1];

xwithdummies=MapIndexed[ReplaceAll[#1,Thread[indexlist[[First[#2],1]]->\
	indexlistunique[[First[#2],1]]]]&,xplist]/.List->Plus;

dummylistunique=Transpose[indexlistunique][[1]];
freelistunique=Union[Flatten[Transpose[indexlistunique][[2]]]];
(*Set IndexQ to True for each of the dummies*)
alldummies=DeleteDuplicates[Flatten[dummylistunique]];
addTFUpValues[alldummies,IndexQ];
(*Set DindexQ to True for each of the dummies*)
If[OptionValue[writeDummyDefs],
	addTFUpValues[alldummies,DIndexQ];
	(* To cancel terms, use non-unique terms for deltaSimplfy etc.*)
	If[Length[freelistunique]>1,Message[FindDummies::generalTensor]];
	ToExpression[exportlist<>"={"<>ToString[dummylistunique]<>","<>ToString[freelistunique]<>"}"];,
(*Note, the exportlist is backwards if writeDummyDefs is off*)
	ToExpression[exportlist<>"="<>ToString[freelistunique]];
];

If[$printToCanonicalProgress && OptionValue[writeDummyDefs],\
	Print["FindDummies: Dummy indices identified"]];
xwithdummies
];



(* Private function to find the dummies over a single term *)
FindDummiesSingle[expr_,exportlist_]:=Module[{xwsq,vecfreeind,ilist,itally,itallydum,itallyfree},
xwsq=expr/.obj_^2:>obj**obj;
ilist=Cases[xwsq,_?IndexQ,-1];
(*Form a tally of each index*)
itally=Tally[ilist];

(*Check for correct number of dummy occurences*)
If[Max[itally[[All,2]]]>2||Count[{expr},f_[i__]^n_/;n>=3]>0,\
	Message[FindDummies::dummies,itally];\
	Print[xwsq];CleanIndexList[exportlist];Abort[]];
(*Remove Tallies from the result*)
itallydum=Select[itally,#[[2]]==2&];
itallyfree=Select[itally,#[[2]]==1&];
itally={If[Length[itallydum]!= 0,Transpose[Sort[itallydum]][[1]],{}],\
	If[Length[itallyfree]!= 0,Transpose[Sort[itallyfree]][[1]],{}]};

itally
];


ClearAll[expandForFindDummies];
General::distExpand="large expression contains powers, using Expand";
(*Expand can be very slow, similar to for deltaSimplify. Use a separate function
here just to have slightly different checks*)
(*Uses FixedPoint around whole function since sometimes doesn't finish on one pass (if head is Times)*)
expandForFindDummies[expr_]:=Module[{},
FixedPoint[
	If[!(Head[#]===Plus||Head[#]===Times),Expand[#],\
	(*If expression contains objects^2 that contain +, Distribute won't work*)
		If[Position[Cases[#,obj_^n_,Infinity],Plus]=!={},\
			Message[expandForFindDummies::distExpand];Expand[#],
			If[Head[#]===Plus,
				FixedPoint[(Plus@@Distribute/@List@@#)&,#],
				Distribute[#]]
	]]&,expr]
]




(*
(*Seems to work without needing this, keep just in case I find a problem*)

(*This was in FindDummiesSingle*)
(*Want to be able to do things like curl[b]+b, without any indices*)
(*This is slow, only run if vector functions have just been run. Reset after ToCanonical*)
If[$findDumsAllowVecObjects,
(*See if the term matches scalar*(vector without index), if it does, add placeholder index to be recognized later*)
	vecfreeind=xwsq/.mult_ obj_/;(IObjectRank[obj]===1&&!MatchQ[obj,_[___]])->True;
	If[Length[ilist]=!=1&&(TrueQ@vecfreeind||IObjectRank[xwsq]===1),\
		AppendTo[ilist,freeDummyPlaceholder]];
];

(*This was from FindDummies, just after FindDummiesSingle call*)
(*If terms do not match accros sum, try running it with $findDumsAllowVecObjects to true*)
(*This is necessary so that ToCanonical works with Map with vector objects, won't make any difference
for non-vector sums without mistakes*)
If[Not[SameQ@@(Transpose[indexlist][[2]])]&&!$findDumsAllowVecObjects,Print["Running second FindDummiesSingle"];
	$findDumsAllowVecObjects=True;indexlist=Map[(FindDummiesSingle[#,exportlist]&),xplist]];
(*If any of the individual terms contain freeDummyPlaceholder, run function to fix this*)
If[$findDumsAllowVecObjects,
	If[MemberQ[Flatten[indexlist],freeDummyPlaceholder],
		{xplist,indexlist}=fixObjWithoutIndex[xplist,indexlist]]
];

(*Extra function*)
ClearAll[fixObjWithoutIndex];
(*Returns xplist and indexlist with freeDummyPlaceholder replaced by correct index*)
(*Needed for nice vector input*)
fixObjWithoutIndex[xplistf_,indexlistf_]:=Module[{ilist,freeind,fDPpos},
	ilist=indexlistf[[All,2,1]];
(*If there are other objects, choose the free index from these*)
	freeind=If[SameQ@@#,If[Length[#]==0,$vecFreeI,#[[1]]],
		Message[FindDummies::indmatch]]&@Cases[ilist,Except[freeDummyPlaceholder]];
	fDPpos=Flatten@Position[ilist,freeDummyPlaceholder];
(*Function to put the index on the vector term*)
	addFreeInd=Function[xterm,If[Head[#]===List,Times[#[[1]],#[[2]][freeind]],#[freeind]]&@\
		((#/.mult_ obj_/;(IObjectRank[obj]===1&&!MatchQ[obj,_[___]]):>{mult ,obj})&@xterm)];
(*Replace the relevant bits of xlist and change placeholder to freeind*)
	{ReplacePart[xplistf,Thread[fDPpos->addFreeInd/@xplistf[[fDPpos]]]],
		indexlistf/.freeDummyPlaceholder->freeind}
];



(*OLD FUNCTION - now using Mathematica canonicalizer *)
(*Function to reorder dummy indices (so equal terms will be recognized by Mathematica).
Also resets $IndexList *)
dummyReorder[expr_,exportlist_:"$IndexList"]:=
	Module[{xplist,fulldlist,inindorder,outindorder,reorderrules,xsimpfor,xsimp,
		randoutindorder,randreorderrules,xsimprand},
xplist=plusToList[expr];
(* Form rule for set for reorder of the indices *)
fulldlist=Cases[{#},_?DIndexQ,Infinity]&/@xplist;
inindorder=DeleteDuplicates/@fulldlist;
If[Or@@MapThread[Length[#1]<Length[#2]/2&,{inindorder,fulldlist}],\
	Print["More than two instances of a dummy have appeared"];];
outindorder=dumUnique[Length[#]]&/@inindorder;
removeTFUpValues[Union[Flatten[inindorder]],DIndexQ];
addTFUpValues[Union[Flatten[outindorder]],IndexQ];
addTFUpValues[Union[Flatten[outindorder]],DIndexQ];
reorderrules=MapThread[MapThread[Rule,{#1,#2}]&,{inindorder,outindorder}];
(*Update $IndexList*)
ToExpression[exportlist<>"[[1]]="<>ToString[outindorder]];
(*Reorder dummies in each term *)
xsimpfor=Plus@@(MapThread[ReplaceAll[#1,#2]&,{xplist,reorderrules}]\
	/.leviC[inds__]:> Signature[{inds}]leviC[Sort[{inds}]/.List->Sequence]);

(* Now do in reverse to catch symmetry and anti-symmetry*)
xsimpfor=plusToList[xsimpfor];
inindorder=DeleteDuplicates[Reverse[Cases[{#},_?DIndexQ,Infinity]]]&/@ xsimpfor;
outindorder=dumUnique[Length[#]]&/@inindorder;
reorderrules=MapThread[MapThread[Rule,{#1,#2}]&,{inindorder,outindorder}];
(*Update $IndexList*)
ToExpression[exportlist<>"[[1]]="<>ToString[outindorder]];
(*Reorder dummies in each term *)
xsimp=Plus@@(MapThread[ReplaceAll[#1,#2]&,{xsimpfor,reorderrules}]\
	/.leviC[inds__]:> Signature[{inds}]leviC[Sort[{inds}]/.List->Sequence]);
xsimp
];

*)


(* Check if expression has n free indices *)
nFIndexQ[obj_,n_Integer]:=Module[{nFIndexQlist},
	FindDummies[obj,ToString[nFIndexQlist],writeDummyDefs->False];
	Length[nFIndexQlist]===n
];



(*Return the free index of expression as a list*)
freeIndexFind[expr_]:=Module[{ilist},
FindDummies[expr,ToString[ilist],writeDummyDefs->False];
ilist]


(* ::Subsection::Closed:: *)
(*Replace Indexed Object*)


(*Use standard assignment, will work provided dummies aren't reset
(by running SymmetryCanonicalize) in between calls to expr *)
setIndEqual::indmatch="`1` supplied with incorrect number of indices!";

SetAttributes[setIndEqual,HoldFirst];
setIndEqual[lhs_,expr_/;(Head[expr]=!=List)]:=Module[{flistname},
	ClearAll[lhs];
	If[expr===0,lhs[___]=0;, \
			lhs[ind___?(And@@IndexQ/@{##}&)]:=\
		FindDummies[expr,ToString[flistname],writeDummyDefs->False]/.\
			If[Length[flistname]==Length[{ind}],Thread[flistname:>{ind}],\
				Message[setIndEqual::indmatch,lhs];Abort[]]
		]
];

Notation[
	ParsedBoxWrapper[RowBox[{"lhs_", " ", OverscriptBox["=", "ind"], "expr_"}]] \
	\[DoubleLongLeftRightArrow] ParsedBoxWrapper[RowBox[{"setIndEqual", "[", RowBox[{"lhs_", ",", "expr_"}], "]"}]]];
AddInputAlias[ParsedBoxWrapper[OverscriptBox["=", "ind"]],"i=",EvaluationNotebook[]];



(* ::Subsection::Closed:: *)
(*Canonicalizing with respect to tensor symmetries*)


(*This function canonicalizes each term in a sum using the new Mathematica 9 function
TensorReduce. It is performing the same task as dummyReorder but in a much more 
consistent way and should be guaranteed to put each term in a canonical form.
It is a bit slower though, so could be helpful to run dummyReorder sometimes (e.g. 
after deltaSimplify) as a first step to remove most of the terms.*)



ClearAll[singleSymmetryCanon];
SymmetryCanonicalize::newterms="Warning, new terms generated during summation, be worried!";
(*Execute for each term seperately, executes at the same speed*)
SymmetryCanonicalize[expr_,exportlist_:"$IndexList"]:=Module[{xplist,xsimp,dlistoutsimp},
(*Used to store dummies across entire sum*)
dlistin={};dlistout={};
xplist=expandForFindDummies[expr];
xplist=plusToList[xplist];
$progressIter=1;$totalterms=Length@xplist;
If[$printToCanonicalProgress,
	Print["SymmetryCanonicalize: progress..."];
	Monitor[xsimp=Map[singleSymmetryCanon,xplist];,\
		ProgressIndicator[Dynamic[$progressIter/$totalterms]]];
	Print[];,
	xsimp=Map[singleSymmetryCanon,xplist]];

(*Update dummy definitions and $IndexList, have to do this after simplifying expr*)
xsimp=plusToList[Plus@@xsimp];
dlistoutsimp=Union[Cases[#,Alternatives@@DeleteDuplicates[Flatten[dlistout]],-1]]&/@xsimp;
removeTFUpValues[DeleteDuplicates[dlistin],DIndexQ];
addTFUpValues[DeleteDuplicates[Flatten[dlistoutsimp]],IndexQ];
addTFUpValues[DeleteDuplicates[Flatten[dlistoutsimp]],DIndexQ];
If[Nor[Complement[dlistoutsimp,dlistout]=!={},Complement[dlistoutsimp,dlistout]=!={{}}],
	Message[SymmetryCanonicalize::newterms]];
ToExpression[exportlist<>"[[1]]="<>ToString[dlistoutsimp]];
Plus@@xsimp
]

(*Canonicalizes a single term not involving sums*)
singleSymmetryCanon[xsing_]:=
Module[{freeindices,tenlist,multfac,it,dlist,num2dummies,xsout,maxderiv,duminit,indpairs},
	freeindices=Cases[xsing,_?(DIndexQ[#]=!=True &&IndexQ[#] &),-1];
(*If the expression contains derivatives that haven't been initialized, do this*)
	maxderiv=Max[Cases[xsing,Xd[inds___]|Vd[inds___]:>Length[{inds}],Infinity]];
	If[maxderiv>$MaxStoredTensorDeriv,IncreaseStoredTensorDeriv[maxderiv]];
	(*Find dummies and sort into list for TensorContract*)
	xsout=TensorProduct@@timesToList[xsing/.obj_^2:>obj\[TensorProduct]obj];
	duminit=Cases[xsout,_?IndexQ,-1];
	indpairs=Select[DeleteDuplicates[Flatten[Position[duminit,#]]&/@duminit],Length[#]==2 &];
	xsout=xsout/._?IndexQ->d$m;
(*Store old dummy indices for removal*)
	dlistin=Join[dlistin,Complement[duminit,freeindices]];
(*Put tensor in canonical form*)
	tenlist=TensorContract[xsout,indpairs]//TensorReduce;
(*Put back into index notation*)
(*Trick to deal with TensorContract^n*)
	tenlist=tenlist/.obj_TensorContract^n_:>Nest[NonCommutativeMultiply[#,obj]&,obj,n-1];
(*Some objects are factored out of TensorContract*)
	multfac=Times@@DeleteCases[timesToList[tenlist],_TensorContract|_NonCommutativeMultiply,Infinity];
(*Convert into a list and number dummy appearances*)
	tenlist=Cases[timesToList[tenlist],TensorContract[iobj_,dlist_List]_.:>{iobj,dlist},Infinity];
	tenlist=(it=1;ReplaceAll[#,d$m:>{d$m,it++}])&/@tenlist;
(*Generate new dummy indices*)
	dlist=dumUnique[Length@Flatten[tenlist[[All,2]],1]];it=1;
	num2dummies=Map[Alternatives@@#1->dlist[[it++]]&,\
		(#/.num_Integer:>{d$m,num})]&/@tenlist[[All,2]];
	AppendTo[dlistout,dlist];
(*Put into a single product*)
	tenlist=multfac Times@@(MapThread[#1/.#2&,{tenlist[[All,1]],num2dummies}]/.TensorProduct->Times);
(*For the progress bar, should make sure this doesn't stuff up parallization, if there is any*)
$progressIter++;
(*Also need dlist for $IndexList*)
(*Sometimes free index is factored out of TensorContract, use two replaceAll's*)
	it=1;
	tenlist/.{d$m,_Integer}:>freeindices[[it++]]/.d$m:>freeindices[[it++]]
]





(* ::Subsection:: *)
(*Partial derivative rules*)


PDExpand::disall="Disallowed object `1` in partial derivative expansion";
ClearAll[expandPDrules,contractPDrules,productPDrules];
ClearAll[ExpandPDexpr,ContractPDexpr,ProductPDexpr];
(*Partial Derivative expansion functions *)
expandPDrules:={PD[obj_,Xd[ind1_,ind2__],Vd[ind3_,ind4__]]:>\
		PD[PD[obj,Xd[ind1],Vd[ind3]],Xd[ind2],Vd[ind4]],
	PD[obj_,Xd[indi___],Vd[indk_,indl__]]:> PD[PD[obj,Xd[indi],Vd[indk]],Xd[],Vd[indl]],
	PD[obj_,Xd[indi_,indj__],Vd[indk___]]:> PD[PD[obj,Xd[indi],Vd[indk]],Xd[indj],Vd[]],
	PD[obj_,Xd[indi_],Vd[indk_]]:> PD[PD[obj,Xd[indi],Vd[]],Xd[],Vd[indk]]};
contractPDrules:=PD[PD[obj_,Xd[indi___],Vd[indk___]],Xd[indj___],Vd[indl___]]:>\
	PD[obj,Xd[indi,indj],Vd[indk,indl]];
productPDrules:={PD[objf_ objg_,Xd[indi_],Vd[]]:>\
		objf PD[objg,Xd[indi],Vd[]]+objg PD[objf,Xd[indi],Vd[]],
	PD[objf_ objg_,Xd[],Vd[indi_]]:>objf PD[objg,Xd[],Vd[indi]]+objg PD[objf,Xd[],Vd[indi]],
	PD[objf_^n_,Xd[indi_],Vd[]]:>If[n!=2&&(IObjectRank[objf]=!=0),Message[PDExpand::disall,objf^n],\
		n objf^(n-1) PD[objf,Xd[indi],Vd[]]],
	PD[objf_^n_,Xd[],Vd[indi_]]:>If[n!=2&&(IObjectRank[objf]=!=0),Message[PDExpand::disall,objf^n],\
		n objf^(n-1) PD[objf,Xd[],Vd[indi]]]
	};

ExpandPDexpr[expr_]:=expr//.expandPDrules;
ContractPDexpr[expr_]:=expr//.contractPDrules;
ProductPDexpr[expr_]:=expr//.productPDrules;
PDExpand[expr_]:=Module[{xexp},
	xexp=ContractPDexpr[ProductPDexpr[ExpandPDexpr[expr]]];
	If[$printToCanonicalProgress,Print["PDExpand: completed expansion of partial derivatives"]];
	xexp
]


(* ::Subsection:: *)
(*Applying rules as part of ToCanonical*)


$UserRuleList={};
(*rules should be a non-nested list of rules*)
defineUserRules::norule="Input must be a list of Rule or RuleDelay objects";
defineUserRules::dummies="Odd number of instances of new dummy index on RHS of rules: `1`";
defineUserRules[rules_List]:=Module[{patinds,ruleRHS,indsRHS,dummies,newrules},
	(*Find inds used in the patterns themselves*)
	patinds=DeleteDuplicates[HoldPattern[#][[1,1]]&/@Cases[rules,obj_Pattern:>obj,Infinity]];
(*Check it is actually a rule*)
	ruleRHS=If[Head[#]===Rule||Head[#]===RuleDelayed,#[[2]]/.obj_^2:>obj**obj,\
		Message[defineUserRules::norule,#];Abort[]]&/@rules;
(*Find indices on the RHS and delete those from the rule*)	
	indsRHS=Cases[ruleRHS,_?IndexQ,Infinity];
	dummies=Cases[Tally[indsRHS],{ind_,n_}/;!MemberQ[Flatten[patinds],ind]];
(*Check for an even number, could let something bad through, but annoying to do a more thorough check*)
	dummies=If[TrueQ@(And@@(EvenQ/@dummies[[All,2]])),dummies[[All,1]],\
		Message[defineUserRules::dummies,dummies];Abort[]];
(*Change squares into NCM*)
	newrules=MapThread[(#1/.obj_?IObjectQ^n_/;n>1:>NonCommutativeMultiply@@conA[obj,n]):>#2 &,{rules[[All,1]],ruleRHS}];
	newrules=#/.Thread[dummies->(Hold[createInd[#]]&/@dummies)]&/@newrules;
	$UserRuleList=newrules;
];
	

applyUserRules[expr_]:=Module[{out},
(* Clear index definitions so all start at zero*)	
	ClearAll[createIndPrevious];
	out=((expr/.obj_?IObjectQ^n_/;n>1:>NonCommutativeMultiply@@conA[obj,n])//.ReleaseHold[$UserRuleList])\
		/.NonCommutativeMultiply->Times;
out
];

(*Creates a Unique index each time the rule is evaluated. After two calls to createInd for a 
particular variable the value is incremented*)
ClearAll[createInd,createIndPrevious];
createInd[ind_]:=Module[{prev,newind},
(*createIndPrevious constains info about the previous call, since each index must be called twice*)
	prev=If[Head[createIndPrevious[ind]]===createIndPrevious,Null,createIndPrevious[ind]];
	If[prev===Null,newind=dumUnique[1,++$maxDummyVal];createIndPrevious[ind]=newind,
		newind=prev;createIndPrevious[ind]=Null];
	addTFUpValues[newind,IndexQ];
	addTFUpValues[newind,DIndexQ];
	newind[[1]]];



(* ::Subsection::Closed:: *)
(*Unit vector simplification*)


$UnitVecRules={(obj_?UnitVecQ[ind_])^2->1,obj_?UnitVecQ[ind_] PD[obj_?UnitVecQ[ind_],Xd[ind2_],Vd[]]->0,\
	obj_?UnitVecQ[ind_] PD[obj_?UnitVecQ[ind_],Xd[],Vd[ind2_]]->0};
$UnitVecRuleOrder=1;

setUnitVec[obj_?IObjectQ,order_Integer:2]:=Module[{},
obj/:UnitVecQ[obj]=True;
If[Length[$FullIndexList]=!=0,
	Print[obj[$FullIndexList[[1]]]," set as unit vector"],
	Print[ToString[obj]<>"[i]  set as unit vector"]];
If[order>$UnitVecRuleOrder,
	($UnitVecRules=Join[$UnitVecRules,createUnitVecRule[#]])&/@Range[$UnitVecRuleOrder+1,order];
	If[order!=2,Print["Unit vector rules added up to order ",order]]];
];

(*Set up and apply unit vector rules*)
(*Just add first 3 partial derivaives for now*)
(*UPDATE this to generate more rules!*)
unitVecReduce[expr_]:=Module[{xsimp},
xsimp=expr//.$UnitVecRules;
If[$printToCanonicalProgress,Print["unitVecReduce: rules applied to unit vectors"]];
xsimp
];

(*Creates unit vector rules (above order 1) in the initialization of setUnitVecs*)
ClearAll[createUnitVecRule]
createUnitVecRule[n_/;n>=2]:=Module[{dinds,allder,ruleas0},
	dinds=ToExpression["VEST`Private`ind"<>ToString[#]]&/@Range[n];
	allder={Take[dinds,#],Take[dinds,-(n-#)]}&/@Range[0,n];
	ruleas0=PDExpand[PD[obj[ind]^2,Xd@@#[[1]],Vd@@#[[2]]]/2]&/@allder;
	(#[[-1]]/.Join[{obj->obj_?UnitVecQ,ind->ind_},Thread[dinds->(Pattern[#,_]&/@dinds)]])\
		:>Evaluate[-Most[#]]&/@ruleas0
]


(* ::Subsection::Closed:: *)
(*Subscript[\[CurlyEpsilon], ijk] simplification*)


ClearAll[levi2rules];
(* Simplify expressions involving Subscript[\[CurlyEpsilon], ijk] *)

(*More general version, use full expansion of an arbitrary leviC product*)
levi2rules:={leviC[indi_?IndexQ,indj_?IndexQ,indk_?IndexQ]leviC[indl_?IndexQ,indm_?IndexQ,indn_?IndexQ]:>\
	delta[indi,indl](delta[indj,indm]delta[indk,indn]-delta[indj,indn]delta[indk,indm])-\
	delta[indi,indm](delta[indj,indl]delta[indk,indn]-delta[indj,indn]delta[indk,indl])+\
	delta[indi,indn](delta[indj,indl]delta[indk,indm]-delta[indj,indm]delta[indk,indl]),
	leviC[indi_?IndexQ,indj_?IndexQ,indk_?IndexQ]^2->6}
(* Use non-commutative multiply to account for leviC^2 so that expansion is always in the order they appear*)
(* This is important in simplifyVertorForm to ensure interesting expansions of leviC products*)
leviCSimplify[expr_]:=expr//.levi2rules;





(* ::Subsection::Closed:: *)
(*Subscript[\[Delta], ij] simplification*)


ClearAll[expandForDelta,delSimprules];

delSimprules={delta[ind1_?DIndexQ,ind2_]obj_[inds1___,ind1_,inds2___]:>\
		obj[inds1,ind2,inds2],
	delta[ind1_?DIndexQ,ind2_]PD[obj_[inds1___,ind1_,inds2___],Xd[indsx___],Vd[indsv___]]:>\
		PD[obj[inds1,ind2,inds2],Xd[indsx],Vd[indsv]],
	delta[ind1_?DIndexQ,ind2_]PD[obj_[inds1___],Xd[indsx1___,ind1_,indsx2___],Vd[indsv___]]:> \
		PD[obj[inds1],Xd[indsx1,ind2,indsx2],Vd[indsv]],
	delta[ind1_?DIndexQ,ind2_]PD[obj_[inds1___],Xd[indsx___],Vd[indsv1___,ind1_,indsv2___]]:> \
		PD[obj[inds1],Xd[indsx],Vd[indsv1,ind2,indsv2]],
	delta[ind1_,ind2_]^2->3};

(* Most obvious way turns out to be the fastest anyway*)
deltaSimplify[expr_]:=Module[{xsimp},
xsimp=expandForDelta[expr]//.delSimprules;
If[$printToCanonicalProgress,Print["deltaSimplify: completed \!\(\*SubscriptBox[\(\[CurlyEpsilon]\), \(ijk\)]\) and \!\(\*SubscriptBox[\(\[Delta]\), \(ij\)]\) expansions "]];
xsimp
]

(*Expand can be very very slow in expanding after levi2rules for some reason, maybe the 
very large number of variables? Expression is always expanded after FindDummies so know
that expansion is of the form expr(delta(delta delta-delta delta)...)+... So can use 
Distribute instead and this is much much faster. Check that there is nothing squared first.*)
expandForDelta::FPiter="more than 3 FixedPoint iterations required, be worried";

expandForDelta[expr_]:=Module[{fpit,xexp},
If[Head[expr]=!=Plus || Length[expr]<10,xexp=Expand[expr],
(*If expression contains objects^2 that contain +, Distribute won't work*)
	If[Position[Cases[expr,obj_^n_,Infinity],Plus]!={},xexp=Expand[expr];\
		Message[expandForDelta::distExpand],
	fpit=0;xexp=FixedPoint[(fpit++;Plus@@Distribute/@List@@#)&,expr];\
		If[fpit>3,Message[expandForDelta::FPiter]];
	]];
xexp
]

	
	



(* ::Subsection::Closed:: *)
(*ToCanonical function*)


(* FULL SIMPLIFICATION - concatenation of the other functions *)
Options[ToCanonical]={printProgress->False,exportlist->"$IndexList"};
$printToCanonicalProgress=False;

ToCanonical[expr_/;Head[expr]=!=List,OptionsPattern[]]:=Module[{xsimp,dumlist},
	If[TrueQ[OptionValue[printProgress]],$printToCanonicalProgress=True,\
		$printToCanonicalProgress=False];
	xsimp=printDummies[SymmetryCanonicalize[
		unitVecReduce[applyUserRules[deltaSimplify[leviCSimplify[\
		FindDummies[PDExpand[expr],OptionValue[exportlist]]]]]],\
		OptionValue[exportlist]],OptionValue[exportlist]];
	$printToCanonicalProgress=False;
	(*Clear d$ variables if max value gets large to save memory. If working with
	VERY large expressions (more than 100 dummies in a term!), may need to increase this*)
	If[$maxDummyVal>$clearInternalDummiesAt,clearHighDummyDefs[100];$maxDummyVal=101];
	xsimp
]



(* ::Subsection::Closed:: *)
(*Printing dummy indices*)


(*Print dummy variables in a more easily readable format using inital DefIndices data*)
$DummyColorList={Blue,Red,Black,Purple,Brown,Gray,Green,Orange,Magenta,Pink};
$PrintNiceDummies=True;


printDummies[expr_,exportlist_:"$IndexList"]:=
	Module[{indlist,dummymax,printindlist,numcolors,colorit},
indlist=ToExpression[exportlist];
(* Find largest dummy *)
dummymax=Max@@Length/@indlist[[1]];
(*available indices*)
printindlist=Complement[$FullIndexList,indlist[[2]]];
printindlist=If[Length[printindlist]>dummymax,Take[printindlist,dummymax]\
	,Join[printindlist,Drop[dumUnique[dummymax],Length[printindlist]]]];
(*For each of the set of dummies change the print display using var$Print*)
(*setDummyPrint sets print display for a set number of times (2) then deletes it*)
colorit=1;
numcolors=Length[$DummyColorList];
setDummyPrint[vsym_,letsym_]:=Module[{vstr,vstrP,let},
	vsym/:SymbolName[vsym]:=If[$PrintNiceDummies,ToString[letsym],ToString[vsym]];
	vsym/:indColor[vsym]=$DummyColorList[[Mod[colorit,numcolors,1]]];
	colorit=colorit+1;
	];
MapThread[setDummyPrint,{dumUnique[dummymax],printindlist}];
expr
];





(* ::Subsection::Closed:: *)
(*Print as a vector expression*)


(*Function to print indexed expression as a vector expression*)
(*Only works with x derivatives up to first order at the moment*)
vectorForm[expr_]:=Module[{exprwrules,freeind,replev},
	(*Sub in rules, gives expression with free index still present*)
	exprwrules=expandForFindDummies[expr]//.vecsubrules;
	freeind=freeIndexFind[exprwrules];
(*Remove the free index. Only replace at level 2 since otherwise could end up
replacing bits in PD*)
	If[Head[exprwrules]===Times,replev={0,1},replev={0,2}];
	DisplayForm@(Replace[exprwrules,obj_@@freeind:>objPrintAs[obj],replev]\
		/.ind_?IndexQ:>ToExpression@SymbolName[ind])
]



(*Works out the permutations of a rule involving leviC tensors and prints out 6 rules that cover all possibilities*)
ClearAll[leviRuleCreate];
leviRuleCreate[rule_RuleDelayed]:=Module[{LHSlist,indperms,rhsmult},
	LHSlist=rule[[1]]/.lbit_leviC mul_:>{lbit,mul};
	indperms=Permutations[{##}]&@@(LHSlist[[1]]);
	rhsmult=(Signature/@indperms)/Signature[indperms[[1]]];
	Thread[Thread[(leviC@@@indperms) LHSlist[[2]]]:>Evaluate[rule[[2]] rhsmult]]
]

(* Set default for objPrintAs *)
objPrintAs[obj_RowBox]:=obj;

(* Rules used for substitutions*)
ClearAll[vecsubrules1,vecsubrules2,vecsubrules3,vecsubrules];
vecsubrules1={PD[o1_[],Xd[i1_],Vd[]]:>prinGobj[o1][i1],o2_[i1_?IndexQ]PD[o1_[i1_],Xd[i2_],Vd[]]:>ncDot[prinGobj[o1],\
	prinobj[o2]][i2],o2_[i2_?IndexQ]PD[o1_[i1_],Xd[i2_],Vd[]]:>ncDot[prinobj[o2],prinGobj[o1]][i1],\
	PD[o1_[i1_],Xd[i1_],Vd[]]:>prinDiv[o1],PD[o1_[i1_],Xd[i2_],Vd[]]^2:>prinGdd[o1],\
	PD[o1_[i1_],Xd[i2_],Vd[]]PD[o1_[i2_],Xd[i1_],Vd[]]:>prinGddT[o1]};
vecsubrules2=Join[leviRuleCreate[leviC[i1_,i2_,i3_]o1_[i2_?IndexQ]o2_[i3_?IndexQ]:>prinCP[o1,o2][i1]],\
	leviRuleCreate[leviC[i1_,i2_,i3_]PD[o1_[i3_?IndexQ],Xd[i2_?IndexQ],Vd[]]:>prinCurl[o1][i1]]];
vecsubrules3={o1_[i1_?IndexQ]^2:>prinobj[o1]^2,o1_[i1_?IndexQ]o2_[i1_?IndexQ]:>ncDot[prinobj[o1],prinobj[o2]]};
vecsubrules=Join[vecsubrules1,vecsubrules2,vecsubrules3];

(*Vector printing functions, uses objPrintAs which is essentially just printas as input for DefObject*)
ClearAll[prinobj,prinGobj,ncDot,prinDiv,prinGdd,prinGddT,prinCP,prinCurl];
prinobj[obj_]:=RowBox[{objPrintAs[obj]}]
prinGobj[obj_]:=RowBox[{Style["\[Del]",Bold],objPrintAs[obj]}]
ncDot[sb1_,sb2_]:=RowBox[{"(",sb1,"\[CenterDot]",sb2,")"}]
prinDiv[obj_]:=RowBox[{"(",Style["\[Del]",Bold],"\[CenterDot]",objPrintAs[obj],")"}]
prinGdd[obj_]:=RowBox[{"(",Style["\[Del]",Bold],objPrintAs[obj],":",Style["\[Del]",Bold],objPrintAs[obj],")"}]
prinGddT[obj_]:=RowBox[{"(",Style["\[Del]",Bold],objPrintAs[obj],":",Style["\[Del]",Bold],objPrintAs[obj],"\[Transpose]",")"}]
prinCP[o1_,o2_]:=RowBox[{"(",objPrintAs[o1],"\[Times]",objPrintAs[o2],")"}]
prinCurl[obj_]:=RowBox[{Style["\[Del]",Bold],"\[Times]",objPrintAs[obj]}]


(* ::Subsection::Closed:: *)
(*Inputing vector expressions*)


(*Functions to input in vector notation*)
(*All work by simply outputing an indexed expression with unique indices that can be 
used as a standard input expression*)

(*Specifies default free index of vector expressions, won't matter for scalar expressions*)
$vecFreeI=Global`n;


cross[o1_?(rankn[#,{1}]&),o2_?(rankn[#,{1}]&),ind_:False]:=Module[{ni,o1sub,o2sub},
(*New indices to use*)
ni=dumUnique[2,($maxDummyVal+=2)-1];ni=Join[{$vecFreeI},ni];
addTFUpValues[ni,IndexQ];
(*Both individual objects and expressions work with setIndEqual since FindDummies adds in a free index if necessary*)
setIndEqual[o1sub,o1];
setIndEqual[o2sub,o2];
If[ind===False,leviC@@ni o1sub[ni[[2]]]o2sub[ni[[3]]],leviC@@ni o1sub[ni[[2]]]o2sub[ni[[3]]]/.ni[[1]]->ind]
]

dot::ambig="Non derivative rank 2 tensor `1` supplied, result may be transposed";
dot::scalar="Index supplied to scalar object";
dot::twotens="Cannot dot two rank 2 tensors";
dot[o1_?(rankn[#,{1,2}]&),o2_?(rankn[#,{1,2}]&),ind_:False]:=Module[
	{ni,tensdot,o1sub,o2sub,r2inds,outobj},
(*don't necessarily need these indices, but get them anyway*)
ni=dumUnique[1,($maxDummyVal+=1)];ni=Join[{$vecFreeI},ni];
addTFUpValues[ni,IndexQ];
(*Get "dotting index" from each object*)
tensdot={False,False};(*If False, corresponding object is rank 1*)
If[rankn[o1,{1}],setIndEqual[o1sub,o1],
(*Find indices of rank 2 object, this is the correct way around for dotting gradients*)
	If[Head[o1]=!=PD,Message[dot::ambig,o1]];\
	tensdot[[1]]=True;setIndEqual[o1sub,o1]];
If[rankn[o2,{1}],setIndEqual[o2sub,o2],
(*Find indices of rank 2 object, this is the correct way around for dotting gradients*)
	If[Head[o2]=!=PD,Message[dot::ambig,o2]];
	tensdot[[2]]=True;setIndEqual[o2sub,o2]];
(*Dot together the two objects*)
Which[tensdot==={False,False},
	outobj=o1sub[ni[[2]]]o2sub[ni[[2]]];If[ind=!=False,Message[dot::scalar]],
tensdot==={True,False},
	outobj=o1sub[ni[[2]],ni[[1]]]o2sub[ni[[2]]];If[ind=!=False,outobj=outobj/.ni[[1]]->ind],
tensdot==={False,True},
	outobj=o1sub[ni[[2]]]o2sub[ni[[1]],ni[[2]]];If[ind=!=False,outobj=outobj/.ni[[1]]->ind],
tensdot==={True,True},
	Message[dot::twotens];Abort[]];
outobj
]

div[o1_?(rankn[#,{1}]&)]:=Module[{ni,o1sub},
ni=dumUnique[1,($maxDummyVal+=1)];
addTFUpValues[ni,IndexQ];
setIndEqual[o1sub,o1];
PD[o1sub[ni[[1]]],Xd[ni[[1]]],Vd[]]
]

curl[o1_?(rankn[#,{1}]&),ind_:False]:=Module[{ni,o1sub},
ni=dumUnique[2,($maxDummyVal+=2)-1];ni=Join[{$vecFreeI},ni];
addTFUpValues[ni,IndexQ];
setIndEqual[o1sub,o1];
If[ind===False,leviC@@ni PD[o1sub[ni[[3]]],Xd[ni[[2]]],Vd[]],
	leviC@@ni PD[o1sub[ni[[3]]],Xd[ni[[2]]],Vd[]]/.ni[[1]]->ind]
]

grad::indnum="Number of supplied indices is not equal to rank+1";
grad[o1_?(rankn[#,{0,1}]&),ind_:False]:=Module[{ni,o1sub,indi},
ni=dumUnique[1,($maxDummyVal+=1)];ni=Join[{$vecFreeI},ni];
addTFUpValues[ni,IndexQ];
(*Make sure ind is a list and check its length matches the supplied object*)
If[Head[ind]===Symbol&&ind=!=False,indi={ind},indi=ind];
If[indi=!=False&&((rankn[o1,{0}]&&Length[indi]!=1)||(rankn[o1,{1}]&&Length[indi]!=2)),
	Message[grad::indnum];indi=False];
(*Form the gradient*)
setIndEqual[o1sub,o1];
If[rankn[o1,{0}],
	If[indi===False,PD[o1sub[],Xd[ni[[1]]],Vd[]],PD[o1sub[],Xd@@indi,Vd[]]],
(*$vecFreeI is put onto the grad index for derivative tensor*)
	If[indi===False,PD[o1sub[ni[[2]]],Xd[ni[[1]]],Vd[]],PD[o1sub[indi[[1]]],Xd[indi[[2]]],Vd[]]]
	]
]


gradv::indnum="Number of supplied indices is not equal to rank+1";
gradv[o1_?(rankn[#,{0,1}]&),ind_:False]:=Module[{ni,o1sub,indi},
ni=dumUnique[1,($maxDummyVal+=1)];ni=Join[{$vecFreeI},ni];
addTFUpValues[ni,IndexQ];
(*Make sure ind is a list and check its length matches the supplied object*)
If[Head[ind]===Symbol&&ind=!=False,indi={ind},indi=ind];
If[indi=!=False&&((rankn[o1,{0}]&&Length[indi]!=1)||(rankn[o1,{1}]&&Length[indi]!=2)),
	Message[grad::indnum];indi=False];
(*Form the gradient*)
setIndEqual[o1sub,o1];
If[rankn[o1,{0}],
	If[indi===False,PD[o1sub[],Xd[],Vd[ni[[1]]]],PD[o1sub[],Xd[],Vd@@indi]],
(*$vecFreeI is put onto the grad index for derivative tensor*)
	If[indi===False,PD[o1sub[ni[[2]]],Xd[],Vd[ni[[1]]]],PD[o1sub[indi[[1]]],Xd[],Vd[indi[[2]]]]]
	]
]


(*Tests if object or expression has any of given ranks*)
(* Also set $findDumsAllowVecObjects to True, as this is necessary of nFIndexQ to work correctly*)
ClearAll[rankn];
rankn[obj_,n_List]:=obj===0||Or@@(TrueQ@nFIndexQ[obj,#]&/@n)



(*NOTATION FOR CROSS, DOT and GRAD*)
AddInputAlias[ParsedBoxWrapper[OverscriptBox["\[Times]","\[RightVector]"]],"v*",EvaluationNotebook[]];
AddInputAlias[ParsedBoxWrapper[OverscriptBox[StyleBox["\[Bullet]",FontColor->Black],"\[RightVector]"]],"v.",EvaluationNotebook[]];
AddInputAlias[ParsedBoxWrapper[OverscriptBox["\[Del]","\[RightVector]"]],"vg",EvaluationNotebook[]];
AddInputAlias[ParsedBoxWrapper[SubscriptBox[OverscriptBox["\[PartialD]","\[RightVector]"],"v"]],"vgv",EvaluationNotebook[]];

Notation[ParsedBoxWrapper[RowBox[{"o1_",OverscriptBox["\[Times]","\[RightVector]"],"o2_"}]] \[DoubleLongRightArrow] \
	ParsedBoxWrapper[RowBox[{"cross", "[", RowBox[{"o1_", ",", "o2_"}], "]"}]]];
Notation[ParsedBoxWrapper[RowBox[{"o1_",OverscriptBox[StyleBox["\[Bullet]",FontColor->Black],"\[RightVector]"],"o2_"}]]\[DoubleLongRightArrow]\
	ParsedBoxWrapper[RowBox[{"dot", "[", RowBox[{"o1_", ",", "o2_"}], "]"}]]];
Notation[ParsedBoxWrapper[RowBox[{OverscriptBox["\[Del]","\[RightVector]"],"o1_"}]]\[DoubleLongRightArrow] \
	ParsedBoxWrapper[RowBox[{"grad", "[", "o1_", "]"}]]];
Notation[ParsedBoxWrapper[RowBox[{SubscriptBox[OverscriptBox["\[PartialD]","\[RightVector]"],"v"],"o1_"}]]\[DoubleLongRightArrow] \
	ParsedBoxWrapper[RowBox[{"gradv", "[", "o1_", "]"}]]];



(* ::Subsection:: *)
(*Equivalent forms of an indexed expression*)


(*Apply Subscript[\[CurlyEpsilon], ijk] Subscript[\[CurlyEpsilon], ljk]=2Subscript[\[Delta], il ]to different indices of a term in an expression to obatin equivalent
forms. 
This is necessary for more complex vector identities, for instance
(a.b\[Cross]c)d+(b.c\[Cross]d)a+(c.d\[Cross]a)b+(d.a\[Cross]b)c=0*)

(*Iterates with simplifyVectorForm until expression non longer changes and returns the
result*)

Options[FullSimplifyVectorForm]={printSubs->False,maxPatternsToSub->100000,
		randomFinalForm->0,unitVec->False,numLeviPairs->1};

FullSimplifyVectorForm[expr_,OptionsPattern[]]:=
Module[{iterprint,itexpr,finalans,randitM},
$simplifyVecFormFinsished=False;
$sVFRandomIt=0;(*Iterator for random steps if it gets stuck*)
randitM=OptionValue[randomFinalForm];(* Just to reduce typing*)
(*To remember the forms from the previous iteration to save computation on large simplifications*)
$SVFpreviousStepRules={{},{}};

iterprint=1;
itexpr=expr;
Catch[While[True,
itexpr=simplifyVectorForm[itexpr,OptionValue[unitVec],{printSubs->OptionValue[printSubs],\
		maxPatternsToSub->OptionValue[maxPatternsToSub],randomReturn->0=!=randitM,\
		numLeviPairs->OptionValue[numLeviPairs]}];
If[itexpr===0||($simplifyVecFormFinsished&&(randitM<$sVFRandomIt||randitM==0)),\
	Throw[itexpr]];
Print[" "];
Print["Completed iteration ", iterprint++];
Print[" "];
]]
]



ClearAll[simplifyVectorForm];
(* Applies equivalentVectorForms to each term and substitutes in each rule, returning
the shortest expression. Not guaranteed to give shortest expression, since this could require
substituting two rules at once*)
(* maxPatterns is maximum number of patterns to try subbing in, if larger just subs in each 
seperately *)
Options[simplifyVectorForm]={printSubs->False,maxPatternsToSub->100000,\
		randomReturn->False,numLeviPairs->1};
(*THIS WILL NOT WORK PROPERLY IF SCOPE OF WHOLE PACKAGE IS CHANGED TO PUBLIC!!*)
SVFIndexListName="VEST`Private`SVFIndexList";

simplifyVectorForm[expr_,unitvec_:False,OptionsPattern[]]:=
Module[{initexpr,xplist,lengthfun,fullrulelist,rulecombs,subexprlens,shortestexpr},
initexpr=ToCanonical[expr,exportlist->SVFIndexListName];
xplist=plusToList[initexpr];
(*Calculate all equivalent forms for each term*)
$progressIterSimp=0;$totaltermsSimp=Length[xplist];
Print["Calculating equivalent forms of each term: progress..."];
Monitor[fullrulelist=equivalentVectorForms[xplist[[#]],unitvec,OptionValue[numLeviPairs]]\
		&/@Range[Length[xplist]];,ProgressIndicator[Dynamic[$progressIterSimp/$totaltermsSimp]]];
$SVFpreviousStepRules={xplist,fullrulelist};

rulecombs=createRuleSubList[fullrulelist,OptionValue[maxPatternsToSub]];
Print["Trying ", Length[rulecombs]," substitution combinations"];
If[OptionValue[printSubs],Print[fullrulelist];Print["Rule combinations to try: ", rulecombs]];
(*Returns the number of terms in an expression (with no brackets)*)
lengthfun=If[Head[#]===Plus,Length[#],If[#===0,0,1]]&;
(*Save storing everything in memory by storing the length of each term*)
subexprlens=MapIndexed[{First[#2],lengthfun@(initexpr/.Flatten[Extract[fullrulelist,#1]])}&,rulecombs];

(*Choose the shortest (not efficient, but really doesn't matter) *)
shortestexpr=Sort[subexprlens,#1[[2]]<#2[[2]]&];
If[shortestexpr==={},$simplifyVecFormFinsished=True;
	Print["Substitutions yielded no expression shorter than input"];initexpr,
(*If other forms of some objects have been found *)
shortestexpr=First[shortestexpr];
Which[shortestexpr[[2]]>Length[xplist],$simplifyVecFormFinsished=True;\
		Print["Substitutions yielded no expression shorter than input"];Print[];initexpr,\
	shortestexpr[[2]]==Length[xplist],
		$simplifyVecFormFinsished=True;$sVFRandomIt++;\
			Print["Substitutions yielded no expression shorter than input"];Print[];\
			If[OptionValue[randomReturn],Print["Trying random choice of same length"];\
				initexpr/.Flatten[Extract[fullrulelist,rulecombs\
					[[RandomChoice[Select[subexprlens,#[[2]]==Length[xplist]&]][[1]]]]]],\
			initexpr],\
	shortestexpr[[2]]<Length[xplist],
		If[shortestexpr[[2]]===0,$simplifyVecFormFinsished=True;0,\
			Print["Expression reduced from ", Length[xplist], " to ", shortestexpr[[2]]," terms"];\
			initexpr/.Flatten[Extract[fullrulelist,rulecombs[[shortestexpr[[1]]]]]]]]
]
]


(*Function takes a single term in an expression and returns the equivalent 
forms of it *)
ClearAll[equivalentVectorForms];
equivalentVectorForms::singleterm="single indexed term required!";
equivalentVectorForms::unitvec="Specified vector is not defined as unit vector!";
(*Including a unitvec also checks what you get multiplying each term by vec^2*)
equivalentVectorForms[xsing_,unitvec_,numlevis_]:=
	Module[{posInPrevious,orig,origNC,levidums,unitdums,unitmult,dummylist,dumlist2add,\
		levitosub,repdumlist,exprlevipairs,leviorders,equivFormlist},
(*Add support for multiple terms in a separate function*)
Catch[
If[Position[xsing,Plus,Heads->True]!={},
	Message[equivalentVectorForms::singleterm];Abort[]];
If[Cases[xsing,leviC[___],Infinity]==={}&&numlevis==1,Throw[{}]];
(*If this term has already been calculated, use this*)
posInPrevious=Position[$SVFpreviousStepRules[[1]],xsing];
If[Length[posInPrevious]!=0,
	Throw[Extract[$SVFpreviousStepRules[[2]],posInPrevious[[1]]]]];

orig=ToCanonical[xsing,exportlist->SVFIndexListName];
(*Extra indices necessary for leviC and unit vectors*)
levidums=Unique/@conA["lc$", 3numlevis];
If[unitvec=!=False,unitdums=Unique/@conA["un$",numlevis];
		unitmult=Times@@(unitvec[#]unitvec[#]&/@unitdums);
		If[UnitVecQ[unitvec]=!=True,Message[equivalentVectorForms::unitvec];Abort[]],
	unitmult=1;unitdums={}];
addTFUpValues[Join[levidums,unitdums],IndexQ];
origNC=(orig unitmult)/.Times->NonCommutativeMultiply/.obj_^2:>obj**obj;
(*List of index positions and names*)
dummylist={Position[origNC,_?IndexQ,-1],Cases[origNC,_?IndexQ,-1]};
dumlist2add=Join[Flatten[SVFIndexList[[1]]],levidums,unitdums];
(*List of leviC tensors to sub in*)
levitosub=leviC@@@Partition[levidums,3];
(* Possible ways to sub in \[Delta][i,j]*)
repdumlist=Subsets[Range[Length[dummylist[[2]]]],{numlevis}];
(*lists of; expr with leviC, all the leviCs to sub in*)
exprlevipairs=createLeviList[origNC,dummylist,levitosub,#]&/@repdumlist;
(*Possible orders in which to substitute leviC*)
leviorders=generateLeviOrder[Length[exprlevipairs[[1,2]]]];
(*Find all equivalent forms*)
equivFormlist=Flatten[Outer[orderedLeviExpand[#1[[1]],#1[[2]],#2,dumlist2add]&,\
	exprlevipairs,leviorders,1]];
equivFormlist=Rest@DeleteDuplicates[Join[{orig},equivFormlist]];
$progressIterSimp++;
CleanIndexList[SVFIndexListName];
(*Return forms as a rule for easy replacement*)
Throw[{Rule[orig,#]}&/@ equivFormlist];
]]



(* ::Subsubsection:: *)
(*Subfunctions*)


(* Expand leviC in levilist in the order given by orderlist*)
ClearAll[orderedLeviExpand];

orderedLeviExpand[expr_,levilist_,orderlist_,dums2add_]:=Module[{uncanexpr},
	addTFUpValues[dums2add,DIndexQ];
	uncanexpr=1/2^Floor[Length[levilist]/2]Fold[deltaSimplify[(#1 Times@@levilist[[#2]])\
		//.levi2rules]&,expr,orderlist];
	SymmetryCanonicalize[applyUserRules[unitVecReduce[uncanexpr]],SVFIndexListName]
]


(* Generate all possible ordering given number of leviC*)
(* Applying sort, not entirely convinced this catches everything *)
ClearAll[generateLeviOrder];

generateLeviOrder[llen_]:=Module[{numpairs,fulllist},
	numpairs=Floor[llen/2];
	fulllist=Select[Permutations[Subsets[Range[llen],{2}],{numpairs}],\
	Length[DeleteDuplicates[Flatten[#]]]==Length[Flatten[#]]&];
(*Rest is to remove the (1,2)(3,4)*)
	fulllist=Rest@DeleteDuplicates[Sort/@fulllist];
(*Add in the extra number *)
	If[OddQ[llen],Join[#,{Complement[Range[llen],Flatten[#]]}]&/@fulllist,fulllist]
]

ClearAll[createLeviList];
(*Generate expr and leviC list with dummies substituted in repdumnum*)
createLeviList[expr_,dummylist_,levitosub_,repdumnum_]:=Module[{exprindrep,levilist},
(*Replace parts of expression with leviC dummies*)
	exprindrep=ReplacePart[expr,Thread[dummylist[[1,repdumnum]]->levitosub[[All,1]]]];
(*Create list of leviC*)
	levilist=Flatten@MapThread[{#1,ReplacePart[#1,1->#2]}&,{levitosub,dummylist[[2,repdumnum]]}];
(*If expr contains leviC, remove and add to list*)
	levilist=Join[levilist,Cases[exprindrep,leviC[___],Infinity]];
	exprindrep=DeleteCases[exprindrep,leviC[___],Infinity];
	{exprindrep/.NonCommutativeMultiply->Times,levilist}
]



ClearAll[equivalentVectorFormsingle];	
(* Takes a single index position and returns forms obatined by replacing this*)
equivalentVectorFormsingle[dnum_,fullindpos_,xsingns_,dummiestokeep_]:=
Module[{reppos,repdum,newexpr1,newexpr2},
	reppos=fullindpos[[dnum]];
	repdum=Extract[xsingns,reppos];
(*There are two choices for which leviC to multiply out, these can lead to different 
results*)
	newexpr1=1/2(ReplacePart[xsingns,reppos->lc$n]/.NonCommutativeMultiply->Times)\
		(leviC@@lcinds);
	newexpr2=1/2(ReplacePart[xsingns,reppos->lc$n]/.NonCommutativeMultiply->Times)\
		(leviC@@(lcinds/.lc$n->repdum));
	newexpr1=(leviCSimplify[newexpr1])leviC@@(lcinds/.lc$n->repdum);
	newexpr2=(leviCSimplify[newexpr2])(leviC@@lcinds);
(* Need to add Dindex upvalues for all dummies in the initial term in case they are 
eliminated and deleted by SymmetryCanonicalize*)
	addTFUpValues[dummiestokeep,DIndexQ];
	newexpr1=SymmetryCanonicalize[unitVecReduce[applyUserRules[\
		deltaSimplify[leviCSimplify[newexpr1]]]]];
	addTFUpValues[dummiestokeep,DIndexQ];
	newexpr2=SymmetryCanonicalize[unitVecReduce[applyUserRules[\
		deltaSimplify[leviCSimplify[newexpr2]]]]];
(* Add to the list *)
	{newexpr1,newexpr2}
];


(* ::Subsubsection::Closed:: *)
(*Creating rule combinations to try*)


ClearAll[createRuleSubList];

(*Given the rule list, prints out all combinations of rules that could possibly 
join together and make a shorter expression. This is done by looking for matching 
terms between the rules and choosing subsets only of the groups with terms that match*)
createRuleSubList::head="Head of term is not Times";

createRuleSubList[fullrulelist_,maxsize_Integer]:=
Module[{allterms,chooseSame,matchlist,collectGroups,combgroups},
allterms=Drop[#,-1]&/@Position[fullrulelist,Times,Infinity];
(*Function returns a list of the rules that contain the element as the one at position patpos*)
(*Do not include terms on LHS of rule (position list is length 4)*)
chooseSame[patpos_]:=Block[{patt,matchpos,fullmatches},
(*Remove numeric factor for the pattern*)
	patt=Optional[Pattern[Null,_?NumericQ],1]If[Head[#]===Times,If[NumericQ[#[[1]]],Rest[#],#],\
			Message[createRuleSubList::head]]&@Extract[fullrulelist,patpos];
	matchpos=Select[Position[fullrulelist,patt,Infinity],Length[#]==5&];
	DeleteCases[Take[#,2]&/@matchpos,i_List/;(i[[1]]===patpos[[1]]&&i=!=patpos[[{1,2}]])]];
matchlist=DeleteDuplicates[chooseSame/@Select[allterms,Length[#]==5&]];

(*Collect the matching elements into large groups*)
collectGroups[mlist_,elnum_Integer]:=Block[{joinelfun,sameel},
	joinelfun=Join[{Union[Join@@Part[#1,#2]]},Part[#1,Complement[Range[Length[#1]],#2]]]&;
	sameel=Flatten[Most/@Position[mlist,Flatten[mlist,1][[elnum]],Infinity]];
	joinelfun[mlist,sameel]];
combgroups=Catch[Do[If[i<=Length[Flatten[matchlist,1]],matchlist=collectGroups[matchlist,i],\
	Throw[matchlist]],{i,Length[Flatten[matchlist,1]]+1}]];

Print[Length[combgroups]," subgroups of lengths ", Length/@combgroups];
(*Find subsets of these groups, i.e. rules that should be plugged in at the same time*)
(* If the number would be larger than maxsize, print out each individual *)
If[2^(Max@@(Length/@combgroups))>maxsize,\
	Print["maxPatternsToSub exceded, substituting individual terms for some combinations"]];\
Select[Flatten[If[2^Length[#]<maxsize,Subsets[#],List/@#]&/@combgroups,1],Length[#]>0&]
]


(* ::Subsection:: *)
(*Equality Checking*)


(*Very complicated indexed expressions are regularly generated. Often equality of these 
can manifest itself in a very non-trivial way, essentially with large vector identities. 
To check equality of two very general expressions we expand out the tensor products
explicity. Function can be very slow if there are large numbers of dummies*)

Options[CheckTensorZero]={userFormSpec->{},simplifyFunction->Simplify,subsetSearch->False};

CheckTensorZero[expr_,OptionsPattern[]]:=
	Module[{xlist,expReprules,finalans,finallist},
xlist=plusToList[ToCanonical[expr]];
(*Without running ToCanonical, there is a possibility that dummies won't be correctly recognized*)
expReprules=generateExplicitRepRules[xlist,OptionValue[userFormSpec]];
(*Add expressions together each time TensorProduct is calculated, may reduce memory requirements?*)
(* Will want to change this if function is parallelized*)
If[OptionValue[subsetSearch]===False,
	finalans=Fold[#1+Expand[indexToSumContract[#2,expReprules]]&,0,xlist];
	finalans=OptionValue[simplifyFunction]@finalans;,
(* Option to search inside terms for combinations that could give 0 *)
	finallist=indexToSumContract[#1,expReprules]&/@xlist;
	calcSubsetSums[Normal[finallist],OptionValue[subsetSearch],xlist];
	finalans=Plus@@finallist];

If[And@@(SameQ[#,0]||SameQ[#,True]&/@Flatten[{finalans}]),True,\
	printCTZfulloutput[finalans]]
]



ClearAll[indexToSumContract];
(*Put indexed expression into explicit Sum. This is necessary because TensorProduct explicitly 
calculates out the product and with more than 7 or 8 or dummies the computer runs out of memory*)
indexToSumContract[xsing_,expReprules_]:=
Module[{xtlist,partexpr,dumlist,outsum,sumlen,nothingiter},
xtlist=timesToList[xsing]/.obj_^2:>Sequence[obj,obj];
(*Part evaluates, so just briefly turn off message, using Hold has its own problems*)
Off[General::pspec];
partexpr=Times@@(Part[#/._?IndexQ->d$m,Sequence@@Cases[#,_?IndexQ,-1]]&/@xtlist);
partexpr=partexpr/.ind_/;(IndexQ[ind]===True && DIndexQ[ind]=!=True)->All;
dumlist=DeleteDuplicates[Cases[partexpr,_?DIndexQ,-1]];
(*Deal with no dummy index case*)
sumlen=3;If[dumlist==={},sumlen=1;dumlist={nothingiter}];
(*Evaluate sum explicitly*)
outsum=Sum[Evaluate[partexpr/.expReprules],Evaluate[Sequence@@({#,sumlen}&/@dumlist)]];
On[General::pspec];
outsum
]


ClearAll[generateExplicitRepRules];

$covecsTensZero={{$x1,$x2,$x3},{$v1,$v2,$v3},{$x1,$x2,$x3,$v1,$v2,$v3}};
(* Create a list of replacement rules based on functional dependence of objects and user
specification*)
generateExplicitRepRules[xlist_,userFormSpec_]:=
	Module[{covecs,conBnotN,allobjs,indivobjs,objcodeps,explicitvecs,pdrulelist,fullrulelist},
covecs=$covecsTensZero;
(*Find objects from the $Assumptions list *)
conBnotN=(ConstQ[#]&&!NumericQ[#])&;
allobjs=DeleteDuplicates@Cases[xlist/._?IndexQ->d$m,obj_?IObjectQ[___]|PD[___]|obj_?conBnotN,Infinity];
indivobjs=DeleteCases[Cases[allobjs,obj_[___]|obj_?conBnotN],PD[__]|leviC[__]|delta[__]|v[__]|x[__]];
(*Dependence on co-ordinates is found from the number of partial derivatives in $Assumptions*)
objcodeps={Length[Position[$Assumptions,PD[#,Xd[__],Vd[]]]],\
	Length[Position[$Assumptions,PD[#,Xd[],Vd[__]]]]}&/@indivobjs;

(*Form vectors for base objects*)
explicitvecs=MapThread[generateExplicitVecs[#1,#2,covecs]&,{indivobjs,objcodeps}];
(*Differentiate vectors and create list of substitutions*)
pdrulelist=MapThread[generateDerivativeTensors[#1,#2,#3,covecs]&,\
	{indivobjs,objcodeps,explicitvecs}];
fullrulelist=Flatten[{MapThread[Rule[#1,#2]&,{indivobjs,explicitvecs}],pdrulelist}];
fullrulelist=Join[fullrulelist,{leviC[d$m,d$m,d$m]->Normal@LeviCivitaTensor[3],\
	delta[d$m,d$m]->IdentityMatrix[3],v[d$m]->covecs[[2]],x[d$m]->covecs[[1]]}];
fullrulelist=If[userFormSpec=!={},userFormSpecification[fullrulelist,userFormSpec],fullrulelist]
]


(* ::Subsubsection:: *)
(*Searching through all subsets of list for parts that are zero*)


ClearAll[calcSubsetSums];
calcSubsetSums::symbolic="List must be numeric!";
calcSubsetSums::meanval="Mean value of list is less tan \!\(\*SuperscriptBox[\(10\), \(-6\)]\), best choose larger values";

calcSubsetSums[numlist_,subsetorder_,xlist_]:=
Module[{numsubsets,subsorder,sequencelen,findSSlimited,subiter,cutoffval,poss0sums},
(*Put subsetorder into a more easily used form*)
Print["Values of terms in sum:"]
Print[numlist];
If[Head[subsetorder]=!=List, subsorder={1,subsetorder},\
	subsorder=If[Length[subsetorder]==1,{subsetorder[[1]],subsetorder[[1]]},subsetorder]];
If[Not[And@@NumericQ/@Flatten[numlist]],
	Message[calcSubsetSums::symbolic];Abort[]];
(*Calculate total number of subsets*)
numsubsets=Sum[Length[numlist]!/subit!/(Length[numlist]-subit)!,{subit,Range@@subsorder}];
Print["Trying ", numsubsets," substitution combinations"];
cutoffval=Mean[Abs[Flatten[numlist]]]/10^7;
If[cutoffval*10^7 <10^(-6),Message[calcSubsetSums::meanval]];

sequencelen=1000;(*Break up search into bits so memory is not an issue. Empirical tests showed ~1000
combinations at once to be the fastest*)
(*Treat scalar and vector cases seperately*)
poss0sums={};
If[Head[numlist[[1]]]=!=List,
	findSSlimited[startval_]:=Block[{stopval},
		stopval=If[startval+sequencelen<=numsubsets,startval+sequencelen,numsubsets];
		poss0sums=Join[poss0sums,\
			Position[Abs[Total/@Subsets[numlist,subsorder,{startval,stopval}]],i_/;i<=cutoffval]\
				+startval-1];
		stopval],
	findSSlimited[startval_]:=Block[{stopval},
		stopval=If[startval+sequencelen<=numsubsets,startval+sequencelen,numsubsets];
		poss0sums=Join[poss0sums,\
			Position[Total/@Abs[Total/@Subsets[numlist,subsorder,{startval,stopval}]],i_/;i<=cutoffval]\
				+startval-1];
		stopval]\
];

(* Break up in While loop *)
subiter=1;
Monitor[While[subiter<numsubsets,\
	subiter=findSSlimited[subiter]],ProgressIndicator[Dynamic[subiter/numsubsets]]];
(* Print out possible zeros *)
poss0sums=Flatten[Subsets[Range[Length[numlist]],subsorder,{#}],1]&/@Flatten[poss0sums];
If[poss0sums==={},Print["No zero subset groups found"];,
	Print["Possible zero combinations in sum:"];
	Print[Plus@@Evaluate[xlist[[#]]]," = ",Total[numlist[[#]]]]&/@poss0sums;
	Print["Combinations output to $CTZsubsetzeros"];
	ToExpression["$CTZsubsetzeros="<>ToString[Plus@@@(xlist[[#]]&/@poss0sums),InputForm]];
];	
]


(* ::Subsubsection:: *)
(*Subfunctions*)


ClearAll[generateExplicitVecs];
(*Generates explicit vectors given the functional dependence and the object name*)
(*Generates explicit vectors given the functional dependence and the object name*)
generateExplicitVecs[vecobj_,vecdep_,covecs_]:=Module[{objsize,vecname,covec,inds4it,vecexplic},
objsize=If[TrueQ@ConstQ[vecobj],0,Length[vecobj/.obj_[inds___]:>{inds}]];
vecname=If[TrueQ@ConstQ[vecobj],vecobj,vecobj/.obj_[inds___]:>obj];
Which[(#>0&/@vecdep)==={False,False},
	covec={},
(#>0&/@vecdep)==={True,False},
	covec=covecs[[1]],
(#>0&/@vecdep)==={False,True},
	covec=covecs[[2]],
(#>0&/@vecdep)==={True,True},
	covec=covecs[[3]]];
If[objsize==0,If[TrueQ@ConstQ[vecobj],vecexplic=vecname,vecexplic=vecname@@covec]];
inds4it=Unique/@conA["i",objsize];
If[objsize>=1,vecexplic=Table[(vecname@@inds4it)@@covec,\
	Evaluate[Sequence@@Inner[List,inds4it,conA[3,objsize],List]]]];
(*Not the fastest for expansion but nice and simple*)
If[objectSymmetry[vecname]=!=False,vecexplic=Normal@Symmetrize[vecexplic,objectSymmetry[vecname]]];
If[UnitVecQ[vecname]===True,\
	vecexplic[[3]]=Sqrt[1-vecexplic[[1]]^2-vecexplic[[2]]^2]];
vecexplic]


ClearAll[generateDerivativeTensors];
(*Generates partial derivative tensors for a given vector and object up to specified order*)
generateDerivativeTensors[vecobj_,vecdep_,explicitvec_,covecs_]:=
	Module[{derivlist,PDrulefunc},
(*List of all the necessary tensors to form as orders in x and v*)
derivlist=Rest@Select[Tuples[Range[0,Max@@vecdep],2],\
	(#[[1]]<=vecdep[[1]]&&#[[2]]<=vecdep[[2]]&&#[[1]]+#[[2]]<=Max@@vecdep)&];
(*Function to generate the rule from the generic PD[] tensor to explicit arrays*)
PDrulefunc=PD[vecobj,Xd@@conA[d$m,#2[[1]]],Vd@@conA[d$m,#2[[2]]]]->\
	D[#1,Sequence@@conA[{covecs[[1]]},#2[[1]]],Sequence@@conA[{covecs[[2]]},#2[[2]]]]&;
MapThread[PDrulefunc,{conA[explicitvec,Length[derivlist]],derivlist}]
]


ClearAll[userFormSpecification];
(*Given the full replacement rule list, applies user specified functions to chosen elements*)
userFormSpecification::ruleform="User defined rules must be specified as a function!";
userFormSpecification::recobj="Unrecognized object `1`";

userFormSpecification[rlist_,subrules_]:=
	Module[{indpos,subrwd$m,aruleparts,arules,notinexpr},
If[Not[And@@(Head[#[[2]]]===Function&/@subrules)],\
	Message[userFormSpecification::ruleform];Abort[]];
If[Not@MatchQ[#[[1]],obj_?IObjectQ[___]|PD[__]|obj_?ConstQ],
	Message[userFormSpecification::recobj,#[[1]]];Abort[]]&/@subrules;
indpos=Position[#[[1]],_?IndexQ,-1]&/@subrules;
subrwd$m=MapThread[
	If[#2==={}&&IObjectRank[#1[[1]]]=!=0&&!TrueQ@ConstQ[#1[[1]]],\
		Message[userFormSpecification::recobj,#1[[1]]];Abort[],
		{ReplacePart[#1[[1]],#2->d$m],#1[[2]]}]&,
	{subrules,indpos}];

aruleparts=Flatten[Position[rlist,Rule[#[[1]],_]|RuleDelayed[#[[1]],_]]]&/@subrwd$m;
arules=Extract[rlist,#]&/@aruleparts;
(*Check for user rules that aren't in the expression*)
notinexpr=Position[arules,{}];
If[Length[notinexpr]!=0,Print["userFormSpec: Warning, ",Extract[subrules[[All,1]],\
		notinexpr]," not found in expression"];
	arules=Delete[arules,notinexpr];aruleparts=Delete[aruleparts,notinexpr];\
	subrwd$m=Delete[subrwd$m,notinexpr]];

ReplacePart[rlist,Thread[aruleparts->\
	MapThread[ReplacePart[#1,2->#2[[2]]@#1[[2]]]&,{arules,subrwd$m}]]]
]


(*Just removes the arguments of the functions so everything is more readable*)
printCTZfulloutput[expr_]:=Module[{printnicerules},
	printnicerules=Join[{obj_?IObjectQ[ind__][___]:>Subscript[obj,FromDigits[{ind}]],\
		Derivative[ders__][obj_?IObjectQ[ind__]][___]:>Derivative[ders][Subscript[obj,FromDigits[{ind}]]],
		obj_?IObjectQ[___]:>obj,Derivative[ders__][obj_?IObjectQ][___]:>Derivative[ders][obj]},
		Thread[$covecsTensZero[[3]]->
			{Subscript[x, 1],Subscript[x, 2],Subscript[x, 3],
				Subscript[v, 1],Subscript[v, 2],Subscript[v, 3]}]
	];
	expr/.printnicerules
]


(*Probably best to change this so you specify object to change with any index you 
want*)
(*First check that object is in b[i] or \[Mu][] format, then check that has indices and 
change indices into d$m. After this program will run exactly the same, just somewhat
more user friendly input*)


End[]
EndPackage[]
