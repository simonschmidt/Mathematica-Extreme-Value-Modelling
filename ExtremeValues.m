(* ::Package:: *)

BeginPackage["ExtremeValues`"]


If[Not@ValueQ[EVModelFit::usage],
EVModelFit::usage=
"EVModelFit[data,Options] Provides a FittedEVModel"]

If[Not@ValueQ[FittedEVModel::usage],
FittedEVModel::usage=
"FittedEVModel[...] represents the fitted model from EVModelFit"]

(* Unprotect[...] *)
Unprotect[EVModelFit];

Begin["`Private`"]
ClearAll[EVModelFit,FittedEVModel];

(* Display as FittedEVModel[<data length>] *)
MakeBoxes[FittedEVModel[s_],_] ^:= InterpretationBox[
  RowBox[{"FittedEVModel","[","<",#,">","]"}],
    FittedEVModel[s]]&@ToBoxes[Length[FittedEVModel[s]["Data"]]];

FittedEVModel[s_]["Properties"]:={
  "Data","Paramaters","Distribution","ObservedInformation",
  "ProbabilityPlot","QuantilePlot","ReturnLevelPlot",
  "V","VarZp","ConfidenceInterval","LogLikelihood",
  "DevianceConfidenceInterval"
};


(* Constructor *)
(* TODO:
  Block maxima: Do the blocking
  Threshold:    Implement *)
EVModelFit[data_,OptionsPattern[]] :=(With[{s=Unique["EVModel"]},
  FittedEVModel[s];
  FittedEVModel[s]["Data"]=data;

  FittedEVModel[s]["paramNames"]={\[Mu]$,\[Sigma]$,\[Xi]$};
  FittedEVModel[s]["nParams"]=3;
  FittedEVModel[s]["Parameters"]=
    FindDistributionParameters[data,MaxStableDistribution[\[Mu]$,\[Sigma]$,\[Xi]$]];

  FittedEVModel[s]["Distribution"]=
    MaxStableDistribution@@(FittedEVModel[s]["paramNames"]/.FittedEVModel[s]["Parameters"]);

  (* Return model *)
  FittedEVModel[s]
  ])

FittedEVModel[s_]["ProbabilityPlot",opts:OptionsPattern[ProbabilityPlot]]:=
  ProbabilityPlot[
    FittedEVModel[s]["Distribution"],FittedEVModel[s]["Data"],
    PlotLabel->"Probability Plot",
    FrameLabel->{"Empirical","Model"},
    opts];

FittedEVModel[s_]["QuantilePlot",opts:OptionsPattern[QuantilePlot]]:=
  QuantilePlot[
    FittedEVModel[s]["Distribution"],FittedEVModel[s]["Data"],
    PlotLabel->"Quantile Plot",
    FrameLabel->{"Model","Empirical"}
    opts];

(* TODO: Verify empirical expression (+1?), custom conf interval and period *)
FittedEVModel[s_]["ReturnLevelPlot",opts:OptionsPattern[LogLinearPlot]]:=
  LogLinearPlot[{
    Quantile[FittedEVModel[s]["Distribution"], 1 - 1/p],
    Quantile[FittedEVModel[s]["Distribution"], 1 - 1/p]+1.96Sqrt[FittedEVModel[s]["VarZp",1/p]],
    Quantile[FittedEVModel[s]["Distribution"], 1 - 1/p]-1.96Sqrt[FittedEVModel[s]["VarZp",1/p]]
    }, {p, 1, 1000}, 
	PlotStyle->{Automatic,Dashed,Dashed},
    (* Odd LogLinearPlot behaviour, coordinates do not correspond to axes ticks!
       so x-axis must be Log'ed *)
    Epilog -> Point@MapIndexed[
      {Log@((Length[FittedEVModel[s]["Data"]] + 1)/(Length[FittedEVModel[s]["Data"]] - First@#2)), #} &
      ,Most@Sort[FittedEVModel[s]["Data"]]],
    PlotLabel->"Return Level PLot",
    Frame->True,
    FrameLabel->{"Return Period","Return Level"}
    opts];
  
(* TODO: Clip epilog *)
FittedEVModel[s_]["ProbabilityDensityPlot"]:=
  Histogram[FittedEVModel[s]["Data"], "Wand", "PDF",
   Epilog -> {
     First@Plot[
       PDF[FittedEVModel[s]["Distribution"], x],
       {x, Min@FittedEVModel[s]["Data"], Max[FittedEVModel[s]["Data"]]},
       PlotRange -> All],
     Point[{#, 0} & /@ FittedEVModel[s]["Data"]]},
   PlotLabel->"Density Plot",
   Frame->True,
   FrameLabel->{"z","f(z)"}];

FittedEVModel[s_]["PlotDiagnostics"]:=
  GraphicsGrid[{
    {FittedEVModel[s]["ProbabilityPlot"],FittedEVModel[s]["QuantilePlot"]},
    {FittedEVModel[s]["ReturnLevelPlot"],FittedEVModel[s]["ProbabilityDensityPlot"]}},
    ImageSize->Large];

FittedEVModel[s_]["deltaVar",f_]:=Module[{
  df=D[f @@ (FittedEVModel[s]["paramNames"]), {FittedEVModel[s]["paramNames"]}]},
  df.FittedEVModel[s]["V"].df/.FittedEVModel[s]["Parameters"]
 ];

FittedEVModel[s_]["ivarZp"]:=FittedEVModel[s]["ivarZp"]=Function[p,Evaluate@PiecewiseExpand@
  FittedEVModel[s]["deltaVar", 
    Function[
     Evaluate@FittedEVModel[s]["paramNames"], 
     Evaluate@Quantile[Head[FittedEVModel[s]["Distribution"]] @@ FittedEVModel[s]["paramNames"], 1-p]
  ]]]

(* TODO: More stable method *)
FittedEVModel[s_]["VarZp",p_]:=FittedEVModel[s]["ivarZp"][p];

plusminus[v_,d_]:={v-d,v+d};

Clear@clusters;
clusters[list_,u_:0,r_:1]:=Module[{
  runs=Split[Thread[list>u]]
  },
  (* ugh *)
  Split[Flatten[Position[Flatten[runs/.{a___,c1:{True..},f:{Repeated[False,{0,r}]},c2:{True..},b___}:>{a,Join[c1,f/.False->True,c2],b}],True]],#2==#1+1&]
];

HessianH[f_, x_List?VectorQ] := D[f, {x, 2}];

FittedEVModel[s_]["ObservedInformation"]:=FittedEVModel[s]["ObservedInformation"]=Module[{
  hess},
  hess[x_]=HessianH[
    -LogLikelihood[Head[FittedEVModel[s]["Distribution"]]@@FittedEVModel[s]["paramNames"],{x}],
    FittedEVModel[s]["paramNames"]]/.FittedEVModel[s]["Parameters"];
  Total[hess/@FittedEVModel[s]["Data"]]
]

FittedEVModel[s_]["V"]:=FittedEVModel[s]["V"]=Inverse[FittedEVModel[s]["ObservedInformation"]];

FittedEVModel[s_]["ConfidenceInterval",\[Alpha]_:0.05]:=Table[FittedEVModel[s]["Parameters"][[i,1]]->
    Sort@plusminus[
     FittedEVModel[s]["Parameters"][[i,2]],
     Quantile[NormalDistribution[],1-\[Alpha]/2]Sqrt[FittedEVModel[s]["V"][[i,i]]]
    ],{i,1,FittedEVModel[s]["nParams"]}];

FittedEVModel[s_]["DevianceConfidenceInterval"]:=Message[FittedEVModel::argm,"DevianceConfidenceInterval",0,2];
FittedEVModel[s_]["DevianceConfidenceInterval",p_Integer,\[Alpha]_:0.05]:=
  FittedEVModel[s]["DevianceConfidenceInterval",FittedEVModel[s]["paramNames"][[p]],\[Alpha]];

FittedEVModel[s_]["LogLikelihood"]:=FittedEVModel[s]["LogLikelihood"]=
  LogLikelihood[FittedEVModel[s]["Distribution"],FittedEVModel[s]["Data"]];


FittedEVModel[s_]["DevianceConfidenceInterval",p_Symbol,\[Alpha]_:0.05]:=Module[{
  var,
  maxl = FittedEVModel[s]["LogLikelihood"],
  bestp = p/.FittedEVModel[s]["Parameters"],
  start,end,
  target,
  b
 },
  (* Solving 2(maxl-l(var)) = Subscript[\[Chi], \[Alpha]] *)
  target = maxl-Quantile[ChiSquareDistribution[1],1-\[Alpha]]/2;

  (* Initial guess, yes .5 *)
  {start,end} = p/.FittedEVModel[s]["ConfidenceInterval",.5];

  (* Dist[a,b,c] -> Dist[var,b,c] etc *)
  d[var_]=LogLikelihood[Head[FittedEVModel[s]["Distribution"]]@@FittedEVModel[s]["paramNames"]/.
       (FittedEVModel[s]["Parameters"]/.(p->_):>(p->var))
     ,FittedEVModel[s]["Data"]];

  (* PrincipalAxis since methods that need to compute derivative are slow
     TODO: analytic starting point and upper bound that is
           within theoretical maximum *)
  ReplaceList[var,Flatten@{
   Last@FindMinimum[
    Min[Abs[d[var]-target],maxl-target],
    {var,(9bestp+start)/10,-Infinity,bestp},
    Method->"PrincipalAxis"],
   Last@FindMinimum[
     Min[Abs[d[var]-target],maxl-target],
     {var,(9bestp+end)/10,bestp,Infinity},
     Method->"PrincipalAxis"]}]
]

End[]
(* Protect[...] *)
Protect[EVModelFit];

EndPackage[]






