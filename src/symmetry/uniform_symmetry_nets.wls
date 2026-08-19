(* uniform_symmetry_nets.wls -- 2026-08-19
   Maximum-symmetry face-rotation nets of a convex uniform 4-polytope.
   One polytope and one cell-orbit per invocation:

     wolframscript -file uniform_symmetry_nets.wls <symbol> <orbitIndex> <nSamples>

   Cell orbits are sorted by size, so orbitIndex 1 is the smallest.
   For each root it computes MAX |H| (the largest subgroup of Stab(r)
   admitting an H-invariant BFS tree), samples that many H-invariant
   trees, and samples the same number of UNCONSTRAINED random BFS trees
   as a control.  Both are checked for self-intersection by SAT.

   NOTE: GroupElements blows the kernel out of memory at degree 2640
   even for a group of order 120.  The stabilizer is therefore expanded
   from GroupGenerators by hand, on integer arrays. *)

(* one polytope, one cell-orbit per invocation; permutations as integer arrays *)
baseDir = "/Users/yoshino/Library/CloudStorage/Dropbox/260324Peeling4D/";
Module[{src, cut},
  src = Import[baseDir <> "src/bfsnet/face_rotation_net_all4D_v2.m", "Text"];
  cut = First@First@StringPosition[src, "Print[\"Face-rotation BFS net check v2"];
  ToExpression[StringTake[src, {1, cut - 1}]]];
Get[baseDir <> "src/bfsnet/bfs_sat.m"]; Get[baseDir <> "data/wythoff_gen_all.mx"];
SeedRandom[5];
sym = $ScriptCommandLine[[2]];
orbIdx = ToExpression[$ScriptCommandLine[[3]]];
nSamp = ToExpression[$ScriptCommandLine[[4]]];

rec = SelectFirst[w4All, #["sym"] === sym &]; data = rec["data"];
versP = N[data[[1]], MachinePrecision]; facesP = data[[4]]; cellsP = data[[6]];
nCP = Length[cellsP]; versCenteredP = (# - Mean[versP]) & /@ versP;
adjListP = Module[{f2c},
  f2c = GroupBy[Flatten@MapIndexed[Function[{cf,k},{# -> First[k]} & /@ cf], cellsP], First -> Last];
  Union /@ Lookup[GroupBy[Join[#, Reverse /@ #] &@ Select[Values[f2c], Length[#] === 2 &],
      First -> Last], Range[nCP], {}]];
adjSetP = Association@Flatten[Table[Sort[{i,j}]->True,{i,nCP},{j,adjListP[[i]]}],1];
satPrepCells[];
gG = Graph[Range[nCP], Flatten[Table[UndirectedEdge[i,j],{i,nCP},{j,Select[adjListP[[i]],#>i&]}]]];
aut = GraphAutomorphismGroup[gG];
orbs = SortBy[GroupOrbits[aut, Range[nCP]], Length];
r = Min[orbs[[orbIdx]]];
(* GroupElements explodes at degree 2640 even for a group of order 120,
   so expand the stabilizer from its generators on integer arrays. *)
grow[gens_] := Module[{S = DeleteDuplicates[Prepend[gens, Range[nCP]]], new},
  While[True,
   new = DeleteDuplicates[Join[S, Flatten[Table[a[[b]], {a, S}, {b, gens}], 1]]];
   If[Length[new] === Length[S], Break[]]; S = new];
  S];
stabA = grow[PermutationList[#, nCP] & /@ GroupGenerators[GroupStabilizer[aut, {r}]]];
Clear[aut, gG];
Print[sym, "  cells=", nCP, "  root=", r, "  orbit size=", Length[orbs[[orbIdx]]],
      "  |Stab(r)|=", Length[stabA]];

lay = <|r->0|>; q={r};
While[q=!={}, v=First[q]; q=Rest[q];
 Do[If[!KeyExistsQ[lay,u], lay[u]=lay[v]+1; AppendTo[q,u]],{u,Sort@adjListP[[v]]}]];
cand = Table[If[v===r,{},Select[adjListP[[v]], lay[#]==lay[v]-1 &]],{v,nCP}];

(* subgroup closure on integer arrays *)
closeA[gens_] := grow[gens];
feasA[H_] := AllTrue[Complement[Range[nCP],{r}], Function[v,
   With[{st = Select[H, #[[v]] === v &]},
     AnyTrue[cand[[v]], Function[c, AllTrue[st, #[[c]] === c &]]]]]];
S1 = Select[stabA, feasA[closeA[{#}]] &];
frontier = DeleteDuplicates[closeA[{#}] & /@ S1]; maximal = {};
While[frontier =!= {}, next = {};
 Do[grown = False;
    Do[K2 = closeA[Join[K, {h}]];
       If[Length[K2] > Length[K] && feasA[K2], grown = True; AppendTo[next, Sort[K2]]], {h, S1}];
    If[! grown, AppendTo[maximal, K]], {K, frontier}];
 frontier = DeleteDuplicates[next]];
H = If[maximal === {}, {Range[nCP]}, First@MaximalBy[DeleteDuplicates[maximal], Length]];
orbits = DeleteDuplicates[Table[Union[Table[h[[v]], {h, H}]], {v, nCP}]];
reps = DeleteCases[First /@ orbits, r];
choices = Table[With[{st = Select[H, #[[v]] === v &]},
    Select[cand[[v]], Function[c, AllTrue[st, #[[c]] === c &]]]], {v, reps}];
Print["  MAX |H| = ", Length[H], "   H-orbits = ", Length[orbits]];

buildPar[sel_] := Module[{pp = <||>},
  MapThread[Function[{v,c}, Do[pp[h[[v]]] = h[[c]], {h, H}]], {reps, sel}]; pp];
unfoldGiven[pm_] := Module[{tfs = <|r -> idAff4|>, order},
  order = SortBy[Complement[Range[nCP],{r}], lay[#]&];
  Do[tfs[v] = unfoldTFP[pm[v], v, tfs[pm[v]]], {v, order}];
  Table[Module[{tf = tfs[i]}, Take[#,3] & /@ (versP[[satCellIdx[[i]]]] . Transpose[tf[[1]]] +
     Threaded[tf[[2]]])], {i, nCP}]];
check[u_] := Module[{bMin,bMax,cent,radius,cd,ne},
  bMin=(Min/@Transpose[#])&/@u; bMax=(Max/@Transpose[#])&/@u; cent=Mean/@u;
  radius=2 Max[MapThread[Function[{pts,c},Max[Norm[#-c]&/@pts]],{u,cent}]];
  cd=satBBoxPairs[cent,bMin,bMax,radius];
  If[Length[cd]===0, Return[0]];
  ne=Association@Table[c->satNormalsEdges[u[[c]],satFaceTri[[c]],satEdgePair[[c]]],{c,Union@@cd}];
  Length@Select[cd, Function[p, satOverlapQ[u[[p[[1]]]],ne[p[[1]]][[1]],ne[p[[1]]][[2]],
              u[[p[[2]]]],ne[p[[2]]][[1]],ne[p[[2]]][[2]]]]]];
versP = satFaceUp[r, versCenteredP];
ovS = Table[check[unfoldGiven[buildPar[RandomChoice /@ choices]]], {nSamp}];
Print["  symmetric trees : VALID ", Count[ovS,0], "/", nSamp,
      "   min ", Min[ovS], "   median ", Median[ovS]];
ovR = Table[check[unfoldGiven[Association@Table[v -> RandomChoice[cand[[v]]],
       {v, Complement[Range[nCP],{r}]}]]], {nSamp}];
Print["  random    trees : VALID ", Count[ovR,0], "/", nSamp,
      "   min ", Min[ovR], "   median ", Median[ovR]];
