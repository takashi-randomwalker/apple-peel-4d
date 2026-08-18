(* =====================================================================
   order_induced_trees.m
   How many BFS trees can a traversal actually reach?

   A traversal processes each layer in some order and gives v the first
   of its candidates to be dequeued, so the trees a traversal can
   produce are exactly those of the form

       p(v) = the minimum of opts(v) under some total order  <

   Such a total order exists iff the constraints it must satisfy --
   p(v) < q for every other candidate q of v -- can be met at once,
   i.e. iff the digraph with an edge p(v) -> q for each such pair is
   acyclic.  All candidates of v lie in one layer, so the digraph
   splits by layer and the test is cheap.

   We enumerate every BFS tree of a root and count how many pass.
   ===================================================================== *)

baseDir = "/Users/yoshino/Library/CloudStorage/Dropbox/260324Peeling4D/";
SetDirectory[baseDir];

layersAndOptions[adj_, nC_, r_] := Module[{dist, q, head, p},
  dist = ConstantArray[-1, nC]; dist[[r]] = 0;
  q = {r}; head = 1;
  While[head <= Length[q],
   p = q[[head]]; head++;
   Do[If[dist[[nb]] < 0, dist[[nb]] = dist[[p]] + 1; AppendTo[q, nb]],
    {nb, adj[[p]]}]];
  {dist, Table[If[i === r, {},
     Select[adj[[i]], dist[[#]] === dist[[i]] - 1 &]], {i, nC}]}];

(* is this parent assignment induced by some total order? *)
orderInducedQ[opts_, rest_, choice_] := Module[{edges},
  edges = Join @@ Table[
    With[{v = rest[[k]], p = choice[[k]]},
     {p -> #} & /@ DeleteCases[opts[[v]], p]],
    {k, Length[rest]}];
  edges = Flatten[edges];
  If[edges === {}, True, AcyclicGraphQ[Graph[edges]]]];

run[file_, label_, roots_] := Module[
  {raw, cellsF, nC, adj, out},
  raw = Get[baseDir <> "_4DData/" <> file];
  cellsF = raw[[6]]; nC = Length[cellsF];
  adj = Table[Select[Range[nC], Function[j, j =!= i &&
     Length[Intersection[cellsF[[i]], cellsF[[j]]]] >= 1]], {i, nC}];
  Print["\n", label, "  (", nC, " cells)"];
  out = Table[
    Module[{dist, opts, rest, trees, nTot, nOK, t, shared, byLayer},
     {dist, opts} = layersAndOptions[adj, nC, r];
     rest = Complement[Range[nC], {r}];
     nTot = Times @@ (Length /@ opts[[rest]]);
     byLayer = GroupBy[rest, dist[[#]] &];
     shared = Total@Table[
        Length@Select[Subsets[byLayer[k], {2}],
          Length[Intersection[opts[[#[[1]]]], opts[[#[[2]]]]]] >= 2 &],
        {k, Keys[byLayer]}];
     trees = Tuples[opts[[rest]]];
     {t, nOK} = AbsoluteTiming@Count[trees, c_ /; orderInducedQ[opts, rest, c]];
     Print["  root ", StringPadLeft[ToString[r], 3], ":  ",
       StringPadLeft[ToString[nOK], 7], " / ", nTot,
       "  reachable by a traversal  (",
       NumberForm[100. nOK/nTot, {4, 2}], "%)",
       "   sibling pairs sharing >=2 parents: ", shared,
       "   ", Round[t, 0.1], " s"];
     {r, nTot, nOK, shared}],
    {r, roots}];
  <|"label" -> label, "cells" -> nC, "out" -> out|>];

Print[StringRepeat["=", 78]];
Print["BFS trees reachable by a traversal, out of all trees the"];
Print["layer condition allows"];
Print[StringRepeat["=", 78]];

oiAll = {
  run["f8.m",  "8-cell",  Range[8]],
  run["f16.m", "16-cell", {1, 2, 3}],
  run["f24.m", "24-cell", {1, 2, 3}]};

DumpSave[baseDir <> "order_induced_trees.mx", oiAll];
Print["\nsaved to order_induced_trees.mx"];
