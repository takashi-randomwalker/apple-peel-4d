(* symmetry_spanning_trees.wls -- 2026-08-19
   Same question as symmetry_ceiling.wls but for GENERAL rooted spanning
   trees: parents may be any neighbour, not only one layer up.  Greedy
   attachment of whole H-orbits is exact, because attachability is
   monotone in the attached set.

   Usage: wolframscript -file symmetry_spanning_trees.wls '{"f16.m"}'

   Result: identical to the BFS ceiling for all six polytopes
   (24, 8, 1, 6, 10, 1).  Dropping the shortest-path condition buys
   nothing -- the obstruction is not the BFS constraint. *)

(* Max |H| admitting an H-invariant SPANNING tree rooted at r
   (parents may be any neighbour, not only one layer up).
   Greedy on orbits is exact: attachability is monotone in A. *)
baseDir = "/Users/yoshino/Library/CloudStorage/Dropbox/260324Peeling4D/";
img[pm_, x_] := PermutationReplace[x, pm];
files = ToExpression[$ScriptCommandLine[[2]]];
Do[
 raw = Get[baseDir <> "_4DData/" <> f]; cellsP = raw[[6]]; n = Length[cellsP];
 adj = Table[Select[Range[n], Function[j, j =!= i &&
      Length[Intersection[cellsP[[i]], cellsP[[j]]]] >= 1]], {i, n}];
 g = Graph[Range[n], Flatten[Table[UndirectedEdge[i, j], {i, n},
      {j, Select[adj[[i]], # > i &]}]]];
 stabR = GroupElements[GroupStabilizer[GraphAutomorphismGroup[g], {1}]];
 r = 1;
 lay = <|r -> 0|>; q = {r};
 While[q =!= {}, v = First[q]; q = Rest[q];
  Do[If[! KeyExistsQ[lay, u], lay[u] = lay[v] + 1; AppendTo[q, u]], {u, Sort@adj[[v]]}]];
 candBFS = Table[If[v === r, {}, Select[adj[[v]], lay[#] == lay[v] - 1 &]], {v, n}];
 close[gens_] := Sort@GroupElements[PermutationGroup[gens]];

 (* BFS-tree feasibility, for comparison *)
 feasBFS[H_] := AllTrue[Complement[Range[n], {r}], Function[v,
    With[{st = Select[H, img[#, v] === v &]},
      AnyTrue[candBFS[[v]], Function[c, AllTrue[st, img[#, c] === c &]]]]]];

 (* spanning-tree feasibility: greedy attach of whole orbits *)
 feasSPAN[H_] := Module[{A, par, stab, grew, orbOf, todo},
   stab = Association@Table[v -> Select[H, img[#, v] === v &], {v, n}];
   orbOf = Association@Table[v -> Union[Table[img[h, v], {h, H}]], {v, n}];
   A = <|r -> True|>; par = <||>; todo = Complement[Range[n], {r}];
   While[True,
    grew = False;
    Do[If[! KeyExistsQ[A, v],
       Module[{c = SelectFirst[adj[[v]], Function[u,
            KeyExistsQ[A, u] && AllTrue[stab[v], img[#, u] === u &]]]},
        If[! MissingQ[c],
          Do[par[img[h, v]] = img[h, c]; A[img[h, v]] = True, {h, H}];
          grew = True]]], {v, todo}];
    If[! grew, Break[]]];
   {Length[A] === n, par}];

 climb[feas_] := Module[{S1, frontier, maximal = {}, next, grown, K2},
   S1 = Select[stabR, feas[close[{#}]][[1]] &];
   frontier = DeleteDuplicates[close[{#}] & /@ S1];
   While[frontier =!= {}, next = {};
    Do[grown = False;
       Do[K2 = close[Join[K, {h}]];
          If[Length[K2] > Length[K] && feas[K2][[1]], grown = True; AppendTo[next, K2]], {h, S1}];
       If[! grown, AppendTo[maximal, K]], {K, frontier}];
    frontier = DeleteDuplicates[next]];
   DeleteDuplicates[maximal]];

 mB = climb[{feasBFS[#], 0} &];
 mS = climb[feasSPAN];
 Print[StringPadRight[f, 8], " cells=", StringPadLeft[ToString[n], 4],
   "  |Stab(r)|=", StringPadLeft[ToString[Length[stabR]], 5],
   "   MAX|H| BFS trees = ", If[mB === {}, 1, Max[Length /@ mB]],
   "     MAX|H| spanning trees = ", If[mS === {}, 1, Max[Length /@ mS]],
   "   (maximal orders ", Sort[Length /@ mS], ")"];
 If[MemberQ[{"f16.m", "f600.m"}, f] && mS =!= {},
   Module[{Hb = First@MaximalBy[mS, Length], res},
    res = feasSPAN[Hb];
    DumpSave[baseDir <> "spansym_" <> StringDrop[StringDrop[f, 1], -2] <> ".mx",
     {Hb, res, adj, lay, n, r}];
    Print["         saved an H-invariant spanning tree, |H| = ", Length[Hb],
          ", element orders ", Sort@Tally[PermutationOrder /@ Hb]]]],
 {f, files}];
