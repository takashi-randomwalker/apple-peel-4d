(* symmetry_fixed_column.wls -- 2026-08-19
   Why the 16-cell and 600-cell admit no symmetric net.

   Lemma: if a rooted spanning tree is H-invariant then the parent of an
   H-fixed cell is itself H-fixed, so Fix(H) must induce a CONNECTED
   subgraph of the cell graph containing the root.

   This script prints, for one representative of each element order in
   Stab(r), the size of Fix(h), how it is spread over the BFS layers,
   and how many components it has.

   Result: for the 16-cell and the 600-cell every non-identity element
   has Fix(h) disconnected, so MAX|H| = 1 for any spanning tree.  For
   the 120-cell the order-5 element has Fix = {1,2,2,2,2,1} per layer in
   ONE component -- the connected 10-cell column along the 5-fold axis
   that the sector construction runs on. *)

baseDir = "/Users/yoshino/Library/CloudStorage/Dropbox/260324Peeling4D/";
img[pm_, x_] := PermutationReplace[x, pm];
Do[
 raw = Get[baseDir <> "_4DData/" <> f]; cellsP = raw[[6]]; n = Length[cellsP];
 adj = Table[Select[Range[n], Function[j, j =!= i &&
      Length[Intersection[cellsP[[i]], cellsP[[j]]]] >= 1]], {i, n}];
 gg = Graph[Range[n], Flatten[Table[UndirectedEdge[i, j], {i, n},
      {j, Select[adj[[i]], # > i &]}]]];
 stabR = GroupElements[GroupStabilizer[GraphAutomorphismGroup[gg], {1}]];
 lay = <|1 -> 0|>; q = {1};
 While[q =!= {}, v = First[q]; q = Rest[q];
  Do[If[! KeyExistsQ[lay, u], lay[u] = lay[v] + 1; AppendTo[q, u]], {u, Sort@adj[[v]]}]];
 nl = Max[Values[lay]];
 Print["\n", f, "   layer sizes ", Table[Count[Values[lay], d], {d, 0, nl}]];
 Do[Module[{h, fix, prof},
    h = SelectFirst[Rest[stabR], PermutationOrder[#] === ord &];
    If[! MissingQ[h],
      fix = Select[Range[n], img[h, #] === # &];
      prof = Table[Count[lay[#] & /@ fix, d], {d, 0, nl}];
      Print["   element order ", ord, " : |Fix| = ", StringPadLeft[ToString[Length[fix]], 3],
            "   Fix per layer ", prof,
            "   components ", Length[ConnectedComponents[Subgraph[gg, fix]]],
            If[MemberQ[Rest[prof[[1 ;; -2]]], 0], "   <-- GAP in a middle layer", ""]]]],
  {ord, Union[PermutationOrder /@ Rest[stabR]]}],
 {f, {"f16.m", "f600.m", "f120.m"}}];
