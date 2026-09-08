(* check_dhr_chain_120cell.wls
   Is the overlapping chain of nine dodecahedra of Devadoss, Harvey and
   Richter (arXiv:2608.30285) a path in some BFS tree of the 120-cell?

   Cells of the 120-cell are vertices of the 600-cell, so the chain is a
   path in the 600-cell's edge graph, built here from the unit icosians in
   the same coordinates the paper uses.  Their Table 3 lists the nine
   vertices; the seventh came out of the PDF with its last two entries
   transposed, and is repaired here by the fact that points 6 and 8 have
   exactly one common neighbour.

   Answer: no, for none of the 120 roots.
*)

(* Can the DHR overlapping chains lie in a BFS tree?
   A path in a tree rooted at r has depths differing by exactly 1 along it,
   and no interior local MAXIMUM (a cell with both neighbours shallower would
   need two parents).  So the depth sequence must be V-shaped.  That is
   necessary and sufficient: any such assignment extends to a spanning BFS
   tree by giving every other cell any admissible parent. *)
phi = GoldenRatio; psi = 1/GoldenRatio;
v8  = Flatten[Table[s UnitVector[4,i], {i,4},{s,{1,-1}}],1];
v16 = Tuples[{1/2,-1/2},4];
even = Select[Permutations[Range[4]], Signature[#]==1&];
v96 = DeleteDuplicates@Flatten[Table[Permute[{s1 phi/2, s2 1/2, s3 psi/2, 0}, p],
   {p,even},{s1,{1,-1}},{s2,{1,-1}},{s3,{1,-1}}],3];
V = N[Join[v8,v16,v96]];
c = N[phi/2];
adj = Table[Select[Range[120], #!=i && Abs[V[[i]].V[[#]] - c] < 10^-8 &], {i,120}];

dists[a_, r_] := Module[{d = ConstantArray[-1, Length[a]], q = {r}, p},
  d[[r]] = 0;
  While[q =!= {}, p = First[q]; q = Rest[q];
    Do[If[d[[nb]] < 0, d[[nb]] = d[[p]]+1; AppendTo[q,nb]], {nb, a[[p]]}]]; d];

(* V-shaped: strictly down to a minimum, then strictly up *)
pathOK[dd_] := Module[{k = Length[dd], st},
  st = Table[dd[[i+1]] - dd[[i]], {i, k-1}];
  If[Union[Abs[st]] =!= {1}, Return[False]];
  (* no interior local max: no  +1 followed by -1 *)
  !MemberQ[Partition[st,2,1], {1,-1}]];

chain = {1, 33, 9, 113, 99, 114, 97, 41, 49};
Print["=== 120-cell (cells = vertices of the 600-cell) ==="];
Print["graph: ", Length[V], " nodes, degrees ", Union[Length/@adj]];
Print["layers from node 1: ", Module[{d=dists[adj,1]}, Table[Count[d,k],{k,0,Max[d]}]]];
Print["chain is a path: ", And @@ Table[MemberQ[adj[[chain[[i]]]], chain[[i+1]]], {i,8}]];
ok = Select[Range[120], pathOK[dists[adj,#][[chain]]] &];
Print["roots for which the chain fits in a BFS tree: ", Length[ok], " of 120  ", ok];
Do[Print["   root ", r, "  depths along the chain: ", dists[adj,r][[chain]]], {r, Take[ok, UpTo[3]]}];
If[ok === {},
  Print["   depth sequences from a few roots: "];
  Do[Print["   root ", r, ": ", dists[adj,r][[chain]]], {r, {1,2,33,60}}]];
