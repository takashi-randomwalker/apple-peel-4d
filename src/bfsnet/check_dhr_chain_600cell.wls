(* check_dhr_chain_600cell.wls
   The same question for the 600-cell, answered without identifying their
   chain.  Theirs is eight tetrahedra sharing a vertex, so it lies in one
   vertex figure, whose 20 cells form a dodecahedral graph; this unfolds
   every 8-cell path in that graph and asks which overlap and which fit in
   a BFS tree.

   Answer: 120 of the 1440 overlap, 600 fit in a BFS tree, and no path does
   both -- so whichever one they exhibit, it is not part of a BFS tree.
*)

(* Unfold every 8-cell chain inside one vertex figure of the 600-cell and
   look for overlap.  A chain that (a) overlaps and (b) fits in a BFS tree
   would mean the 600-cell has an INVALID BFS tree. *)
base = "/Users/yoshino/Library/CloudStorage/Dropbox/260324Peeling4D/";
raw = Get[base <> "_4DData/f600.m"];
vers0 = N@raw[[1]]; faces = raw[[4]]; cells = raw[[6]]; nC = Length[cells];
f2c = Merge[Flatten[Table[# -> i & /@ cells[[i]], {i, nC}]], Identity];
adj = Table[{}, {nC}];
Do[If[Length[f2c[k]] == 2,
   adj[[f2c[k][[1]]]] = Append[adj[[f2c[k][[1]]]], f2c[k][[2]]];
   adj[[f2c[k][[2]]]] = Append[adj[[f2c[k][[2]]]], f2c[k][[1]]]], {k, Keys[f2c]}];
adj = Union /@ adj;
cellVerts = Table[Union @@ faces[[cells[[i]]]], {i, nC}];

composeAff[{A2_, t2_}, {A1_, t1_}] := {A2 . A1, A2 . t1 + t2};
applyAff[{A_, t_}, x_] := A . x + t;
idAff4 = {IdentityMatrix[4], {0., 0., 0., 0.}};
rot4[e1_, e2_, th_] := IdentityMatrix[4] + (Cos[th]-1)(Outer[Times,e1,e1]+Outer[Times,e2,e2]) +
   Sin[th](Outer[Times,e2,e1] - Outer[Times,e1,e2]);
unfoldEdge[p_, c_, pTF_, vv_] := Module[
  {fI, fv, fvT, cf, rel, u1, u2, v, nP, nC2, cosT, perp, e1, e2, th},
  fI = First[Intersection[cells[[p]], cells[[c]]]];
  fv = vv[[faces[[fI]]]]; fvT = applyAff[pTF, #] & /@ fv; cf = Mean[fvT];
  rel = (# - cf) & /@ fvT; u1 = Normalize[rel[[1]]];
  u2 = Normalize[rel[[2]] - (rel[[2]].u1) u1];
  v = applyAff[pTF, Mean[vv[[cellVerts[[p]]]]]] - cf;
  nP = Normalize[v - (v.u1) u1 - (v.u2) u2];
  v = applyAff[pTF, Mean[vv[[cellVerts[[c]]]]]] - cf;
  nC2 = Normalize[v - (v.u1) u1 - (v.u2) u2];
  cosT = nC2 . (-nP); e1 = nC2; perp = -nP - cosT nC2;
  e2 = If[Norm[perp] > 10^-8, Normalize[perp], Normalize@First@NullSpace[{u1,u2,nC2}]];
  th = ArcCos[Clip[cosT, {-1,1}]];
  composeAff[{rot4[e1,e2,th], cf - rot4[e1,e2,th].cf}, pTF]];

faceUp[root_] := Module[{loc, cc, wHat = {0.,0.,0.,1.}, ang, R},
  loc = (# - Mean[vers0]) & /@ vers0;
  cc = Mean[loc[[cellVerts[[root]]]]];
  Which[Abs[(cc/Norm[cc]).wHat - 1] < 10^-5, loc,
        Abs[(cc/Norm[cc]).wHat + 1] < 10^-5, -loc,
        True, ang = ArcCos[Clip[cc.wHat/Norm[cc], {-1.,1.}]];
              R = RotationMatrix[ang, {cc/Norm[cc], wHat}]; R.# & /@ loc]];

chainPts[path_] := Module[{vv = faceUp[path[[1]]], tf = idAff4, out = {}},
  AppendTo[out, (applyAff[tf, #] & /@ vv[[cellVerts[[path[[1]]]]]])[[All, 1;;3]]];
  Do[tf = unfoldEdge[path[[i-1]], path[[i]], tf, vv];
     AppendTo[out, (applyAff[tf, #] & /@ vv[[cellVerts[[path[[i]]]]]])[[All, 1;;3]]],
     {i, 2, Length[path]}];
  out];

boxes[pts_] := {Min /@ Transpose[pts], Max /@ Transpose[pts]};
overlapPairs[path_] := Module[{P = chainPts[path], bb, res = {}},
  bb = boxes /@ P;
  Do[If[j > i+1 && And @@ Table[bb[[i,1,k]] < bb[[j,2,k]] - 10^-7 &&
                                bb[[j,1,k]] < bb[[i,2,k]] - 10^-7, {k,3}],
      If[RegionMeasure[RegionIntersection[ConvexHullMesh[P[[i]]], ConvexHullMesh[P[[j]]]]] > 10^-8,
         AppendTo[res, {i,j}]]], {i, Length[P]}, {j, Length[P]}];
  res];

(* vertex figure and its 8-cell paths, deduplicated by direction *)
fig = Select[Range[nC], MemberQ[cellVerts[[#]], 1] &];
figAdj = Association@Table[c -> Intersection[adj[[c]], fig], {c, fig}];
paths = {};
grow[p_] := If[Length[p] == 8, AppendTo[paths, p],
  Do[If[!MemberQ[p, nb], grow[Append[p, nb]]], {nb, figAdj[Last[p]]}]];
Do[grow[{s}], {s, fig}];
paths = DeleteDuplicates[paths, #1 === Reverse[#2] &];
Print["distinct 8-cell paths in one vertex figure: ", Length[paths]];

dists[r_] := Module[{d = ConstantArray[-1, nC], q = {r}, p},
  d[[r]] = 0; While[q =!= {}, p = First[q]; q = Rest[q];
    Do[If[d[[nb]] < 0, d[[nb]] = d[[p]]+1; AppendTo[q,nb]], {nb, adj[[p]]}]]; d];
D600 = Table[dists[r], {r, nC}];
pathOK[dd_] := Module[{st = Differences[dd]},
  Union[Abs[st]] === {1} && !MemberQ[Partition[st,2,1], {1,-1}]];
bfsable = Select[paths, Function[p, AnyTrue[Range[nC], pathOK[D600[[#, p]]] &]]];
Print["of these, fit in some BFS tree: ", Length[bfsable]];

{t, ov} = AbsoluteTiming[Table[overlapPairs[p], {p, paths}]];
Print["unfolded all ", Length[paths], " chains in ", Round[t,0.1], " s"];
bad = Select[Range[Length[paths]], ov[[#]] =!= {} &];
Print["chains whose unfolding overlaps: ", Length[bad], " of ", Length[paths]];
If[bad =!= {},
  Print["  example: path ", paths[[First[bad]]], "  overlapping pairs ", ov[[First[bad]]]]];
badBfs = Select[bad, MemberQ[bfsable, paths[[#]]] &];
Print[">>> chains that BOTH overlap AND fit in a BFS tree: ", Length[badBfs]];
If[badBfs =!= {},
  Print["    path ", paths[[First[badBfs]]], "  pairs ", ov[[First[badBfs]]]];
  Print["    a root that works: ",
     SelectFirst[Range[nC], pathOK[D600[[#, paths[[First[badBfs]]]]]] &]]];
