(* face_rotation_net_600cell.m
   Face-rotation net check for the 600-cell.
   Each tetrahedral cell is unfolded by rotating it around the shared
   triangular face into the common 3D hyperplane (BFS spanning tree).
   Self-contained: does not depend on prior session state.

   600-cell:  120 vertices, 1200 triangular faces, 600 tetrahedral cells
              each cell adjacent to 4 others (face-sharing)
*)

baseDir = "/Users/yoshino/Library/CloudStorage/Dropbox/260324Peeling4D/";

(* ===== 1. Load raw data ===== *)
Print["[1] Loading 600-cell data..."];
{tLoad, raw600} = AbsoluteTiming[Get[baseDir <> "_4DData/f600.m"]];
vers600  = raw600[[1]];
faces600 = raw600[[4]];
cells600 = raw600[[6]];
nC600    = Length[cells600];
Print["    done in ", tLoad, " s"];
Print["    vertices=", Length[vers600],
      "  faces=", Length[faces600],
      "  cells=", nC600,
      "  (expect 120, 1200, 600)"];

(* ===== 2. Cell geometry helpers ===== *)
cVertsIdx600[i_] := Union @@ faces600[[ cells600[[i]] ]];
cVerts600[i_]    := vers600[[ cVertsIdx600[i] ]];
cCentroid600[i_] := Mean @ cVerts600[i];

(* ===== 3. Centered vertices (root-independent base) ===== *)
versCentered600 = (# - Mean[N[vers600]]) & /@ N[vers600];

(* ===== 4. Cell adjacency ===== *)
Print["[2] Building cell adjacency..."];
{tAdj, adjList600} = AbsoluteTiming[
  Table[
    Select[Range[nC600],
      Function[j, j =!= i &&
        Length[Intersection[cells600[[i]], cells600[[j]]]] >= 1]],
    {i, nC600}]];
Print["    done in ", tAdj, " s"];
Print["    degree min=", Min[Length /@ adjList600],
      "  max=", Max[Length /@ adjList600],
      "  (expect 4 for all cells)"];

adjSet600 = Association @ Flatten[
  Table[Sort[{i,j}] -> True, {i, nC600}, {j, adjList600[[i]]}], 1];
adjQ600[i_,j_] := KeyExistsQ[adjSet600, Sort[{i,j}]];

(* Vertices of the shared triangular face between cells i and j *)
sharedFaceVerts600[i_,j_] :=
  vers600[[ faces600[[ First[Intersection[cells600[[i]], cells600[[j]]]] ]] ]];

(* ===== 5. 4D affine transform utilities ===== *)
composeAff600[{A2_,t2_},{A1_,t1_}] := {A2 . A1, A2 . t1 + t2};
applyAff600[{A_,t_},x_]            := A . x + t;
idAff600 = {IdentityMatrix[4], {0.,0.,0.,0.}};

rot4DInPlane600[e1_,e2_,theta_] :=
  IdentityMatrix[4] +
  (Cos[theta]-1)(Outer[Times,e1,e1] + Outer[Times,e2,e2]) +
  Sin[theta](Outer[Times,e2,e1] - Outer[Times,e1,e2]);

affRotAroundPt600[cf_,e1_,e2_,theta_] :=
  With[{R = rot4DInPlane600[N@e1, N@e2, N@theta]}, {R, cf - R . cf}];

(* ===== 6. Per-edge unfolding transform ===== *)
(* Identical algorithm to 120-cell version; shared face is now a triangle *)
unfoldTF600[pIdx_,cIdx_,pTF_] := Module[
  {fv, fvT, cf, u1, u2, rel, nPar, nChd, v, perp, e1, e2, cosT, theta},
  fv   = N @ sharedFaceVerts600[pIdx, cIdx];
  fvT  = applyAff600[pTF, #] & /@ fv;
  cf   = Mean[fvT];
  rel  = (# - cf) & /@ fvT;
  u1   = Normalize[rel[[1]]];
  u2   = Normalize[rel[[2]] - (rel[[2]] . u1) u1];
  v    = applyAff600[pTF, N @ cCentroid600[pIdx]] - cf;
  nPar = Normalize[v - (v . u1) u1 - (v . u2) u2];
  v    = applyAff600[pTF, N @ cCentroid600[cIdx]] - cf;
  nChd = Normalize[v - (v . u1) u1 - (v . u2) u2];
  cosT  = nChd . (-nPar);
  e1    = nChd;
  perp  = -nPar - cosT nChd;
  e2    = If[Norm[perp] > 10^-8,
             Normalize[perp],
             Normalize @ First @ NullSpace[{u1, u2, nChd}]];
  theta = ArcCos[Clip[cosT, {-1,1}]];
  composeAff600[affRotAroundPt600[cf, e1, e2, theta], pTF]
];

(* ===== 7. BFS unfolding ===== *)
bfsUnfold600[root_] := Module[
  {q = {root}, vis = {root},
   tfs = Association[root -> idAff600], par, nb},
  While[q =!= {},
    par = First[q]; q = Rest[q];
    Do[nb = neigh;
       If[!MemberQ[vis, nb],
         AppendTo[vis, nb]; AppendTo[q, nb];
         tfs[nb] = unfoldTF600[par, nb, tfs[par]]],
       {neigh, adjList600[[par]]}]];
  tfs
];

(* ===== 8. Face-up rotation per root ===== *)
faceUpForRoot600[r_, locV_] := Module[
  {ccIdx, cc, wHat = {0.,0.,0.,1.}, angle, R},
  ccIdx = Union @@ faces600[[ cells600[[r]] ]];
  cc    = Mean[locV[[ ccIdx ]]];
  Which[
    Abs[(cc/Norm[cc]) . wHat - 1] < 0.0001, locV,
    Abs[(cc/Norm[cc]) . wHat + 1] < 0.0001, -locV,
    True,
      angle = ArcCos[Clip[cc . wHat / Norm[cc], {-1.,1.}]];
      R     = RotationMatrix[angle, {cc/Norm[cc], wHat}];
      R . # & /@ locV
  ]
];

(* ===== 9. Full net check for one root ===== *)
checkRoot600[r_] := Module[
  {tBFS, tfsR, unf4D, wSpread, unf3D,
   bMinR, bMaxR, bboxList, nBBox,
   tHulls, hulls, tHits, trueHits},

  vers600 = faceUpForRoot600[r, versCentered600];
  {tBFS, tfsR} = AbsoluteTiming[bfsUnfold600[r]];

  unf4D   = Table[applyAff600[tfsR[i], #] & /@ N[cVerts600[i]], {i, nC600}];
  wSpread = (Max[#] - Min[#])& @ Flatten[unf4D[[All, All, 4]]];
  unf3D   = unf4D[[All, All, 1;;3]];

  bMinR    = (Min /@ Transpose[#]) & /@ unf3D;
  bMaxR    = (Max /@ Transpose[#]) & /@ unf3D;
  bboxList = Reap[
    Do[If[!adjQ600[i,j] &&
          And @@ Thread[bMinR[[i]] <= bMaxR[[j]]] &&
          And @@ Thread[bMinR[[j]] <= bMaxR[[i]]],
         Sow[{i,j}]],
       {i, nC600}, {j, i+1, nC600}]][[2]];
  bboxList = If[bboxList === {}, {}, First[bboxList]];
  nBBox    = Length[bboxList];

  If[nBBox == 0,
    Print["  root ", r, ": w-spread=", N[wSpread],
          "  bbox=0  *** VALID ***  (BFS ", tBFS, " s)"];
    Return[{r, 0, 0, True}]
  ];

  {tHulls, hulls} = AbsoluteTiming[
    Table[ConvexHullMesh[unf3D[[i]]], {i, nC600}]];
  {tHits, trueHits} = AbsoluteTiming[
    Select[bboxList, Function[p,
      With[{i = p[[1]], j = p[[2]]},
        AnyTrue[unf3D[[j]], RegionMember[hulls[[i]], #]&] ||
        AnyTrue[unf3D[[i]], RegionMember[hulls[[j]], #]&]]]]];

  Print["  root ", r, ": w-spread=", N[wSpread],
        "  bbox=", nBBox, "  true=", Length[trueHits],
        If[Length[trueHits] == 0, "  *** VALID ***", "  (INVALID)"],
        "  (BFS ", tBFS, " s, hulls ", tHulls, " s, hits ", tHits, " s)"];
  {r, nBBox, Length[trueHits], Length[trueHits] == 0}
];

(* ===== 10. Quick check: roots 1-5 with full output ===== *)
Print["\n[A] Detailed check for roots 1-5..."];
{tSample, sampleResults} = AbsoluteTiming[
  Table[checkRoot600[r], {r, 1, 5}]];
Print["    Done in ", tSample, " s"];

validSample = Select[sampleResults, #[[4]]&][[All,1]];
Print["    Valid among roots 1-5: ", validSample];

(* ===== 11. Full check: all 600 roots (bbox proxy) ===== *)
(* If roots 1-5 are all valid, run all 600 with bbox-only proxy for speed *)
Print["\n[B] Full check: all 600 roots (bbox proxy for speed)..."];
Print["    (roots where bbox=0 are certified VALID; others get full test)"];
{tAll, allResults} = AbsoluteTiming[
  Table[
    Module[{tfsR, unf3D, bMinR, bMaxR, nBBox = 0},
      vers600 = faceUpForRoot600[r, versCentered600];
      tfsR    = bfsUnfold600[r];
      unf3D   = Table[applyAff600[tfsR[i], #] & /@ N[cVerts600[i]],
                      {i, nC600}][[All, All, 1;;3]];
      bMinR   = (Min /@ Transpose[#]) & /@ unf3D;
      bMaxR   = (Max /@ Transpose[#]) & /@ unf3D;
      Do[If[!adjQ600[i,j] &&
            And @@ Thread[bMinR[[i]] <= bMaxR[[j]]] &&
            And @@ Thread[bMinR[[j]] <= bMaxR[[i]]],
           nBBox++],
         {i, nC600}, {j, i+1, nC600}];
      {r, nBBox}],
    {r, 1, nC600}]];
Print["    Done in ", tAll, " s"];

zeroBBox  = Select[allResults, #[[2]] == 0 &];
nonzero   = Select[allResults, #[[2]] > 0  &];
Print["    Roots with bbox=0 (certainly VALID): ", Length[zeroBBox], " / ", nC600];
Print["    Roots with bbox>0 (need full test):  ", Length[nonzero]];
If[Length[nonzero] > 0,
  Print["    bbox counts for bbox>0 roots: ",
        MinMax[nonzero[[All, 2]]]]];

(* ===== 12. Full intersection test for bbox>0 roots ===== *)
If[Length[nonzero] > 0,
  Print["\n[C] Full intersection test for ", Length[nonzero], " roots with bbox>0..."];
  {tFull, fullResults} = AbsoluteTiming[
    Table[checkRoot600[r[[1]]], {r, nonzero}]];
  Print["    Done in ", tFull, " s"];
  validFull  = Select[fullResults, #[[4]]&][[All,1]];
  invalidFull = Select[fullResults, !#[[4]]&][[All,1]];
  Print["    Valid: ", Length[validFull], "  Invalid: ", Length[invalidFull]];
  If[Length[invalidFull] > 0,
    Print["    First few invalid roots: ",
          Take[invalidFull, Min[5, Length[invalidFull]]]]],

  Print["\n[C] All roots have bbox=0 -- skipping full test"]
];

(* ===== 13. Final summary ===== *)
totalValid = Length[zeroBBox] +
  If[Length[nonzero] > 0, Length[validFull], 0];
Print["\n========================================"];
Print["600-cell face-rotation BFS net summary"];
Print["  Total roots checked: ", nC600];
Print["  Valid (no self-intersection): ", totalValid, " / ", nC600];
Print["  Result: ",
  If[totalValid == nC600,
    "*** ALL VALID -- 600-cell admits face-rotation nets from every root ***",
    If[totalValid == 0,
      "*** ALL INVALID -- face-rotation BFS cannot net the 600-cell ***",
      "Mixed: " <> ToString[totalValid] <> " valid, " <>
        ToString[nC600 - totalValid] <> " invalid"]]];
Print["========================================"];
