(* face_rotation_net_120cell.m
   Experiment: does the 120-cell have valid face-rotation nets?
   A "face-rotation net" unfolds each dodecahedral cell by rotating it
   around the shared pentagonal face into the common 3D hyperplane.
   Method: BFS spanning tree from a root cell; compose 4D affine transforms;
           then check for self-intersection in the resulting 3D layout.
*)

baseDir = "/Users/yoshino/Library/CloudStorage/Dropbox/260324Peeling4D/";

(* ===== 1. Load data ===== *)
Print["[1] Loading 120-cell data..."];
(* f120.m returns a ≥6-element list; indices confirmed in run4DGlobal_120cell.m *)
Module[{raw = Get[baseDir <> "_4DData/f120.m"]},
  vers  = raw[[1]];
  faces = raw[[4]];
  cells = raw[[6]]];
nC = Length[cells];  (* 120 *)
Print["    cells=", nC, "  faces=", Length[faces], "  vertices=", Length[vers]];

(* Cell geometry helpers *)
cVertsIdx[i_] := Union @@ faces[[ cells[[i]] ]];
cVerts[i_]    := vers[[ cVertsIdx[i] ]];
cCentroid[i_] := Mean @ cVerts[i];

(* ===== 2. Cell adjacency ===== *)
Print["[2] Building cell adjacency (degree should be 12 for all cells)..."];
adjList = Table[
  Select[Range[nC],
    Function[j, j =!= i && Length[Intersection[cells[[i]], cells[[j]]]] >= 1]],
  {i, nC}];
{Print["    degree min=", Min[Length /@ adjList],
      "  max=", Max[Length /@ adjList]]};

(* Fast O(1) adjacency test *)
adjSet = Association @ Flatten[
  Table[Sort[{i, j}] -> True, {i, nC}, {j, adjList[[i]]}], 1];
adjQ[i_, j_] := KeyExistsQ[adjSet, Sort[{i, j}]];

(* Vertices of the shared pentagonal face between cells i and j *)
sharedPentVerts[i_, j_] :=
  vers[[ faces[[ First[Intersection[cells[[i]], cells[[j]]]] ]] ]];

(* ===== 2b. Face-up rotation =====
   Same approach as peeling4Df4.m: RotationMatrix[angle, {u,v}] is more
   robust than the manual Rodrigues formula for near-degenerate cases. *)
Print["[2b] Applying face-up rotation (root = cell 1)..."];
Module[{locV, cc1Idx, cc1, wHat = {0., 0., 0., 1.}, angle, R},
  locV   = (# - Mean[N[vers]]) & /@ N[vers];           (* center at origin, numerical *)
  cc1Idx = Union @@ faces[[ cells[[1]] ]];              (* vertex indices of cell 1 *)
  cc1    = Mean[ locV[[ cc1Idx ]] ];                    (* centroid of cell 1 *)
  Print["    Cell 1 centroid (pre-rotation): ", cc1, "  (expect large nonzero)"];
  locV   = Which[
    Abs[(cc1 / Norm[cc1]) . wHat - 1] < 0.0001,        (* already aligned to +w *)
      locV,
    Abs[(cc1 / Norm[cc1]) . wHat + 1] < 0.0001,        (* anti-aligned: flip sign *)
      -locV,
    True,
      angle = ArcCos[Clip[cc1 . wHat / Norm[cc1], {-1., 1.}]];
      R     = RotationMatrix[angle, {cc1 / Norm[cc1], wHat}];
      R . # & /@ locV
  ];
  vers = locV;
];
Print["    Cell 1 centroid after face-up: ",
      N @ Mean[ vers[[ Union @@ faces[[ cells[[1]] ]] ]] ],
      "  (expect {~0,~0,~0,r>0})"];

(* ===== 3. 4D affine transform utilities ===== *)
(* Represent an affine map as {matrix 4x4, translation 4-vector} *)
composeAff[{A2_, t2_}, {A1_, t1_}] := {A2 . A1, A2 . t1 + t2};
applyAff[{A_, t_}, x_]             := A . x + t;
idAff = {IdentityMatrix[4], {0., 0., 0., 0.}};

(* 4D rotation in the 2D plane span{e1, e2} by angle theta.
   Convention: e1 -> cos(theta) e1 + sin(theta) e2 *)
rot4DInPlane[e1_, e2_, theta_] :=
  IdentityMatrix[4] +
  (Cos[theta] - 1)(Outer[Times, e1, e1] + Outer[Times, e2, e2]) +
  Sin[theta](Outer[Times, e2, e1] - Outer[Times, e1, e2]);

(* Affine rotation around a fixed point cf *)
affRotAroundPt[cf_, e1_, e2_, theta_] :=
  With[{R = rot4DInPlane[N@e1, N@e2, N@theta]}, {R, cf - R . cf}];

(* ===== 4. Per-edge unfolding transform ===== *)
(*
  Given:
    pIdx      - parent cell index (already unfolded by pTF)
    cIdx      - child cell index (to be unfolded)
    pTF       - parent's cumulative 4D affine transform

  Returns: child's cumulative transform = R_local . pTF
  where R_local rotates the child across the shared pentagon into the
  parent's 3D hyperplane.

  The rotation maps the child's inward normal nC to -nP (parent's inward
  normal reflected), bringing the child's hyperplane flush with the parent's.
  Rotation angle = pi - dihedral_angle = arccos(1/sqrt(5)) ~ 63.43 deg.
*)
unfoldTF[pIdx_, cIdx_, pTF_] := Module[
  {fv, fvT, cf, u1, u2, rel, nP, nC, v, perp, e1, e2, cosT, theta},

  (* Shared pentagon in the parent's unfolded frame *)
  fv  = N @ sharedPentVerts[pIdx, cIdx];
  fvT = applyAff[pTF, #] & /@ fv;
  cf  = Mean[fvT];

  (* Orthonormal basis of the pentagon plane in 4D *)
  rel = (# - cf) & /@ fvT;
  u1  = Normalize[rel[[1]]];
  u2  = Normalize[rel[[2]] - (rel[[2]] . u1) u1];

  (* Normals from face toward each cell, projected onto face complement *)
  v  = applyAff[pTF, N @ cCentroid[pIdx]] - cf;
  nP = Normalize[v - (v . u1) u1 - (v . u2) u2];

  v  = applyAff[pTF, N @ cCentroid[cIdx]] - cf;
  nC = Normalize[v - (v . u1) u1 - (v . u2) u2];

  (* Rotation that maps nC -> -nP (unfold child into parent hyperplane) *)
  cosT = nC . (-nP);
  e1   = nC;
  perp = -nP - cosT nC;
  e2   = If[Norm[perp] > 10^-8,
            Normalize[perp],
            Normalize @ First @ NullSpace[{u1, u2, nC}]];  (* fallback *)
  theta = ArcCos[Clip[cosT, {-1, 1}]];

  composeAff[affRotAroundPt[cf, e1, e2, theta], pTF]
];

(* ===== 5. BFS unfolding from a root cell ===== *)
bfsUnfold[root_] := Module[
  {q = {root}, vis = {root}, tfs = Association[root -> idAff], par, nb},
  While[q =!= {},
    par = First[q];  q = Rest[q];
    Do[nb = neigh;
       If[!MemberQ[vis, nb],
         AppendTo[vis, nb];  AppendTo[q, nb];
         tfs[nb] = unfoldTF[par, nb, tfs[par]]],
       {neigh, adjList[[par]]}]];
  tfs
];

(* ===== 6. Analyze BFS net from root = 1 ===== *)
Print["[5] BFS unfolding from cell 1..."];
{t1, tfs} = AbsoluteTiming[bfsUnfold[1]];
Print["    Transforms computed: ", Length[tfs], "  time: ", t1, " s"];

Print["[6] Extracting unfolded 4D coordinates..."];
unf4D = Table[applyAff[tfs[i], #] & /@ N[cVerts[i]], {i, nC}];

(* Sanity check 1: w-spread across all unfolded cells should be ~0 *)
wVals = Flatten[unf4D[[All, All, 4]]];
Print["    w-coordinate spread: ", (Max[wVals] - Min[wVals]) // N,
      "  (should be ~0 if all cells lie in one 3D hyperplane)"];

(* Sanity check 2: shared face between cells 1 and its first neighbor
   should coincide after transformation *)
nb1  = First[adjList[[1]]];
sf1  = applyAff[tfs[1],   #] & /@ N[sharedPentVerts[1, nb1]];
sf2  = applyAff[tfs[nb1], #] & /@ N[sharedPentVerts[1, nb1]];
Print["    Shared face mismatch (cells 1 & ", nb1, "): ",
      Max[Norm /@ (Sort[sf1] - Sort[sf2])] // N,
      "  (should be ~0)"];

(* 3D layout: drop the w coordinate (all cells share the same w after unfolding) *)
unf3D = unf4D[[All, All, 1;;3]];

(* ===== 7. Bounding-box overlap detection ===== *)
Print["[7] Bounding-box collision check..."];
bMin = (Min /@ Transpose[#]) & /@ unf3D;
bMax = (Max /@ Transpose[#]) & /@ unf3D;

bboxOK[i_, j_] :=
  And @@ Thread[bMin[[i]] <= bMax[[j]]] &&
  And @@ Thread[bMin[[j]] <= bMax[[i]]];

bboxPairs = {};
{t2, dummy$} = AbsoluteTiming[
  Do[If[!adjQ[i, j] && bboxOK[i, j], AppendTo[bboxPairs, {i, j}]],
     {i, nC}, {j, i + 1, nC}]];
Print["    Non-adjacent bbox-overlapping pairs: ", Length[bboxPairs],
      "  (time: ", t2, " s)"];

(* ===== 8. Detailed intersection test ===== *)
If[Length[bboxPairs] == 0,

  Print["*** No bbox overlaps: BFS net from cell 1 is certainly VALID! ***"],

  (* Build convex hull meshes once *)
  Print["[8] Building ", nC, " convex hull meshes..."];
  {t3, hulls} = AbsoluteTiming[Table[ConvexHullMesh[unf3D[[i]]], {i, nC}]];
  Print["    Done in ", t3, " s"];

  Print["    Checking ", Length[bboxPairs], " candidate pairs..."];
  {t4, trueHits} = AbsoluteTiming[
    Select[bboxPairs, Function[p,
      With[{i = p[[1]], j = p[[2]]},
        AnyTrue[unf3D[[j]], RegionMember[hulls[[i]], #] &] ||
        AnyTrue[unf3D[[i]], RegionMember[hulls[[j]], #] &]]]]];
  Print["    True intersecting pairs: ", Length[trueHits],
        "  (time: ", t4, " s)"];

  If[Length[trueHits] == 0,
    Print["*** BFS net from cell 1 is VALID (no self-intersection)! ***"],
    Print["    Net has self-intersections."];
    Print["    First few intersecting pairs: ",
          Take[trueHits, Min[5, Length[trueHits]]]]
  ]
];

(* ===== 9. Sample nets from several root cells (bbox proxy only) ===== *)
Print["[9] Sampling nets from root cells {1, 2, 3, 4, 5} (bbox proxy)..."];
sampleRoots = {1, 2, 3, 4, 5};
sampleResults = Table[
  Module[{tfsR, unf3DR, bMinR, bMaxR, nBBox = 0},
    tfsR   = bfsUnfold[root];
    unf3DR = Table[applyAff[tfsR[i], #] & /@ N[cVerts[i]], {i, nC}][[All, All, 1;;3]];
    bMinR  = (Min /@ Transpose[#]) & /@ unf3DR;
    bMaxR  = (Max /@ Transpose[#]) & /@ unf3DR;
    Do[If[!adjQ[i, j] &&
          And @@ Thread[bMinR[[i]] <= bMaxR[[j]]] &&
          And @@ Thread[bMinR[[j]] <= bMaxR[[i]]],
         nBBox++],
       {i, nC}, {j, i + 1, nC}];
    {root, nBBox}],
  {root, sampleRoots}];
Print["    root -> bbox_overlaps: ", sampleResults];

(* ===== 10. Visualize the first net ===== *)
Print["[10] Generating 3D visualization of net from cell 1..."];
(* Rebuild hulls if not already built *)
If[!NameQ["hulls"],
  hulls = Table[ConvexHullMesh[unf3D[[i]]], {i, nC}]];

netViz = Graphics3D[
  Table[{Opacity[0.55], ColorData["Rainbow"][i/nC],
    GraphicsComplex[
      unf3D[[i]],
      Polygon @ MeshCells[hulls[[i]], 2][[All, 1]]]},
    {i, nC}],
  Boxed -> False, Lighting -> "Neutral",
  PlotLabel -> "120-cell face-rotation net (BFS from cell 1)"];

Export[baseDir <> "face_rotation_net_120cell.png", netViz, ImageSize -> 1000];
Print["    Saved: face_rotation_net_120cell.png"];

Print["Done."];
