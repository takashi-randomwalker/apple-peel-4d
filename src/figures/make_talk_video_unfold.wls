#!/usr/bin/env wolframscript
(* make_talk_video_unfold.wls
   The 8-cell opening from the polytope into its net, for slide 3 of the
   JCDCG^3 talk -- the slide whose own speaker note says the audience has
   to picture the rotation.

   HOW THE PARTIAL UNFOLDING IS DEFINED
   For a tree edge (p,c) the child is rotated about the shared square by
   the angle theta that brings it into the parent's hyperplane.  At time
   t we rotate by t*theta instead, recursively:

       T_r(t) = id ,      T_c(t) = rot(t*theta_e) o T_p(t) .

   The hinge plane is recomputed in the parent's *current* frame, and
   theta is measured against the child in its *unrotated* position, so
   theta does not depend on t.  Every cell therefore stays rigid and
   stays hinged to its parent throughout; t=0 is the polytope itself and
   t=1 is the net of src/figures/make_talk_video_8cell.wls.

   WHY NOT SIMPLY DROP w
   Step 3 of the algorithm drops w.  That is right at t=1, where every
   vertex shares one w, but at t=0 it flattens the six side cells to
   squares and drops the two end cells onto each other, so the polytope
   would read as a single cube and the animation would open out of a
   featureless block.  We instead put an eye on the w axis beyond the
   polytope and project

       (x,y,z,w)  |-->  (x,y,z) * dScal / (eyeW - w) .

   At t=0 this is the picture everyone knows: the near cell as the outer
   cube, the far cell as the inner one, the six side cells as the
   frustums between them.  At t=1 every vertex shares one w, so the map
   is a uniform scaling, and dScal is chosen to make that scaling 1 --
   the last frame is the net at its true size, the same one that
   make_talk_video_8cell.wls renders.  Between the two the cells are
   distorted, as anything drawn from 4D must be.

   The timeline unfolds, holds, refolds and holds, so the video loops
   without a snap.  Refolding is also worth watching: it is the claim of
   Remark "Why overlap must be excluded" that the rotations run
   backwards and reassemble P.
*)

baseDir = "/Users/yoshino/Library/CloudStorage/Dropbox/260324Peeling4D/";

outDir = If[Length[$ScriptCommandLine] >= 2, $ScriptCommandLine[[2]],
            baseDir <> "attic/outputs/frames_unfold"];
(* How far past the polytope the 4D eye sits, as a fraction of the w
   range.  Small = strong perspective and a small inner cell; large =
   nearly orthographic, and the folded state collapses onto itself. *)
persMargin = If[Length[$ScriptCommandLine] >= 3,
                ToExpression[$ScriptCommandLine[[3]]], 0.7];
If[!DirectoryQ[outDir], CreateDirectory[outDir, CreateIntermediateDirectories -> True]];

fps     = 30;
totalT  = 10.0;
nFrames = Round[fps totalT];
(* Rendered oversized because build_talk_video.sh crops the frames to
   what is actually drawn before scaling them down to 600: an explicit
   PlotRange may not be tightened past the geometry without clipping it,
   so the last third of the picture can only be recovered afterwards. *)
imgSize = 900;
bgCol   = GrayLevel[0.13];

(* ── 4D affine utilities (from face_rotation_net_viz_all.m) ─────────────── *)
composeAff[{A2_, t2_}, {A1_, t1_}] := {A2 . A1, A2 . t1 + t2};
applyAff[{A_, t_}, x_]             := A . x + t;
idAff4 = {IdentityMatrix[4], {0., 0., 0., 0.}};

rot4DInPlane[e1_, e2_, theta_] :=
  IdentityMatrix[4] +
  (Cos[theta] - 1)(Outer[Times, e1, e1] + Outer[Times, e2, e2]) +
  Sin[theta](Outer[Times, e2, e1] - Outer[Times, e1, e2]);

affRotAroundPt[cf_, e1_, e2_, theta_] :=
  With[{R = rot4DInPlane[N @ e1, N @ e2, N @ theta]}, {R, cf - R . cf}];

(* Same as unfoldEdge, but turning only the fraction s of the way. *)
unfoldEdgePartial[pIdx_, cIdx_, pTF_, versL_, facesL_, cellsL_, s_] :=
  Module[{fIdx, fv, fvT, cf, rel, u1, u2, v, nP, nC, cosT, perp, e1, e2, theta},
    fIdx = First[Intersection[cellsL[[pIdx]], cellsL[[cIdx]]]];
    fv   = N @ versL[[ facesL[[fIdx]] ]];
    fvT  = applyAff[pTF, #] & /@ fv;
    cf   = Mean[fvT];
    rel  = (# - cf) & /@ fvT;
    u1   = Normalize[rel[[1]]];
    u2   = Normalize[rel[[2]] - (rel[[2]] . u1) u1];
    v    = applyAff[pTF, N @ Mean[versL[[ Union @@ facesL[[ cellsL[[pIdx]] ]] ]]]] - cf;
    nP   = Normalize[v - (v . u1) u1 - (v . u2) u2];
    v    = applyAff[pTF, N @ Mean[versL[[ Union @@ facesL[[ cellsL[[cIdx]] ]] ]]]] - cf;
    nC   = Normalize[v - (v . u1) u1 - (v . u2) u2];
    cosT  = nC . (-nP);
    e1    = nC;
    perp  = -nP - cosT nC;
    e2    = If[Norm[perp] > 10^-8, Normalize[perp],
               Normalize @ First @ NullSpace[{u1, u2, nC}]];
    theta = s ArcCos[Clip[cosT, {-1, 1}]];
    composeAff[affRotAroundPt[cf, e1, e2, theta], pTF]
  ];

(* ── Load the 8-cell, face-up on root cell 1 ────────────────────────────── *)
raw   = Get[baseDir <> "_4DData/f8.m"];
vers  = N @ raw[[1]];
faces = raw[[4]];
cells = raw[[6]];
nCell = Length[cells];

adjList = Table[
  Select[Range[nCell],
    Function[j, j =!= i && Length[Intersection[cells[[i]], cells[[j]]]] >= 1]],
  {i, nCell}];

wHat = {0., 0., 0., 1.};
locV = (# - Mean[vers]) & /@ vers;
cc1  = Mean[locV[[ Union @@ faces[[ cells[[1]] ]] ]]];
vers = Which[
  Abs[(cc1 / Norm[cc1]) . wHat - 1] < 0.0001, locV,
  Abs[(cc1 / Norm[cc1]) . wHat + 1] < 0.0001, -locV,
  True,
    Module[{angle, R},
      angle = ArcCos[Clip[cc1 . wHat / Norm[cc1], {-1., 1.}]];
      R     = RotationMatrix[angle, {cc1 / Norm[cc1], wHat}];
      R . # & /@ locV]];

cellVerts = Table[Union @@ faces[[ cells[[i]] ]], {i, nCell}];

(* BFS tree: parents and layers, from root cell 1. *)
parent   = Association[1 -> None];
bfsLayer = Association[1 -> 0];
order    = {1};
bfsQ     = {1};
While[bfsQ =!= {},
  Module[{p = First[bfsQ]},
    bfsQ = Rest[bfsQ];
    Do[If[!KeyExistsQ[parent, nb],
        parent[nb]   = p;
        bfsLayer[nb] = bfsLayer[p] + 1;
        AppendTo[order, nb];
        AppendTo[bfsQ, nb]],
      {nb, adjList[[p]]}]]];
maxLayer = Max[Values[bfsLayer]];
Print["cells=", nCell, "  BFS layers 0..", maxLayer,
      "  sizes=", Values[Counts[Values[bfsLayer]]]];

(* Cell coordinates in 4D, given one unfolding fraction per BFS layer:
   the edge into a cell of layer d turns by sv[[d]] of its angle.  A
   single number for all layers would open the whole tree at once; a
   vector lets each layer wait for the one above it, which is the point
   of the staging in the timeline below. *)
coords4[sv_List] := Module[{tfs = Association[1 -> idAff4]},
  Do[If[c =!= 1,
      tfs[c] = unfoldEdgePartial[parent[c], c, tfs[parent[c]],
                                 vers, faces, cells, sv[[ bfsLayer[c] ]]]],
    {c, order}];
  Table[applyAff[tfs[i], #] & /@ vers[[ cellVerts[[i]] ]], {i, nCell}]];

allOpen = ConstantArray[1., maxLayer];

(* ── Timeline ──────────────────────────────────────────────────────────── *)
(* The layers open one after another, with a beat of stillness between
   them: the six neighbours of the root swing out and stop, and only
   then does the last cell -- the one two steps from the root -- follow.
   That pause is the whole reason the film is staged rather than opened
   all at once; it is what makes the BFS layers visible. *)
smooth[u_] := u^2 (3 - 2 u);

tUnf0 = 0.5;   tUnf1 = 6.0;    (* the opening occupies this window   *)
pauseD = 0.6;                  (* the beat between two layers        *)
tOpen1 = 7.8;  tFold1 = 9.4;   (* hold open, then close in one go    *)
stageD = (tUnf1 - tUnf0 - pauseD (maxLayer - 1)) / maxLayer;
stageStart[d_] := tUnf0 + (d - 1) (stageD + pauseD);

sVec[tt_] := Which[
  tt < tOpen1,
    Table[smooth[Clip[(tt - stageStart[d]) / stageD, {0, 1}]], {d, maxLayer}],
  tt < tFold1,
    ConstantArray[1. - smooth[(tt - tOpen1) / (tFold1 - tOpen1)], maxLayer],
  True, ConstantArray[0., maxLayer]];

(* One scalar for the camera and for fading the perspective out: the
   staging must not make either of them stutter. *)
camP[tt_] := Which[
  tt < tOpen1, smooth[Clip[(tt - tUnf0) / (tUnf1 - tUnf0), {0, 1}]],
  tt < tFold1, 1. - smooth[(tt - tOpen1) / (tFold1 - tOpen1)],
  True, 0.];

times  = Table[k / fps, {k, 0, nFrames - 1}];
svecs  = sVec /@ times;
cams   = camP /@ times;

(* ── Projection: perspective at the start, drop-w at the end ───────────── *)
netPts = coords4[allOpen];
wSpread = Max[Flatten[netPts[[All, All, 4]]]] - Min[Flatten[netPts[[All, All, 4]]]];
Print["w-spread of the net = ", wSpread];

(* The eye sits on the $w$ axis beyond every vertex the animation ever
   reaches; margin is in units of the whole w range. *)
wAll  = Flatten[Table[coords4[sVec[tt]][[All, All, 4]], {tt, 0., totalT, 0.05}]];
wLo   = Min[wAll]; wHi = Max[wAll];
wNet  = Mean[Flatten[netPts[[All, All, 4]]]];
eyeW  = wHi + persMargin (wHi - wLo);
dScal = eyeW - wNet;   (* so that the net comes out at its true size *)
Print["w in [", wLo, ", ", wHi, "]   net at w=", wNet,
      "   eye at w=", eyeW, "   inner/outer = ",
      (eyeW - wHi) / (eyeW - wLo)];

(* The perspective is faded out as the net opens.  Kept at full strength
   throughout, the cells that swing down to w = -2.5 would recede from
   the eye and shrink by a third on the way, so the net would appear to
   contract before it spreads.  At s=1 the weight is zero and the map is
   the plain drop-w of Step 3. *)
proj[pts_, s_] := With[{a = 1. - s},
  (((1 - a) + a dScal / (eyeW - #[[4]])) Take[#, 3]) & /@ pts];

(* Triangulation is computed once, on the net, where no cell is
   degenerate.  The cells are rigid in 4D, so the same vertex indices
   describe each cell's boundary at every t; only the projection moves. *)
tris = Table[
  MeshCells[ConvexHullMesh[netPts[[i, All, 1 ;; 3]]], 2][[All, 1]], {i, nCell}];

allPts = Table[proj[#, cams[[k]]] & /@ coords4[svecs[[k]]], {k, nFrames}];

(* The camera pulls back as the net opens.  With one range for the whole
   film the tesseract would occupy a fifth of the frame and be
   unreadable; with a range refitted every frame the growth would not be
   felt at all.  So the centre and the half-width are interpolated
   linearly in the unfolding parameter, then widened frame by frame
   wherever that would clip -- the max of two continuous functions, so
   the camera never jumps. *)
(* An explicit PlotRange CLIPS, so the box must contain the frame; it is
   fitted per axis rather than as a cube, because a cube wastes a fifth
   of the picture on empty corners.  The box is interpolated from the
   folded frame to the whole film and then, axis by axis, widened to
   whatever the frame actually needs -- again a max of continuous
   functions, so the camera never jumps. *)
(* Takes anything whose leaves are 3-vectors: one frame, or all of them. *)
bbox[pts_] := Module[{f = Partition[Flatten[N @ pts], 3]},
  Transpose[{Min /@ Transpose[f], Max /@ Transpose[f]}]];

gBox  = bbox[allPts];
fBox  = bbox[allPts[[1]]];
pad   = 0.02 Max[gBox[[All, 2]] - gBox[[All, 1]]];
zoom0 = 1.35;   (* room left around the folded polytope *)

grow[b_, f_] := Module[{c = Mean /@ b},
  Table[{c[[j]] + f (b[[j, 1]] - c[[j]]),
         c[[j]] + f (b[[j, 2]] - c[[j]])}, {j, 3}]];

fBox0 = grow[fBox, zoom0];
boxes = Table[
  Module[{lin = fBox0 + (gBox - fBox0) cams[[k]], own = bbox[allPts[[k]]]},
    Table[{Min[lin[[j, 1]], own[[j, 1]] - pad],
           Max[lin[[j, 2]], own[[j, 2]] + pad]}, {j, 3}]],
  {k, nFrames}];
Print["box at the start ", First[boxes]];
Print["box at the end   ", boxes[[Round[7 fps]]]];

plotRange[k_] := boxes[[k]];

frame[k_] := Graphics3D[
  Table[{Opacity[0.72], EdgeForm[{GrayLevel[0.3], Thin}],
         GrayLevel[0.3 + 0.6 bfsLayer[i] / maxLayer],
         GraphicsComplex[allPts[[k + 1, i]], Polygon @ tris[[i]]]},
        {i, nCell}],
  Boxed           -> False,
  Lighting        -> "Neutral",
  ImageSize       -> {imgSize, imgSize},
  PlotRange       -> plotRange[k + 1],
  PlotRangePadding-> 0,
  ViewPoint       -> {2.4, -2.0, 1.8},
  ViewVertical    -> {0, 0, 1},
  (* No SphericalRegion here, unlike the turntable: it fits the sphere
     around the plot range, which for a cube costs a factor sqrt(3) and
     leaves the net occupying half the frame.  It is only needed when
     the camera moves, and here it does not; the explicit PlotRange
     already keeps the framing from jumping. *)
  Background      -> bgCol];

Do[
  Export[FileNameJoin[{outDir, StringTemplate["f``.png"][IntegerString[k, 10, 4]]}],
         frame[k], ImageSize -> {imgSize, imgSize}];
  If[Mod[k, 50] == 0, Print["  frame ", k, "/", nFrames]],
  {k, 0, nFrames - 1}];

Print["Wrote ", nFrames, " frames to ", outDir];
