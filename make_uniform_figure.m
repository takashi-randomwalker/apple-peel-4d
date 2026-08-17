(* =====================================================================
   make_uniform_figure.m
   Provisional figure for Section 4 (The Uniform Case).

   Two panels, chosen to show the two ways Theorem 1 fails:
     (a) runcinated 24-cell t_{0,3}{3,4,3}, 240 cells, at its BEST root
         (3 intersecting pairs, the fewest of any of its 240 roots).
         No root works for this polytope.
     (b) runcitruncated 120-cell t_{0,1,3}{5,3,3}, 2640 cells, at a
         FAILING root.  2570 of its roots are fine; the failures all
         have exactly one intersecting pair, so the panel is a close-up
         of that single pair.

   Intersecting cells are drawn solid; the rest of the net is drawn
   translucent so the failure is visible inside it.
   ===================================================================== *)

baseDir = "/Users/yoshino/Library/CloudStorage/Dropbox/260324Peeling4D/";
SetDirectory[baseDir];

Get[baseDir <> "wythoff_gen_all.mx"];
Get[baseDir <> "wythoff_bfs_allroots.mx"];
Get[baseDir <> "h4_bfs_allroots.mx"];
Module[{src, cut},
  src = Import[baseDir <> "face_rotation_net_all4D_v2.m", "Text"];
  cut = First@First@StringPosition[src, "Print[\"Face-rotation BFS net check v2"];
  ToExpression[StringTake[src, {1, cut - 1}]]];
Get[baseDir <> "bfs_sat.m"];

setup[sym_] := Module[{rec, data, faceToCells, pairs},
  rec = SelectFirst[w4All, #["sym"] === sym &];
  data = rec["data"];
  versP  = N[data[[1]], MachinePrecision];
  facesP = data[[4]];
  cellsP = data[[6]];
  nCP    = Length[cellsP];
  versCenteredP = (# - Mean[versP]) & /@ versP;
  faceToCells = GroupBy[
    Flatten@MapIndexed[Function[{cf, k}, {# -> First[k]} & /@ cf], cellsP],
    First -> Last];
  pairs = Select[Values[faceToCells], Length[#] === 2 &];
  adjListP = Union /@ Lookup[
    GroupBy[Join[pairs, Reverse /@ pairs], First -> Last], Range[nCP], {}];
  adjSetP = Association@Flatten[
    Table[Sort[{i, j}] -> True, {i, nCP}, {j, adjListP[[i]]}], 1];
  satPrepCells[]];

(* unfolded coordinates and the intersecting pairs, for one root *)
netAtRoot[r_] := Module[{tfsR, unf3D, cent, radius, bMin, bMax, cand, ne, hits},
  versP = satFaceUp[r, versCenteredP];
  tfsR  = satBfsUnfold[r];
  unf3D = Table[
    Module[{tf = tfsR[i]},
     Take[#, 3] & /@ (versP[[satCellIdx[[i]]]] . Transpose[tf[[1]]] +
                      Threaded[tf[[2]]])], {i, nCP}];
  bMin = (Min /@ Transpose[#]) & /@ unf3D;
  bMax = (Max /@ Transpose[#]) & /@ unf3D;
  cent = Mean /@ unf3D;
  radius = 2 Max[MapThread[Function[{p, c}, Max[Norm[# - c] & /@ p]],
                           {unf3D, cent}]];
  cand = satBBoxPairs[cent, bMin, bMax, radius];
  ne = Association@Table[
     c -> satNormalsEdges[unf3D[[c]], satFaceTri[[c]], satEdgePair[[c]]],
     {c, Union @@ cand}];
  hits = Select[cand, Function[p,
     satOverlapQ[unf3D[[p[[1]]]], ne[p[[1]]][[1]], ne[p[[1]]][[2]],
                 unf3D[[p[[2]]]], ne[p[[2]]][[1]], ne[p[[2]]][[2]]]]];
  {unf3D, hits}];

(* polygons of cell i in unfolded coordinates *)
cellPolys[unf3D_, i_] := Module[{loc},
  loc = AssociationThread[satCellIdx[[i]] -> Range[Length[satCellIdx[[i]]]]];
  Polygon[unf3D[[i]][[Lookup[loc, facesP[[#]]]]]] & /@ cellsP[[i]]];

(* The two members of an intersecting pair get different colours, so the
   interpenetration is legible; everything else is faint context.

   Two styles.  "black" matches Figure 1 (260612AllFaceRotation.pdf),
   which is white-on-black, so the context cells have to be light rather
   than dark; the highlight colours are brightened to survive the dark
   ground.  "white" is the same figure for a white page. *)
style["white"] = <|
  "bg" -> White,
  "ctxFace" -> GrayLevel[0.85, 0.05], "ctxEdge" -> GrayLevel[0.55, 0.35],
  "hotEdge" -> GrayLevel[0.1],
  "hotA" -> RGBColor[0.82, 0.12, 0.12, 0.92],
  "hotB" -> RGBColor[0.13, 0.35, 0.75, 0.92]|>;
(* "black" is measured off Figure 1: the background there is
   GrayLevel[0.13], not pure black, the cells are opaque light gray, and
   no edges are drawn -- faces are separated by shading alone.  We keep
   the no-edge treatment and the background, and trade full opacity for
   0.25 because this figure has to show cells buried inside the net.
   Drawn edges were tried and rejected: they accumulate into a
   wireframe and lose Figure 1's shaded look. *)
style["black"] = <|
  "bg" -> GrayLevel[0.13],
  "ctxGray" -> 0.92, "ctxEdge" -> None, "hotEdge" -> None,
  "hotA" -> RGBColor[1.00, 0.28, 0.24, 0.95],
  "hotB" -> RGBColor[0.35, 0.62, 1.00, 0.95]|>;

(* Context opacity is per panel, not per style: it has to fall as the
   number of overlapping cells rises.  0.25 suits panel (b), a close-up
   of 36 cells; at panel (a)'s 240 cells the same value saturates the
   lower half into a single white mass, and 0.13 restores it. *)
opacityFor = <|"a" -> 0.13, "b" -> 0.25|>;

render[unf3D_, hits_, context_, sty_, op_, opts___] :=
 Module[{hotSet, s = style[sty]},
  hotSet = Union @@ hits;
  Graphics3D[{
    If[s["ctxEdge"] === None, EdgeForm[None],
     EdgeForm[{Thickness[0.0008], s["ctxEdge"]}]],
    FaceForm[If[sty === "black", GrayLevel[s["ctxGray"], op], s["ctxFace"]]],
    Table[cellPolys[unf3D, i], {i, Complement[context, hotSet]}],
    If[s["hotEdge"] === None, EdgeForm[None],
     EdgeForm[{Thickness[0.0035], s["hotEdge"]}]],
    Table[{FaceForm[s["hotA"]], cellPolys[unf3D, p[[1]]],
           FaceForm[s["hotB"]], cellPolys[unf3D, p[[2]]]}, {p, hits}]},
   Boxed -> False, Lighting -> "Neutral", ImageSize -> 520,
   Background -> s["bg"], opts]];

(* ---- panel (a): runcinated 24-cell, best root ------------------- *)
rec24 = SelectFirst[w4BFSAll, #["sym"] === "x3o4o3x" &];
best24 = First@SortBy[rec24["results"], #[[3]] &];
Print["panel (a): x3o4o3x  root ", best24[[1]], "  overlaps ", best24[[3]],
      "  (min over its ", Length[rec24["results"]], " roots)"];
setup["x3o4o3x"];
{u24, h24} = netAtRoot[best24[[1]]];
Print["  intersecting pairs: ", h24];
Do[Export[baseDir <> "fig_uniform_a" <> If[sty === "black", "_black", ""] <> ".png",
   render[u24, h24, Range[nCP], sty, opacityFor["a"], ViewPoint -> {2.2, -2.4, 1.6}],
   ImageResolution -> 300], {sty, {"white", "black"}}];

(* ---- panel (b): runcitruncated 120-cell, a failing root --------- *)
recRT = SelectFirst[h4All, #["sym"] === "x5x3o3x" &];
badRT = First@Select[recRT["results"], ! #[[4]] &];
Print["panel (b): x5x3o3x  root ", badRT[[1]], "  overlaps ", badRT[[3]],
      "  (", recRT["valid"], " of ", recRT["cells"], " roots are valid)"];
setup["x5x3o3x"];
{uRT, hRT} = netAtRoot[badRT[[1]]];
Print["  intersecting pair: ", hRT];

(* Close-up: the offending pair plus the cells around it.  Cells are
   picked by centroid distance, which keeps the neighborhood connected,
   and the view range is left Automatic.

   An explicit PlotRange was tried first and was wrong: a net is
   irregular, so a cube centered on the pair always slices whatever
   crosses its faces, and the cells appeared cut off.  Fixing it by
   keeping only the cells that fit entirely inside the box is worse
   still -- the neighborhood loses cells from its middle and reads as
   floating debris rather than as part of a net.  Widening the box to
   contain every picked cell wastes most of the frame, because one
   outlying cell sets the scale.  With no PlotRange the view simply
   fits what is drawn: nothing is clipped and no margin is wasted. *)
Module[{hot, ctr, cellR, near},
  hot = Union @@ hRT;
  ctr = Mean[Join @@ uRT[[hot]]];
  cellR = Median[Table[Max[Norm[# - Mean[uRT[[i]]]] & /@ uRT[[i]]],
                       {i, nCP}]];
  near = Select[Range[nCP], Norm[Mean[uRT[[#]]] - ctr] < 3 cellR &];
  Print["  median cell radius ", Round[cellR, 0.01],
        "   cells drawn: ", Length[near]];
  Do[Export[baseDir <> "fig_uniform_b" <> If[sty === "black", "_black", ""] <> ".png",
    render[uRT, hRT, near, sty, opacityFor["b"],
     ViewPoint -> {1.8, -2.6, 1.4}],
    ImageResolution -> 300], {sty, {"white", "black"}}]];

(* Panel (b) comes out about 17% taller than panel (a), so the two sit
   slightly unevenly side by side.  This is not padding and cannot be
   cropped away: the apparently empty band at the bottom of (b) holds
   translucent cells that are faint but present.  Measured with a
   content mask, both renders are already tight. *)

Print["\nwrote fig_uniform_a.png and fig_uniform_b.png"];
