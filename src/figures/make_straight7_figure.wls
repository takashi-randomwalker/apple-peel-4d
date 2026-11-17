(* make_straight7_figure.wls -- Figure 1(b) of bfs_nets_regular.tex.

   A straight 7-cycle in a net of the 120-cell, rooted at cell 0 (= L).
   data/figures/net120_flip.json holds two BFS nets that differ only in the
   parent of cell 13 (a): in one it hangs from 96 (p_a), straight through
   the chain 9 -> 96 -> 13 (x -> p_a -> a), and meets cell 21 (b); in the
   other it hangs from 6, along the bent route 0 -> 1 -> 6, and is clear of b.
   "cells" are the second net's placements, "movedBad" is a in the first.
   The data were extracted from the 2026-11-15 interactive viewer.

   Drawn: L as a gray outline; x, p_a, a yellow; b blue; the bent route
   1, 6 as green outlines and a' (a in the second net) solid green.  The
   view is the viewer's default direction turned 20 degrees about the line
   through the centroids of x and a, which sends b back.

   wolframscript -file src/figures/make_straight7_figure.wls   (about 1 min)
*)
base = ParentDirectory[DirectoryName[$InputFileName], 2];
net = Import[FileNameJoin[{base, "data", "figures", "net120_flip.json"}], "RawJSON"];
cell[i_] := net["cells"][[i + 1]];                 (* 0-based cell index *)
pts[c_] := Partition[c["v"], 3];
cen[c_] := Mean[pts[c]];
tri[c_] := Polygon[pts[c][[# + 1]] & /@ Partition[c["t"], 3]];
(* the edges of a dodecahedron are its shortest vertex pairs *)
edges[p_] := Module[{d = Min[EuclideanDistance @@@ Subsets[p, {2}]]},
   Line /@ Select[Subsets[p, {2}], Abs[EuclideanDistance @@ # - d] < 10^-3 &]];
solid[c_, col_, op_] := {{EdgeForm[None], FaceForm[Directive[col, Opacity[op]]], tri[c]},
   {Darker[col, 0.5], AbsoluteThickness[1.1], edges[pts[c]]}};
frame[c_, col_] := {col, AbsoluteThickness[1.6], edges[pts[c]]};
yellow = RGBColor["#f2c14e"]; blue = RGBColor["#6ea8ff"]; green = RGBColor["#4fc7a4"];
label[t_, p_, col_: Black] := Text[Style[t, col, Bold, FontSize -> 30,
    FontFamily -> "Times New Roman", FontSlant -> Italic], p];

{L, x, pa, b, aGood} = cell /@ {0, 9, 96, 21, 13};
aBad = net["movedBad"];
axis = Normalize[cen[aBad] - cen[x]];
view = 3 RotationMatrix[20 Degree, axis] . Normalize[{-0.55, 0.45, 0.7}];

g = Graphics3D[{
    frame[L, GrayLevel[0.78]], frame[cell[#], Lighter[green, 0.15]] & /@ {1, 6},
    solid[x, yellow, 0.6], solid[pa, yellow, 0.6], solid[b, blue, 0.8],
    solid[aGood, green, 0.7], solid[aBad, yellow, 0.75],
    label["x", cen[x]], label[Subscript["p", "a"], cen[pa]], label["b", cen[b]],
    label["a", cen[aBad]],
    label[Row[{"a", Style["\[Prime]", FontFamily -> "Arial", FontWeight -> Plain, FontSlant -> Plain]}], cen[aGood]],
    label["L", cen[L] + 0.55 Normalize[cen[L] - cen[x]], GrayLevel[0.85]]},
   Boxed -> False, Background -> GrayLevel[0.13], Lighting -> "Neutral",
   ViewPoint -> view, ViewVertical -> {0, 1, 0}, ImageSize -> 700];

(* trim: drop the rasterizer's 1-pixel frame, cut the background, pad evenly *)
img = ImagePad[Rasterize[g, ImageResolution -> 200], -6];
img = ImagePad[ImagePad[img, -BorderDimensions[img, 0.03]], 40, GrayLevel[0.13]];
Export[FileNameJoin[{base, "fig_straight7_net.png"}], img];
Print["fig_straight7_net.png ", ImageDimensions[img]];
