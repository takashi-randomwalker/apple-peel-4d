(* =====================================================================
   make_talk_figures.m
   Slide images for the JCDCG^3 talk: candidate "valid net" pictures to
   pair with the two failure panels already made for the paper.

   Candidates are chosen for legibility rather than size -- a net of
   2640 small cells reads as a cloud, while a few hundred large cells
   reads as a structure.
   ===================================================================== *)

baseDir = "/Users/yoshino/Library/CloudStorage/Dropbox/260324Peeling4D/";
SetDirectory[baseDir];

Get[baseDir <> "wythoff_gen_all.mx"];
Get[baseDir <> "uniform4D_special.m"];
Module[{src, cut},
  src = Import[baseDir <> "face_rotation_net_all4D_v2.m", "Text"];
  cut = First@First@StringPosition[src, "Print[\"Face-rotation BFS net check v2"];
  ToExpression[StringTake[src, {1, cut - 1}]]];
Get[baseDir <> "bfs_sat.m"];

setup[rec_] := Module[{data, faceToCells, pairs},
  data   = rec["data"];
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

netAt[r_] := Module[{tfsR},
  versP = satFaceUp[r, versCenteredP];
  tfsR  = satBfsUnfold[r];
  Table[Module[{tf = tfsR[i]},
     Take[#, 3] & /@ (versP[[satCellIdx[[i]]]] . Transpose[tf[[1]]] +
                      Threaded[tf[[2]]])], {i, nCP}]];

cellPolys[unf_, i_] := Module[{loc},
  loc = AssociationThread[satCellIdx[[i]] -> Range[Length[satCellIdx[[i]]]]];
  Polygon[unf[[i]][[Lookup[loc, facesP[[#]]]]]] & /@ cellsP[[i]]];

(* colour by BFS depth from the root, so the layer structure reads *)
depths[root_] := Module[{d = <|root -> 0|>, q = {root}, head = 1, p},
  While[head <= Length[q],
   p = q[[head]]; head++;
   Do[If[! KeyExistsQ[d, nb], d[nb] = d[p] + 1; AppendTo[q, nb]],
    {nb, adjListP[[p]]}]];
  Lookup[d, Range[nCP], 0]];

renderNet[unf_, root_, file_, size_] := Module[{dp, mx},
  dp = depths[root]; mx = Max[dp];
  Export[baseDir <> file,
   Graphics3D[
    Table[{EdgeForm[{Thickness[0.0006], GrayLevel[0.25, 0.55]}],
      FaceForm[Blend[{RGBColor[0.16, 0.36, 0.72], RGBColor[0.92, 0.85, 0.25],
                      RGBColor[0.80, 0.22, 0.16]}, dp[[i]]/mx]],
      cellPolys[unf, i]}, {i, nCP}],
    Boxed -> False, Lighting -> "Neutral", ImageSize -> size,
    ViewPoint -> {2.4, -2.2, 1.5}],
   ImageResolution -> 200];
  {file, FileByteCount[baseDir <> file]/1024.^2, mx}];

cands = {
  {"snub 24-cell",            "s3s4o3o", sp4All},
  {"grand antiprism",         "gap",     sp4All},
  {"truncated 120-cell",      "x5x3o3o", w4All},
  {"bitruncated 120/600-cell","o5x3x3o", w4All}};

Print[StringPadRight["candidate", 26], StringPadLeft["cells", 6],
      StringPadLeft["depth", 7], StringPadLeft["MB", 8]];
Do[Module[{nm = c[[1]], sym = c[[2]], src = c[[3]], rec, unf, r},
   rec = SelectFirst[src, #["sym"] === sym &];
   setup[rec];
   unf = netAt[1];
   r = renderNet[unf, 1, "talk_" <> StringReplace[sym, " " -> "_"] <> ".png", 900];
   Print[StringPadRight[nm, 26], StringPadLeft[ToString[nCP], 6],
         StringPadLeft[ToString[r[[3]]], 7],
         StringPadLeft[ToString[Round[r[[2]], 0.01]], 8]]],
 {c, cands}];

Print["\ncolour = BFS depth from the root (blue = root, red = last layer)"];
