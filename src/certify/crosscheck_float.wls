(* crosscheck_float.wls
   Unfold the trees written by crosscheck_float.py with the paper's
   floating-point hinge rotations (unfoldTFP) and compare the squared
   distances between cell centroids with the exact ones.

   A case without "parent" is the paper's own tree T_r: it is unfolded by
   bfsUnfoldP with adjListP built exactly as in runPolytopeV2, so a match
   also shows that the Python code rebuilt the same tree.

   wolframscript -file crosscheck_float.wls <json> *)

baseDir = "/Users/yoshino/Library/CloudStorage/Dropbox/260324Peeling4D/";
Module[{src, cut},
  src = Import[baseDir <> "src/bfsnet/face_rotation_net_all4D_v2.m", "Text"];
  cut = First@First@StringPosition[src, "Print[\"Face-rotation BFS net check v2"];
  ToExpression[StringTake[src, {1, cut - 1}]]];

data = Import[Last[$ScriptCommandLine], "RawJSON"];
raw = Get[baseDir <> "_4DData/f" <> data["polytope"] <> ".m"];
versP = N[raw[[1]]]; facesP = raw[[4]]; cellsP = raw[[6]]; nCP = Length[cellsP];
cellV = Table[Union @@ facesP[[cellsP[[i]]]], {i, nCP}];
adjListP = Table[Select[Range[nCP],
    Function[j, j =!= i && Length[Intersection[cellsP[[i]], cellsP[[j]]]] >= 1]], {i, nCP}];

worst = 0.;
Do[
  Module[{r = c["root"], parent, dist, order, tfs, cen, d2},
   If[KeyExistsQ[c, "parent"],
    parent = c["parent"];
    (* depth of each cell along the given parents *)
    dist = Table[Module[{v = i, d = 0}, While[v =!= r, v = parent[[v]]; d++]; d], {i, nCP}];
    order = SortBy[DeleteCases[Range[nCP], r], dist[[#]] &];
    tfs = <|r -> idAff4|>;
    Do[tfs[v] = unfoldTFP[parent[[v]], v, tfs[parent[[v]]]], {v, order}],
    tfs = bfsUnfoldP[r]];
   cen = Table[Mean[applyAff[tfs[i], #] & /@ versP[[cellV[[i]]]]], {i, nCP}];
   d2 = Outer[(#1 - #2).(#1 - #2) &, cen, cen, 1];
   worst = Max[worst, Max@Abs[d2 - c["d2"]]]],
  {c, data["cases"]}];

Print[data["polytope"], "-cell: ", Length[data["cases"]],
  If[KeyExistsQ[First[data["cases"]], "parent"], " trees", " trees T_r via bfsUnfoldP"],
  ", max |d2_float - d2_exact| = ", ToString[worst, CForm]];
