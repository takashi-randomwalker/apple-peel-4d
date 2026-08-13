(* =====================================================================
   run_wythoff_bfs.m
   Feed the Wythoff-generated uniform 4-polytopes (uniform4D_wythoff.m)
   into the face-rotation BFS net checker (face_rotation_net_all4D_v2.m).

   Phase A : all roots for the regular members that also exist as
             _4DData/f*.m files.  Reproducing ALL VALID there validates
             the whole data path (generator -> BFS -> volume test).
   Phase B : root 1 only, for every generated polytope, smallest first.

   Only the definitions of face_rotation_net_all4D_v2.m are loaded; its
   main loop (which re-runs the 120- and 600-cell) is cut off.
   ===================================================================== *)

baseDir = "/Users/yoshino/Library/CloudStorage/Dropbox/260324Peeling4D/";
SetDirectory[baseDir];

(* --- 1. generated polytopes ------------------------------------- *)
Get[baseDir <> "uniform4D_wythoff.m"];

(* keep one representative per combinatorial signature *)
w4Distinct = DeleteDuplicatesBy[w4All, {#["counts"], #["hist"]} &];
w4Distinct = SortBy[w4Distinct, #["counts"][[4]] &];
Print["\n", Length[w4Distinct], " distinct polytopes to feed"];

(* --- 2. BFS machinery, definitions only -------------------------- *)
Module[{src, cut},
  src = Import[baseDir <> "face_rotation_net_all4D_v2.m", "Text"];
  cut = First@First@StringPosition[src, "Print[\"Face-rotation BFS net check v2"];
  ToExpression[StringTake[src, {1, cut - 1}]]];
Print["BFS definitions loaded (epsVol = ", epsVol, ")\n"];

(* --- 3. driver ---------------------------------------------------- *)
runPolytopeData[name_, data_, roots_, cap_] := Module[
  {res, nValid, t, nTried},
  versP  = N[data[[1]], MachinePrecision];
  facesP = data[[4]];
  cellsP = data[[6]];
  nCP    = Length[cellsP];
  versCenteredP = (# - Mean[versP]) & /@ versP;

  adjListP = Table[
    Select[Range[nCP],
      Function[j, j =!= i &&
        Length[Intersection[cellsP[[i]], cellsP[[j]]]] >= 1]],
    {i, nCP}];
  adjSetP = Association@Flatten[
    Table[Sort[{i, j}] -> True, {i, nCP}, {j, adjListP[[i]]}], 1];

  {t, res} = AbsoluteTiming@TimeConstrained[
     Table[checkRootV2[r], {r, roots}], cap, $Aborted];

  If[res === $Aborted,
   Print[StringPadRight[name, 10],
     StringPadLeft[ToString[nCP], 6], " cells  ",
     StringPadLeft[ToString[Length[roots]], 5], " roots   TIMEOUT after ",
     Round[t], "s"];
   Return[<|"name" -> name, "cells" -> nCP, "status" -> "timeout"|>]];

  nTried = Length[res];
  nValid = Count[res, {_, _, _, True}];
  Print[StringPadRight[name, 10],
   StringPadLeft[ToString[nCP], 6], " cells  ",
   StringPadLeft[ToString[nTried], 5], " roots  valid ",
   StringPadLeft[ToString[nValid], 5], "/", nTried, "  ",
   StringPadRight[
     Which[nValid === nTried, "ALL VALID", nValid === 0, "ALL INVALID",
           True, "MIXED"], 12],
   "deg ", Min[Length /@ adjListP], "-", Max[Length /@ adjListP],
   "   ", Round[t, 0.1], "s"];
  <|"name" -> name, "cells" -> nCP, "tried" -> nTried,
    "valid" -> nValid, "status" -> "done"|>];

(* --- 4. Phase A : regular members, all roots ---------------------- *)
Print[StringRepeat["=", 78]];
Print["Phase A : regular polytopes from the generator, ALL roots"];
Print["          (expected: ALL VALID, matching the paper)"];
Print[StringRepeat["=", 78]];

phaseA = Table[
  Module[{r = w4Find[p[[2]]]},
   runPolytopeData[p[[1]], r["data"], Range[r["counts"][[4]]], 900]],
  {p, {{"5-cell", "x3o3o3o"}, {"8-cell", "x4o3o3o"},
       {"16-cell", "o4o3o3x"}, {"24-cell", "x3o4o3o"}}}];

(* --- 5. Phase B : every generated polytope, root 1 ---------------- *)
Print["\n", StringRepeat["=", 78]];
Print["Phase B : all generated uniform polytopes, root 1 only"];
Print[StringRepeat["=", 78]];

phaseB = Table[
  runPolytopeData[r["sym"], r["data"], {1}, 900],
  {r, w4Distinct}];

(* --- 6. summary --------------------------------------------------- *)
Print["\n", StringRepeat["=", 78]];
Print["SUMMARY (Phase B, root 1)"];
Print[StringRepeat["=", 78]];
Print["  done    : ", Count[phaseB, KeyValuePattern["status" -> "done"]],
      " / ", Length[phaseB]];
Print["  valid   : ", Count[phaseB, KeyValuePattern[{"status" -> "done",
                                                     "valid" -> 1}]]];
Print["  invalid : ", Count[phaseB, KeyValuePattern[{"status" -> "done",
                                                     "valid" -> 0}]]];
Print["  timeout : ", Count[phaseB, KeyValuePattern["status" -> "timeout"]]];
Print["\n  invalid at root 1: ",
  Select[phaseB, #["status"] === "done" && #["valid"] === 0 &][[All, "name"]]];
