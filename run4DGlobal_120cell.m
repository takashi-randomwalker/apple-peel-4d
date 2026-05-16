(* ================================================================
   run4DGlobal_120cell.m
   ================================================================
   120胞体を現行 peeling4Df4.m で計算し，
   ans4DGlobal_v4.mx に追記して上書き保存する．
   ================================================================ *)

scriptDir = If[ValueQ[$InputFileName],
               DirectoryName[$InputFileName],
               NotebookDirectory[]];
Get[scriptDir <> "peeling4Df4.m"];
dataDir = scriptDir <> "_4DData/";

(* v4 を読み込む *)
v4File = scriptDir <> "ans4DGlobal_v4.mx";
If[FileExistsQ[v4File],
  Get[v4File]; Print["Loaded ans4DGlobal_v4.mx"],
  Print["ERROR: ans4DGlobal_v4.mx not found."]; Abort[]
];

(* 120胞体データ *)
raw   = Get[dataDir <> "f120.m"];
vers  = raw[[1]]; faces = raw[[4]]; cells = raw[[6]];
nc    = Length[cells];

Print[""];
Print[StringRepeat["=", 60]];
Print["  Apple-Peel 4D: 120-cell (", nc, " cells)"];
Print[StringRepeat["=", 60]];
Print["  (Previous: Possible -- 513/1440)"];

t0 = AbsoluteTime[];

ans = Table[
  If[Mod[top, 10] == 0, Print["  top = ", top, " / ", nc]];
  peeling4Df4[N[vers], faces, cells, top],
  {top, nc}
];

dt         = Round[AbsoluteTime[] - t0, 0.1];
allBools   = Flatten[Map[Last, ans, {2}]];
trueCount  = Count[allBools, True];
falseCount = Count[allBools, False];
uniqOrders = DeleteDuplicates[
  Select[Flatten[ans, 1], Last[#]&][[All, 2]]];
classif    = Which[falseCount == 0, "Perfect",
                   trueCount  == 0, "Impossible",
                   True,            "Possible"];

Print[""];
Print["  ", dt, " s  |  ", classif,
      "  (True: ", trueCount,
      ", False: ", falseCount,
      ", Unique: ", Length[uniqOrders], ")"];

results4DGlobal[120] = <|
  "name"           -> "120-cell",
  "n"              -> 120,
  "classification" -> classif,
  "trueCount"      -> trueCount,
  "falseCount"     -> falseCount,
  "uniqueCount"    -> Length[uniqOrders],
  "results"        -> ans
|>;

(* v4 に上書き保存 *)
DumpSave[v4File, results4DGlobal];
Print["Saved: ", v4File];
Print["Done."];
