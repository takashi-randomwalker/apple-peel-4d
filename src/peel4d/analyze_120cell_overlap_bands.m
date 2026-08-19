(* ================================================================
   analyze_120cell_overlap_bands.m
   120胞体展開図の重なりが同一帯内 vs 異帯間かを解析
   C1=3 の複数 C2 について帯ペア統計を集計する
   ================================================================ *)

scriptDir = If[ValueQ[$InputFileName], DirectoryName[$InputFileName], NotebookDirectory[]];
Get[scriptDir <> "unfold3DExport.m"];

dataDir = scriptDir <> "_4DData/";
raw   = Get[dataDir <> "f120.m"];
vers  = N[raw[[1]]]; faces = raw[[4]]; cells = raw[[6]];
nc    = Length[cells];

Get[scriptDir <> "ans4DGlobal_v4.mx"];

cellVerts[ci_] := Union @@ faces[[ cells[[ci]] ]];

(* ----------------------------------------------------------------
   C1=3 の全成功ペアを取得（order 重複除去）
   ---------------------------------------------------------------- *)
fixedC1 = 3;
allRes = results4DGlobal[120]["results"][[fixedC1]];
succRes = DeleteDuplicatesBy[Select[allRes, Last[#] &], #[[2]] &];
Print["C1=3 の成功ペア数（order unique）: ", Length[succRes]];

(* ----------------------------------------------------------------
   Face-up 後のセル重心 w を代表ケースから取得しバンド割り当てを作成
   （全ペアで localVers は同一方向に face-up されているので1件で代表可）
   ---------------------------------------------------------------- *)
repLV = succRes[[1, 1]];  (* localVers（face-up 済み頂点座標） *)
cC    = Mean[repLV[[ cellVerts[#] ]]] & /@ Range[nc];
wAll  = cC[[#, 4]] & /@ Range[nc];

(* w 値を 0.01 刻みで丸めてバンドグループを作成 *)
wRounded   = Round[wAll, 0.01];
wLevels    = Sort[DeleteDuplicates[wRounded], Greater];
bandOf     = Association @ Table[
  c -> Position[wLevels, wRounded[[c]]][[1, 1]], {c, nc}];

Print["帯数: ", Length[wLevels]];
Print["帯構造 (帯番号: w値 → セル数):"];
Do[
  Print["  帯", b, ": w=", wLevels[[b]],
        "  (", Length[Select[Range[nc], bandOf[#] == b &]], " cells)"],
  {b, Length[wLevels]}];

(* ----------------------------------------------------------------
   cellLocal3Ds を一度だけ計算
   ---------------------------------------------------------------- *)
Print["\nPre-computing cellLocal3Ds ..."];
t0 = AbsoluteTime[];
cellLocal3Ds = Table[unfoldCellLocal3D[ci, vers, cells, faces], {ci, nc}];
Print["  Done: ", Round[AbsoluteTime[] - t0, 0.1], " s"];

(* ----------------------------------------------------------------
   各ペアについて重なりペアを帯分類
   ---------------------------------------------------------------- *)
Print["\n== 各 C2 の重なり帯分析 =="];

allSame  = 0; allCross = 0;
crossTallyAll = <||>;

Do[
  res = succRes[[idx]];
  ord = res[[2]];
  c2  = ord[[2]];

  embs     = unfoldTo3DFast[cells, faces, ord, cellLocal3Ds];
  ovPairs  = checkOverlaps[ord, cells, faces, embs];
  nOv      = Length[ovPairs];

  If[nOv == 0,
    Print["C2=", c2, "  重なりなし"],
    (* 帯分類 *)
    same  = Count[ovPairs, p_ /; bandOf[ord[[p[[1]]]]] == bandOf[ord[[p[[2]]]]]];
    cross = nOv - same;
    allSame  += same;
    allCross += cross;

    (* 帯ペアの Tally *)
    crossPairs = Select[ovPairs, bandOf[ord[[#[[1]]]]] != bandOf[ord[[#[[2]]]]] &];
    bpTally = Tally[Sort[{bandOf[ord[[#[[1]]]]], bandOf[ord[[#[[2]]]]]}] & /@ crossPairs];
    Do[
      key = bpTally[[i, 1]];
      crossTallyAll[key] = Lookup[crossTallyAll, key, 0] + bpTally[[i, 2]],
      {i, Length[bpTally]}];

    Print["C2=", c2,
          "  total=", nOv,
          "  同一帯=", same, " (", Round[100. same/nOv, 1], "%)",
          "  異帯=",   cross, " (", Round[100. cross/nOv, 1], "%)"]
  ],
  {idx, Length[succRes]}];

(* ----------------------------------------------------------------
   集計サマリ
   ---------------------------------------------------------------- *)
total = allSame + allCross;
Print["\n== 全 C2 合計 =="];
Print["同一帯内の重なり: ", allSame,  " ペア (", Round[100. allSame/total,  0.1], "%)"];
Print["異なる帯間の重なり: ", allCross, " ペア (", Round[100. allCross/total, 0.1], "%)"];

Print["\n異帯重なりの帯ペア組み合わせ（上位15）:"];
sortedCP = SortBy[Normal[crossTallyAll], -#[[2]] &];
Do[
  Print["  帯", sortedCP[[i, 1, 1]], " ↔ 帯", sortedCP[[i, 1, 2]],
        "  : ", sortedCP[[i, 2]], " ペア"],
  {i, Min[15, Length[sortedCP]]}];
