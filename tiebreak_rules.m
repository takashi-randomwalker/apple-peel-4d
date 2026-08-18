(* =====================================================================
   tiebreak_rules.m
   Why does index order do so well on the two MIXED uniform polytopes,
   when uniformly random BFS trees almost never work?

   Measured fact (idx2 experiment): the cell index is not arbitrary.
   The generator emits cells in the order its group closure discovers
   the group elements, i.e. by word length from the identity, so index
   correlates with distance from the prototype cell and, after face-up,
   anti-correlates with the w-coordinate: corr(index, w) = -0.32 for
   these two polytopes and -0.95 for the 120-cell.  Choosing the
   smallest-index parent therefore leans toward the parent with larger
   w, which is what the apple-peel RZ rule does deliberately.

   Two hypotheses:
     A (direction)  what helps is the w-bias.  Then an explicit max-w
                    rule should do at least as well, and max-index or
                    min-w should do badly.
     B (coherence)  what helps is that all parents are chosen by one
                    global order, unlike uniform random trees where each
                    cell picks independently.  Then a random global
                    order should also do well.

   Five rules on a common sample of roots separate these.
   ===================================================================== *)

baseDir = "/Users/yoshino/Library/CloudStorage/Dropbox/260324Peeling4D/";
SetDirectory[baseDir];
logFile = baseDir <> "tiebreak_rules_progress.txt";
SeedRandom[20260817];

nRoots = 300;   (* common sample of roots per polytope *)

Get[baseDir <> "wythoff_gen_all.mx"];
Get[baseDir <> "h4_bfs_allroots.mx"];
Module[{src, cut},
  src = Import[baseDir <> "face_rotation_net_all4D_v2.m", "Text"];
  cut = First@First@StringPosition[src, "Print[\"Face-rotation BFS net check v2"];
  ToExpression[StringTake[src, {1, cut - 1}]]];
Get[baseDir <> "bfs_sat.m"];

say[s__] := (Print[s];
  PutAppend[OutputForm[StringJoin[ToString /@ {s}]], logFile]);

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

layersAndOptions[r_] := Module[{dist, q, head, p},
  dist = ConstantArray[-1, nCP]; dist[[r]] = 0;
  q = {r}; head = 1;
  While[head <= Length[q],
   p = q[[head]]; head++;
   Do[If[dist[[nb]] < 0, dist[[nb]] = dist[[p]] + 1; AppendTo[q, nb]],
    {nb, adjListP[[p]]}]];
  {dist, Table[If[i === r, {},
     Select[adjListP[[i]], dist[[#]] === dist[[i]] - 1 &]], {i, nCP}]}];

(* Parent-selection rules.  cw = cell centroid w after face-up.

   Note on "queue": this is the rule the rest of the paper uses, and it
   is NOT the same as min-index.  With a queue, the parent of v is
   whichever neighbor reaches v first, which depends on the order the
   queue was built in; min-index picks the smallest-indexed neighbor one
   layer closer, regardless of when it was reached.  Comparing the two
   is the point of including both.

   Note on ties: MinimalBy takes ONE criterion.  Passing a list of them
   silently compares unevaluated expressions -- an earlier version did
   that and the w rules were quietly reduced to index order.  Ties are
   broken here by adding a small index term to the key instead. *)
pickParents[rule_, opts_, rest_, cw_, r_] := Module[{perm, vis, q, head, par, p},
  Switch[rule,
   "min-index", Table[Min[opts[[v]]], {v, rest}],
   "max-index", Table[Max[opts[[v]]], {v, rest}],
   "max-w",  Table[First@MinimalBy[opts[[v]], -cw[[#]] + 10^-9 # &], {v, rest}],
   "min-w",  Table[First@MinimalBy[opts[[v]],  cw[[#]] + 10^-9 # &], {v, rest}],
   "rand-order",
     perm = Ordering@RandomSample[Range[nCP]];
     Table[First@MinimalBy[opts[[v]], perm[[#]] &], {v, rest}],
   "queue",
     par = ConstantArray[0, nCP];
     vis = <|r -> True|>; q = {r}; head = 1;
     While[head <= Length[q],
      p = q[[head]]; head++;
      Do[If[! KeyExistsQ[vis, nb],
         vis[nb] = True; AppendTo[q, nb]; par[[nb]] = p],
       {nb, adjListP[[p]]}]];
     par[[rest]]]];

validQ[r_, parent_, dist_] := Module[
  {tfsR, unf3D, bMin, bMax, cent, radius, cand, ne, hits, treeSet, order},
  order = SortBy[Complement[Range[nCP], {r}], dist[[#]] &];
  tfsR = <|r -> idAff4|>;
  Do[tfsR[v] = unfoldTFP[parent[[v]], v, tfsR[parent[[v]]]], {v, order}];
  unf3D = Table[
    Module[{tf = tfsR[i]},
     Take[#, 3] & /@ (versP[[satCellIdx[[i]]]] . Transpose[tf[[1]]] +
                      Threaded[tf[[2]]])], {i, nCP}];
  treeSet = Association@Table[Sort[{v, parent[[v]]}] -> True, {v, order}];
  bMin = (Min /@ Transpose[#]) & /@ unf3D;
  bMax = (Max /@ Transpose[#]) & /@ unf3D;
  cent = Mean /@ unf3D;
  radius = 2 Max[MapThread[Function[{p, c}, Max[Norm[# - c] & /@ p]],
                           {unf3D, cent}]];
  cand = satBBoxPairs[cent, bMin, bMax, radius];
  cand = Select[cand, ! KeyExistsQ[treeSet, #] &];
  If[Length[cand] === 0, Return[True]];
  ne = Association@Table[
     c -> satNormalsEdges[unf3D[[c]], satFaceTri[[c]], satEdgePair[[c]]],
     {c, Union @@ cand}];
  hits = Select[cand, Function[p,
     satOverlapQ[unf3D[[p[[1]]]], ne[p[[1]]][[1]], ne[p[[1]]][[2]],
                 unf3D[[p[[2]]]], ne[p[[2]]][[1]], ne[p[[2]]][[2]]]]];
  Length[hits] === 0];

rules = {"queue", "min-index", "max-index", "max-w", "min-w", "rand-order"};

runPolytope[sym_, label_] := Module[{roots, res, t, full},
  setup[sym];
  full = SelectFirst[h4All, #["sym"] === sym &];
  roots = Sort@RandomSample[Range[nCP], nRoots];
  say["\n", StringRepeat["-", 70]];
  say[label, "  ", sym, "   ", nCP, " cells"];
  say["  full index-order result over all roots: ",
      full["valid"], "/", full["cells"],
      " (", Round[100. full["valid"]/full["cells"], 0.1], "%)"];
  say["  sample of ", nRoots, " roots:"];
  {t, res} = AbsoluteTiming@Association@Table[
    Module[{ok},
     ok = Total@Table[
        Module[{dist, opts, rest, par, parent},
         {dist, opts} = layersAndOptions[r];
         rest = Complement[Range[nCP], {r}];
         versP = satFaceUp[r, versCenteredP];
         par = pickParents[rule, opts, rest,
                Table[Mean[versP[[satCellIdx[[i]]]]][[4]], {i, nCP}], r];
         parent = ConstantArray[0, nCP];
         Do[parent[[rest[[k]]]] = par[[k]], {k, Length[rest]}];
         If[validQ[r, parent, dist], 1, 0]],
        {r, roots}];
     say["    ", StringPadRight[rule, 12], StringPadLeft[ToString[ok], 5],
         "/", nRoots, "  (", Round[100. ok/nRoots, 0.1], "%)"];
     rule -> ok],
    {rule, rules}];
  say["  time ", Round[t, 0.1], " s"];
  <|"sym" -> sym, "label" -> label, "roots" -> roots, "res" -> res|>];

say[StringRepeat["=", 70]];
say["Which property of index order makes it work?"];
say["direction (w-bias) vs coherence (one global order)"];
say[StringRepeat["=", 70]];

tbRules = {
  runPolytope["x5o3x3x", "runcitruncated 600-cell"],
  runPolytope["x5x3o3x", "runcitruncated 120-cell"]};

DumpSave[baseDir <> "tiebreak_rules.mx", tbRules];
say["\nsaved to tiebreak_rules.mx"];
