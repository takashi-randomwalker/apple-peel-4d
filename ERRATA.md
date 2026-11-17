# Errata (2026-11-16)

Several results on **face-rotation BFS nets** that this repository published
were wrong. This file says what was wrong, why, and what is now known. The
apple-peel results (`STLs_v4/`, `src/peel4d/`) are not affected: their
intersection test (`src/peel4d/unfold3DExport.m`) checks every pair of cells
and has neither bug.

## The two bugs

Both were in the floating-point intersection test of the face-rotation BFS
nets, and they act in opposite directions.

1. **Adjacent cells were exempt.** Every pair of cells sharing a face in the
   polytope was skipped (`src/bfsnet/bfs_sat.m`, and the volume test in
   `src/bfsnet/face_rotation_net_all4D_v2.m`). Only pairs joined by an edge of
   the tree are glued in the net; the others can overlap, and in the 120-cell
   and in many uniform polytopes they do. *Effect: real overlaps were missed.*
2. **Edge directions were dropped.** The separating-axis test deduplicated edge
   directions by a sign taken from a zero component, so every edge with a zero
   coordinate collapsed onto one key and all but one were discarded
   (`bfs_sat.m`). Missing edge-cross-edge axes made some disjoint pairs look
   overlapping. *Effect: a few spurious overlaps.*

Both are fixed (commit `37524b3`). The fixed test reproduces, tree by tree,
the overlap counts of an exact certifier (`src/certify/`) on all 4,448 BFS
trees of the 120-cell and 1,200 of the 600-cell that were checked.

## What was wrong, and what is true

**Regular 4-polytopes.** Now certified in exact arithmetic (rational, or
$\mathbb{Q}(\sqrt5)$), with every pair of cells checked; see
[`src/certify/README.md`](src/certify/README.md).

| claim as published | status |
|---|---|
| The tree $\mathcal{T}_r$ gives a net from every root, for all six | **true** (exact) |
| Every BFS tree of the 5-, 8-, 16- and 24-cell gives a net | **true** (exact) |
| 20 random BFS trees per root of the 120- and 600-cell all give nets ("1,132,661 trees, none invalid") | **false for the 120-cell**: 2,238 of its 2,400 sampled trees overlap. True for the 600-cell. |
| Regular polytopes are insensitive to the choice of BFS tree | **false for the 120-cell** |
| (new) Every BFS tree of the 600-cell gives a net | **true** (exact, all pairs of root paths) |
| (new) A BFS tree of the 120-cell gives a net iff it contains no same-layer 7-cycle around a vertex whose chord branch passes straight through its middle cell | **true** (exact, all pairs of root paths) |

**Convex uniform 4-polytopes.** Recomputed with the fixed test, for
$\mathcal{T}_r$ from every root (`data/uniform_recheck/`). Floating point; the
exact certifier does not yet handle non-regular cells.

| | published | corrected |
|---|---:|---:|
| valid roots (of 24,487) | 21,584 | **9,872** |
| polytopes valid from every root | 61 | **51** |
| from some roots only | 2 | **5** |
| from no root | 1 | **8** |

The 13 polytopes that changed, and the 16 self-intersecting files among the
published nets, are listed in [`uniform_nets/README.md`](uniform_nets/README.md).
No root changed from invalid to valid.

**Tie-breaking rules and random trees, recomputed (2026-11-17).** The
max-index and min-index rules from every root, and 10 random BFS trees from
root 1, for all 64 polytopes
([`data/uniform_recheck/rules_compare.txt`](data/uniform_recheck/rules_compare.txt)):

| rule | from every root | some roots | no root | valid roots |
|---|---:|---:|---:|---:|
| queue ($\mathcal{T}_r$) | 51 | 5 | 8 | 9,872 |
| max-index | 49 | 10 | 5 | 10,901 |
| min-index | 49 | 9 | 6 | 10,297 |

48 polytopes are valid from every root under all three rules; 38 of them
also give a net in all 10 random trees, and 20 polytopes in none of them.
The examples published before change as follows (valid roots, queue /
max-index / min-index):

| polytope | published | corrected |
|---|---|---|
| grand antiprism | 320 / 320 / 0 | 320 / 256 / 0 |
| runcinated 24-cell `x3o4o3x` | 0 / 0 / 0 | 0 / 26 / 26 |
| runcitruncated 600-cell `x5o3x3x` | 47 / 2,635 / 2,627 | 0 / 720 / 0 |
| runcitruncated 120-cell `x5x3o3x` | 2,570 / 2,640 / 2,151 | 0 / 0 / 110 |

So "the grand antiprism is valid from all 320 roots under two rules and from
none under the third" is false: two rules disagree on it entirely, the third
gives 256. The runcinated 24-cell is not invalid under every rule.

**Not yet recomputed, and to be treated as void:** the symmetric-tree
experiments (e.g. "70 of 600 symmetric trees of the runcinated 24-cell are
valid").

## Where the wrong claims appeared

- `face_rotation_net_paper.tex` / `.pdf`, the draft submitted to JoCG
  (tag `jocg-submitted-2026-10-12`). **Removed from the repository** until it is
  revised; it remains in the git history.
- `jcdcg3_slides.tex` / `.pdf`, the JCDCG³ 2026 talk (slides 8 and 9 and the
  backup slides). Kept as the record of the talk. The JCDCG³ 2026 abstract
  itself claims only that $\mathcal{T}_r$ gives a net from every root,
  which is true.
- `uniform_nets/manifest.csv`, `uniform_nets/allroots.csv` and
  `uniform_nets/README.md`. **Corrected.** The STL and OFF nets themselves were
  right and are unchanged; only which of them self-intersect has changed.
