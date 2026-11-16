# Exact certification of the face-rotation BFS nets

Everything here is exact: rational arithmetic for the 5-, 8-, 16- and
24-cell and arithmetic in Q(√5) for the 120- and 600-cell.  No floating
point is used to decide anything, so no tolerance (such as the volume
threshold of the original Mathematica check) is involved.  Plain Python 3,
standard library only.

## What is proved

1. **Every BFS spanning tree** of the 5-, 8-, 16- and 24-cell, from every
   root, gives a net whose cells have pairwise disjoint interiors
   (5 + 48 + 331,776 + 786,432 trees).
2. For all six regular 4-polytopes, the tree T_r of the paper (queue BFS,
   neighbours in ascending order) gives such a net from every root.

3. **Given BFS trees** (`--trees=`), each checked to be a BFS tree first:

   | trees | invalid | overlapping pairs | time |
   |---|---:|---:|---:|
   | 120-cell, 20 random per root (2,400) | **2,238** | 6,585 | 3.7 h |
   | 120-cell, D_5-invariant from root 0 (2,048) | **1,280** | 15,360 | 3.9 h |
   | 600-cell, 20 random per root (12,000) | **0** | 0 | 11.2 h |

   These are the random trees of `src/tiebreak/tiebreak_regular.m`, which
   the floating-point check reported all valid: it exempted every pair of
   cells adjacent in P, and in the 120-cell the overlaps are between such
   pairs.  So most BFS trees of the 120-cell do not give a net, while T_r
   does from every root.

4. **Every BFS tree of the 120- and 600-cell** (`all_bfs.py`).  From root 0
   (all roots, by cell-transitivity) the 388,464 pairs of root paths that
   occur together in some BFS tree are all decided; exactly 240 overlap.
   They are all of one type: a and b adjacent in P and in the same layer,
   their tree paths closing a 7-cycle L,x,pa,a,b,pb,y of G(P) with
   {a,b,pa,pb} the four cells around a vertex.  Exactly one of the chords
   pa~y, pb~x exists; the cycle overlaps iff the branch carrying it
   (x→pa→a when pa~y) passes straight through its middle cell (the two
   faces are opposite).  The other branch is always bent.  Hence a BFS
   tree of the 120-cell gives a net iff it contains no such straight
   7-cycle; even cycles, cycles of length 3, 5, 9 and pairs not adjacent
   in P never overlap.  `straight_loops_120.py` counts these cycles.
   **600-cell**: the 1,757,416 pairs of its 1,889 root paths are all
   disjoint, so every BFS tree of the 600-cell gives a net.

5. **Where BFS overlaps can occur** (`path_classes.py`, observations, not
   claims).  Every pair of cells is placed by its tree path alone, so the
   pairs fall into classes under the symmetry group.  Of these classes
   exactly one overlaps (120-cell), and it is across a *cut ridge* (two
   cells adjacent in P, not joined in the tree).  The holonomy of the
   shortest cycle of a cut ridge is a
   rotation by the angle defect of an edge (24-cell 31.59°, 120-cell
   10.30°, 600-cell 7.36°; the 16-cell has no odd cycles, its shortest
   is 4, 77.88°).  The one overlapping class (120-cell) is a twist of 36°
   about the cut ridge's normal: one vertex of b lies 0.019 inside a.
   In the known non-BFS counterexamples the overlapping cells are not
   adjacent (Seamons and DHR 24-cell, DHR 120-cell) except in the DHR
   600-cell chain (`known_counterexamples.py`,
   `src/bfsnet/check_dhr_chain_120cell.wls`).

Pairs joined by a tree edge are checked too, so the exemption in the
paper's definition of a valid net is not needed.

## Files

| file | role |
|---|---|
| `certify_bfs_exact.py` | the certifier (claim 1 for 5/8/16/24, claim 2 with `--tree` or for 120/600); `--export=` writes a certificate |
| `qsqrt5.py` | exact arithmetic in Q(√5); `python3 -I qsqrt5.py` runs its self-test |
| `verify_certificate.py` | checks a certificate without importing any of the above |
| `test_verifier_mutations.py` | six kinds of tampering; the verifier must reject each |
| `second_checker.py` | an independent re-verification of both claims; shares no code |
| `negative_control.py` | drops the BFS restriction and finds overlapping chains, with a rational point interior to both cells |
| `export_random_trees.wls` | redraws the random trees of `tiebreak_regular.m` (same seed) as parent arrays |
| `export_sector120_trees.wls` | writes the 2,048 D_5-invariant trees of `src/symmetry/sector120_nets.wls` |
| `all_bfs.py` | claim 4: certifies all root-path pairs of the 120- or 600-cell; `classify` sorts them by cycle type |
| `path_classes.py` | the root-path pairs up to symmetry (16/24-cell 36 classes, 120-cell 2,876, 600-cell 37,169) and the holonomy of every cut-ridge class |
| `known_counterexamples.py` | in the overlapping unfoldings of Seamons and of Devadoss–Harvey–Richter, are the overlapping cells adjacent in P? |
| `tree_observations_120.py` | 120-cell: bent-preferring trees are valid; straight branches that overlap miss the root; tree shape does not decide validity |
| `straight_loops_120.py` | counts the straight 7-cycles of a 120-cell tree, which by claim 4 is its number of overlapping pairs |
| `crosscheck_float.py` / `.wls` | the exact placements agree with the paper's floating-point unfolding (`unfoldTFP`, `bfsUnfoldP`) |

Data and logs are in `data/certify/`: `cert_{5,8,16,24}cell.json.gz`
(certificates, 0.2 MB at most), the given trees
`random_trees_{120,600}.json.gz` and `sector_trees_120.json.gz`
(gzip JSON: per root a list of parent arrays, 0-based, −1 at the root; the
sector file has one entry, root 0), `all_bfs_120_overlaps.json.gz` (the 240
overlapping root-path pairs), and the run logs.  `validate_no_exemption.log`
records the float pipeline of `src/bfsnet/bfs_sat.m`, after two fixes,
reproducing the exact overlap count tree by tree.

## How it works

**Unfolding without trigonometry.**  A regular cell is determined by one of
its faces together with the side on which it lies.  The certifier places
the child as the isometry that fixes the shared face and sends the child's
centroid to the reflection of the parent's centroid in the face plane;
the second checker instead uses the symmetry of the polytope that swaps
the two cells and reflects the parent's placement vertex by vertex.  Both
use only field operations, so coordinates stay in Q or Q(√5).

**All BFS trees at once.**  A cell's placement depends only on its path to
the root.  Two root paths occur together in some BFS tree iff they assign
the same parent to every cell they share; certifying every such pair
covers every tree.  The second checker does not use this: it enumerates
the trees one by one, and arrives at the same number of distinct pairs.

**Disjointness.**  The certifier uses separating axes: the interiors of two
convex polytopes A, B are disjoint iff some facet normal u of A − B has
max_A u·x ≤ min_B u·x, and those normals are face normals of A or B or
cross products of an edge of each.  A separating u is the certificate;
if none works the interiors provably meet.  The second checker uses
circumscribed balls and, for the remaining pairs, enumerates the vertices
of A ∩ B exactly: the interiors are disjoint iff they span fewer than
three dimensions.

## Running

```
python3 -I certify_bfs_exact.py 24                 # claim 1, all roots (17 s)
python3 -I certify_bfs_exact.py 600                # claim 2, all roots (31 min, 8 processes)
python3 -I certify_bfs_exact.py 24 --export=cert.json.gz
python3 -I verify_certificate.py cert.json.gz
python3 -I certify_bfs_exact.py 120 --trees=../../data/certify/random_trees_120.json.gz --jobs=7
python3 -I certify_bfs_exact.py 120 0 --trees=../../data/certify/sector_trees_120.json.gz
python3 -I straight_loops_120.py ../../data/certify/random_trees_120.json.gz   # 1 s
python3 -I all_bfs.py 120                          # claim 4 (53 min, 7 processes)
python3 -I all_bfs.py 120 classify
python3 -I all_bfs.py 600                          # 34 min, 7 processes
python3 -I path_classes.py 120 --all               # 1 min; 600: 10 min
python3 -I known_counterexamples.py                # about 10 min
python3 -I tree_observations_120.py                # 6 s
python3 -I second_checker.py 600                   # about 5.5 h with 8 processes (120-cell: 72 min)
```
