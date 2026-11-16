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
| `straight_loops_120.py` | predicts the 120-cell overlaps without geometry; matches the exact count on all 4,448 trees above (observed, not proved) |
| `crosscheck_float.py` / `.wls` | the exact placements agree with the paper's floating-point unfolding (`unfoldTFP`, `bfsUnfoldP`) |

Data and logs are in `data/certify/`: `cert_{5,8,16,24}cell.json.gz`
(certificates, 0.2 MB at most), the given trees
`random_trees_{120,600}.json.gz` and `sector_trees_120.json.gz`
(gzip JSON: per root a list of parent arrays, 0-based, −1 at the root; the
sector file has one entry, root 0), and the run logs.

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
python3 -I second_checker.py 600                   # about 5.5 h with 8 processes (120-cell: 72 min)
```
