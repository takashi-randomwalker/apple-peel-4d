"""straight_loops_120.py -- predict the overlaps of a 120-cell BFS net
without any geometry.

The count equals the number of overlapping cell pairs for every BFS tree:
all_bfs.py certifies every pair of root paths and finds that exactly
these straight 7-cycles overlap (2026-11-16).  It was first observed on the
2,400 random and 2,048 D_5-invariant trees; it is 0 for T_r from every root.

Usage: python3 -I straight_loops_120.py TREES.json.gz [roots...]
       prints a JSON list of counts, one per tree, in file order."""
# The rule (no coordinates are used):
# For each pair (a,b) of adjacent cells in one layer, with parents pa != pb and
# {a,b,pa,pb} the four cells around a vertex, follow both root paths to their
# meeting cell.  If the loop has 7 cells L,x,pa,a,b,pb,y, orient it so that
# pa shares a face with y.  Count the pair if x -> pa -> a passes straight
# through pa: the faces pa∩x and pa∩a are opposite (no edge joins them).
import gzip, json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # -I drops the script dir
from certify_bfs_exact import Polytope
P = Polytope("120")
V = [set(s) for s in P.cellV]; ADJ = [set(a) for a in P.adj]
def straight(x, m, q):
    F1, F2 = V[x] & V[m], V[q] & V[m]
    return not any((u in F1 and w in F2) or (w in F1 and u in F2) for u, w in P.cellE[m])
def count(r, parent):
    dist, _ = P.layers(r)
    def path(v):
        out = [v]
        while out[-1] != r: out.append(parent[out[-1]])
        return out[::-1]
    n = 0
    for a in range(P.nC):
        for b in P.adj[a]:
            if b <= a or dist[a] != dist[b]: continue
            pa, pb = parent[a], parent[b]
            if pa == pb or pb not in ADJ[pa] or b not in ADJ[pa] or a not in ADJ[pb]: continue
            A, B = path(a), path(b)
            if len(A) < 4 or A[-4] != B[-4] or A[-3] == B[-3]: continue   # 7-cell loop
            x, y = A[-3], B[-3]
            if pa in ADJ[y]: n += straight(x, pa, a)
            else:
                assert pb in ADJ[x]
                n += straight(y, pb, b)
    return n
f, roots = sys.argv[1], [int(z) for z in sys.argv[2:]]
trees = json.load(gzip.open(f, "rt")); roots = roots or range(len(trees))
print(json.dumps([count(r, p) for r in roots for p in trees[r]]))
