"""locality_balls.py -- which root-path pairs are separated by distance alone.

For every class of path_classes.py except tree edges, compare, exactly,
|c_u - c_v|^2 with (2R)^2 (R the circumradius of a cell): if it is larger
the circumscribed balls are disjoint and u, v cannot overlap, whatever else
is true.  The question is whether this holds beyond some length of the tree
path u ... L ... v, so that only short, local configurations need a finer
argument.  Distances are relative to that of two glued cells.

  python3 -I locality_balls.py 24|120|600
"""
import collections, io, os, sys
from contextlib import redirect_stdout
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)                     # -I drops the script dir
NAME = sys.argv[1]
sys.argv = [sys.argv[0], NAME]
with redirect_stdout(io.StringIO()):
    import path_classes as PC
import certify_bfs_exact as C

Pt, U, ADJ = PC.Pt, PC.U, PC.ADJ
d2 = lambda x, y: C.dot(C.sub(x, y), C.sub(x, y))
cent = lambda path: C.centroid(list(U.place(path).values()))
r0 = PC.r
touch = d2(cent((r0,)), cent((r0, Pt.adj[r0][0])))
ball = 4 * max(d2(x, cent((r0,))) for x in U.place((r0,)).values()) / touch

by = collections.defaultdict(lambda: [0, 0, None, 0, 0])   # (kind, length) -> [classes, balls, min far, overlapping, pairs]
for cls, (cnt, verdict, (p, q)) in PC.classes.items():
    u, v = p[-1], q[-1]
    k = 0
    while k < min(len(p), len(q)) and p[k] == q[k]:
        k += 1
    if v in ADJ[u] and k == min(len(p), len(q)):
        continue                                         # tree edge
    kind = "cut ridge" if v in ADJ[u] else "not adjacent"
    e = by[(kind, (len(p) - k) + (len(q) - k) + 1)]
    far = d2(cent(p), cent(q)) / touch
    e[0] += 1
    e[1] += far > ball
    e[2] = far if e[2] is None or far < e[2] else e[2]
    e[3] += True in verdict
    e[4] += cnt

print("%s-cell, root 0: (2R)^2 = %.4f (glued distance^2 = 1)" % (NAME, float(ball)))
print("kind           path   classes   by balls   not   min |c_u-c_v|^2   overlapping")
for kind in ("cut ridge", "not adjacent"):
    for L in sorted(l for k_, l in by if k_ == kind):
        n, b, m, o, _ = by[(kind, L)]
        print("%-13s  %4d   %7d   %8d   %3d   %10.4f        %d" % (kind, L, n, b, n - b, float(m), o))
fail = [L for (k_, L), e in by.items() if e[1] < e[0]]
print("longest tree path with a class not separated by balls: %s cells" % (max(fail) if fail else "-"))
print("classes not separated by balls: %d of %d"
      % (sum(e[0] - e[1] for e in by.values()), sum(e[0] for e in by.values())))
