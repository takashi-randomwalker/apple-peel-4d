"""locality_axes.py -- the pairs the ball test leaves, tried on the centroid axis.

For every class of path_classes.py that locality_balls.py does not separate
by circumscribed balls, project both cells, exactly, onto the line through
their centroids (Euclidean, in R^4): if the two intervals do not overlap,
the plane normal to c_v - c_u separates them -- still an argument from the
centroids alone, but one that does not treat a cell as a ball.  What is
left is tried on the face planes of u and v and then edge x edge planes.

  python3 -I locality_axes.py 600
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


def method(u, v, p, q):
    A4, B4 = list(U.place(p).values()), list(U.place(q).values())
    if C.separates(C.sub(C.centroid(B4), C.centroid(A4)), A4, B4):
        return "centroid axis"
    A, B = U.chart3(p), U.chart3(q)
    gA, gB = C.cell_geometry(Pt, u, A), C.cell_geometry(Pt, v, B)
    if any(C.separates(n, gA[0], gB[0]) for n in gA[1] + gB[1]):
        return "face plane"
    return "edge x edge" if C.certificate(gA, gB) is not None else "OVERLAP"


METHODS = ["centroid axis", "face plane", "edge x edge", "OVERLAP"]
by = collections.defaultdict(collections.Counter)
for cls, (cnt, verdict, (p, q)) in PC.classes.items():
    u, v = p[-1], q[-1]
    k = 0
    while k < min(len(p), len(q)) and p[k] == q[k]:
        k += 1
    adj = v in ADJ[u]
    if adj and k == min(len(p), len(q)):
        continue                                         # tree edge
    if d2(cent(p), cent(q)) / touch > ball:
        continue
    by[("cut ridge" if adj else "not adjacent", (len(p) - k) + (len(q) - k) + 1)][method(u, v, p, q)] += 1

print("%s-cell, root 0: classes not separated by balls, by the first plane that separates them"
      % NAME)
print("kind           path   " + "  ".join("%13s" % m for m in METHODS))
tot = collections.Counter()
for key in sorted(by):
    tot += by[key]
    print("%-13s  %4d   " % key + "  ".join("%13d" % by[key][m] for m in METHODS))
print("total                " + "  ".join("%13d" % tot[m] for m in METHODS))
