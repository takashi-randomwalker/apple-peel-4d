"""close_cut_ridges.py -- the cut-ridge classes that the ball test leaves open.

locality_balls.py separates almost every root-path pair by distance alone
(disjoint circumscribed balls).  This script lists, for one polytope, the
cut-ridge classes it leaves -- b, sharing ridge F with a, placed by the tree
closer to a than 2R -- with:

  * the fundamental cycle a ... L ... b (a): layers of a, b, its length, and
    the turn at every cell of the cycle -- 's' straight through (the two
    faces are opposite), 'a' through two faces sharing an edge, 'b' between
    (faces sharing no vertex, not opposite),
  * the holonomy (float, from path_classes.py): rotation angle, centroid
    distance, and the displacement of b along F's outward normal and across,
    all relative to the distance of a glued pair,
  * exactly: whether the plane of F separates a from b (and the least signed
    distance of a vertex of b beyond it, relative to the glued centroid),
    otherwise the simplest separating plane, or OVERLAP.

  python3 -I close_cut_ridges.py 120
"""
import io, os, sys
from contextlib import redirect_stdout
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)                     # -I drops the script dir
NAME = sys.argv[1]
sys.argv = [sys.argv[0], NAME]
with redirect_stdout(io.StringIO()):
    import path_classes as PC
import certify_bfs_exact as C

Pt, U, ADJ, dist = PC.Pt, PC.U, PC.ADJ, PC.dist
V = [set(s) for s in Pt.cellV]
d2 = lambda x, y: C.dot(C.sub(x, y), C.sub(x, y))
cent = lambda path: C.centroid(list(U.place(path).values()))
r0 = PC.r
touch = d2(cent((r0,)), cent((r0, Pt.adj[r0][0])))
ball = 4 * max(d2(x, cent((r0,))) for x in U.place((r0,)).values()) / touch


def turn(x, m, q):
    F1, F2 = V[x] & V[m], V[q] & V[m]
    if F1 & F2:
        return "a"
    linked = any((s in F1 and t in F2) or (t in F1 and s in F2) for s, t in Pt.cellE[m])
    return "b" if linked else "s"


def separation(a, b, p, q):
    A, B, G = U.chart3(p), U.chart3(q), U.chart3(p + (b,))
    fv = Pt.faces[Pt.shared[(a, b)]]
    i0, i1, i2 = C.independent([A[v] for v in fv], 2)
    P0 = A[fv[i0]]
    n = C.cross(C.sub(A[fv[i1]], P0), C.sub(A[fv[i2]], P0))
    if C.dot(n, C.sub(C.centroid(list(A.values())), P0)) > 0:
        n = [-x for x in n]
    side = lambda x: C.dot(n, C.sub(x, P0))
    margin = min(side(x) for x in B.values()) / side(C.centroid(list(G.values())))
    if margin >= 0:
        return "plane of F", margin
    gA, gB = C.cell_geometry(Pt, a, A), C.cell_geometry(Pt, b, B)
    if any(C.separates(w, gA[0], gB[0]) for w in gA[1]):
        return "face of a", margin
    if any(C.separates(w, gA[0], gB[0]) for w in gB[1]):
        return "face of b", margin
    if C.separates(C.sub(C.centroid(gB[0]), C.centroid(gA[0])), gA[0], gB[0]):
        return "centroids", margin
    return ("edge x edge" if C.certificate(gA, gB) is not None else "OVERLAP"), margin


rows = []
for cls, (cnt, verdict, (p, q)) in PC.classes.items():
    a, b = p[-1], q[-1]
    k = 0
    while k < min(len(p), len(q)) and p[k] == q[k]:
        k += 1
    if b not in ADJ[a] or k == min(len(p), len(q)):
        continue
    if d2(cent(p), cent(q)) / touch > ball:
        continue
    cyc = list(p[k - 1:][::-1]) + list(q[k:])          # a ... L ... b
    turns = "".join(turn(cyc[i - 1], cyc[i], cyc[(i + 1) % len(cyc)]) for i in range(len(cyc)))
    ang, dt, nrm, acr = PC.holonomy(p, q)
    how, margin = separation(a, b, p, q)
    rows.append((len(cyc), dist[a], dist[b], turns, ang, dt, nrm, acr, how, float(margin), cnt,
                 True in verdict))

print("%s-cell, root 0: %d cut-ridge classes closer than the ball bound (2R)^2 = %.4f"
      % (NAME, len(rows), float(ball)))
print("turns: cycle a ... L ... b, one letter per cell from a (s straight, a sharp, b bent)")
print("cycle layers turns        rot(deg) dist  normal  across  separated by   F-margin  pairs  overlap")
for L, la, lb, t, ang, dt, nrm, acr, how, m, cnt, ov in sorted(rows):
    print("%3d   %d,%d   %-12s %7.2f  %5.3f  %+6.3f  %6.3f  %-13s  %+7.3f  %5d  %s"
          % (L, la, lb, t, ang, dt, nrm, acr, how, m, cnt, ov))
