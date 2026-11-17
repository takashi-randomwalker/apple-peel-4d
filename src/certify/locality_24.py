"""locality_24.py -- the 26 classes of 24-cell cell pairs that share no ridge.

Part 1 of the case analysis (Section 7 of the paper): cells u, v that are
not adjacent in P never overlap in a BFS net.  For each class of
path_classes.py with u, v not adjacent, this reports, exactly in Q:

  * what u and v share in P (a vertex, an edge, nothing) and the length of
    their tree path,
  * |c_u - c_v|^2 over the squared distance of two glued cells, and whether
    it exceeds (2R)^2 (circumscribed balls disjoint),
  * otherwise the simplest separating plane found: a face plane of u or v,
    the plane normal to c_v - c_u, or an edge-cross-edge plane.

  python3 -I locality_24.py
"""
import collections, io, os, sys
from contextlib import redirect_stdout
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)                     # -I drops the script dir
sys.argv = [sys.argv[0], "24"]
with redirect_stdout(io.StringIO()):
    import path_classes as PC
import certify_bfs_exact as C

Pt, U, ADJ, dist = PC.Pt, PC.U, PC.ADJ, PC.dist
V = [set(s) for s in Pt.cellV]
d2 = lambda x, y: C.dot(C.sub(x, y), C.sub(x, y))

# squared distance of two glued cells and (2R)^2, from the root and a neighbour
r0, n0 = PC.r, Pt.adj[PC.r][0]
c0 = C.centroid(list(U.place((r0,)).values()))
touch = d2(c0, C.centroid(list(U.place((r0, n0)).values())))
ball = 4 * max(d2(x, c0) for x in U.place((r0,)).values()) / touch


def plane(u, v, p, q):
    """The simplest separating plane in the root chart, by kind."""
    A, B = U.chart3(p), U.chart3(q)
    gA, gB = C.cell_geometry(Pt, u, A), C.cell_geometry(Pt, v, B)
    if any(C.separates(n, gA[0], gB[0]) for n in gA[1]):
        return "face of u"
    if any(C.separates(n, gA[0], gB[0]) for n in gB[1]):
        return "face of v"
    if C.separates(C.sub(C.centroid(gB[0]), C.centroid(gA[0])), gA[0], gB[0]):
        return "centroids"
    return "edge x edge" if C.certificate(gA, gB) is not None else "NONE"


rows = []
for cls, (cnt, verdict, (p, q)) in PC.classes.items():
    u, v = p[-1], q[-1]
    if v in ADJ[u]:
        continue
    k = 0
    while k < min(len(p), len(q)) and p[k] == q[k]:
        k += 1
    share = len(V[u] & V[v])
    rel = {0: "nothing", 1: "a vertex", 2: "an edge"}[share]
    far = d2(*(C.centroid(list(U.place(x).values())) for x in (p, q))) / touch
    how = "balls" if far > ball else plane(u, v, p, q)
    rows.append((rel, (len(p) - k) + (len(q) - k) + 1, far, how, cnt))

print("24-cell, root 0: %d classes of pairs not adjacent in P; (2R)^2 = %s (glued distance^2 = 1)"
      % (len(rows), ball))
print("share in P   tree path   |c_u-c_v|^2   separated by     pairs")
for rel, L, far, how, cnt in sorted(rows, key=lambda t: (t[0], t[1], t[2])):
    print("%-10s   %5d       %-8s      %-12s   %5d" % (rel, L, far, how, cnt))
print("by method:", dict(collections.Counter(t[3] for t in rows)))
print("all separated:", all(t[3] != "NONE" for t in rows))
