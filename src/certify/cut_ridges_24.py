"""cut_ridges_24.py -- the 9 cut-ridge classes of the 24-cell, in exact arithmetic.

For a cut ridge F shared by cells a, b (adjacent in P, not joined in the
tree), compare b as the tree places it with b glued to a across F.  For
each class (path_classes.py) this reports, exactly in Q:

  * the cycle: layers of a, b and the length of the fundamental cycle,
  * the holonomy (glued b -> tree b): det, cos of the rotation angle
    (trace of the linear part; basis-independent), and the angle,
  * whether the plane of F separates a from the tree-placed b, with the
    least signed distance of a vertex of b from that plane, relative to
    the distance of the glued b's centroid (> 0: strictly beyond F).

  python3 -I cut_ridges_24.py
"""
import io, math, os, sys
from contextlib import redirect_stdout
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)                     # -I drops the script dir
sys.argv = [sys.argv[0], "24"]
with redirect_stdout(io.StringIO()):
    import path_classes as PC               # classes, Unfolder, symmetry group
import certify_bfs_exact as C

Pt, U, dist, ADJ = PC.Pt, PC.U, PC.dist, PC.ADJ


def det3(M):
    return (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1])
            - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
            + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]))


def inv3(M):
    d = det3(M)
    cof = [[(M[(j + 1) % 3][(i + 1) % 3] * M[(j + 2) % 3][(i + 2) % 3]
             - M[(j + 1) % 3][(i + 2) % 3] * M[(j + 2) % 3][(i + 1) % 3]) for j in range(3)]
           for i in range(3)]
    return [[cof[i][j] / d for j in range(3)] for i in range(3)]


def mul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(3)) for j in range(3)] for i in range(3)]


rows = []
for cls, (cnt, verdict, (p, q)) in PC.classes.items():
    a, b = p[-1], q[-1]
    k = 0
    while k < min(len(p), len(q)) and p[k] == q[k]:
        k += 1
    if b not in ADJ[a] or k == min(len(p), len(q)):
        continue                                       # not a cut ridge
    A, B, G = U.chart3(p), U.chart3(q), U.chart3(p + (b,))
    W = sorted(B)
    # linear part of the holonomy G -> B (chart coordinates, column vectors)
    base = W[0]
    for tri in __import__("itertools").combinations(W[1:], 3):
        Dg = [[G[w][i] - G[base][i] for w in tri] for i in range(3)]
        if det3(Dg) != 0:
            break
    Db = [[B[w][i] - B[base][i] for w in tri] for i in range(3)]
    M = mul(Db, inv3(Dg))
    cos = (M[0][0] + M[1][1] + M[2][2] - 1) / 2
    # plane of F in a's chart, oriented away from a
    fv = Pt.faces[Pt.shared[(a, b)]]
    i0, i1, i2 = C.independent([A[v] for v in fv], 2)
    P0 = A[fv[i0]]
    n = C.cross(C.sub(A[fv[i1]], P0), C.sub(A[fv[i2]], P0))
    ca = C.centroid(list(A.values()))
    if C.dot(n, C.sub(ca, P0)) > 0:
        n = [-x for x in n]
    side = lambda x: C.dot(n, C.sub(x, P0))
    ref = side(C.centroid(list(G.values())))
    margin = min(side(B[w]) for w in W) / ref
    length = (len(p) - k) + (len(q) - k) + 1
    # balls: |c_a - c_b|^2 against (2R)^2, both relative to the glued distance^2 (exact, in R^4)
    Ap, Bp, Gp = U.place(p), U.place(q), U.place(p + (b,))
    pa, pb, pg = (C.centroid(list(X.values())) for X in (Ap, Bp, Gp))
    d2 = lambda x, y: C.dot(C.sub(x, y), C.sub(x, y))
    touch = d2(pa, pg)
    far, ball = d2(pa, pb) / touch, 4 * max(d2(x, pa) for x in Ap.values()) / touch
    rows.append((dist[a], dist[b], length, det3(M), cos, margin, far, ball, cnt))

print("24-cell, root 0: %d cut-ridge classes" % len(rows))
print("(far = |c_a-c_b|^2, ball = (2R)^2, both over |c_a-c_glued|^2; far > ball: circumballs disjoint)")
print("layers cycle det cos(rot)  angle   F-plane sep.  min vertex beyond F    far     ball  balls sep.  pairs")
for la, lb, L, d, cs, m, far, ball, cnt in sorted(rows):
    print("%d,%d    %2d   %+d  %-8s %7.2f   %-5s         %-8s (%+7.3f)  %-6s  %-4s  %-5s      %4d"
          % (la, lb, L, d, cs, math.degrees(math.acos(float(cs))), m >= 0, m, float(m),
             far, ball, far > ball, cnt))
print("every class separated by the plane of F or by circumballs:",
      all(m >= 0 or far > ball for _, _, _, _, _, m, far, ball, _ in rows))
