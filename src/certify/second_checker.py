"""
second_checker.py
An independent re-verification of the exact certification.  It shares no
code with certify_bfs_exact.py, qsqrt5.py or verify_certificate.py, and
each step is done by a different method:

  numbers     its own Q(sqrt5) type (class K);
  data        its own parser of _4DData/ (the 5-cell, stored with
              GoldenRatio, is replaced by conv(e_1..e_5) in R^5);
  placement   by symmetry, not by solving for an affine map.  The
              reflection of R^d in the hyperplane through the centre of
              P and the shared face f is a symmetry of the regular
              polytope that swaps the two cells on f (checked: it must
              map the parent's vertices exactly onto the child's).  It
              gives the vertex correspondence w -> sigma(w), and the
              unfolded child is the parent's placement reflected in the
              plane of f inside the root's 3-flat, vertex by vertex;
  disjointness  no separating axes.  Two cells whose circumscribed balls
              meet in at most a point are disjoint.  For the rest, the
              vertices of the intersection polytope A cap B are
              enumerated exactly, and the interiors are disjoint iff
              those vertices span less than 3 dimensions;
  all trees   for the 5/8/16/24-cell, every BFS tree is enumerated
              explicitly and every pair of its cells is looked up; no
              lemma about compatible paths is used.  For the 120/600-cell
              only the paper's tree T_r (queue BFS, ascending neighbours).

Usage: python3 -I second_checker.py <5|8|16|24|120|600> [--jobs=N] [root ...]
"""

import itertools
import json
import os
import re
import sys
import time
from fractions import Fraction

BASE = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))


# ------------------------------------------------------------ Q(sqrt 5)

class K:
    """a + b*sqrt(5), a and b rational."""
    __slots__ = ("a", "b")

    def __init__(self, a, b=0):
        self.a, self.b = Fraction(a), Fraction(b)

    def __add__(s, o):
        o = o if isinstance(o, K) else K(o)
        return K(s.a + o.a, s.b + o.b)

    def __sub__(s, o):
        o = o if isinstance(o, K) else K(o)
        return K(s.a - o.a, s.b - o.b)

    def __mul__(s, o):
        o = o if isinstance(o, K) else K(o)
        return K(s.a * o.a + 5 * s.b * o.b, s.a * o.b + s.b * o.a)

    def __truediv__(s, o):
        o = o if isinstance(o, K) else K(o)
        n = o.a * o.a - 5 * o.b * o.b
        return s * K(o.a / n, -o.b / n)

    __radd__ = __add__
    __rmul__ = __mul__

    def __rsub__(s, o):
        return K(o) - s

    def __neg__(s):
        return K(-s.a, -s.b)

    def __eq__(s, o):
        o = o if isinstance(o, K) else K(o)
        return s.a == o.a and s.b == o.b

    def __hash__(s):
        return hash((s.a, s.b))

    def positive(s):
        # a + b sqrt5 > 0
        if s.b == 0:
            return s.a > 0
        if s.b > 0:
            return s.a >= 0 or 5 * s.b * s.b > s.a * s.a
        return s.a > 0 and s.a * s.a > 5 * s.b * s.b

    def __gt__(s, o): return (s - o).positive()
    def __lt__(s, o): return (K(o) - s if not isinstance(o, K) else o - s).positive()
    def __ge__(s, o): return not s < o
    def __le__(s, o): return not s > o


PHI = K(Fraction(1, 2), Fraction(1, 2))


# --------------------------------------------------------------- data

def load(name):
    if name == "5":
        n = 5
        V = [[Fraction(int(i == j)) for j in range(n)] for i in range(n)]
        F = [list(t) for t in itertools.combinations(range(n), 3)]
        C = [[k for k, f in enumerate(F) if set(f) <= set(c)] for c in itertools.combinations(range(n), 4)]
        return V, F, C
    txt = open(os.path.join(BASE, "_4DData", "f%s.m" % name)).read()
    txt = re.sub(r"(-?)(Sqrt\[5\]|GoldenRatio(?:\^\(-?\d+\)|\^\d+)?)", r'"\1\2"', txt)
    raw = json.loads(txt.replace("{", "[").replace("}", "]"))
    sym = {"Sqrt[5]": K(0, 1), "GoldenRatio": PHI, "GoldenRatio^2": PHI * PHI,
           "GoldenRatio^(-1)": PHI - 1, "GoldenRatio^(-2)": (PHI - 1) * (PHI - 1)}

    def num(x):
        if isinstance(x, int):
            return K(x)
        neg = x.startswith("-")
        v = sym[x.lstrip("-")]
        return -v if neg else v

    if any(isinstance(x, str) for v in raw[0] for x in v):
        V = [[num(x) for x in v] for v in raw[0]]
    else:
        V = [[Fraction(x) for x in v] for v in raw[0]]
    return V, [[i - 1 for i in f] for f in raw[3]], [[i - 1 for i in c] for c in raw[5]]


# ------------------------------------------------------ linear algebra

def dot(a, b):
    s = a[0] * b[0]
    for x, y in zip(a[1:], b[1:]):
        s = s + x * y
    return s


def sub(a, b):
    return [x - y for x, y in zip(a, b)]


def nullspace(rows, d):
    """Basis of {x in F^d : r.x = 0 for every row r}, by reduced row echelon form."""
    m = [list(r) for r in rows]
    piv = []
    rk = 0
    for c in range(d):
        p = next((i for i in range(rk, len(m)) if m[i][c] != 0), None)
        if p is None:
            continue
        m[rk], m[p] = m[p], m[rk]
        inv = m[rk][c]
        m[rk] = [x / inv for x in m[rk]]
        for i in range(len(m)):
            if i != rk and m[i][c] != 0:
                f = m[i][c]
                m[i] = [x - f * y for x, y in zip(m[i], m[rk])]
        piv.append(c)
        rk += 1
    zero, one = m[0][0] * 0 if m else 0, None
    basis = []
    for free in (c for c in range(d) if c not in piv):
        x = [zero] * d
        x[free] = zero + 1
        for i, c in enumerate(piv):
            x[c] = zero - m[i][free]
        basis.append(x)
    return basis


def rank(rows):
    if not rows:
        return 0
    return len(rows[0]) - len(nullspace(rows, len(rows[0])))


def reflect(x, point, n):
    """Reflection of x in the hyperplane through point with normal n."""
    t = dot(sub(x, point), n) * 2 / dot(n, n)
    return [xi - t * ni for xi, ni in zip(x, n)]


def det3(m):
    return (m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]))


def solve3(m, rhs):
    d = det3(m)
    out = []
    for k in range(3):
        mk = [row[:k] + [rhs[i]] + row[k + 1:] for i, row in enumerate(m)]
        out.append(det3(mk) / d)
    return out


def cross(a, b):
    return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]


# -------------------------------------------------------- the checker

class Checker:
    def __init__(self, name):
        self.V, self.F, self.C = load(name)
        self.d = len(self.V[0])
        self.nC = len(self.C)
        self.cv = [sorted({v for f in c for v in self.F[f]}) for c in self.C]
        own = {}
        for c, fs in enumerate(self.C):
            for f in fs:
                own.setdefault(f, []).append(c)
        self.face = {}
        self.adj = [[] for _ in range(self.nC)]
        for f, cs in own.items():
            assert len(cs) == 2
            a, b = cs
            self.face[(a, b)] = self.face[(b, a)] = f
            self.adj[a].append(b)
            self.adj[b].append(a)
        for a in self.adj:
            a.sort()
        n = len(self.V)
        self.centre = [sum((v[k] for v in self.V[1:]), self.V[0][k]) / n for k in range(self.d)]
        # directions orthogonal to aff(P) (only for the 5-cell in R^5)
        self.offP = nullspace([sub(v, self.V[0]) for v in self.V[1:]], self.d)
        self.sigma = {}
        self.index = {tuple(v): i for i, v in enumerate(self.V)}

    def swap(self, p, c):
        """Vertex correspondence parent -> child from the symmetry fixing the shared face."""
        if (p, c) not in self.sigma:
            fv = self.F[self.face[(p, c)]]
            rows = [sub(self.V[v], self.centre) for v in fv] + self.offP
            ns = nullspace(rows, self.d)
            assert len(ns) == 1
            m = ns[0]
            mp = {}
            for w in self.cv[p]:
                img = reflect(self.V[w], self.centre, m)
                v = self.index.get(tuple(img))
                assert v is not None and v in self.cv[c], "reflection is not a symmetry swapping the cells"
                mp[w] = v
            assert sorted(mp.values()) == self.cv[c]
            for v in fv:
                assert mp[v] == v
            self.sigma[(p, c)] = mp
        return self.sigma[(p, c)]

    def setup_root(self, r):
        self.r = r
        rv = [self.V[v] for v in self.cv[r]]
        dirs = [sub(x, rv[0]) for x in rv[1:]]
        self.flatN = nullspace(dirs, self.d)          # normals of the root 3-flat
        assert len(self.flatN) == self.d - 3
        # three coordinates on which the flat projects bijectively
        b = [dirs[i] for i in range(len(dirs))]
        basis = []
        for x in b:
            if rank(basis + [x]) > len(basis):
                basis.append(x)
        self.coords = next(c for c in itertools.combinations(range(self.d), 3)
                           if det3([[x[k] for k in c] for x in basis]) != 0)
        self.memo = {(r,): {v: self.V[v] for v in self.cv[r]}}

    def place(self, path):
        if path in self.memo:
            return self.memo[path]
        par = self.place(path[:-1])
        p, c = path[-2], path[-1]
        fv = self.F[self.face[(p, c)]]
        F = [par[v] for v in fv]
        E = [sub(x, F[0]) for x in F[1:]]
        n = nullspace(E + self.flatN, self.d)
        assert len(n) == 1
        n = n[0]
        s = self.swap(p, c)
        out = {s[w]: reflect(x, F[0], n) for w, x in par.items()}
        for v in fv:
            assert out[v] == par[v]
        self.memo[path] = out
        return out

    # --------------------------------------------------- disjointness
    def solid(self, path):
        """Projected 3D vertices and inward half-spaces n.x >= b of a placed cell."""
        pl = self.place(path)
        c = path[-1]
        X = {v: [x[k] for k in self.coords] for v, x in pl.items()}
        cen = [sum((X[v][k] for v in self.cv[c][1:]), X[self.cv[c][0]][k]) / len(self.cv[c]) for k in range(3)]
        hs = []
        for f in self.C[c]:
            P = [X[v] for v in self.F[f]]
            n = next(w for a, b in itertools.combinations(P[1:], 2)
                     for w in [cross(sub(a, P[0]), sub(b, P[0]))] if any(t != 0 for t in w))
            b = dot(n, P[0])
            if dot(n, cen) < b:
                n, b = [-t for t in n], -b
            hs.append((n, b))
        return list(X.values()), hs

    def interiors_disjoint(self, pa, pb):
        A, hA = self.solid(pa)
        B, hB = self.solid(pb)
        hs = hA + hB
        pts = set()
        for i, j, k in itertools.combinations(range(len(hs)), 3):
            m = [hs[i][0], hs[j][0], hs[k][0]]
            if det3(m) == 0:
                continue
            x = solve3(m, [hs[i][1], hs[j][1], hs[k][1]])
            if all(dot(n, x) >= b for n, b in hs):
                pts.add(tuple(x))
        pts = list(pts)
        if len(pts) < 4:
            return True
        return rank([sub(list(q), list(pts[0])) for q in pts[1:]]) < 3

    def centroid4(self, path):
        pl = self.place(path)
        xs = list(pl.values())
        return [sum((x[k] for x in xs[1:]), xs[0][k]) / len(xs) for k in range(self.d)]

    def radius2(self):
        c = self.cv[0]
        cen = [sum((self.V[v][k] for v in c[1:]), self.V[c[0]][k]) / len(c) for k in range(self.d)]
        r2 = [dot(sub(self.V[v], cen), sub(self.V[v], cen)) for v in c]
        assert all(x == r2[0] for x in r2)
        return r2[0]

    def far(self, ca, cb, R2):
        """Circumscribed balls meet in at most a point: |ca-cb|^2 >= 4 R^2."""
        d = sub(ca, cb)
        return dot(d, d) >= 4 * R2


# ------------------------------------------------------------- drivers

def bfs(ck, r):
    dist = [-1] * ck.nC
    dist[r] = 0
    q = [r]
    for p in q:
        for nb in ck.adj[p]:
            if dist[nb] < 0:
                dist[nb] = dist[p] + 1
                q.append(nb)
    return dist, q


def all_trees(ck, r):
    """Every BFS tree rooted at r, every pair of its cells."""
    ck.setup_root(r)
    R2 = ck.radius2()
    dist, order = bfs(ck, r)
    rest = [v for v in order if v != r]
    opts = [[p for p in ck.adj[v] if dist[p] == dist[v] - 1] for v in rest]
    pid, verdict = {}, {}
    cen = {}

    def id_of(path):
        if path not in pid:
            pid[path] = len(pid)
            cen[pid[path]] = (path, ck.centroid4(path))
        return pid[path]

    trees = 0
    n_far = n_enum = 0
    for choice in itertools.product(*opts):
        trees += 1
        path = {r: (r,)}
        for v, p in zip(rest, choice):         # rest is in BFS order
            path[v] = path[p] + (v,)
        ids = [id_of(path[v]) for v in range(ck.nC)]
        for i, j in itertools.combinations(ids, 2):
            key = (i, j) if i < j else (j, i)
            if key not in verdict:
                (pa, ca), (pb, cb) = cen[key[0]], cen[key[1]]
                if ck.far(ca, cb, R2):
                    verdict[key] = True
                    n_far += 1
                else:
                    verdict[key] = ck.interiors_disjoint(pa, pb)
                    n_enum += 1
            if not verdict[key]:
                return {"root": r, "trees": trees, "bad": key, "far": n_far, "enum": n_enum}
    return {"root": r, "trees": trees, "bad": None, "far": n_far, "enum": n_enum}


def queue_tree(ck, r):
    """The paper's T_r: every pair of its cells."""
    ck.setup_root(r)
    R2 = ck.radius2()
    parent = {r: None}
    q = [r]
    for p in q:
        for nb in ck.adj[p]:
            if nb not in parent:
                parent[nb] = p
                q.append(nb)
    path = {r: (r,)}
    for v in q[1:]:
        path[v] = path[parent[v]] + (v,)
    cen = [ck.centroid4(path[v]) for v in range(ck.nC)]
    # sweep on centroid coordinate 0: a gap >= 2R there already makes the balls disjoint
    order = sorted(range(ck.nC), key=lambda v: cen[v][0])
    R2x4 = 4 * R2
    active, n_far, n_enum, n_gap = [], 0, 0, 0
    bad = None
    for b in order:
        keep = []
        for a in active:
            g = cen[b][0] - cen[a][0]
            if g * g >= R2x4:            # g >= 0 by the sort; this and every later b
                continue
            keep.append(a)
        active = keep
        for a in active:
            if ck.far(cen[a], cen[b], R2):
                n_far += 1
            else:
                n_enum += 1
                if not ck.interiors_disjoint(path[a], path[b]) and bad is None:
                    bad = (a, b)
        active.append(b)
    total = ck.nC * (ck.nC - 1) // 2
    return {"root": r, "trees": 1, "bad": bad, "far": n_far + (total - n_far - n_enum), "enum": n_enum}


_ck = {}


def work(args):
    name, r = args
    if name not in _ck:
        _ck[name] = Checker(name)
    t = time.time()
    res = (queue_tree if name in ("120", "600") else all_trees)(_ck[name], r)
    res["time"] = time.time() - t
    return res


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    name = args[0]
    jobs = next((int(a[7:]) for a in sys.argv if a.startswith("--jobs=")), os.cpu_count())
    ck = Checker(name)
    roots = [int(x) for x in args[1:]] or list(range(ck.nC))
    mode = "tree T_r" if name in ("120", "600") else "every BFS tree"
    print("second checker: %s-cell, %s, %d roots" % (name, mode, len(roots)), flush=True)
    from concurrent.futures import ProcessPoolExecutor
    t0 = time.time()
    trees = far = enum = 0
    bad = []
    with ProcessPoolExecutor(jobs) as ex:
        for i, res in enumerate(ex.map(work, [(name, r) for r in roots])):
            trees += res["trees"]
            far += res["far"]
            enum += res["enum"]
            if res["bad"]:
                bad.append(res)
            if i < 2 or (i + 1) % max(1, len(roots) // 10) == 0 or res["bad"]:
                print("  root %3d  trees %7d  pairs: far %8d, enumerated %6d  %s  (%.1fs)"
                      % (res["root"], res["trees"], res["far"], res["enum"],
                         "OVERLAP %s" % (res["bad"],) if res["bad"] else "ok", res["time"]), flush=True)
    print("SECOND CHECKER  %s-cell  roots %d  trees %d  far %d  enumerated %d  overlapping roots %d  (%.1fs)"
          % (name, len(roots), trees, far, enum, len(bad), time.time() - t0))


if __name__ == "__main__":
    main()
