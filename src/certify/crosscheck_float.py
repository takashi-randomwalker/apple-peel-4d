"""
crosscheck_float.py
Does the exact reflection construction reproduce the paper's unfolding?

The paper unfolds with floating-point hinge rotations (unfoldTFP in
src/bfsnet/face_rotation_net_all4D_v2.m).  Here we pick random BFS trees,
place them exactly, and write the parent arrays together with the exact
squared distances between placed cell centroids.  crosscheck_float.wls
unfolds the same trees with unfoldTFP and compares.  Distances are used
because they do not see the face-up rotation of Step 1, which the exact
construction omits.

With "--queue r1 r2 ..." in place of the tree count, the trees are T_r
(Polytope.queue_tree) and the .wls side calls bfsUnfoldP itself, so the
comparison also checks that the two codes build the same tree.

Usage: python3 -I crosscheck_float.py <name> <trees per root> <out.json>
       python3 -I crosscheck_float.py <name> --queue <root> ... <out.json>
"""

import json
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # -I drops the script dir
from certify_bfs_exact import Polytope, Unfolder, centroid, dot, sub


def main():
    name, out = sys.argv[1], sys.argv[-1]
    P = Polytope(name)
    rng = random.Random(20261114)
    if sys.argv[2] == "--queue":
        jobs = [(r, P.queue_tree(r), True) for r in map(int, sys.argv[3:-1])]
    else:
        jobs = []
        for r in range(P.nC):
            dist, opts = P.layers(r)
            for _ in range(int(sys.argv[2])):
                jobs.append((r, [None if v == r else rng.choice(opts[v]) for v in range(P.nC)], False))
    cases = []
    for r, parent, queue in jobs:
        U = Unfolder(P, r)

        def path(v):
            return (v,) if v == r else path(parent[v]) + (v,)

        cen = [centroid(list(U.place(path(v)).values())) for v in range(P.nC)]
        d2 = [[float(dot(sub(a, b), sub(a, b))) for b in cen] for a in cen]
        case = {"root": r + 1, "d2": d2}
        if not queue:
            case["parent"] = [0 if p is None else p + 1 for p in parent]
        cases.append(case)
    with open(out, "w") as fh:
        json.dump({"polytope": name, "cases": cases}, fh)
    print("wrote %d trees" % len(cases))


if __name__ == "__main__":
    main()
