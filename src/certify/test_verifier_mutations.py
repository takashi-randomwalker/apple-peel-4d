"""
test_verifier_mutations.py
A verifier that accepts everything proves nothing.  Tamper with a valid
certificate in six ways and require verify_certificate.py to reject each.

The subtle one is "mirror": the placement is replaced by its reflection in
the shared face, which keeps conditions (a)-(c) and breaks only (d).

Usage: python3 -I test_verifier_mutations.py <cert.json.gz>
"""

import copy
import gzip
import json
import os
import subprocess
import sys
import tempfile
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # -I drops the script dir
from certify_bfs_exact import coords_in, dot, scl, sub

HERE = os.path.dirname(os.path.abspath(__file__))


def verify(doc):
    with tempfile.NamedTemporaryFile(suffix=".json.gz", delete=False) as tmp:
        path = tmp.name
    with gzip.open(path, "wt") as fh:
        json.dump(doc, fh)
    out = subprocess.run([sys.executable, "-I", os.path.join(HERE, "verify_certificate.py"), path],
                         capture_output=True, text=True).stdout.strip()
    os.unlink(path)
    return out


def mirror(doc):
    """Reflect one non-root placement (and nothing else) in its shared face."""
    rc = doc["roots"][0]
    faces, cells = doc["faces"], doc["cells"]
    pl = {tuple(p): m for p, m in rc["placements"]}
    path = next(p for p in pl if len(p) == 2)
    par, m = pl[path[:-1]], pl[path]
    f = (set(cells[path[0]]) & set(cells[path[1]])).pop()
    F = [[Fraction(x) for x in m[str(v)]] for v in faces[f]]
    root = [[Fraction(x) for x in v] for v in pl[(rc["root"],)].values()]
    flat = [sub(x, root[0]) for x in root[1:]]
    e = [sub(F[1], F[0]), sub(F[2], F[0])]
    n = next(res for b in flat for _, res in [coords_in(e, b)] if any(res))
    for k, x in m.items():
        y = [Fraction(t) for t in x]
        y2 = sub(y, scl(2 * dot(sub(y, F[0]), n) / dot(n, n), n))
        m[k] = [str(t) for t in y2]
    for i, (p, _) in enumerate(rc["placements"]):
        if tuple(p) == path:
            rc["placements"][i][1] = m
    return doc


def main():
    with gzip.open(sys.argv[1], "rt") as fh:
        good = json.load(fh)
    print("untouched      :", verify(good))
    cases = {}

    d = copy.deepcopy(good)
    p, m = d["roots"][0]["placements"][-1]
    k = next(iter(m))
    m[k][0] = str(Fraction(m[k][0]) + Fraction(1, 10**9))
    cases["coordinate"] = d

    cases["mirror"] = mirror(copy.deepcopy(good))

    d = copy.deepcopy(good)
    d["roots"][0]["pairs"].pop()
    cases["drop pair"] = d

    d = copy.deepcopy(good)
    d["roots"][0]["pairs"][0][2] = [str(-Fraction(x)) for x in d["roots"][0]["pairs"][0][2]]
    cases["flip u"] = d

    d = copy.deepcopy(good)
    d["roots"][0]["pairs"][0][2] = ["0"] * len(d["roots"][0]["pairs"][0][2])
    cases["zero u"] = d

    d = copy.deepcopy(good)
    d["vertices"][0][0] = str(Fraction(d["vertices"][0][0]) + 1)
    cases["polytope"] = d

    ok = True
    for name, doc in cases.items():
        out = verify(doc)
        rejected = out.startswith("FAIL")
        ok &= rejected
        print("%-15s: %s  %s" % (name, "rejected" if rejected else "*** ACCEPTED ***", out))
    print("ALL MUTATIONS REJECTED" if ok else "*** SOME MUTATION ACCEPTED ***")


if __name__ == "__main__":
    main()
