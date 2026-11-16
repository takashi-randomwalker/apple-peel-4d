"""compare_recheck.py -- the re-run of 2026-11-16 (recheck_uniform_allroots.wls,
data/uniform_recheck/*.json) against the published per-root results
(uniform_nets/allroots.csv, made with the buggy test).  Writes the
classification change and every polytope whose verdicts or overlap counts
changed.  Output saved as data/uniform_recheck/compare_with_allroots_csv.txt.

  python3 -I compare_recheck.py
"""
import csv,json,os,collections
B="/Users/yoshino/Library/CloudStorage/Dropbox/260324Peeling4D/"
idx={r["file"]:r for r in csv.DictReader(open(B+"_4DData_uniform/index.csv"))}
old=collections.defaultdict(dict)
for r in csv.DictReader(open(B+"uniform_nets/allroots.csv")): old[r["coxeter_symbol"]][int(r["root"])]=(r["verdict"]=="valid",int(r["intersecting_pairs"]))
def cls(n,t): return "ALL VALID" if n==t else ("ALL INVALID" if n==0 else "MIXED")
rows=[]; trans=collections.Counter()
for f,r in idx.items():
    new=json.load(open(B+"data/uniform_recheck/%s.json"%f)); sym=r["coxeter_symbol"]
    o=old[sym.replace("_"," ")]
    t=len(new); assert len(o)==t
    nv=sum(x[2]==0 for x in new); ov=sum(v for v,_ in o.values())
    ch=sum(1 for x in new if (x[2]==0)!=o[x[0]][0])
    lost=sum(1 for x in new if x[2]>0 and o[x[0]][0]); gained=sum(1 for x in new if x[2]==0 and not o[x[0]][0])
    trans[(cls(ov,t),cls(nv,t))]+=1
    rows.append((sym,r["name"],t,ov,nv,lost,gained,sum(q for _,q in o.values()),sum(x[2] for x in new)))
print("roots %d  valid old %d new %d"%(sum(r[2] for r in rows),sum(r[3] for r in rows),sum(r[4] for r in rows)))
print("classification old -> new:")
for k,v in sorted(trans.items()): print("   %-11s -> %-11s %2d"%(k[0],k[1],v))
print("\nchanged polytopes:")
print("%-22s %-26s %5s %7s %7s %5s %6s  %s"%("symbol","name","cells","old ok","new ok","lost","gained","overlap pairs old -> new"))
for r in sorted(rows,key=lambda r:r[2]):
    if r[5] or r[6] or r[7]!=r[8]: print("%-22s %-26s %5d %7d %7d %5d %6d  %d -> %d"%(r[0],r[1][:26],r[2],r[3],r[4],r[5],r[6],r[7],r[8]))
