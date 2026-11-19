"""Red controls for the C_6^13 certificate: planted errors that both checkers must refuse.
  ctl_flipz   : in block S1, the paper's term z = (2,2,2,4,1) (row u) is replaced by the image of (2,2,2,4,2)
                (the paper-check control; breaks Lemma 3.1 for that copy).
  ctl_vzero   : in block V, the term with V-coordinates 001;20021 is replaced by 100;20000 = -(100;10000).
  ctl_append  : the term -sigma(S) is appended (the whole sequence then has sum 0; length 67).
  ctl_oneline : block S2 is glued along ha instead of hb (both copies transfer their h-coefficient to ha).
Each control is written as controls/<name>.json (labelled) and controls/<name>.txt (plain)."""
import json, os, copy
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
cert = json.load(open(os.path.join(ROOT, "certs", "c6_13_len66.json")))
terms = cert["terms"]; r = 13

def c6(b, t):
    return [(3 * b[i] + 4 * t[i]) % 6 for i in range(r)]
def bt(v):
    return [x % 2 for x in v], [x % 3 for x in v]

def write(name, ts):
    json.dump({"terms": ts}, open(os.path.join(ROOT, "controls", name + ".json"), "w"))
    with open(os.path.join(ROOT, "controls", name + ".txt"), "w") as f:
        for t in ts:
            f.write(" ".join(map(str, t["c6"])) + "\n")

paper = json.load(open(os.path.join(ROOT, "seeds", "paper_sequences_2611.24775.json")))
w = [x % 3 for x in paper["T"]["h"]]
ha, hb = [2, 0, 1, 0, 2], [1, 1, 1, 1, 2]
def s_image(t, copyno, hX):
    tb, tt = [x % 2 for x in t], [x % 3 for x in t]
    lam = tt[4]
    vpart = [(lam * c) % 3 for c in hX]
    sp = [(tt[j] - lam * w[j]) % 3 for j in range(5)]
    b = [0, 0, 0] + (tb + [0] * 5 if copyno == 1 else [0] * 5 + tb)
    tr = vpart + (sp[:4] + [0] * 4 if copyno == 1 else [0] * 4 + sp[:4])
    return c6(b, tr)

# flipz
ts = copy.deepcopy(terms)
old = s_image([2, 2, 2, 4, 1], 1, ha); new = s_image([2, 2, 2, 4, 2], 1, ha)
i = next(k for k, t in enumerate(ts) if t["block"] == "S1" and t["c6"] == old); ts[i]["c6"] = new
write("ctl_flipz", ts)
# vzero
ts = copy.deepcopy(terms)
old = c6([0, 0, 1] + [0] * 10, [2, 0, 0, 2, 1] + [0] * 8); new = c6([1, 0, 0] + [0] * 10, [2, 0, 0, 0, 0] + [0] * 8)
i = next(k for k, t in enumerate(ts) if t["block"] == "V" and t["c6"] == old); ts[i]["c6"] = new
write("ctl_vzero", ts)
# append -sigma
ts = copy.deepcopy(terms)
s = [sum(t["c6"][i] for t in terms) % 6 for i in range(r)]
ts.append({"block": "V", "c6": [(-x) % 6 for x in s]})
write("ctl_append", ts)
# oneline: S2 glued along ha
ts = [t for t in copy.deepcopy(terms) if t["block"] != "S2"]
T = []
for i, m in enumerate(paper["T"]["basis_multiplicities"]):
    T += [[1 if j == i else 0 for j in range(5)]] * m
for row in paper["T"]["rows"]:
    z, d = row["z"], row["d"]; T.append(list(z)); T.append([(z[j] + d[j]) % 6 for j in range(5)])
ts += [{"block": "S2", "c6": s_image(t, 2, ha)} for t in T]
write("ctl_oneline", ts)
# sanity: the unmodified S2 block equals the certificate's
assert [t["c6"] for t in terms if t["block"] == "S2"] == [s_image(t, 2, hb) for t in T]
print("controls written")
