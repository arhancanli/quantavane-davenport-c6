"""Build the length-66 sequence over C_6^13 (claimed zero-sumfree, so D(C_6^13) >= 67 > 5*13 + 1).

Ingredients
  V  : a zero-sumfree sequence of length 20 over H_{3,5} = F_2^3 + F_3^5 with two terms of order 3, ha and hb.
       Found by our k-exchange search (src/kx.c) from the paper's U (arXiv 2611.24775, Lemma 3.2): first U' = U minus
       one term (100;10000) plus -sigma(U) = (000;20102) (removing a term of the minimal zero-sum sequence U.(-sigma U)),
       then 4 terms of U' exchanged for 4 others (one of them hb = (000;11112)).
  S  : the paper's sequence T.h over C_6^5 (Lemma 3.1), length 25, zero-sumfree, h = (0,0,2,0,4) of order 3. Used twice.
Gluing (our lemma; the computation below does not depend on it, both checkers test the final sequence directly):
  glue V and S along ha ~ h, then the result and a second S along hb ~ h.  Ternary coordinate 5 of each S copy is
  eliminated: x = (xV, xS) -> (xV + xS_5 * hX, xS - xS_5 * w), w = ternary part of h (w_5 = 1), hX = ha or hb.
Output coordinates
  binary  : V(3), S1(5), S2(5);  ternary : V(5), S1(1..4), S2(1..4);  C_6 coordinate i = CRT(binary_i, ternary_i)
  = (3*b_i + 4*t_i) mod 6.
Writes certs/c6_13_len66.txt (one term per line, 13 integers in 0..5) and certs/c6_13_len66.json (with block labels)."""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

V = """100;10000 100;10000 100;10000 010;01000 010;01000 010;01000 010;01000 011;00100 101;21100 101;21100
001;00010 111;12110 111;12110 101;00001 000;20102 011;11001 001;20021 000;11112 111;01212 001;01111""".split()
# the line above is solution 1 of out/kx_Uprime_k4.txt; checked again below
paper = json.load(open(os.path.join(ROOT, "seeds", "paper_sequences_2611.24775.json")))

def parse(s):
    b, t = s.split(";")
    return [int(c) for c in b], [int(c) for c in t]

Vp = [parse(s) for s in V]
ha, hb = parse("000;20102"), parse("000;11112")
assert Vp.count(ha) == 1 and Vp.count(hb) == 1

# the paper's S = T.h over C_6^5
Tsp = paper["T"]
T = []
for i, m in enumerate(Tsp["basis_multiplicities"]):
    T += [[1 if j == i else 0 for j in range(5)]] * m
for r in Tsp["rows"]:
    z, d = r["z"], r["d"]
    T.append(list(z)); T.append([(z[j] + d[j]) % 6 for j in range(5)])
assert len(T) == 24
h = Tsp["h"]; w = [x % 3 for x in h]; assert [x % 2 for x in h] == [0] * 5 and w[4] == 1

out = []  # (label, binary 13, ternary 13)
for vb, vt in Vp:
    if (vb, vt) in (ha, hb):
        continue
    out.append(("V", vb + [0] * 10, vt + [0] * 8))
for copy, hX in ((1, ha[1]), (2, hb[1])):
    for t in T:
        tb, tt = [x % 2 for x in t], [x % 3 for x in t]
        lam = tt[4]
        vpart = [(lam * c) % 3 for c in hX]
        spart = [(tt[j] - lam * w[j]) % 3 for j in range(5)]
        assert spart[4] == 0
        sbin = [0, 0, 0] + (tb + [0] * 5 if copy == 1 else [0] * 5 + tb)
        ster = vpart + (spart[:4] + [0] * 4 if copy == 1 else [0] * 4 + spart[:4])
        out.append(("S%d" % copy, sbin, ster))
assert len(out) == 66
c6 = [[(3 * b[i] + 4 * t[i]) % 6 for i in range(13)] for _, b, t in out]
os.makedirs(os.path.join(ROOT, "certs"), exist_ok=True)
with open(os.path.join(ROOT, "certs", "c6_13_len66.txt"), "w") as f:
    for v in c6:
        f.write(" ".join(map(str, v)) + "\n")
json.dump({"group": "C_6^13", "length": 66, "claim": "zero-sumfree, hence D(C_6^13) >= 67 = 5*13 + 2",
           "coordinates": "C_6 coordinate i = (3*binary_i + 4*ternary_i) mod 6; binary order V(3),S1(5),S2(5); "
                          "ternary order V(5),S1(1..4),S2(1..4)",
           "terms": [{"block": lab, "c6": v} for (lab, _, _), v in zip(out, c6)],
           "ingredients": {"V_over_H35": V, "ha": "000;20102", "hb": "000;11112",
                           "S_over_C6_5": [h] + T}},
          open(os.path.join(ROOT, "certs", "c6_13_len66.json"), "w"), indent=0)
print("wrote", len(c6), "terms over C_6^13")
