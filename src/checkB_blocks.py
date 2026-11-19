"""Program B (Python + numpy; written separately from Program A, different method).

Decides whether a labelled sequence over C_6^r is zero-sumfree by exact enumeration over blocks, without any lemma:
  * every term x in C_6^r is split into x mod 2 (in F_2^r) and x mod 3 (in F_3^r): 2r coordinates;
  * for each block label, the coordinates where some term of that block is nonzero are computed; a coordinate is
    PRIVATE to a block when no other block touches it; the remaining coordinates are SHARED;
  * a nonempty sub-multiset R = R_1 u ... u R_m (R_i from block i) has sum 0 iff for every block the private part of
    sigma(R_i) is 0 and the shared parts add up to 0;
  * for each block all sub-multisets are enumerated (numpy broadcasting over multiplicity choices); we keep the
    shared parts of those with zero private part, with a flag 'empty / nonempty';
  * finally every combination of one kept shared part per block is tested for sum 0 (via a set lookup on the last
    block), excluding the all-empty combination.
Input: JSON with "terms": [{"block": label, "c6": [r integers]}, ...].
Usage: python3 checkB_blocks.py cert.json"""
import json, sys, itertools
import numpy as np

cert = json.load(open(sys.argv[1]))
terms = cert["terms"]
r = len(terms[0]["c6"])
mods = np.array([2] * r + [3] * r, dtype=np.int64)
def split(v):
    return [x % 2 for x in v] + [x % 3 for x in v]
blocks = {}
for t in terms:
    assert len(t["c6"]) == r and all(0 <= x <= 5 for x in t["c6"])
    blocks.setdefault(t["block"], []).append(tuple(split(t["c6"])))
labels = sorted(blocks)
touch = {b: {j for v in blocks[b] for j in range(2 * r) if v[j]} for b in labels}
private = {b: sorted(touch[b] - set().union(*[touch[c] for c in labels if c != b])) for b in labels}
shared = sorted(set(range(2 * r)) - set().union(*[set(private[b]) for b in labels]))
print("checkB: r=%d, L=%d, blocks %s" % (r, len(terms), {b: len(blocks[b]) for b in labels}))
print("        private coordinates per block %s; shared %d" % ({b: len(private[b]) for b in labels}, len(shared)))

def enumerate_block(vecs):
    """all sub-multisets: returns (sums array [n, 2r], size array [n])"""
    distinct = {}
    for v in vecs:
        distinct[v] = distinct.get(v, 0) + 1
    sums = np.zeros((1, 2 * r), dtype=np.int64); size = np.zeros(1, dtype=np.int64)
    for v, m in distinct.items():
        va = np.array(v, dtype=np.int64)
        ch = np.array([(k * va) % mods for k in range(m + 1)])
        sums = ((sums[:, None, :] + ch[None, :, :]) % mods).reshape(-1, 2 * r)
        size = (size[:, None] + np.arange(m + 1)[None, :]).reshape(-1)
    return sums, size

kept = {}
for b in labels:
    sums, size = enumerate_block(blocks[b])
    pz = np.all(sums[:, private[b]] == 0, axis=1) if private[b] else np.ones(len(sums), bool)
    d = {}
    for row, sz in zip(sums[pz][:, shared], size[pz]):
        key = tuple(int(x) for x in row)
        d.setdefault(key, set()).add(sz > 0)
    kept[b] = d
    print("        block %s: %d sub-multisets enumerated, %d with zero private part, %d distinct shared parts"
          % (b, len(sums), int(pz.sum()), len(d)))

smod = mods[shared]
last = labels[-1]
found = None
for combo in itertools.product(*[list(kept[b].items()) for b in labels[:-1]]):
    acc = np.zeros(len(shared), dtype=np.int64); nonempty_possible = [False]
    flags = []
    for key, fl in combo:
        acc = (acc + np.array(key, dtype=np.int64)) % smod
        flags.append(fl)
    need = tuple(int(x) for x in ((-acc) % smod))
    if need not in kept[last]:
        continue
    fl_last = kept[last][need]
    # some choice of flags (one per block) must contain a nonempty part
    if any(True in fl for fl in flags) or True in fl_last:
        found = (combo, need); break
if found:
    print("checkB: NOT zero-sumfree (a nonempty zero-sum sub-multiset exists)")
    sys.exit(1)
print("checkB: ZERO-SUMFREE (no nonempty sub-multiset of the %d terms has sum 0)" % len(terms))
