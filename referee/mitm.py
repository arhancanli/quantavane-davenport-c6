"""Second, structurally different zero-sum check (referee): meet in the middle on a coordinate split.

Valid for ANY sequence S over (Z/6)^r. Pick a set P of coordinates ("tail", default the last 4) and split
S = A . B, where A = terms vanishing on P (checked here, not assumed), B = the rest. A nonempty zero-sum Z = X . Y
(X in A, Y in B) needs sigma(Y)|_P = 0 (since sigma(X)|_P = 0) and sigma(X) = -sigma(Y) on the other coordinates Q.
So: S has a nonempty zero sum iff
  (i)  A has a nonempty zero sum (DP over (Z/6)^|Q| with 'last term' test), or
  (ii) some NONEMPTY sub-multiset Y of B has sigma(Y)|_P = 0 and -sigma(Y)|_Q in Sigma(A) (empty sum allowed).
Sigma(A) is a boolean array over (Z/6)^|Q| built with numpy.roll; Y runs over all count vectors of B's distinct
elements. When a zero sum is found it is reconstructed from the stored DP layers and re-added with plain integers.

Usage: python3 -I mitm.py FILE [--tail 9,10,11,12]
"""
import itertools
import sys

import numpy as np


def load(fn):
    seq = []
    with open(fn) as f:
        for ln, line in enumerate(f, 1):
            toks = line.split()
            if not toks:
                raise SystemExit(f"blank line {ln}")
            for t in toks:
                if t not in {"0", "1", "2", "3", "4", "5"}:
                    raise SystemExit(f"bad token {t!r} line {ln}")
            seq.append(tuple(int(t) for t in toks))
    r = len(seq[0])
    assert all(len(g) == r for g in seq), "ragged"
    return seq, r


def main():
    fn = sys.argv[1]
    tail = [9, 10, 11, 12]
    if "--tail" in sys.argv:
        tail = [int(x) for x in sys.argv[sys.argv.index("--tail") + 1].split(",")]
    seq, r = load(fn)
    Q = [k for k in range(r) if k not in tail]
    P = tail
    A_idx = [i for i, g in enumerate(seq) if all(g[k] == 0 for k in P)]
    B_idx = [i for i, g in enumerate(seq) if i not in set(A_idx)]
    print(f"file {fn}: length {len(seq)}, rank {r}; tail P={P}; |A|={len(A_idx)} (terms vanishing on P), |B|={len(B_idx)}")
    zero_terms = [i + 1 for i, g in enumerate(seq) if all(v == 0 for v in g)]
    print(f"zero terms: {zero_terms}")
    q = len(Q)
    shape = (6,) * q
    layers = [np.zeros(shape, dtype=bool)]
    layers[0][(0,) * q] = True
    found = None
    for j, i in enumerate(A_idx):
        a = tuple(seq[i][k] for k in Q)
        prev = layers[-1]
        neg = tuple((-v) % 6 for v in a)
        if prev[neg] and found is None:
            found = ("A", j, neg)
        layers.append(prev | np.roll(prev, shift=a, axis=tuple(range(q))))
    SigA = layers[-1]
    print(f"|Sigma(A)| incl. empty = {int(SigA.sum())} of 6^{q} = {6 ** q}")

    def reconstruct(target, upto):
        """indices (into A_idx) of a sub-multiset of the first `upto` A-terms with sum target on Q"""
        v = tuple(target)
        chosen = []
        for j in range(upto, 0, -1):
            if layers[j - 1][v]:
                continue
            a = tuple(seq[A_idx[j - 1]][k] for k in Q)
            v = tuple((x - y) % 6 for x, y in zip(v, a))
            chosen.append(j - 1)
            assert layers[j - 1][v]
        assert v == (0,) * q
        return chosen

    witness = None
    if found is not None:
        _, j, neg = found
        X = reconstruct(neg, j) + [j]
        witness = [A_idx[t] for t in X]
        print("case (i): A has a nonempty zero sum")
    else:
        print("case (i): A is zero-sumfree")
        # enumerate sub-multisets of B
        distinct = sorted(set(seq[i] for i in B_idx))
        mult = [sum(1 for i in B_idx if seq[i] == d) for d in distinct]
        total = int(np.prod([m + 1 for m in mult]))
        print(f"B: {len(distinct)} distinct elements, multiplicities {mult}; sub-multisets to scan: {total}")
        Dm = np.array(distinct, dtype=np.int64)  # (nd, r)
        counts = np.array(list(itertools.product(*[range(m + 1) for m in mult])), dtype=np.int64)  # (total, nd)
        sig = (counts @ Dm) % 6  # (total, r)
        nonempty = counts.sum(axis=1) > 0
        okP = np.all(sig[:, P] == 0, axis=1) & nonempty
        cand = np.nonzero(okP)[0]
        print(f"nonempty Y with sigma(Y) = 0 on P: {len(cand)}")
        negQ = (-sig[cand][:, Q]) % 6
        hits = SigA[tuple(negQ[:, t] for t in range(q))]
        nh = int(hits.sum())
        print(f"of those, -sigma(Y)|Q in Sigma(A): {nh}")
        if nh:
            c = cand[np.nonzero(hits)[0][0]]
            Y = []
            for d, cnt in zip(distinct, counts[c]):
                pool = [i for i in B_idx if seq[i] == d]
                Y += pool[: int(cnt)]
            X = reconstruct(tuple(int(v) for v in negQ[np.nonzero(hits)[0][0]]), len(A_idx))
            witness = [A_idx[t] for t in X] + Y
    if witness is None:
        print("VERDICT (mitm): ZERO-SUMFREE")
        return 0
    s = [0] * r
    for i in witness:
        s = [(x + y) % 6 for x, y in zip(s, seq[i])]
    assert len(witness) == len(set(witness)) and len(witness) > 0
    print(f"explicit zero-sum subsequence: {len(witness)} terms, lines {sorted(i + 1 for i in witness)}; sum = {s}")
    assert all(v == 0 for v in s)
    print("VERDICT (mitm): NOT zero-sumfree")
    return 1


if __name__ == "__main__":
    sys.exit(main())
