"""Referee's red controls: tampered copies of the witness, each with an explicit zero-sum checked here by plain
integer arithmetic. Both checkers must refuse every one.
Run: python3 -I make_controls.py WITNESS OUTDIR
"""
import os
import random
import sys


def load(fn):
    with open(fn) as f:
        return [[int(t) for t in line.split()] for line in f if line.strip()]


def save(fn, seq):
    with open(fn, "w") as f:
        for g in seq:
            f.write(" ".join(map(str, g)) + "\n")


def neg(g):
    return [(-v) % 6 for v in g]


def add(*gs):
    s = [0] * len(gs[0])
    for g in gs:
        s = [(x + y) % 6 for x, y in zip(s, g)]
    return s


def check(name, seq, zs_lines):
    s = add(*[seq[i - 1] for i in zs_lines])
    nz = [i for i in zs_lines if any(seq[i - 1])]
    assert all(v == 0 for v in s), (name, s)
    assert len(nz) == len(zs_lines), (name, "a term of the zero sum is 0")
    print(f"{name}: length {len(seq)}; explicit zero-sum on lines {zs_lines} (sum {s}); all its terms nonzero")


def main():
    w = load(sys.argv[1])
    out = sys.argv[2]
    os.makedirs(out, exist_ok=True)
    # R1: change ONE coordinate of ONE term: line 9 (5,4,1,0..) -> (5,4,0,0..); then lines 1,4,5,9 sum to 0
    c = [list(g) for g in w]
    assert c[8][:3] == [5, 4, 1]
    c[8][2] = 0
    check("R1_onecoord", c, [1, 4, 5, 9])
    save(os.path.join(out, "R1_onecoord.txt"), c)
    # R2: replace line 66 by minus the sum of five others from different blocks (lines 1, 20, 40, 50, 60)
    c = [list(g) for g in w]
    c[65] = neg(add(w[0], w[19], w[39], w[49], w[59]))
    check("R2_neg5", c, [1, 20, 40, 50, 60, 66])
    save(os.path.join(out, "R2_neg5.txt"), c)
    # R3: duplicate a term's inverse: line 66 replaced by -(line 35)
    c = [list(g) for g in w]
    c[65] = neg(w[34])
    check("R3_inv_replace", c, [35, 66])
    save(os.path.join(out, "R3_inv_replace.txt"), c)
    # R4: append the inverse of line 17 (length 67)
    c = [list(g) for g in w] + [neg(w[16])]
    check("R4_inv_append", c, [17, 67])
    save(os.path.join(out, "R4_inv_append.txt"), c)
    # R5: line 66 replaced by -sigma(lines 1..65): the whole sequence is a zero sum (hardest: needs all 66 terms)
    c = [list(g) for g in w]
    c[65] = neg(add(*w[:65]))
    check("R5_negsum65", c, list(range(1, 67)))
    save(os.path.join(out, "R5_negsum65.txt"), c)
    # R6: append -sigma(S) (length 67): the whole sequence is a zero sum
    c = [list(g) for g in w] + [neg(add(*w))]
    check("R6_append_negsum", c, list(range(1, 68)))
    save(os.path.join(out, "R6_append_negsum.txt"), c)
    # R7: multiplicity control: three more copies of e1 (lines 1-3 are e1) -> e1^6 is a zero sum
    c = [list(g) for g in w] + [list(w[0])] * 3
    check("R7_e1_six", c, [1, 2, 3, 67, 68, 69])
    save(os.path.join(out, "R7_e1_six.txt"), c)
    # I1..I5: informational random single-coordinate changes (no explicit zero sum known in advance)
    rng = random.Random(1119)
    for t in range(1, 6):
        c = [list(g) for g in w]
        i = rng.randrange(66)
        k = rng.randrange(13)
        old = c[i][k]
        c[i][k] = rng.choice([v for v in range(6) if v != old])
        print(f"I{t}_random: line {i + 1} coord {k + 1}: {old} -> {c[i][k]}")
        save(os.path.join(out, f"I{t}_random.txt"), c)


if __name__ == "__main__":
    main()
