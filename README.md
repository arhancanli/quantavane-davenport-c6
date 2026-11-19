# D(C₆¹³) ≥ 67: the van Emde Boas–Kruyswijk conjecture fails for C₆ʳ from rank 13

**Status: verified by two independent programs and an internal referee; not yet peer reviewed.**

Author: Arhan Canli (Quantavane, arhan@quantavane.com). On the use of AI, see
[Authorship and the role of AI](#authorship-and-the-role-of-ai).

## The result

For a finite abelian group G, the Davenport constant D(G) is the least ℓ such that every sequence of ℓ elements of G
has a nonempty subsequence with sum 0. Write d(G) = D(G) − 1 for the largest length of a zero-sumfree sequence. In
1969 van Emde Boas and Kruyswijk conjectured that D(C_n^r) = r(n − 1) + 1 for all n and r. For n = 6 this says
D(C₆ʳ) = 5r + 1.

**Theorem.** D(C₆¹³) ≥ 67 = 5·13 + 2. Consequently D(C₆ʳ) > 5r + 1 for every r ≥ 13, and
D(C₆ʳ) ≥ 5r + ⌊r/13⌋ + 1 for every r.

The witness is the explicit sequence of 66 elements of C₆¹³ in `certs/c6_13_len66.txt` (one element per line, 13
coordinates in 0..5). It has no nonempty zero-sum subsequence, and every nonzero element of C₆¹³ is the sum of one of
its subsequences, so it cannot be lengthened. The consequences follow by direct sums: for r = 13q + s, take q copies
of the witness and s further coordinates carrying e_i⁵ each.

## The earlier work this builds on

Everything here is built on W. Gao, A. Geroldinger, G. Wang and X. Zeng, *A disproof of the van Emde Boas–Kruyswijk
conjecture on the Davenport constant*, [arXiv:2611.24775](https://arxiv.org/abs/2611.24775) (v1 submitted 18 November
2026), which proved D(C₆ʳ) ≥ 5r + ⌊r/23⌋ + 1, so D(C₆ʳ) > 5r + 1 for r ≥ 23. Their sequences U (Lemma 3.2, length 20
over F₂³ ⊕ F₃⁵) and T·h (Lemma 3.1, length 25 over C₆⁵, h = (0,0,2,0,4)) are the building blocks of our witness, and
their gluing along a subgroup of order three is the idea we extend. The sequences are copied from their paper into
`seeds/paper_sequences_2611.24775.json`. The conjecture is due to P. van Emde Boas and D. Kruyswijk (1969); see the paper above for its history and the
full references.

## How the witness is built

Notation as in the paper: H_{a,b} = F₂ᵃ ⊕ F₃ᵇ, so C₆ʳ = H_{r,r}.

1. **Gluing lemma.** Let S₁ over G₁ and S₂ over G₂ be zero-sumfree, and let h₁ ∈ S₁, h₂ ∈ S₂ have order 3. Then
   (S₁h₁⁻¹)(S₂h₂⁻¹), mapped into (G₁ ⊕ G₂)/⟨(h₁, −h₂)⟩, is zero-sumfree. *Proof.* A zero sum R₁R₂ in the quotient
   means σ(R₁) = j·h₁ and σ(R₂) = −j·h₂ for some j ∈ {0,1,2}. If R₁ is nonempty, j = 0 gives a zero sum in S₁ and
   j = 2 gives the zero sum R₁·h₁ in S₁; so j = 1, and then R₂ empty gives −h₂ = 0, while R₂ nonempty gives the zero
   sum R₂·h₂ in S₂. If R₁ is empty, j = 0 and R₂ is a zero sum in S₂. ∎ The paper's Lemma 3.1 with a shared subgroup
   of order three is the case S₁ = S₂.
2. **An order-3 term in U.** U is zero-sumfree, so U·(−σ(U)) is a minimal zero-sum sequence, and deleting any one of
   its terms leaves a zero-sumfree sequence. Here −σ(U) = (000;20102) has order 3, so U′ = U with (100;10000) replaced
   by (000;20102) is zero-sumfree of length 20 over H_{3,5} with one order-3 term.
3. **A second order-3 term.** An exchange search (`src/kx.c`: delete k terms and list every set of k new terms whose
   nonempty subsums avoid −σ of the rest) found at k = 4 nine zero-sumfree length-20 sequences over H_{3,5} with two
   order-3 terms (`results/kx_Uprime_k4.txt`). We use the first, V, with ha = (000;20102) and hb = (000;11112).
4. **Gluing.** Glue V to a copy of T·h along ha, then to a second copy along hb. The length is 20 + 25 + 25 − 4 = 66,
   over H_{3+5+5, 5+5+5−2} = H_{13,13} = C₆¹³, and 66 = 5·13 + 1. `src/build_c6_13.py` does this and writes the
   certificate.

The certificate does not depend on the lemma: every check below tests the final 66 terms directly.

## Checking it: one command

    sh check.sh              # about 3 minutes, 3.3 GB of memory; needs cc and Python 3 with numpy
    sh check.sh --referee    # also runs the referee's programs and controls

`check.sh` checks the file hashes; rebuilds the 66 terms from the paper's sequences and compares them byte for byte;
runs Program A and Program B on the certificate; and runs both programs on four red controls, each of which must be
refused. It ends with `ALL CHECKS PASSED` or stops at the first failure.

- **Program A** (`src/checkA_fulldp.c`, C): a subset-sum computation over the whole group C₆¹³ = F₂¹³ × F₃¹³
  (13,060,694,016 elements, as bitsets). A sequence is zero-sumfree exactly when, for every i, −g_i is not a sum of a
  subsequence of g₁, …, g_{i−1} (the empty sum included). Result: zero-sumfree, and all 13,060,694,015 nonzero
  elements are subsums.
- **Program B** (`src/checkB_blocks.py`, Python with numpy, written separately with a different method): splits the
  coordinates into those private to one block (V, first T·h, second T·h) and the shared ones, enumerates every
  sub-multiset of every block (307,200 + 307,200 + 23,040), keeps the shared parts of those whose private part
  vanishes, and tests every combination. No lemma is used. Result: zero-sumfree.
- **Red controls** (`controls/`, made by `src/make_controls.py`): a term changed inside a copy of T·h, a term of V
  replaced by the negative of another, −σ appended (length 67), and both copies glued along the same term. Both
  programs refuse all four.
- **Internal referee** (`referee/`): a fresh referee wrote two programs from the definition before reading ours:
  `zsf_dp.c` (a DP over F₃¹³ cells holding F₂¹³ bitsets) and `mitm.py` (meet in the middle on a coordinate split).
  Both accept the witness and refuse twelve further controls (`referee/controls/`). Its logs are in `results/`.

Recorded outputs from our runs are in `results/` (19 November 2026).

The few files edited for this public copy (status sentences, local paths) are listed with their hashes in
`PUBLIC_COPY_CHANGES.md`.

## Not yet done

- Outside review: none yet. The authors of arXiv:2611.24775 were told of this result before this repository was made
  public.
- Lean 4: not yet; a kernel-checked version is planned and its status will be added here.
- Paper: a short paper is in preparation.
- The least r with D(C₆ʳ) > 5r + 1 is open. Our searches for smaller ranks found nothing; they are evidence only.

## Layout

| path | contents |
| --- | --- |
| `certs/` | the witness (`.txt`; `.json` with block labels V, S1, S2) and its SHA-256 |
| `seeds/` | U, T and h as printed in arXiv:2611.24775 |
| `src/` | `build_c6_13.py` (construction), `checkA_fulldp.c` (Program A), `checkB_blocks.py` (Program B), `make_controls.py`, `kx.c` (the exchange search) |
| `controls/` | our four red controls |
| `referee/` | the referee's programs and twelve controls |
| `results/` | recorded outputs |
| `tools/sha256check.py` | hash check (standard library) |

## Authorship and the role of AI

Author and maintainer: Arhan Canli (Quantavane, https://quantavane.com), who takes responsibility for the results.

The searches, the construction, the checkers and this text were written and carried out by AI agents (Claude, by
Anthropic) working for Quantavane under Arhan Canli's direction. Programs described as written separately share no
code but were written with the same AI system. No mathematician outside Quantavane has reviewed this repository; the
witness can be checked by anyone with the programs here or with their own.

## Licence and citation

The witness, controls, outputs and text are under CC BY 4.0 (`LICENSE`); the code is under the MIT licence
(`LICENSE-CODE`). Please cite as: Arhan Canli, *D(C₆¹³) ≥ 67: the van Emde Boas–Kruyswijk conjecture fails for C₆ʳ
from rank 13*, https://github.com/arhancanli/quantavane-davenport-c6 (2026), together with Gao, Geroldinger, Wang and
Zeng, arXiv:2611.24775.
