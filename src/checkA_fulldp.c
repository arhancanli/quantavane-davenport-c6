/* checkA_fulldp: Program A.  Decides whether a sequence over C_6^r is zero-sumfree by a subset-sum DP over the
   WHOLE group, written as F_2^r x F_3^r (x_i mod 2, x_i mod 3; CRT).  The reachable set is a bitset: one block of 2^r
   bits for each ternary vector (base-3 index).  Adding a term g: R' = R | (R + g); before that, if -g is already in R
   (R contains the empty sum 0) then some nonempty subsequence has sum 0.
   Input: a text file, one term per line, r integers in 0..5.  Prints the verdict and |Sigma(S)|.
   Usage: checkA_fulldp r file           (r <= 13 keeps memory at 2 x 1.6 GB)                                       */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static int r, NBW; static long NT;
static uint64_t perm_word(uint64_t x, int gl) {
  if (gl & 1) x = ((x & 0x5555555555555555ULL) << 1) | ((x >> 1) & 0x5555555555555555ULL);
  if (gl & 2) x = ((x & 0x3333333333333333ULL) << 2) | ((x >> 2) & 0x3333333333333333ULL);
  if (gl & 4) x = ((x & 0x0F0F0F0F0F0F0F0FULL) << 4) | ((x >> 4) & 0x0F0F0F0F0F0F0F0FULL);
  if (gl & 8) x = ((x & 0x00FF00FF00FF00FFULL) << 8) | ((x >> 8) & 0x00FF00FF00FF00FFULL);
  if (gl & 16) x = ((x & 0x0000FFFF0000FFFFULL) << 16) | ((x >> 16) & 0x0000FFFF0000FFFFULL);
  if (gl & 32) x = (x << 32) | (x >> 32);
  return x;
}
static long tern_add(long t, const int *gd) { /* base-3 digitwise addition */
  long res = 0, m = 1;
  for (int i = 0; i < r; i++) { int d = (int)(t % 3); t /= 3; res += ((d + gd[i]) % 3) * m; m *= 3; }
  return res;
}
static int getbit(const uint64_t *R, long t, int b) { return (R[t * NBW + (b >> 6)] >> (b & 63)) & 1; }

int main(int argc, char **argv) {
  if (argc < 3) { fprintf(stderr, "usage: checkA_fulldp r file\n"); return 2; }
  r = atoi(argv[1]);
  if (r < 6 || r > 13) { fprintf(stderr, "this program needs 6 <= r <= 13\n"); return 2; }
  NBW = (1 << r) / 64; NT = 1; for (int i = 0; i < r; i++) NT *= 3;
  FILE *f = fopen(argv[2], "r"); if (!f) { perror("input"); return 2; }
  int seq[512][16], L = 0;
  while (L < 512) { int ok = 1; for (int i = 0; i < r; i++) if (fscanf(f, "%d", &seq[L][i]) != 1) { ok = 0; break; } if (!ok) break;
    for (int i = 0; i < r; i++) if (seq[L][i] < 0 || seq[L][i] > 5) { fprintf(stderr, "bad entry\n"); return 2; } L++; }
  fclose(f);
  size_t words = (size_t)NT * NBW;
  uint64_t *R = calloc(words, 8), *Q = malloc(words * 8);
  if (!R || !Q) { fprintf(stderr, "out of memory\n"); return 2; }
  R[0] = 1; /* empty sum */
  long *shift = malloc(sizeof(long) * NT);
  int zero_hit = -1;
  for (int k = 0; k < L; k++) {
    int gb = 0, gd[16], nd[16];
    for (int i = 0; i < r; i++) { gb |= (seq[k][i] % 2) << i; gd[i] = seq[k][i] % 3; nd[i] = (3 - gd[i]) % 3; }
    if (gb == 0) { int z = 1; for (int i = 0; i < r; i++) if (gd[i]) z = 0; if (z) { zero_hit = k; break; } }
    /* is -g (binary gb, ternary nd) already reachable? */
    long tn = 0, m = 1; for (int i = 0; i < r; i++) { tn += nd[i] * m; m *= 3; }
    if (getbit(R, tn, gb)) { zero_hit = k; break; }
    for (long t = 0; t < NT; t++) shift[t] = tern_add(t, gd);
    memcpy(Q, R, words * 8);
    int gl = gb & 63, gh = gb >> 6;
    for (long t = 0; t < NT; t++) {
      const uint64_t *src = R + t * NBW; uint64_t *dst = Q + shift[t] * NBW;
      for (int w = 0; w < NBW; w++) dst[w ^ gh] |= perm_word(src[w], gl);
    }
    uint64_t *tmp = R; R = Q; Q = tmp;
    fprintf(stderr, "term %d/%d added\n", k + 1, L);
  }
  if (zero_hit >= 0) { printf("checkA: r=%d L=%d NOT zero-sumfree (a nonempty zero sum appears when term %d is added)\n", r, L, zero_hit + 1); return 1; }
  long long cnt = 0; for (size_t i = 0; i < words; i++) cnt += __builtin_popcountll(R[i]);
  long long order = (long long)NT << r;
  printf("checkA: r=%d L=%d ZERO-SUMFREE; distinct nonempty subsums %lld of %lld group elements (%lld nonzero elements missed)\n",
         r, L, cnt - 1, order, order - cnt);
  return 0;
}
