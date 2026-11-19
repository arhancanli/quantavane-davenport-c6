/* kx: k-exchange neighbourhood of a zero-sumfree sequence over H_{a,b} = F_2^a + F_3^b (search tool only).
   For every way to delete k terms (sub-multisets, frozen prefix kept) from S, let R be the rest and
   A(R) = { g != 0 : -g not in Sigma(R) u {0} }.  A sequence X can be appended to R keeping it zero-sumfree iff every
   nonempty subsum of X lies in A(R).  We enumerate all multisets X of size j with that property and report those
   for which R.X has at least 'need3' terms of order 3 (and at least 'need2' of order 2).
   Usage: kx a b k j need3 need2 frozen maxprint seqfile
   Element format: binary digits ';' ternary digits. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static int A, B, NB, NT, N, L, K, J, need3, need2, frozen; long maxprint, printed = 0, nsol = 0, nsub = 0;
static int *tsub, *tadd, *negt;
static int seq[256], del[256];
static unsigned char *reach, *tmpr;
static inline int addg(int x, int y) { return ((x & (NB - 1)) ^ (y & (NB - 1))) | (tadd[(x >> A) * NT + (y >> A)] << A); }
static int ord(int e) { int eb = e & (NB - 1), et = e >> A; if (!e) return 1; if (!et) return 2; if (!eb) return 3; return 6; }
static void printel(FILE *fo, int e) {
  int eb = e & (NB - 1), et = e >> A;
  for (int i = 0; i < A; i++) fprintf(fo, "%d", (eb >> i) & 1);
  fprintf(fo, ";");
  for (int i = 0; i < B; i++) { fprintf(fo, "%d", et % 3); et /= 3; }
}
static int parseel(const char *s) {
  int eb = 0, et = 0, i = 0, m = 1;
  for (; i < A; i++) eb |= (s[i] - '0') << i;
  i++;
  for (int j = 0; j < B; j++) { et += (s[i + j] - '0') * m; m *= 3; }
  return eb | (et << A);
}
static unsigned char *allowed; static int *alist, na;
static int X[64], sums[1 << 16], nsums;
static int c3R, c2R;
static void bt(int depth, int start) {
  if (depth == J) {
    int c3 = c3R, c2 = c2R;
    for (int i = 0; i < J; i++) { if (ord(X[i]) == 3) c3++; if (ord(X[i]) == 2) c2++; }
    if (c3 < need3 || c2 < need2) return;
    nsol++;
    if (printed < maxprint) {
      printed++;
      printf("SOL del=");
      for (int i = 0; i < L; i++) if (del[i]) { printel(stdout, seq[i]); printf(","); }
      printf(" add=");
      for (int i = 0; i < J; i++) { printel(stdout, X[i]); printf(","); }
      printf(" o3=%d o2=%d seq=", c3, c2);
      for (int i = 0; i < L; i++) if (!del[i]) { printel(stdout, seq[i]); printf(" "); }
      for (int i = 0; i < J; i++) { printel(stdout, X[i]); printf(i + 1 < J ? " " : "\n"); }
      fflush(stdout);
    }
    return;
  }
  for (int ai = start; ai < na; ai++) {
    int g = alist[ai], ok = 1, base = nsums;
    for (int s = 0; s < base; s++) { int t = addg(sums[s], g); if (!allowed[t]) { ok = 0; break; } }
    if (!ok) continue;
    for (int s = 0; s < base; s++) sums[nsums++] = addg(sums[s], g);
    X[depth] = g;
    bt(depth + 1, ai);
    nsums = base;
  }
}
static void process(void) {
  nsub++;
  memset(reach, 0, N); reach[0] = 1;
  c3R = c2R = 0;
  for (int i = 0; i < L; i++) {
    if (del[i]) continue;
    int g = seq[i]; if (ord(g) == 3) c3R++; if (ord(g) == 2) c2R++;
    memcpy(tmpr, reach, N);
    for (int x = 0; x < N; x++) if (tmpr[x]) reach[addg(x, g)] = 1;
  }
  na = 0;
  for (int g = 1; g < N; g++) { allowed[g] = !reach[negt[g]]; if (allowed[g]) alist[na++] = g; }
  allowed[0] = 0;
  /* sums holds subsums of X including the empty one (0); 0 must stay out of the nonempty part: handled since
     every new sum is checked against allowed[] and allowed[0] = 0 */
  nsums = 0; sums[nsums++] = 0;
  /* the empty sum 0 is in the list; adding g gives g, which must be allowed */
  bt(0, 0);
}
/* enumerate sub-multisets of size K among positions >= frozen; seq sorted so equal values are adjacent */
static void choose(int pos, int left) {
  if (left == 0) { process(); return; }
  if (pos >= L) return;
  /* group of equal values starting at pos */
  int end = pos; while (end < L && seq[end] == seq[pos]) end++;
  int cnt = end - pos;
  for (int t = 0; t <= cnt && t <= left; t++) {
    for (int i = 0; i < cnt; i++) del[pos + i] = i < t;
    choose(end, left - t);
  }
  for (int i = 0; i < cnt; i++) del[pos + i] = 0;
}
static int cmpint(const void *x, const void *y) { return *(const int *)x - *(const int *)y; }
int main(int argc, char **argv) {
  if (argc < 10) { fprintf(stderr, "usage: kx a b k j need3 need2 frozen maxprint seqfile\n"); return 2; }
  A = atoi(argv[1]); B = atoi(argv[2]); K = atoi(argv[3]); J = atoi(argv[4]); need3 = atoi(argv[5]); need2 = atoi(argv[6]);
  frozen = atoi(argv[7]); maxprint = atol(argv[8]);
  NB = 1 << A; NT = 1; for (int i = 0; i < B; i++) NT *= 3; N = NB * NT;
  tsub = malloc(sizeof(int) * NT * NT); tadd = malloc(sizeof(int) * NT * NT); negt = malloc(sizeof(int) * N);
  for (int x = 0; x < NT; x++) for (int y = 0; y < NT; y++) {
    int r = 0, s = 0, m = 1, xx = x, yy = y;
    for (int i = 0; i < B; i++) { r += (((xx % 3) - (yy % 3) + 3) % 3) * m; s += (((xx % 3) + (yy % 3)) % 3) * m; xx /= 3; yy /= 3; m *= 3; }
    tsub[x * NT + y] = r; tadd[x * NT + y] = s; }
  for (int e = 0; e < N; e++) negt[e] = (e & (NB - 1)) | (tsub[e >> A] << A);
  FILE *fi = fopen(argv[9], "r"); char buf[256]; L = 0;
  if (!fi) { perror("seq"); return 2; }
  while (fscanf(fi, "%255s", buf) == 1) seq[L++] = parseel(buf);
  fclose(fi);
  qsort(seq + frozen, L - frozen, sizeof(int), cmpint);
  reach = malloc(N); tmpr = malloc(N); allowed = calloc(N, 1); alist = malloc(sizeof(int) * N);
  /* sanity: S zero-sumfree */
  memset(del, 0, sizeof del);
  memset(reach, 0, N); reach[0] = 1; int bad = 0;
  for (int i = 0; i < L; i++) { memcpy(tmpr, reach, N); for (int x = 0; x < N; x++) if (tmpr[x]) { int y = addg(x, seq[i]); if (y == 0) bad = 1; reach[y] = 1; } }
  if (bad) { printf("INPUT NOT ZERO-SUMFREE\n"); return 3; }
  printf("START a=%d b=%d L=%d k=%d j=%d need3=%d need2=%d frozen=%d\n", A, B, L, K, J, need3, need2, frozen); fflush(stdout);
  choose(frozen, K);
  printf("DONE subsets=%ld solutions=%ld\n", nsub, nsol);
  return 0;
}
