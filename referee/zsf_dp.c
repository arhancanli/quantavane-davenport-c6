/*
 * zsf_dp.c -- referee's own zero-sum checker for a sequence over C_6^r (r <= 13).
 *
 * Written from the definition only. A sequence S = g_1 ... g_L over G = (Z/6)^r is zero-sumfree iff no nonempty
 * sub-multiset sums to 0. Every zero-sum subsequence has a last term g_i (in input order) and the rest lies in
 * g_1..g_{i-1}; so S is zero-sumfree iff for every i, -g_i is NOT in Sigma_{i-1} := { sigma(X) : X subseq of
 * g_1..g_{i-1}, X may be empty }. (A zero term g_i = 0 is caught because 0 = empty sum is in Sigma_{i-1}.)
 *
 * Method (not a flat bitset over G): use the CRT isomorphism Z/6 -> Z/2 x Z/3, x -> (x mod 2, x mod 3), so
 * G = F_2^r x F_3^r. The state is indexed by the F_3^r part t (3^r cells); each cell holds the set of F_2^r parts u
 * such that (u,t) is in Sigma, as a 2^r-bit bitset. Adding g = (b, s): Sigma <- Sigma u (Sigma + g), i.e.
 *     new[t][u] = old[t][u] | old[t - s][u xor b].
 * Done IN PLACE along the orbits {x, x+s, x+2s} of translation by s (orbits of size 3 if s != 0, size 1 if s = 0),
 * with the xor-translation u -> u xor b applied as a word permutation (high bits of b) and an in-word bit permutation
 * by swaps (low 6 bits of b). Up to 4 threads split the orbit representatives.
 *
 * Usage: zsf_dp FILE [--full] [--threads N]
 *   FILE: one term per line, r entries per line, each a single digit 0..5, separated by spaces.
 *   --full: do not stop at the first zero-sum; continue and report |Sigma(S)| (distinct subsums incl. empty).
 * Output: validation facts, then either "ZERO-SUM: first appears when adding term i (1-based)" or "ZERO-SUMFREE",
 * and (if zero-sumfree, or with --full) the count of distinct subsums.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <pthread.h>
#include <time.h>

#define MAXR 13
#define MAXL 4096

static int L = 0, R = -1;
static int seq[MAXL][MAXR];

static uint64_t *arr;       /* 3^R cells x W words */
static uint32_t *perm;      /* perm[x] = x + s (digitwise mod 3) */
static size_t N3, W;
static uint32_t cur_bhi;    /* high part of b: word xor */
static uint32_t cur_blo;    /* low 6 bits of b: in-word xor */
static int cur_s_zero;
static int cur_k;           /* a coordinate with s_k != 0 */
static size_t pow3[MAXR + 1];

static const uint64_t MASKS[6] = {
    0x5555555555555555ULL, 0x3333333333333333ULL, 0x0F0F0F0F0F0F0F0FULL,
    0x00FF00FF00FF00FFULL, 0x0000FFFF0000FFFFULL, 0x00000000FFFFFFFFULL};

static inline uint64_t bitxor_perm(uint64_t x, uint32_t blo) {
    /* returns y with y[j] = x[j xor blo] */
    for (int k = 0; k < 6; k++) {
        if (blo & (1u << k)) {
            int s = 1 << k;
            x = ((x & MASKS[k]) << s) | ((x >> s) & MASKS[k]);
        }
    }
    return x;
}

/* dst = sh(src): dst[w] = bitxor_perm(src[w xor bhi], blo) */
static inline void shift_into(uint64_t *dst, const uint64_t *src) {
    for (size_t w = 0; w < W; w++) dst[w] = bitxor_perm(src[w ^ cur_bhi], cur_blo);
}

typedef struct { size_t lo, hi; } range_t;

static void *worker(void *p) {
    range_t *rg = (range_t *)p;
    uint64_t *tA = malloc(W * 8), *tB = malloc(W * 8), *tC = malloc(W * 8);
    if (!tA || !tB || !tC) { fprintf(stderr, "malloc failed\n"); exit(2); }
    for (size_t x = rg->lo; x < rg->hi; x++) {
        if (cur_s_zero) {
            uint64_t *A = arr + x * W;
            shift_into(tA, A);
            for (size_t w = 0; w < W; w++) A[w] |= tA[w];
        } else {
            if ((x / pow3[cur_k]) % 3 != 0) continue;   /* orbit representative: digit k == 0 */
            size_t y = perm[x], z = perm[y];
            uint64_t *A = arr + x * W, *B = arr + y * W, *C = arr + z * W;
            shift_into(tA, A); shift_into(tB, B); shift_into(tC, C);
            /* new[y] = B | sh(A); new[z] = C | sh(B); new[x] = A | sh(C) */
            for (size_t w = 0; w < W; w++) { B[w] |= tA[w]; C[w] |= tB[w]; A[w] |= tC[w]; }
        }
    }
    free(tA); free(tB); free(tC);
    return NULL;
}

static int parse_file(const char *fn) {
    FILE *f = fopen(fn, "r");
    if (!f) { perror(fn); return -1; }
    char line[1 << 16];
    int lineno = 0;
    while (fgets(line, sizeof line, f)) {
        lineno++;
        size_t n = strlen(line);
        if (n == sizeof line - 1 && line[n - 1] != '\n') { fprintf(stderr, "line %d too long\n", lineno); return -1; }
        int vals[64], nv = 0;
        char *p = line;
        int blank = 1;
        while (*p) {
            while (*p == ' ' || *p == '\t' || *p == '\r' || *p == '\n') p++;
            if (!*p) break;
            blank = 0;
            char *q = p;
            while (*q && *q != ' ' && *q != '\t' && *q != '\r' && *q != '\n') q++;
            if (q - p != 1 || p[0] < '0' || p[0] > '5') {
                fprintf(stderr, "line %d: bad token '%.*s' (must be a single digit 0..5)\n", lineno, (int)(q - p), p);
                return -1;
            }
            if (nv >= 64) { fprintf(stderr, "line %d: too many entries\n", lineno); return -1; }
            vals[nv++] = p[0] - '0';
            p = q;
        }
        if (blank) { fprintf(stderr, "line %d: blank line (rejected)\n", lineno); return -1; }
        if (R < 0) R = nv;
        if (nv != R) { fprintf(stderr, "line %d: %d entries, expected %d\n", lineno, nv, R); return -1; }
        if (R > MAXR || R < 1) { fprintf(stderr, "rank %d unsupported\n", R); return -1; }
        if (L >= MAXL) { fprintf(stderr, "too many terms\n"); return -1; }
        for (int i = 0; i < R; i++) seq[L][i] = vals[i];
        L++;
    }
    fclose(f);
    return 0;
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s FILE [--full] [--threads N]\n", argv[0]); return 2; }
    int full = 0, nthreads = 1;
    for (int i = 2; i < argc; i++) {
        if (!strcmp(argv[i], "--full")) full = 1;
        else if (!strcmp(argv[i], "--threads") && i + 1 < argc) nthreads = atoi(argv[++i]);
        else { fprintf(stderr, "unknown arg %s\n", argv[i]); return 2; }
    }
    if (nthreads < 1) nthreads = 1;
    if (nthreads > 4) nthreads = 4;
    if (parse_file(argv[1])) { printf("INPUT REJECTED\n"); return 2; }

    printf("file: %s\nterms (length): %d\nrank r: %d (all entries single digits in 0..5: yes)\n", argv[1], L, R);
    int nzero = 0;
    for (int i = 0; i < L; i++) {
        int z = 1;
        for (int k = 0; k < R; k++) if (seq[i][k]) z = 0;
        if (z) { nzero++; printf("term %d is the zero element\n", i + 1); }
    }
    printf("zero terms: %d\n", nzero);
    /* distinct elements and multiplicities */
    {
        int distinct = 0, maxmult = 0;
        for (int i = 0; i < L; i++) {
            int seen = 0, mult = 0;
            for (int j = 0; j < L; j++) {
                if (!memcmp(seq[i], seq[j], sizeof(int) * R)) { mult++; if (j < i) seen = 1; }
            }
            if (!seen) { distinct++; if (mult > maxmult) maxmult = mult; }
        }
        printf("distinct elements: %d, max multiplicity: %d\n", distinct, maxmult);
    }

    pow3[0] = 1;
    for (int i = 1; i <= R; i++) pow3[i] = pow3[i - 1] * 3;
    N3 = pow3[R];
    size_t nbits2 = (size_t)1 << R;
    W = nbits2 >= 64 ? nbits2 / 64 : 1;
    double gb = (double)N3 * W * 8 / 1e9;
    printf("state: 3^%d = %zu cells x %zu words = %.3f GB\n", R, N3, W, gb);
    if (gb > 4.0) { fprintf(stderr, "refusing > 4 GB\n"); return 2; }
    arr = calloc(N3 * W, 8);
    perm = malloc(N3 * sizeof(uint32_t));
    if (!arr || !perm) { fprintf(stderr, "alloc failed\n"); return 2; }
    arr[0] = 1;   /* empty sum: t = 0, u = 0 */

    int first_zero = 0;
    time_t t0 = time(NULL);
    for (int i = 0; i < L; i++) {
        uint32_t b = 0;
        int s[MAXR];
        size_t negs_idx = 0;
        cur_s_zero = 1; cur_k = -1;
        for (int k = 0; k < R; k++) {
            int v = seq[i][k];
            if (v % 2) b |= 1u << k;
            s[k] = v % 3;
            if (s[k]) { cur_s_zero = 0; if (cur_k < 0) cur_k = k; }
            negs_idx += (size_t)((3 - s[k]) % 3) * pow3[k];
        }
        /* -g = (b, -s); test membership of -g in Sigma_{i-1} */
        int hit = (arr[negs_idx * W + (b >> 6)] >> (b & 63)) & 1;
        if (hit && !first_zero) {
            first_zero = i + 1;
            printf("ZERO-SUM: first appears when adding term %d (1-based)\n", i + 1);
            fflush(stdout);
            if (!full) { printf("VERDICT: NOT zero-sumfree\n"); return 1; }
        }
        cur_bhi = b >> 6; cur_blo = b & 63;
        if (!cur_s_zero) {
            for (size_t x = 0; x < N3; x++) {
                size_t y = 0, xx = x;
                for (int k = 0; k < R; k++) { int d = xx % 3; xx /= 3; y += (size_t)((d + s[k]) % 3) * pow3[k]; }
                perm[x] = (uint32_t)y;
            }
        }
        pthread_t th[4];
        range_t rg[4];
        size_t chunk = (N3 + nthreads - 1) / nthreads;
        for (int tI = 0; tI < nthreads; tI++) {
            rg[tI].lo = tI * chunk;
            rg[tI].hi = (tI + 1) * chunk < N3 ? (tI + 1) * chunk : N3;
            pthread_create(&th[tI], NULL, worker, &rg[tI]);
        }
        for (int tI = 0; tI < nthreads; tI++) pthread_join(th[tI], NULL);
        if (R >= 10) { fprintf(stderr, "term %d done (%lds)\n", i + 1, (long)(time(NULL) - t0)); }
    }
    /* count distinct subsums (including the empty sum 0) */
    uint64_t mask = nbits2 >= 64 ? ~0ULL : ((1ULL << nbits2) - 1);
    unsigned long long cnt = 0;
    for (size_t j = 0; j < N3 * W; j++) cnt += __builtin_popcountll(arr[j] & mask);
    unsigned long long order = (unsigned long long)N3 * nbits2;
    printf("group order 6^%d = %llu\n", R, order);
    printf("distinct subsums incl. empty sum: %llu\n", cnt);
    if (!first_zero) {
        printf("distinct NONEMPTY subsums (all nonzero, since zero-sumfree): %llu\n", cnt - 1);
        printf("every nonzero element is a nonempty subsum: %s\n", cnt == order ? "yes" : "no");
        printf("VERDICT: ZERO-SUMFREE (length %d over C_6^%d)\n", L, R);
        return 0;
    }
    printf("VERDICT: NOT zero-sumfree (first zero-sum at term %d)\n", first_zero);
    return 1;
}
