// This file was adapted from output produced by ChatGPT 5.6 Sol (high)

#include "gf16_vandermonde_transpose.h"
#include "gf16.h"

#include <stdint.h>
#include <assert.h>

#define FIELD_BITS 16
#define FIELD_ORDER (1 << FIELD_BITS)

typedef struct {
    unsigned count;
    uint32_t exponent[FIELD_BITS];
} subspace_shape_t;

static gf16_t beta[FIELD_BITS];
static gf16_t node[FIELD_ORDER];
static subspace_shape_t shape[FIELD_BITS];

/*
 * Full-field transpose Vandermonde example over
 *
 *   GF(2^16) = GF(2)[z] / (z^16 + z^12 + z^3 + z + 1)
 *
 * polynomial encoded as 0x1100B.  A field element is the 16-bit
 * polynomial-basis representation in z.
 *
 * We compute, for arbitrary a[x], x = 0..65535,
 *
 *   y[k] = sum_{x in GF(2^16)} a[x] * x^k,   k = 0..65535.
 *
 * Equivalently y = V^T a where V[x,k] = x^k.
 *
 * Internally:
 *   V = P^{-1} F C
 *
 * where
 *   C : monomial coefficients -> Cantor novel-basis coefficients,
 *   F : additive FFT, output in Cantor-coordinate point order,
 *   P : permutation between polynomial-basis field labels and Cantor order.
 *
 * Hence
 *   V^T = C^T F^T P.
 *
 * The monomial<->novel conversion is recursive and multiplication-free.
 * No 65536-entry "lower[]" table is used.
 */

/* Absolute trace GF(2^16) -> GF(2). */
static gf16_t gf16_trace(gf16_t a)
{
    gf16_t t = 0;
    gf16_t x = a;

    for (unsigned i = 0; i < FIELD_BITS; ++i) {
        t ^= x;
        x = gf16_sqr(x);
    }

    assert(t == 0 || t == 1);
    return t;
}

/* ------------------------------------------------------------------------- */
/* Cantor basis generation                                                   */
/* ------------------------------------------------------------------------- */

/*
 * Let A(x) = x^2 + x.  A Cantor chain satisfies
 *
 *   beta[0] = 1,
 *   A(beta[i]) = beta[i-1].
 *
 * For extension degree 16, pick any absolute-trace-one beta[15], then
 * repeatedly apply A downward.  We choose the numerically smallest trace-one
 * element so the result is deterministic and reproducible.
 */
static void generate_cantor_basis()
{
    gf16_t top = 0;

    for (uint32_t x = 1; x < FIELD_ORDER; ++x) {
        if (gf16_trace((gf16_t)x) == 1) {
            top = (gf16_t)x;
            break;
        }
    }

    assert(top != 0);

    beta[FIELD_BITS - 1] = top;
    for (unsigned i = FIELD_BITS - 1; i > 0; --i)
        beta[i - 1] = gf16_sqr(beta[i]) ^ beta[i];

    assert(beta[0] == 1);

    for (unsigned i = 1; i < FIELD_BITS; ++i)
        assert((gf16_sqr(beta[i]) ^ beta[i]) == beta[i - 1]);
}

static void build_cantor_permutation()
{
    node[0] = 0;
    for (uint32_t j = 1; j < FIELD_ORDER; ++j) {
        const uint32_t lb = j & (0u - j);       /* lowest set bit */
        const unsigned b = (unsigned)__builtin_ctz(lb);
        node[j] = node[j ^ lb] ^ beta[b];
    }
}

/* ------------------------------------------------------------------------- */
/* Subspace polynomials                                                      */
/* ------------------------------------------------------------------------- */

/*
 * With this Cantor basis, the vanishing polynomial of
 * V_i = span(beta[0],...,beta[i-1]) is
 *
 *   s_i(X) = A^i(X),  A(X)=X^2+X.
 *
 * Since A = Frobenius + identity,
 *
 *   s_i(X) = sum_{j=0}^i C(i,j) X^(2^j)   over GF(2).
 *
 * C(i,j) is odd iff j is a bit-submask of i (Lucas's theorem).
 * Thus s_i is sparse.  We only need its lower exponents, excluding its
 * leading X^(2^i) term.
 */

static void build_subspace_shapes()
{
    for (unsigned i = 0; i < FIELD_BITS; ++i) {
        shape[i].count = 0;
        for (unsigned j = 0; j < i; ++j) {
            if ((j & ~i) == 0) {
                shape[i].exponent[shape[i].count++] = 1u << j;
            }
        }
    }
}

/* Evaluate s_i(x) = A^i(x), where A(x)=x^2+x. */
static gf16_t subspace_poly_eval(unsigned i, gf16_t x)
{
    for (unsigned r = 0; r < i; ++r)
        x = gf16_sqr(x) ^ x;
    return x;
}

/* ------------------------------------------------------------------------- */
/* Recursive monomial -> novel basis conversion C                            */
/* ------------------------------------------------------------------------- */

#if GF16_VANDERMONDE_TESTS_INCLUDED

/*
 * At a node of dimension m, h=2^(m-1), divide
 *
 *   f(X) = f0(X) + s_{m-1}(X) f1(X),
 *
 * with deg(f0), deg(f1) < h.
 *
 * Because s_{m-1} is monic and sparse, polynomial long division is only XOR.
 * The upper half of a[] becomes f1, the lower half becomes f0.  We then
 * recurse on both halves.
 */
static void monomial_to_novel_rec(gf16_t *a, unsigned m)
{
    if (m == 0)
        return;

    const unsigned i = m - 1;
    const uint32_t h = 1u << i;

    /* Long division by s_i, from highest term downward. */
    for (uint32_t t = h; t-- > 0;) {
        const gf16_t c = a[h + t];

        /* c may be zero, but avoiding the branch is often competitive. */
        if (c != 0) {
            for (unsigned r = 0; r < shape[i].count; ++r)
                a[t + shape[i].exponent[r]] ^= c;
        }
    }

    monomial_to_novel_rec(a,     m - 1);
    monomial_to_novel_rec(a + h, m - 1);
}

#endif

/*
 * Transpose C^T.
 *
 * Forward division uses operations
 *
 *   a[t+e] ^= a[h+t]
 *
 * in descending t order.  Transposition reverses operation order and swaps
 * source/destination:
 *
 *   a[h+t] ^= a[t+e]
 *
 * in ascending t order.  Because forward C recurses after division, C^T
 * recurses first and applies the transposed division afterward.
 */
static void monomial_to_novel_transpose_rec(gf16_t *a, unsigned m)
{
    if (m == 0)
        return;

    const unsigned i = m - 1;
    const uint32_t h = 1u << i;

    monomial_to_novel_transpose_rec(a,     m - 1);
    monomial_to_novel_transpose_rec(a + h, m - 1);

    for (uint32_t t = 0; t < h; ++t) {
        gf16_t c = a[h + t];
        for (unsigned r = 0; r < shape[i].count; ++r)
            c ^= a[t + shape[i].exponent[r]];
        a[h + t] = c;
    }
}

/* ------------------------------------------------------------------------- */
/* Additive FFT F and its transpose F^T                                      */
/* ------------------------------------------------------------------------- */

#if GF16_VANDERMONDE_TESTS_INCLUDED

/*
 * Input is in novel basis for dimension m.
 * Output is evaluations at alpha + V_m in Cantor-coordinate order.
 *
 * Split f = f0 + s_i f1, i=m-1.  On alpha+V_i, s_i is the constant
 * c=s_i(alpha); on alpha+beta_i+V_i it is c+1 because s_i(beta_i)=1.
 *
 * One-multiply butterfly:
 *
 *   L = f0 + c*f1
 *   R = L  + f1
 */
static void additive_fft(gf16_t *a, unsigned m, gf16_t alpha)
{
    if (m == 0)
        return;

    const unsigned i = m - 1;
    const uint32_t h = 1u << i;
    const gf16_t c = subspace_poly_eval(i, alpha);

    for (uint32_t t = 0; t < h; ++t) {
        const gf16_t f0 = a[t];
        const gf16_t f1 = a[h + t];
        const gf16_t left = f0 ^ gf16_mul(c, f1);
        a[t]     = left;
        a[h + t] = left ^ f1;
    }

    additive_fft(a,     m - 1, alpha);
    additive_fft(a + h, m - 1, alpha ^ beta[i]);
}

#endif

/*
 * Exact transpose of additive_fft().
 *
 * Transpose of
 *   L = f0 + c*f1
 *   R = L + f1
 * is
 *   f0' = L' + R'
 *   f1' = c*(L'+R') + R'.
 */
static void additive_fft_transpose(gf16_t *a, unsigned m, gf16_t alpha)
{
    if (m == 0)
        return;

    const unsigned i = m - 1;
    const uint32_t h = 1u << i;

    /* Reverse the forward recursion first. */
    additive_fft_transpose(a,     m - 1, alpha);
    additive_fft_transpose(a + h, m - 1, alpha ^ beta[i]);

    const gf16_t c = subspace_poly_eval(i, alpha);

    for (uint32_t t = 0; t < h; ++t) {
        const gf16_t L = a[t];
        const gf16_t R = a[h + t];
        const gf16_t f0 = L ^ R;
        const gf16_t f1 = gf16_mul(c, f0) ^ R;
        a[t]     = f0;
        a[h + t] = f1;
    }
}

/* ------------------------------------------------------------------------- */
/* Fast V and V^T in external polynomial-basis field-element ordering        */
/* ------------------------------------------------------------------------- */

/* Compute y[k] = sum_x a[x] x^k. */
void gf16_vandermonde_transpose_multiply(const gf16_t a[FIELD_ORDER], gf16_t y[FIELD_ORDER]) {
    /* P: gather external field-label order into Cantor-coordinate order. */
    for (uint32_t j = 0; j < FIELD_ORDER; ++j)
        y[j] = a[node[j]];

    additive_fft_transpose(y, FIELD_BITS, 0);
    monomial_to_novel_transpose_rec(y, FIELD_BITS);
}

void gf16_vandermonde_transpose_init() {
    generate_cantor_basis();
    build_subspace_shapes();
    build_cantor_permutation();
}

#if GF16_VANDERMONDE_TESTS_INCLUDED

#include "w1rand.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------------- */
/* Reference helpers                                                         */
/* ------------------------------------------------------------------------- */

/* Directly compute one output y[k] in O(FIELD_ORDER log k) field operations. */
static gf16_t vandermonde_transpose_one_naive(const gf16_t a[FIELD_ORDER], uint32_t k)
{
    gf16_t sum = 0;

    for (uint32_t x = 0; x < FIELD_ORDER; ++x)
        sum ^= gf16_mul(a[x], gf16_pow((gf16_t)x, k));

    return sum;
}

/* Faster direct check of selected k: walk x^k by exponentiation per x. */
static void check_selected_outputs(const gf16_t a[FIELD_ORDER],
                                   const gf16_t y[FIELD_ORDER],
                                   const uint32_t *ks,
                                   size_t nks)
{
    for (size_t q = 0; q < nks; ++q) {
        const uint32_t k = ks[q];
        const gf16_t ref = vandermonde_transpose_one_naive(a, k);
        if (ref != y[k]) {
            fprintf(stderr,
                    "FAIL at k=%u: fast=%04X naive=%04X\n",
                    k, y[k], ref);
            exit(1);
        }
        printf("  k=%5u : %04X  [OK]\n", k, y[k]);
    }
}

/* Evaluate monomial polynomial c at every field element. */
static void vandermonde_forward(const gf16_t c[FIELD_ORDER], gf16_t values[FIELD_ORDER]) {
    memcpy(values, c, FIELD_ORDER * sizeof(gf16_t));

    monomial_to_novel_rec(values, FIELD_BITS);
    additive_fft(values, FIELD_BITS, 0);

    /* values currently indexed by Cantor coordinates; scatter to x labels. */
    gf16_t *tmp = malloc(FIELD_ORDER * sizeof(gf16_t));
    if (!tmp) {
        fprintf(stderr, "allocation failed\n");
        exit(2);
    }
    memcpy(tmp, values, FIELD_ORDER * sizeof(gf16_t));

    for (uint32_t j = 0; j < FIELD_ORDER; ++j)
        values[node[j]] = tmp[j];

    free(tmp);
}

/* Field-valued dot product. */
static gf16_t dot_product(const gf16_t *a, const gf16_t *b)
{
    gf16_t s = 0;
    for (uint32_t i = 0; i < FIELD_ORDER; ++i)
        s ^= gf16_mul(a[i], b[i]);
    return s;
}

int gf16_vandermonde_transpose_test() {
    gf16_t *a    = malloc(FIELD_ORDER * sizeof(gf16_t));
    gf16_t *y    = malloc(FIELD_ORDER * sizeof(gf16_t));
    gf16_t *c    = malloc(FIELD_ORDER * sizeof(gf16_t));
    gf16_t *Vc   = malloc(FIELD_ORDER * sizeof(gf16_t));

    if (!a || !y || !c || !Vc) {
        fprintf(stderr, "allocation failed\n");
        return 2;
    }

    printf("Cantor basis for modulus 0x1100B:\n");
    for (unsigned i = 0; i < FIELD_BITS; ++i)
        printf("  beta[%2u] = 0x%04X\n", i, beta[i]);

    /* Verify node[] is really a permutation of all 65536 field elements. */
    uint8_t *seen = calloc(FIELD_ORDER, 1);
    assert(seen);
    for (uint32_t j = 0; j < FIELD_ORDER; ++j) {
        assert(!seen[node[j]]);
        seen[node[j]] = 1;
    }
    free(seen);

    printf("\nTest 1: a[x] = x\n");

    for (uint32_t x = 0; x < FIELD_ORDER; ++x)
        a[x] = (gf16_t)x;

    gf16_vandermonde_transpose_multiply(a, y);

    /*
     * y[k] = sum_x x^(k+1).
     * In GF(q), q=65536, this is nonzero only when (q-1)|(k+1).
     * For k=0..65535, the only such k is 65534, and the value is 1.
     */
    for (uint32_t k = 0; k < FIELD_ORDER; ++k) {
        const gf16_t expected = (k == FIELD_ORDER - 2) ? 1 : 0;
        if (y[k] != expected) {
            fprintf(stderr,
                    "FAIL structured test at k=%u: got=%04X expected=%04X\n",
                    k, y[k], expected);
            return 1;
        }
    }
    printf("  PASS: only y[65534] is 1; all other outputs are 0.\n");

    /*
     * Test 2: arbitrary deterministic vector; compare selected outputs with
     * the literal definition sum_x a[x] x^k.
     */
    printf("\nTest 2: pseudorandom vector, selected direct checks\n");
    uint64_t rand_state = 123456789;
    for (uint32_t x = 0; x < FIELD_ORDER; ++x)
        a[x] = (gf16_t)w1rand(&rand_state);

    gf16_vandermonde_transpose_multiply(a, y);

    static const uint32_t ks[] = {
        0, 1, 2, 3, 7, 15, 16, 31,
        255, 256, 257, 1023, 4096, 32767,
        65534, 65535
    };
    check_selected_outputs(a, y, ks, sizeof(ks) / sizeof(ks[0]));
    printf("  PASS: selected outputs match the direct definition.\n");

    /*
     * Test 3: exact transpose identity
     *
     *   <V^T a, c> = <a, V c>
     *
     * for an independent deterministic vector c.
     */
    printf("\nTest 3: transpose inner-product identity\n");
    for (uint32_t k = 0; k < FIELD_ORDER; ++k)
        c[k] = (gf16_t)w1rand(&rand_state);

    vandermonde_forward(c, Vc);

    const gf16_t lhs = dot_product(y, c);
    const gf16_t rhs = dot_product(a, Vc);

    printf("  <V^T a, c> = %04X\n", lhs);
    printf("  <a, V c>   = %04X\n", rhs);

    if (lhs != rhs) {
        fprintf(stderr, "FAIL: transpose identity does not hold\n");
        return 1;
    }
    printf("  PASS: transpose identity holds.\n");

    free(a);
    free(y);
    free(c);
    free(Vc);

    printf("\nALL TESTS PASSED\n");
    return 0;
}
#endif
