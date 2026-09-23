// Implements basic arithmetic over Galois field GF(2^16).
//
// Field elements are defined by the irereducible polynomial:
//
//  0x1100B = x^16 + x^12 + x^3 + x + 1,
//
// which conveniently has x (2) as a primitive element.

#ifndef GF16_H_INCLUDED
#define GF16_H_INCLUDED

#include <stdint.h>

#define GF16_POLY  0x1100B
#define GF16_BITS  16
#define GF16_ORDER 65536

typedef uint16_t gf16_t;

// Must be called at startup to initialize the lookup tables. Idempotent.
void gf16_init();

extern gf16_t gf16_log_lut[GF16_ORDER];
extern gf16_t gf16_exp_lut[GF16_ORDER * 2];
extern gf16_t gf16_sqr_lut[GF16_ORDER];

static inline gf16_t gf16_add(gf16_t a, gf16_t b) { return a ^ b; }
static inline gf16_t gf16_sqr(gf16_t a) { return gf16_sqr_lut[a]; }

static inline gf16_t gf16_mul(gf16_t a, gf16_t b) {
    return a == 0 || b == 0 ? 0 : gf16_exp_lut[gf16_log_lut[a] + gf16_log_lut[b]];
}

// assumes x > 0
static inline gf16_t gf16_log(gf16_t x) { return gf16_log_lut[x]; }

static inline gf16_t gf16_exp(gf16_t x) { return gf16_exp_lut[x]; }

static inline gf16_t gf16_mul_log(gf16_t a, gf16_t log_b) {
    return a == 0 ? 0 : gf16_exp_lut[gf16_log_lut[a] + log_b];
}

gf16_t gf16_pow(gf16_t a, unsigned e);

#endif //  ndef GF16_H_INCLUDED
