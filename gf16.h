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

typedef uint16_t gf16_t;

// Must be called at startup to initialize the lookup tables. Idempotent.
void gf16_init();

extern gf16_t gf16_log_lut[1 << 16];
extern gf16_t gf16_exp_lut[2 << 16];
extern gf16_t gf16_sqr_lut[1 << 16];

static inline gf16_t gf16_add(gf16_t a, gf16_t b) { return a ^ b; }
static inline gf16_t gf16_sqr(gf16_t a) { return gf16_sqr_lut[a]; }

static inline gf16_t gf16_mul(gf16_t a, gf16_t b) {
    return a == 0 || b == 0 ? 0 : gf16_exp_lut[gf16_log_lut[a] + gf16_log_lut[b]];
}

gf16_t gf16_pow(gf16_t a, unsigned e);

#endif //  ndef GF16_H_INCLUDED
