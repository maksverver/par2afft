#ifndef GF16_H_INCLUDED
#define GF16_H_INCLUDED

#include <stdint.h>

// Irereducible polynomial = 0x1100B = x^16 + x^12 + x^3 + x + 1
// which conveniently has x (2) as a primitive element.
#define GF16_POLY  0x1100B

typedef uint16_t gf16_t;

extern gf16_t gf16_log[1 << 16];
extern gf16_t gf16_exp[2 << 16];

// Must be called at startup to initialize the log/exp tables. Idempotent.
void gf16_init();

static inline gf16_t gf16_mul(gf16_t a, gf16_t b) {
    return a == 0 || b == 0 ? 0 : gf16_exp[gf16_log[a] + gf16_log[b]];
}
static inline gf16_t gf16_sqr(gf16_t a) { return gf16_mul(a, a); }
static inline gf16_t gf16_add(gf16_t a, gf16_t b) { return a ^ b; }

gf16_t gf16_pow(gf16_t a, unsigned e);

#endif //  ndef GF16_H_INCLUDED
