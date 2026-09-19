#include "gf16.h"

gf16_t gf16_log[1 << 16];
gf16_t gf16_exp[2 << 16];

void gf16_init() {
    gf16_log[0] = (uint16_t) GF16_POLY;
    gf16_exp[65535] = 1;

	uint32_t x = 1;
    for (uint32_t i = 0; i < 65535; ++i) {
		gf16_exp[i] = x;
		gf16_log[x] = i;
		x <<= 1;
		if (x & 65536) x ^= GF16_POLY;
	}

    for (unsigned i = 0; i < 65535; ++i) {
        gf16_exp[i + 65535] = gf16_exp[i];
    }
}

gf16_t gf16_pow(gf16_t a, unsigned e) {
    gf16_t res = 1;
    for ( ; e > 0 ; e >>= 1) {
        if (e & 1) res = gf16_mul(res, a);
        a = gf16_sqr(a);
    }
    return res;
}
