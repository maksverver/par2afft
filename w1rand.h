#ifndef W1RAND_H_INCLUDED
#define W1RAND_H_INCLUDED

#include <stdint.h>

static inline uint64_t w1rand(uint64_t *s) {
    const uint64_t c = 0xd07ebc63274654c7ull;
    *s += c;
    __uint128_t t = (__uint128_t) *s * (*s ^ c);
    return (t >> 64) ^ t;
}

#endif // ndef W1RAND_H_INCLUDED
