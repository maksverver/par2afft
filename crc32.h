#ifndef CRC32_H_INCLUDED
#define CRC32_H_INCLUDED

#include <stdint.h>
#include <stddef.h>

extern const uint32_t crc32_table[256];

static inline uint32_t crc32_update(uint32_t crc, const void *data, size_t size) {
    const uint8_t *p = (const uint8_t *) data;
    while (size --> 0) crc = crc32_table[(crc ^ *p++) & 0xFFu] ^ (crc >> 8);
    return crc;
}

static inline uint32_t crc32_init() {
    return 0xFFFFFFFFu;
}

static inline uint32_t crc32_finish(uint32_t crc) {
    return crc ^ 0xFFFFFFFFu;
}

static inline uint32_t crc32(const void *data, size_t size) {
    return crc32_finish(crc32_update(crc32_init(), data, size));
}

#endif  // ndef CRC32_H_INCLUDED
