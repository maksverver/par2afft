#include "crc32.h"

#include <string.h>
#include <stdio.h>

static int pass = 0;
static int fail = 0;

static void crc_test(const void *data, size_t size, uint32_t expected) {
    uint32_t received = crc32(data, size);
    printf("expected (%08x) %s received (%08x)\n", expected, (expected == received ? "==" : "!="), received);
    if (expected == received) ++pass; else ++fail;
}

int main() {
    const char *digits = "123456789";
    crc_test(digits, strlen(digits), 0xcbf43926);

    char slice[100] = "hello world\n";
    crc_test(slice, sizeof(slice), 0xf09d7d02);

    printf("%d pass. %d fail.\n", pass, fail);
    return fail > 0 ? 1 : 0;
}
