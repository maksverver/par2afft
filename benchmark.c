#include "gf16.h"
#include "gf16_vandermonde_transpose.h"
#include "w1rand.h"

#include <assert.h>
#include <stdio.h>
#include <time.h>
#include <inttypes.h>
#include <string.h>

int main() {
    gf16_init();
    gf16_vandermonde_transpose_init();

    // I use this property below to randomly fill the buffer.
    assert(sizeof(uint64_t) % sizeof(gf16_t) == 0);

    int niter = 5000;
    printf("Running benchmark (%d iterations)...\n", niter);
    clock_t clock_begin = clock();
    uint64_t input_bytes = 0;
    uint64_t checksum = 0;
    uint64_t rng_state = 0x123456789abcdef;
    for (int iter = 0; iter < niter; ++iter) {
        // Fill buffer with pseudo-random data
        gf16_t a[1 << 16];
        for (size_t i = 0; i < (1 << 16); i += sizeof(uint64_t)/sizeof(gf16_t)) {
            uint64_t value = w1rand(&rng_state);
            memcpy(&a[i], &value, sizeof(value));
        }
        input_bytes += sizeof(a);

        // Perform fast multiplication.
        gf16_t y[1 << 16];
        gf16_vandermonde_transpose_multiply(a, y);

        // Use the result to make sure nothing gets optimized away
        for (size_t i = 0; i < (1 << 16); ++i) checksum += y[i];
    }
    double time_elapsed = (double) (clock() - clock_begin) / CLOCKS_PER_SEC;
    printf("Time elapsed: %.6f s\n", time_elapsed);
    printf("Throughput: %.3f MiB/s\n", (double) input_bytes / (1 << 20) / time_elapsed);
    printf("Checksum: %" PRIx64 "\n", checksum);
    // const uint64_t expected_checksum = 0x1f4047e81cf;  // for niter == 1000
    const uint64_t expected_checksum = 0x9c3e5d8b098;
    if (checksum != expected_checksum) {
        fprintf(stderr, "WARNING! Checksum does not match expected value: %"PRIx64"\n", expected_checksum);
        return 1;
    }
}
