#include "gf16.h"
#include "gf16_vandermonde_transpose.h"

#include <stdio.h>

int main()
{
    gf16_init();
    gf16_vandermonde_transpose_init();
#if GF16_VANDERMONDE_TESTS_INCLUDED
    int res = gf16_vandermonde_transpose_test();
    if (res != 0) return res;
#else
    fprintf(stderr, "gf16_vandermonde_transpose_test() not compiled in!\n");
    return 1;
#endif
}
