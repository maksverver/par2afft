CFLAGS=-Wall -Wno-deprecated-declarations -O3 -march=native -g -D_FILE_OFFSET_BITS=64 -DGF16_VANDERMONDE_TESTS_INCLUDED
BENCHMARK_LDLIBS=$(LDLIBS) -lprofiler
PAR2AFFT_LDLIBS=$(LDLIBS) -lcrypto

BINS=benchmark par2afft
TESTS=crc32_test gf16_vandermonde_test
BENCHMARK_OBJS=benchmark.o gf16.o gf16_vandermonde_transpose.o
PAR2AFFT_OBJS=par2afft.o crc32.o gf16.o gf16_vandermonde_transpose.o
CRC32_TEST_OBJS=crc32_test.o crc32.o
GF16_VANDERMONDE_TEST_OBJS=gf16_vandermonde_test.o gf16.o gf16_vandermonde_transpose.o

all: $(BINS)

benchmark: $(BENCHMARK_OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(BENCHMARK_OBJS) $(BENCHMARK_LDLIBS)

par2afft: $(PAR2AFFT_OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(PAR2AFFT_OBJS) $(PAR2AFFT_LDLIBS)

crc32_test: $(CRC32_TEST_OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(CRC32_TEST_OBJS) $(LDLIBS)

gf16_vandermonde_test: $(GF16_VANDERMONDE_TEST_OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(GF16_VANDERMONDE_TEST_OBJS) $(LDLIBS)

test: $(TESTS)
	./run-tests.sh $(TESTS)

clean:
	rm -f ./*.o $(BINS) $(TESTS)

.PHONY: all clean test
