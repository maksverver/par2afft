// A minimal implementation of the par2 client for testing purposes.
//
// Uses the Additive FFT to multiply with the transpose Vandermonde matrix in
// O(q log q) time, where q is the size of the finite field, 65536 for Par2.

#include "crc32.h"
#include "gf16.h"
#include "gf16_vandermonde_transpose.h"

#include <openssl/md5.h>

#include <fcntl.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <unistd.h>

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>
#include <inttypes.h>
#include <string.h>

#define CRC32_LENGTH 4

#define MAX_BLOCK_COUNT 32768
#define MAX_FILE_COUNT MAX_BLOCK_COUNT

// Note: buffer size must be a multiple of 4.
static const char creator_string[16] = "par2afft 0.0";

static const char *packet_type_creator  = "PAR 2.0\0Creator\0";
static const char *packet_type_main     = "PAR 2.0\0Main\0\0\0\0";
static const char *packet_type_filedesc = "PAR 2.0\0FileDesc";
static const char *packet_type_ifsc     = "PAR 2.0\0IFSC\0\0\0\0";
static const char *packet_type_recovery = "PAR 2.0\0RecvSlic";

// Just blindly assuming structs are packed here.
// Also we assume we're on a little-endian system.
struct PacketHeader {
    uint8_t  magic[8];
    uint64_t length;   /* must be multiple of 4 */
    uint8_t  md5[MD5_DIGEST_LENGTH];
    uint8_t  recovery_set_id[MD5_DIGEST_LENGTH];
    uint8_t  type[16];
} __attribute__((packed));

struct MainPacketBody {
    uint64_t slice_size;
    uint32_t file_count;
    uint8_t  md5[MAX_FILE_COUNT][MD5_DIGEST_LENGTH];
}  __attribute__((packed));

struct MainPacket {
    struct PacketHeader header;
    struct MainPacketBody body;
} __attribute__((packed));

struct FileDescriptionPacketBody {
    uint8_t  id[MD5_DIGEST_LENGTH];
    uint8_t  md5[MD5_DIGEST_LENGTH];
    uint8_t  md5_16k[MD5_DIGEST_LENGTH];
    uint64_t length;
    char     name[];
} __attribute__((packed));

struct FileDescriptionPacket {
    struct PacketHeader header;
    struct FileDescriptionPacketBody body;
} __attribute__((packed));

struct SliceChecksums {
    uint8_t md5[MD5_DIGEST_LENGTH];
    uint8_t crc32[CRC32_LENGTH];
} __attribute__((packed));

// Input File Slice Checksums
struct IfscPacketBody {
    uint8_t file_id[MD5_DIGEST_LENGTH];
    struct SliceChecksums checksums[];
} __attribute__((packed));

struct IfscPacket {
    struct PacketHeader header;
    struct IfscPacketBody body;
} __attribute__((packed));

struct RecoverySlicePacketBody {
    uint32_t exponent;
    gf16_t elems[];
} __attribute__((packed));

struct RecoverySlicePacket {
    struct PacketHeader header;
    struct RecoverySlicePacketBody body;
} __attribute__((packed));

struct InputFile {
    char    *name;
    uint64_t size;
    int      fd;
    void    *addr;
    uint8_t  md5_16k[MD5_DIGEST_LENGTH];
    uint8_t  id[MD5_DIGEST_LENGTH];  // must be unsigned for compare_input_files()
};

struct InputSlice {
    gf16_t      constant;
    const void *addr;
    size_t      size;
};

static struct InputFile input_files[MAX_FILE_COUNT];

static int input_file_count = 0;

static struct InputSlice input_blocks[MAX_BLOCK_COUNT];

static struct RecoverySlicePacket *recovery_blocks[MAX_BLOCK_COUNT];

static uint8_t recovery_set_id[MD5_DIGEST_LENGTH];


static struct MainPacket main_packet;

static int out_fd = -1;
static void *out_mapped_addr = NULL;
static size_t out_mapped_size = 0;

static long long arg_input_blocks   = -1;
static long long arg_block_size     = -1;
static long long arg_redundancy     = -1;
static long long arg_output_blocks  = -1;

static const uint8_t zeroes[65536];

static void debug_print_hash(const uint8_t md5[MD5_DIGEST_LENGTH]) {
    for (int i = 0; i < MD5_DIGEST_LENGTH; ++i) {
        fprintf(stderr, "%02x", md5[i]);
    }
    fprintf(stderr, "\n");
}

static void crc32_slice(uint8_t crc32[CRC32_LENGTH], const void *data, size_t size, size_t padding) {
    uint32_t crc = crc32_init();
    crc = crc32_update(crc, data, size);
    while (padding > sizeof(zeroes)) {
        crc = crc32_update(crc, zeroes, sizeof(zeroes));
        padding -= sizeof(zeroes);
    }
    if (padding > 0) {
        crc = crc32_update(crc, zeroes, padding);
        padding -= sizeof(zeroes);
    }
    crc = crc32_finish(crc);

    // Assume we're on a little-endian system
    memcpy(crc32, &crc, CRC32_LENGTH);
}

static void md5_slice(uint8_t md5[MD5_DIGEST_LENGTH], const void *data, size_t size, size_t padding) {
    MD5_CTX md5_ctx;
    MD5_Init(&md5_ctx);
    MD5_Update(&md5_ctx, data, size);
    while (padding > sizeof(zeroes)) {
        MD5_Update(&md5_ctx, zeroes, sizeof(zeroes));
        padding -= sizeof(zeroes);
    }
    if (padding > 0) {
        MD5_Update(&md5_ctx, zeroes, padding);
        padding -= sizeof(zeroes);
    }
    MD5_Final(md5, &md5_ctx);
}

static void show_usage() {
    fputs(creator_string, stderr);
    fputs(
        "\n"
        "Usage: par2afft create [options] <output.par2> <input>...\n"
        "\n"
        "Options:\n"
        "    -b<n>   block count (default 2000)\n"
        "    -s<n>   block size\n"
        "    -r<n>   redundancy percentage (default 5)\n"
        "    -c<n>   recovery block count (default 100)\n"
        "Note: space after options is not allowed.\n",
        stderr);
}

static int parse_option(const char *arg) {
    switch (arg[1]) {
    case 'b':
        if (sscanf(&arg[2], "%lld", &arg_input_blocks) == 1) return 0;
        break;
    case 's':
        if (sscanf(&arg[2], "%lld", &arg_block_size) == 1) return 0;
        break;
    case 'r':
        if (sscanf(&arg[2], "%lld", &arg_redundancy) == 1) return 0;
        break;
    case 'c':
        if (sscanf(&arg[2], "%lld", &arg_output_blocks) == 1) return 0;
        break;
    }
    fprintf(stderr, "Invalid option argument: %s\n", arg);
    return 1;
}

static void fill_packet_header(
        struct PacketHeader *header, const char packet_type[16],
        const void *body_data, size_t body_length) {
    assert(body_length % 4 == 0);
    memcpy(header->magic, "PAR2\0PKT", 8);
    header->length = body_length + sizeof(struct PacketHeader);
    memcpy(header->recovery_set_id, recovery_set_id, MD5_DIGEST_LENGTH);
    memcpy(header->type, packet_type, 16);

    MD5_CTX md5_ctx;
    MD5_Init(&md5_ctx);
    // This assumes the body starts immediately after the header.
    MD5_Update(&md5_ctx, &header->recovery_set_id,
            sizeof(struct PacketHeader) - offsetof(struct PacketHeader, recovery_set_id));
    MD5_Update(&md5_ctx, body_data, body_length);
    MD5_Final(header->md5, &md5_ctx);
}

static int write_packet(struct PacketHeader *header, const void *body_data) {
    size_t header_length = sizeof(struct PacketHeader);
    assert(header->length >= header_length);
    if (write(out_fd, header, header_length) != header_length) {
        perror("Failed to write packet header");
        return -1;
    }
    size_t body_length = header->length - header_length;
    if (write(out_fd, body_data, body_length) != body_length) {
        perror("Failed to write packet body");
        return -1;
    }
    return 0;
}

static int write_creator_packet() {
    struct PacketHeader header;
    fill_packet_header(&header, packet_type_creator, creator_string, sizeof(creator_string));
    return write_packet(&header, creator_string);
}

static int write_main_packet() {
    return write_packet(&main_packet.header, &main_packet.body);
}

static int write_file_description_packet(const struct InputFile *input_file) {
    size_t name_len = strlen(input_file->name);
    if (name_len % 4) name_len += 4 - name_len % 4;  // pad to 4 bytes
    struct FileDescriptionPacket *packet =
        calloc(1, sizeof(struct FileDescriptionPacket) + name_len);
    if (packet == NULL) {
        perror("calloc");
        return -1;
    }
    memcpy(packet->body.id, input_file->id, MD5_DIGEST_LENGTH);
    memcpy(packet->body.md5_16k, input_file->md5_16k, MD5_DIGEST_LENGTH);
    packet->body.length = input_file->size;
    strncpy(packet->body.name, input_file->name, name_len);
    MD5_CTX md5_ctx;
    MD5_Init(&md5_ctx);
    MD5_Update(&md5_ctx, input_file->addr, input_file->size);
    MD5_Final(packet->body.md5, &md5_ctx);
    fill_packet_header(&packet->header, packet_type_filedesc,
            &packet->body, sizeof(struct FileDescriptionPacketBody) + name_len);
    int res = write_packet(&packet->header, &packet->body);
    free(packet);
    return res;
}

// ifsc == Input File Slice Checksum
// This writes the IFSC packet and as a side-effect also populates input_slices.
static int write_ifsc_packet(const struct InputFile *input_file, struct InputSlice **slice_ptr) {
    const uint8_t *data = input_file->addr;
    size_t size = input_file->size;
    uint64_t slices = (size + arg_block_size - 1) / arg_block_size;

    size_t payload_size = sizeof(struct SliceChecksums) * slices;
    struct IfscPacket *packet = calloc(1, sizeof(struct IfscPacket) + payload_size);
    if (packet == NULL) {
        perror("calloc");
        return -1;
    }
    memcpy(packet->body.file_id, input_file->id, MD5_DIGEST_LENGTH);
    for (uint64_t i = 0; i < slices; ++i) {
        assert(size > 0);
        size_t slice_size = size < arg_block_size ? size : arg_block_size;
        size_t padding = arg_block_size - slice_size;
        crc32_slice(packet->body.checksums[i].crc32, data, slice_size, padding);
        md5_slice(packet->body.checksums[i].md5, data, slice_size, padding);
        (*slice_ptr)->addr = data;
        (*slice_ptr)->size = slice_size;
        (*slice_ptr)++;
        data += slice_size;
        size -= slice_size;
    }
    assert(size == 0);
    fill_packet_header(&packet->header, packet_type_ifsc,
            &packet->body, sizeof(struct IfscPacketBody) + payload_size);
    int res = write_packet(&packet->header, &packet->body);
    free(packet);
    return res;
}

// Creates the main packet after input files have been opened.
// This also generates the recovery_set_id, which is needed for any other packets.
static void create_main_packet() {
    main_packet.body.slice_size = arg_block_size;
    main_packet.body.file_count = input_file_count;
    for (int i = 0; i < input_file_count; ++i) {
        memcpy(main_packet.body.md5[i], input_files[i].id, MD5_DIGEST_LENGTH);
    }

    // Generate recovery_set_id as the hash of the body of the main packet.
    size_t body_length = sizeof(uint64_t) + sizeof(uint32_t) +
            (size_t) input_file_count * MD5_DIGEST_LENGTH;
    MD5_CTX md5_ctx;
    MD5_Init(&md5_ctx);
    MD5_Update(&md5_ctx, &main_packet.body, body_length);
    MD5_Final(recovery_set_id, &md5_ctx);

    // Now we can use it to fill in the header
    fill_packet_header(&main_packet.header, packet_type_main, &main_packet.body, body_length);
}

static int open_output_file(const char *filename) {
    if (strcmp(filename, "-") == 0) {
        out_fd = fileno(stdout);
        return 0;
    }
    out_fd = open(filename, O_RDWR | O_CREAT | O_EXCL, 0644);
    if (out_fd == -1) {
        if (errno == EEXIST) {
            fprintf(stderr, "Output file already exists: %s\n", filename);
        } else {
            perror("Failed to open output file");
        }
        return -1;
    }
    return 0;
}

static int open_input_file(const char *filename) {
    if (input_file_count >= MAX_FILE_COUNT) {
        fprintf(stderr, "Too many input files!\n");
        return -1;
    }
    int fd = open(filename, O_RDONLY);
    if (fd == -1) {
        fprintf(stderr, "Failed to open input file (%s): %s\n", filename, strerror(errno));
        return -1;
    }

    struct stat st;
    if (fstat(fd, &st) != 0) {
        fprintf(stderr, "Could not stat input file (%s): %s\n", filename, strerror(errno));
        close(fd);
        return -1;
    }
    // Note: size must be 8 bytes since we hash it below
    uint64_t size = st.st_size;
    if (size == 0) {
        fprintf(stderr, "Warning: %s is an empty file!\n", filename);
    }

    void *addr = NULL;
    if (size > 0) {
        addr = mmap(NULL, size, PROT_READ, MAP_PRIVATE, fd, 0);
        if (addr == MAP_FAILED) {
            fprintf(stderr, "Could not mmap input file (%s): %s\n", filename, strerror(errno));
            close(fd);
            return -1;
        }
    }
    char *name = strdup(filename);
    if (name == NULL) {
        perror("strdup");
        munmap(addr, size);
        close(fd);
        return -1;
    }

    struct InputFile *file = &input_files[input_file_count++];
    file->name = name;
    file->size = size;
    file->fd   = fd;
    file->addr = addr;

    /* Calculate MD5 hashes */
    MD5_CTX md5_ctx;
    MD5_Init(&md5_ctx);
    MD5_Update(&md5_ctx, addr, size < 16384 ? size : 16384);
    MD5_Final(file->md5_16k, &md5_ctx);
    MD5_Init(&md5_ctx);
    MD5_Update(&md5_ctx, file->md5_16k, MD5_DIGEST_LENGTH);
    MD5_Update(&md5_ctx, &size, sizeof(size));
    MD5_Update(&md5_ctx, name, strlen(name));
    MD5_Final(file->id, &md5_ctx);
    return 0;
}

static int reserve_output_slices(int slice_count) {
    size_t recovery_packet_size = sizeof(struct RecoverySlicePacket) + arg_block_size;
    assert(recovery_packet_size % 4 == 0);
    off_t current_pos = lseek(out_fd, 0, SEEK_CUR);
    if (current_pos == -1) {
        perror("Failed to determine output position (writing to a pipe?)");
        return -1;
    }
    uint64_t space_needed = recovery_packet_size * slice_count;
    if (posix_fallocate(out_fd, current_pos, space_needed) != 0) {
        // No perror(), because posix_fallocate() does NOT set errno on failure
        fprintf(stderr, "Failed to allocate disk space for recovery slices!\n");
        return -1;
    }
    uint64_t new_size = current_pos + space_needed;
    out_mapped_addr = mmap(NULL, new_size, PROT_READ | PROT_WRITE, MAP_SHARED, out_fd, 0);
    if (out_mapped_addr == MAP_FAILED) {
        perror("Failed to mmap output file");
        return -1;
    }
    out_mapped_size = new_size;

    char *ptr = (char*) out_mapped_addr + current_pos;
    for (int i = 0; i < slice_count; ++i) {
        struct RecoverySlicePacket *packet = (struct RecoverySlicePacket*)ptr;
        ptr += recovery_packet_size;
        recovery_blocks[i] = packet;
        // TODO: support starting from a different starting exponent
        // (but should not allow exponent to exceed 32768/65536?)
        packet->body.exponent = i;
    }
    return 0;
}

static int finalize_output_slices(int slice_count) {
    size_t body_size = sizeof(struct RecoverySlicePacketBody) + arg_block_size;
    assert(body_size % 4 == 0);
    for (int i = 0; i < slice_count; ++i) {
        fill_packet_header(&recovery_blocks[i]->header, packet_type_recovery, &recovery_blocks[i]->body, body_size);
    }
    return 0;
}

static void generate_recovery_data() {
    assert(arg_block_size % 4 == 0 && sizeof(gf16_t) == 2);
    gf16_t a[1 << 16];
    gf16_t y[1 << 16];
    memset(a, 0, sizeof(a));
    int last_progress = isatty(fileno(stderr)) ? -1 : 100;
    for (int j = 0; j < arg_block_size / sizeof(gf16_t); ++j) {
        // Print progress
        int progress = j * 100 / (arg_block_size / sizeof(gf16_t));
        if (progress > last_progress) {
            fprintf(stderr, "%3d%%\r", progress);
            last_progress = progress;
        }
        for (int i = 0; i < arg_input_blocks; ++i) {
            // Load the j-th element from the i-th input block.
            // This is slightly tricky because the end of the file may occur
            // in the middle of a slice.
            const struct InputSlice *s = &input_blocks[i];
            if (sizeof(gf16_t) * (j + 1) <= s->size) {
                a[s->constant] = ((gf16_t*) s->addr)[j];
            } else if (sizeof(gf16_t) * j >= s->size) {
                a[s->constant] = 0;  // zero pad to fill slice
            } else {
                // One byte left in the file, the other zero-padded. Since we
                // use little endian encoding, load just the one byte.
                assert(2*j + 1 == s->size);
                a[s->constant] = ((uint8_t*) s->addr)[2*j];
            }
        }
        gf16_vandermonde_transpose_multiply(a, y);
        for (int i = 0; i < arg_output_blocks; ++i) {
            struct RecoverySlicePacketBody *body = &recovery_blocks[i]->body;
            body->elems[j] = y[body->exponent];
        }
    }
}

static int compare_input_files_by_id(const void *a, const void *b) {
    const struct InputFile *p = a, *q = b;
    for (int i = MD5_DIGEST_LENGTH - 1; i >= 0; --i) {
        int diff = p->id[i] - q->id[i];
        if (diff != 0) return diff;
    }
    return 0;
}

static void sort_input_files_by_id() {
    qsort(input_files, input_file_count, sizeof(struct InputFile),
            compare_input_files_by_id);
}

static void init_input_constants() {
    // From the Par2 spec: input constants are x^n where n is not divisible
    // by 3, 5, 17 or 257.
    int j = 0;
    uint32_t elem = 1;
    for (int n = 0; n < (1 << 16); ++n) {
        if (n%3 != 0 && n%5 != 0 && n%17 != 0 && n%257 != 0) {
            input_blocks[j++].constant = elem;
        }
        elem <<= 1;
		if (elem & 65536) elem ^= GF16_POLY;
    }
    assert(j == MAX_BLOCK_COUNT);
}

static int is_prefix(const char *prefix, const char *full) {
    while (*prefix) {
        int diff = *prefix++ - *full++;
        if (diff != 0) return diff;
    }
    return 0;
}

static uint64_t calculate_max_file_size() {
    uint64_t size = 0;
    for (int i = 0; i < input_file_count; ++i) {
        if (input_files[i].size > size) size = input_files[i].size;
    }
    return size;
}
static uint64_t count_nonempty_files() {
    uint64_t count = 0;
    for (int i = 0; i < input_file_count; ++i) {
        if (input_files[i].size > 0) ++count;
    }
    return count;
}

static uint64_t calculate_block_count(uint64_t block_size) {
    uint64_t count = 0;
    for (int i = 0; i < input_file_count; ++i) {
        count += (input_files[i].size + block_size - 1) / block_size;
    }
    return count;
}

static uint64_t calculate_block_size(uint64_t block_count) {
    assert(block_count >= count_nonempty_files());
    uint64_t lo = 1, hi = 1;
    while (calculate_block_count(4*hi) > block_count) {
        lo = hi + 1;
        hi *= 2;
    }
    while (lo < hi) {
        uint64_t mid = lo + (hi - lo)/2;
        if (calculate_block_count(4*mid) > block_count) {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    return 4*lo;
}

int parse_arguments(int argc, char *argv[]) {
    if (argc < 2) {
        show_usage();
        return -1;
    }

    if (is_prefix(argv[1], "create") != 0) {
        fprintf(stderr, "Command must be \"create\"\n");
        show_usage();
        return -1;
    }

    // Parse option arguments
    int arg_i = 2;
    while (arg_i < argc && argv[arg_i][0] == '-' && argv[arg_i][1] != '\0') {
        if (parse_option(argv[arg_i++]) != 0) return -1;
    }

    // First positional argument: output file
    if (arg_i == argc) {
        fprintf(stderr, "Missing output filename!\n");
        return -1;
    }
    // Defer opening output until we've opened all the inputs and verified the
    // encoding parameters, to avoid creating an empty file.
    const char *output_filename = argv[arg_i++];

    // Remaining positional arguments: input files
    if (arg_i == argc) {
        fprintf(stderr, "Missing input filenames!\n");
        return -1;
    }
    while (arg_i < argc) {
        if (open_input_file(argv[arg_i++]) != 0) {
            return -1;
        }
    }

    // Fill in default option arguments
    // The logic here is pretty hairy... this probably needs more testing!
    int nonempty_files = count_nonempty_files();
    if (nonempty_files == 0) {
        fprintf(stderr, "No nonempty files in recovery set!\n");
        return -1;
    }
    if (arg_input_blocks == -1) {
        if (arg_block_size == -1) {
            if (nonempty_files > 2000) {
                arg_input_blocks = nonempty_files;
                arg_block_size = calculate_max_file_size();
                if (arg_block_size % 4 != 0) arg_block_size += 4 - arg_block_size % 4;
            } else {
                arg_input_blocks = 2000;
                arg_block_size = calculate_block_size(arg_input_blocks);
            }
        } else if (arg_block_size <= 0 || arg_block_size % 4 != 0) {
            fprintf(stderr, "Block size (%lld) must be a multiple of 4!\n", arg_block_size);
            return -1;
        } else {
            arg_input_blocks = calculate_block_count(arg_block_size);
            if (arg_input_blocks > MAX_BLOCK_COUNT) {
                fprintf(stderr, "Too many input blocks (%lld) for block size!\n", arg_input_blocks);
                return -1;
            }
        }
    } else if (arg_block_size != -1) {
        fprintf(stderr, "Cannot set both block size (-s) and block count (-b)\n");
        return -1;
    } else if (arg_input_blocks < 1 || arg_input_blocks > MAX_BLOCK_COUNT) {
        fprintf(stderr, "Invalid input block count: %lld\n", arg_input_blocks);
        return -1;
    } else if (arg_input_blocks < nonempty_files) {
        fprintf(stderr, "Too few blocks (%lld) for number of nonempty files (%d)\n",
                arg_input_blocks, nonempty_files);
        return -1;
    } else {
        arg_block_size = calculate_block_size(arg_input_blocks);
    }

    long long calculated_block_count = calculate_block_count(arg_block_size);
    if (calculated_block_count != arg_input_blocks) {
        fprintf(stderr, "Warning: calculated block count (%lld) differs from requested block count (%lld)\n"
                "This can happen due to rounding when calculating the block size.\n",
                (long long) calculated_block_count, arg_input_blocks);
        arg_input_blocks = calculated_block_count;
    }
    if (calculated_block_count > MAX_BLOCK_COUNT) {
        fprintf(stderr, "Calculated block count (%lld) too high!\n", calculated_block_count);
        return -1;
    }
    if (arg_output_blocks == -1) {
        if (arg_redundancy == -1) {
            arg_redundancy = 5;
        } else if (arg_redundancy < 0 || arg_redundancy > 100) {
            fprintf(stderr, "Invalid redudancy percentage: %lld\n", arg_redundancy);
            return -1;
        }
        arg_output_blocks = (calculated_block_count * arg_redundancy + 50) / 100;
    } else if (arg_redundancy != -1) {
        fprintf(stderr, "Cannot set both redundancy (-r) and recovery block count (-b)\n");
        return -1;
    } else if (arg_output_blocks < 0 || arg_output_blocks > MAX_BLOCK_COUNT) {
        fprintf(stderr, "Invalid recovery block count: %lld\n", arg_input_blocks);
        return -1;
    }

    // Finally create output file.
    if (open_output_file(output_filename) != 0) return -1;

    return 0;
}

int main(int argc, char *argv[]) {
    int exit_status = 0;

    if (parse_arguments(argc, argv) != 0) goto fail;

    sort_input_files_by_id();

    // Print a brief summary of what we're about to do.
    fprintf(stderr, "Data blocks:     %10lld\n", arg_input_blocks);
    fprintf(stderr, "Recovery blocks: %10lld\n", arg_output_blocks);
    fprintf(stderr, "Block size:      %10lld (%.1f %s)\n", arg_block_size,
            ((double) arg_block_size / (arg_block_size < (1<<20) ? (1<<10) : (1<<20))),
            (arg_block_size < (1<<20) ? "KiB" : "MiB"));
    fprintf(stderr, "Redundancy:      %10.3f\n",
        arg_output_blocks > 0 ? (double) arg_input_blocks / arg_output_blocks : 0);

    // Create main packet. Must be done first to calculate the recovery set id.
    create_main_packet();

    if (write_creator_packet() != 0) {
        fprintf(stderr, "Failed to write creator packet\n");
        goto fail;
    }
    if (write_main_packet() != 0) {
        fprintf(stderr, "Failed to write main packet\n");
        goto fail;
    }
    for (int i = 0; i < input_file_count; ++i) {
        if (write_file_description_packet(&input_files[i]) != 0) {
            fprintf(stderr, "Failed to write file descriptor packet\n");
            goto fail;
        }
    }
    struct InputSlice *input_block_ptr = input_blocks;
    for (int i = 0; i < input_file_count; ++i) {
        if (input_files[i].size == 0) continue;
        if (write_ifsc_packet(&input_files[i], &input_block_ptr) != 0) {
            fprintf(stderr, "Failed to write input slice checksum packet\n");
            goto fail;
        }
    }
    assert(input_block_ptr - input_blocks == arg_input_blocks);

    if (arg_output_blocks == 0) goto finish;

    // Now the fun begins! First, initialize the tables we need:
    gf16_init();
    gf16_vandermonde_transpose_init();
    init_input_constants();

    if (reserve_output_slices(arg_output_blocks) != 0) goto fail;
    generate_recovery_data();  // this is where most time is spent
    if (finalize_output_slices(arg_output_blocks) != 0) goto fail;

    exit_status = 0;
    goto finish;

fail:
    exit_status = 1;

finish:

    // TODO: close open files, free memory, etc. (technically not necessary)
    return exit_status;
}
