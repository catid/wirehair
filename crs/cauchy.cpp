#include <iostream>
#include <cassert>
using namespace std;

#include "BitMath.hpp"
#include "MemXOR.hpp"
#include "MemSwap.hpp"
#include "AbyssinianPRNG.hpp"
#include "Clock.hpp"
using namespace cat;

static Clock m_clock;

#define DLOG(x)

void print(const u8 *data, int bytes) {
	int sep = bytes / 8;
	for (int ii = 0; ii < bytes; ++ii) {
		if (ii % sep == 0) {
			cout << ": ";
		}
		cout << (int)data[ii] << " ";
	}
	cout << endl;
}

void print_words(const u64 *row, int words) {
	for (int ii = 0; ii < words; ++ii) {
		for (int jj = 0; jj < 64; ++jj) {
			if (row[ii] & ((u64)1 << jj)) {
				cout << "1";
			} else {
				cout << "0";
			}
		}
	}
	cout << endl;
}

void print_matrix(const u64 *matrix, int word_stride, int rows) {
	cout << "Printing matrix with " << word_stride << " words per row, and " << rows << " rows:" << endl;
	for (int ii = 0; ii < rows; ++ii) {
		print_words(matrix, word_stride);
		matrix += word_stride;
	}
}

/*
 * Cauchy Reed Solomon (CRS) codes [1]
 *
 * For general purpose error correction under ~32 symbols it is either the best option,
 * or it is more flexible (due to patents/other awkwardness) than the alternatives.
 *
 * CRS codes are parameterized primarily by m, k, and w:
 * 	k = Number of original data blocks.
 * 	m = Number of redundant error correction blocks.
 * 	w = Exponent of the binary extension Galois field used.  eg. GF(2^w)
 *
 * The choice of w limits k and m by the relation: k + m <= 2^w
 * So if w = 8, it can generate up to 256 blocks of original + correction data.
 *
 * In practice if you want to send more than 256 blocks of data there are definitely more
 * efficient options than CRS codes that scale much more gracefully, so w = 8 is a
 * flexible choice that does not require incredibly large tables and does not require an
 * irritating data massaging step to fit it into the field.
 *
 * Note that m = 1 is a degenerate case where the best solution is to just XOR all of the k
 * input data blocks together.  So CRS codes are interesting for 1 < m < 32.
 *
 * These codes have been thoroughly explored by Dr. James Plank over the past ~15 years [1].
 * In this time there has not been a lot of work on improving Jerasure [2] to speed up CRS
 * codes for small datasets.
 *
 * For example, all of the existing work on Jerasure is in reference to disk or cloud
 * storage applications where the file pieces are many megabytes in size.  A neglected area
 * of interest is packet error correction codes, where the data is small and the setup time
 * for the codes is critical.
 *
 * Jerasure is designed to be generic, so it has best matrices for m = 2 for all of the
 * values of w that may be of interest.  But it does not attempt to optimize for m > 2, which
 * is a huge optimization that is helpful for packet error correction use cases.
 *
 * Jerasure only tries one generator polynomial for GF(256) instead of exploring all 16
 * of the possible generators to find minimal Cauchy matrices.  6% improvement is possible!
 *
 * Jerasure uses a "matrix improvement" formula to quickly derive an optimal Cauchy matrix
 * modified to reduce the number of ones.  I came up with a better approach that has
 * roughly 30% fewer ones in the resulting matrix, while also initializing faster.
 *
 * Jerasure is not optimized for one value of w.  It may be possible to speed up the codec
 * using w = 7, but a generic implementation that uses w = 7 will not run faster than a
 * specialized implementation that uses w = 8.
 *
 * Jerasure also misses a number of opportunities for optimization in the solver that are
 * incorporated to speed up this codec: Windowed back and forward substitution for the
 * Gaussian elimination solver, and solution shortcuts to avoid full diagonalization.
 *
 * [1] "Optimizing Cauchy Reed-Solomon Codes for Fault-Tolerant Storage Applications" (2005)
 *	http://web.eecs.utk.edu/~plank/plank/papers/CS-05-569.pdf
 * [2] "Jerasure 2.0 A Library in C/C++ Facilitating Erasure Coding for Storage Applications" (2014)
 * 	http://jerasure2.googlecode.com/svn/trunk/jerasure3/documentation/paper.pdf
 */

// GF(256) math tables:
// Generated with optimal polynomial 0x187 = 110000111

static const u16 GFC256_LOG_TABLE[256] = {
512,255,1,99,2,198,100,106,3,205,199,188,101,126,107,42,4,141,206,78,
200,212,189,225,102,221,127,49,108,32,43,243,5,87,142,232,207,172,79,131,
201,217,213,65,190,148,226,180,103,39,222,240,128,177,50,53,109,69,33,18,
44,13,244,56,6,155,88,26,143,121,233,112,208,194,173,168,80,117,132,72,
202,252,218,138,214,84,66,36,191,152,149,249,227,94,181,21,104,97,40,186,
223,76,241,47,129,230,178,63,51,238,54,16,110,24,70,166,34,136,19,247,
45,184,14,61,245,164,57,59,7,158,156,157,89,159,27,8,144,9,122,28,
234,160,113,90,209,29,195,123,174,10,169,145,81,91,118,114,133,161,73,235,
203,124,253,196,219,30,139,210,215,146,85,170,67,11,37,175,192,115,153,119,
150,92,250,82,228,236,95,74,182,162,22,134,105,197,98,254,41,125,187,204,
224,211,77,140,242,31,48,220,130,171,231,86,179,147,64,216,52,176,239,38,
55,12,17,68,111,120,25,154,71,116,167,193,35,83,137,251,20,93,248,151,
46,75,185,96,15,237,62,229,246,135,165,23,58,163,60,183};

static const u8 GFC256_EXP_TABLE[512*2+1] = {
1,2,4,8,16,32,64,128,135,137,149,173,221,61,122,244,111,222,59,118,
236,95,190,251,113,226,67,134,139,145,165,205,29,58,116,232,87,174,219,49,
98,196,15,30,60,120,240,103,206,27,54,108,216,55,110,220,63,126,252,127,
254,123,246,107,214,43,86,172,223,57,114,228,79,158,187,241,101,202,19,38,
76,152,183,233,85,170,211,33,66,132,143,153,181,237,93,186,243,97,194,3,
6,12,24,48,96,192,7,14,28,56,112,224,71,142,155,177,229,77,154,179,
225,69,138,147,161,197,13,26,52,104,208,39,78,156,191,249,117,234,83,166,
203,17,34,68,136,151,169,213,45,90,180,239,89,178,227,65,130,131,129,133,
141,157,189,253,125,250,115,230,75,150,171,209,37,74,148,175,217,53,106,212,
47,94,188,255,121,242,99,198,11,22,44,88,176,231,73,146,163,193,5,10,
20,40,80,160,199,9,18,36,72,144,167,201,21,42,84,168,215,41,82,164,
207,25,50,100,200,23,46,92,184,247,105,210,35,70,140,159,185,245,109,218,
51,102,204,31,62,124,248,119,238,91,182,235,81,162,195,1,2,4,8,16,
32,64,128,135,137,149,173,221,61,122,244,111,222,59,118,236,95,190,251,113,
226,67,134,139,145,165,205,29,58,116,232,87,174,219,49,98,196,15,30,60,
120,240,103,206,27,54,108,216,55,110,220,63,126,252,127,254,123,246,107,214,
43,86,172,223,57,114,228,79,158,187,241,101,202,19,38,76,152,183,233,85,
170,211,33,66,132,143,153,181,237,93,186,243,97,194,3,6,12,24,48,96,
192,7,14,28,56,112,224,71,142,155,177,229,77,154,179,225,69,138,147,161,
197,13,26,52,104,208,39,78,156,191,249,117,234,83,166,203,17,34,68,136,
151,169,213,45,90,180,239,89,178,227,65,130,131,129,133,141,157,189,253,125,
250,115,230,75,150,171,209,37,74,148,175,217,53,106,212,47,94,188,255,121,
242,99,198,11,22,44,88,176,231,73,146,163,193,5,10,20,40,80,160,199,
9,18,36,72,144,167,201,21,42,84,168,215,41,82,164,207,25,50,100,200,
23,46,92,184,247,105,210,35,70,140,159,185,245,109,218,51,102,204,31,62,
124,248,119,238,91,182,235,81,162,195,1,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

static const u8 GFC256_INV_TABLE[256] = {
0,1,195,130,162,126,65,90,81,54,63,172,227,104,45,42,235,155,27,53,
220,30,86,165,178,116,52,18,213,100,21,221,182,75,142,251,206,233,217,161,
110,219,15,44,43,14,145,241,89,215,58,244,26,19,9,80,169,99,50,245,
201,204,173,10,91,6,230,247,71,191,190,68,103,123,183,33,175,83,147,255,
55,8,174,77,196,209,22,164,214,48,7,64,139,157,187,140,239,129,168,57,
29,212,122,72,13,226,202,176,199,222,40,218,151,210,242,132,25,179,185,135,
167,228,102,73,149,153,5,163,238,97,3,194,115,243,184,119,224,248,156,92,
95,186,34,250,240,46,254,78,152,124,211,112,148,125,234,17,138,93,188,236,
216,39,4,127,87,23,229,120,98,56,171,170,11,62,82,76,107,203,24,117,
192,253,32,74,134,118,141,94,158,237,70,69,180,252,131,2,84,208,223,108,
205,60,106,177,61,200,36,232,197,85,113,150,101,28,88,49,160,38,111,41,
20,31,109,198,136,249,105,12,121,166,66,246,207,37,154,16,159,189,128,96,
144,47,114,133,51,59,231,67,137,225,143,35,193,181,146,79};

u8 * CAT_RESTRICT GFC256_MUL_TABLE = 0;
u8 * CAT_RESTRICT GFC256_DIV_TABLE = 0;

static void GFC256Init() {
	if (GFC256_MUL_TABLE) {
		return;
	}

	// Allocate table memory 65KB x 2
	GFC256_MUL_TABLE = (u8 *)malloc(256 * 256 * 2);
	GFC256_DIV_TABLE = GFC256_MUL_TABLE + 256 * 256;

	u8 *m = GFC256_MUL_TABLE, *d = GFC256_DIV_TABLE;

	// Unroll y = 0 subtable
	for (int x = 0; x < 256; ++x) {
		m[x] = d[x] = 0;
	}

	// For each other y value,
	for (int y = 1; y < 256; ++y) {
		// Calculate log(y) for mult and 255 - log(y) for div
		const u8 log_y = GFC256_LOG_TABLE[y];
		const u8 log_yn = 255 - log_y;

		// Next subtable
		m += 256;
		d += 256;

		// Unroll x = 0
		m[0] = 0;
		d[0] = 0;

		// Calculate x * y, x / y
		for (int x = 1; x < 256; ++x) {
			int log_x = GFC256_LOG_TABLE[x];

			m[x] = GFC256_EXP_TABLE[log_x + log_y];
			d[x] = GFC256_EXP_TABLE[log_x + log_yn];
		}
	}
}

// return x * y in GF(256)
// For repeated multiplication by a constant, it is faster to put the constant in y.
static CAT_INLINE u8 GFC256Multiply(u8 x, u8 y) {
	return GFC256_MUL_TABLE[((u32)y << 8) + x];
}

// return x / y in GF(256)
// Memory-access optimized for constant divisors in y.
static CAT_INLINE u8 GFC256Divide(u8 x, u8 y) {
	return GFC256_DIV_TABLE[((u32)y << 8) + x];
}

#include "CauchyTables.inc"

void cauchy_init() {
	GFC256Init();
}

#define CAT_CAUCHY_MATRIX_STACK_SIZE 1024

// Precondition: m > 1
static const u8 *cauchy_matrix(int k, int m, int &stride,
		u8 stack[CAT_CAUCHY_MATRIX_STACK_SIZE], bool &dynamic_memory) {
	dynamic_memory = false;

	switch (m) {
	case 2:
		stride = 254;
		return CAUCHY_MATRIX_2;
	case 3:
		stride = 253;
		return CAUCHY_MATRIX_3;
	case 4:
		stride = 252;
		return CAUCHY_MATRIX_4;
	case 5:
		stride = 251;
		return CAUCHY_MATRIX_5;
	case 6:
		stride = 250;
		return CAUCHY_MATRIX_6;
	}

	u8 *matrix = stack;
	int matrix_size = k * (m - 1);
	if (matrix_size > CAT_CAUCHY_MATRIX_STACK_SIZE) {
		matrix = new u8[matrix_size];
		dynamic_memory = true;
	}

	// Get X[] and Y[] vectors
	const u8 *Y = CAUCHY_MATRIX_Y; // Y[0] = 0
	int n = m - 7; // X[0] = 1
	const u8 *X = CAUCHY_MATRIX_X + n*249 - n*(n + 1)/2;

	//   A B C D E <- X[]
	// F 1 1 1 1 1
	// G a b c d e
	// H f g h i j
	//
	// F = 0, A = 1

	u8 *row = matrix;
	for (int y = 1; y < m; ++y) {
		u8 G = Y[y - 1];

		// Unroll x = 0
		*row++ = GFC256_INV_TABLE[1 ^ G];
		for (int x = 1; x < k; ++x) {
			u8 B = X[x - 1];

			// b = (B + F) / (B + G), F = 0
			*row++ = GFC256Divide(B, B ^ G);
		}
	}
	stride = k;

	return matrix;
}

/*
 * Cauchy encode
 */

bool cauchy_encode(int k, int m, const u8 *data, u8 *recovery_blocks, int block_bytes) {
	// If only one input block,
	if (k <= 1) {
		// Copy it directly to output
		memcpy(recovery_blocks, data, block_bytes);
	} else {
		// XOR all input blocks together
		memxor_add(recovery_blocks, data, data + block_bytes, block_bytes);
		const u8 *in = data + block_bytes;
		for (int x = 2; x < k; ++x) {
			in += block_bytes;
			memxor(recovery_blocks, in, block_bytes);
		}
	}

	// If only one recovery block needed,
	if (m == 1) {
		// We're already done!
		return true;
	}

	// Otherwise there is a restriction on what inputs we can handle
	if ((k + m > 256) || (block_bytes % 8 != 0)) {
		return false;
	}

	cauchy_init();

	// Generate Cauchy matrix
	int stride;
	u8 stack_space[CAT_CAUCHY_MATRIX_STACK_SIZE];
	bool dynamic_matrix;
	const u8 *matrix = cauchy_matrix(k, m, stride, stack_space, dynamic_matrix);

	// The first 8 rows of the bitmatrix are always the same, 8x8 identity
	// matrices all the way across.  So we don't even bother generating those
	// with a bitmatrix.  In fact the initial XOR for m=1 case has already
	// taken care of these bitmatrix rows.

	// Start on the second recovery block
	u8 *out = recovery_blocks + block_bytes;
	const int subbytes = block_bytes / 8;

	// Clear output buffer
	memset(out, 0, block_bytes * (m - 1));

	// For each remaining row to generate,
	const u8 *row = matrix;
	for (int y = 1; y < m; ++y, row += stride, out += block_bytes) {
		const u8 *src = data;

		// For each symbol column,
		const u8 *column = row;
		for (int x = 0; x < k; ++x, ++column, src += block_bytes) {
			u8 slice = column[0];
			DLOG(cout << "ENCODE: Using " << (int)slice << " at " << x << ", " << y << endl;)
			u8 *dest = out;

			// Generate 8x8 submatrix and XOR in bits as needed
			for (int bit_y = 0;; ++bit_y) {
				const u8 *src_x = src;
				for (int bit_x = 0; bit_x < 8; ++bit_x, src_x += subbytes) {
					if (slice & (1 << bit_x)) {
						memxor(dest, src_x, subbytes);
					}
				}

				if (bit_y >= 7) {
					break;
				}

				slice = GFC256Multiply(slice, 2);
				dest += subbytes;
			}
		}
	}

	if (dynamic_matrix) {
		delete []matrix;
	}

	return true;
}

/*
 * Cauchy decode
 */

// Descriptor for received data block
struct Block {
	u8 *data;
	u8 row;
};

// Specialized fast decoder for m = 1
static void cauchy_decode_m1(int k, Block *blocks, int block_bytes) {
	// Find erased row
	Block *erased = blocks;
	for (int ii = 0; ii < k; ++ii, ++erased) {
		if (erased->row >= k) {
			break;
		}
	}

	// XOR all other blocks into the recovery block
	u8 *out = erased->data;
	const u8 *in = 0;
	int original_block_count = k - 1;

	// For each block,
	for (int ii = 0; ii < original_block_count; ++ii) {
		Block *block = blocks + ii;
		if (block != erased) {
			if (!in) {
				in = block->data;
			} else {
				memxor_add(out, in, block->data, block_bytes);
				in = 0;
			}
		}
	}

	// Complete XORs
	if (in) {
		memxor(out, in, block_bytes);
	}
}

// Sort blocks into original and recovery blocks
static void sort_blocks(int k, Block *blocks,
		Block *original[256], int &original_count,
		Block *recovery[256], int &recovery_count, u8 erasures[256]) {
	Block *block = blocks;
	original_count = 0;
	recovery_count = 0;

	// Initialize erasures to zeroes
	for (int ii = 0; ii < k; ++ii) {
		erasures[ii] = 0;
	}

	// For each input block,
	for (int ii = 0; ii < k; ++ii, ++block) {
		int row = block->row;

		// If it is an original block,
		if (row < k) {
			original[original_count++] = block;
			erasures[row] = 1;
		} else {
			recovery[recovery_count++] = block;
		}
	}

	// Identify erasures
	for (int ii = 0, erasure_count = 0; erasure_count < recovery_count; ++ii) {
		if (!erasures[ii]) {
			erasures[erasure_count++] = ii;
		}
	}
}

static void eliminate_original(Block *original[256], int original_count,
							   Block *recovery[256], int recovery_count,
							   const u8 *matrix, int stride, int subbytes) {
	// If no original blocks,
	if (original_count <= 0) {
		// Nothing to do here
		return;
	}

	DLOG(cout << "Eliminating original:" << endl;)

	int row_offset = original_count + recovery_count + 1;

	// For each recovery block,
	for (int ii = 0; ii < recovery_count; ++ii) {
		Block *recovery_block = recovery[ii];
		int matrix_row = recovery_block->row - row_offset;
		const u8 *row = matrix + stride * matrix_row;

		DLOG(cout << "+ From recovery block " << ii << " at row " << matrix_row << ":" << endl;)

		// For each original block,
		for (int jj = 0; jj < original_count; ++jj) {
			Block *original_block = original[jj];
			int original_row = original_block->row;

			DLOG(cout << "++ Eliminating original column " << original_row << endl;)

			// If this matrix element is an 8x8 identity matrix,
			if (matrix_row < 0 || row[original_row] == 1) {
				// XOR whole block at once
				memxor(recovery_block->data, original_block->data, subbytes * 8);
				DLOG(cout << "XOR" << endl;)
			} else {
				// Grab the matrix entry for this row,
				u8 slice = row[original_row];
				u8 *dest = recovery_block->data;

				// XOR in bits set in 8x8 submatrix
				for (int bit_y = 0;; ++bit_y) {
					const u8 *src = original_block->data;

					for (int bit_x = 0; bit_x < 8; ++bit_x, src += subbytes) {
						if (slice & (1 << bit_x)) {
							memxor(dest, src, subbytes);
						}
					}

					// Stop after 8 bits
					if (bit_y >= 7) {
						break;
					}

					// Calculate next slice
					slice = GFC256Multiply(slice, 2);
					dest += subbytes;
				}
			}
		}
	}
}

static u64 *generate_bitmatrix(int k, Block *recovery[256], int recovery_count,
						const u8 *matrix, int stride, const u8 erasures[256],
						int &bitstride) {
	// Allocate the bitmatrix
	int bitrows = recovery_count * 8;
	bitstride = (bitrows + 63) / 64;
	u64 *bitmatrix = new u64[bitstride * bitrows];
	u64 *bitrow = bitmatrix;

	// For each recovery block,
	for (int ii = 0; ii < recovery_count; ++ii) {
		Block *recovery_block = recovery[ii];

		// If first row of matrix,
		int recovery_row = recovery_block->row - k;
		if (recovery_row == 0) {
			// Write 8x8 identity submatrix pattern across each bit row
			u64 pattern = 0x0101010101010101ULL;
			for (int ii = 0; ii < 8; ++ii, pattern <<= 1, bitrow += bitstride) {
				for (int x = 0; x < bitstride; ++x) {
					bitrow[x] = pattern;
				}
			}
		} else {
			DLOG(cout << "For recovery row " << recovery_row << endl;)
			// Otherwise read the elements of the matrix
			const u8 *row = matrix + (recovery_row - 1) * stride;

			// Generate eight 64-bit columns of the bitmatrix at a time
			int remaining = recovery_count;
			const u8 *erasure = erasures;
			while (remaining > 0) {
				// Take up to 8 columns at a time
				int limit = remaining;
				if (limit > 8) {
					limit = 8;
				}
				remaining -= limit;

				u64 w[8];

				// Unroll first loop
				u8 slice = row[*erasure++];
				w[0] = (u64)slice;

				DLOG(cout << "+ Generating 8x8 submatrix from slice=" << (int)slice << endl;)

				for (int ii = 1; ii < 8; ++ii) {
					slice = GFC256Multiply(slice, 2);
					w[ii] = (u64)slice;
				}

				// For each remaining 8 bit slice,
				for (int shift = 8; --limit > 0; shift += 8) {
					slice = row[*erasure++];
					DLOG(cout << "+ Generating 8x8 submatrix from slice=" << (int)slice << endl;)
					w[0] |= (u64)slice << shift;

					for (int ii = 1; ii < 8; ++ii) {
						slice = GFC256Multiply(slice, 2);
						w[ii] |= (u64)slice << shift;
					}
				}

				// Write 64-bit column of bitmatrix
				u64 *out = bitrow;
				for (int ii = 0; ii < 8; ++ii, out += bitstride) {
					out[0] = w[ii];
				}
			}

			bitrow += bitstride * 8;
		}

		// Set the row to what the final recovered row will be
		recovery_block->row = erasures[ii];
	}

	return bitmatrix;
}

static void gaussian_elimination(int rows, Block *recovery[256], u64 *bitmatrix, int bitstride, int subbytes) {
	const int bit_rows = rows * 8;

	u64 mask = 1;

	// For each pivot to find,
	u64 *base = bitmatrix;
	for (int pivot = 0; pivot < bit_rows - 1; ++pivot, mask = CAT_ROL64(mask, 1), base += bitstride) {
		int pivot_word = pivot >> 6;
		u64 *offset = base + pivot_word;
		u64 *row = offset;

		// For each option,
		for (int option = pivot; option < bit_rows; ++option, row += bitstride) {
			// If bit in this row is set,
			if (row[0] & mask) {
				// Prepare to add in data
				u8 *src = recovery[pivot >> 3]->data + (pivot & 7) * subbytes;
				DLOG(cout << "Found pivot " << pivot << endl;)
				DLOG(print_matrix(bitmatrix, bitstride, bit_rows);)

				// If the rows were out of order,
				if (option != pivot) {
					// Reorder data into the right place
					u8 *data = recovery[option >> 3]->data + (option & 7) * subbytes;
					memswap(src, data, subbytes);

					// Reorder matrix rows
					memswap(row, offset, (bitstride - pivot_word) << 3);
				}

				// For each other row,
				u64 *other = row;
				while (++option < bit_rows) {
					other += bitstride;

					// If that row also has the bit set,
					if (other[0] & mask) {
						DLOG(cout << "Eliminating from row " << option << endl;)
						// For each remaining word,
						for (int ii = 0; ii < bitstride - (pivot >> 6); ++ii) {
							other[ii] ^= offset[ii];
						}

						// Add in the data
						u8 *dest = recovery[option >> 3]->data + (option & 7) * subbytes;
						memxor(dest, src, subbytes);
					}
				}

				// Stop here
				break;
			}
		}
	}
}

static void back_substitution(int rows, Block *recovery[256], u64 *bitmatrix, int bitstride, int subbytes) {
	for (int pivot = rows * 8 - 1; pivot > 0; --pivot) {
		const u8 *src = recovery[pivot >> 3]->data + (pivot & 7) * subbytes;
		const u64 *offset = bitmatrix + (pivot >> 6);
		const u64 mask = (u64)1 << (pivot & 63);

		DLOG(cout << "BS pivot " << pivot << endl;)

		for (int other_row = pivot - 1; other_row >= 0; --other_row) {
			if (offset[bitstride * other_row] & mask) {
				DLOG(cout << "+ Backsub to row " << other_row << endl;)
				u8 *dest = recovery[other_row >> 3]->data + (other_row & 7) * subbytes;

				memxor(dest, src, subbytes);
			}
		}
	}
}

bool cauchy_decode(int k, int m, Block *blocks, int block_bytes) {
	// For the special case of one erasure,
	if (m == 1) {
		cauchy_decode_m1(k, blocks, block_bytes);
		return true;
	}

	// Sort blocks into original and recovery
	Block *recovery[256];
	int recovery_count;
	Block *original[256];
	int original_count;
	u8 erasures[256];
	sort_blocks(k, blocks, original, original_count, recovery, recovery_count, erasures);

	DLOG(cout << "Recovery rows(" << recovery_count << "):" << endl;
	for (int ii = 0; ii < recovery_count; ++ii) {
		cout << "+ Element " << ii << " fills in for erased row " << (int)erasures[ii] << " with recovery row " << (int)recovery[ii]->row << endl;
	}
	cout << "Original rows(" << original_count << "):" << endl;
	for (int ii = 0; ii < original_count; ++ii) {
		cout << "+ Element " << ii << " points to original row " << (int)original[ii]->row << endl;
	})

	// If nothing is erased,
	if (recovery_count <= 0) {
		return true;
	}

	// Otherwise there is a restriction on what inputs we can handle
	if ((k + m > 256) || (block_bytes % 8 != 0)) {
		return false;
	}

	// The Cauchy matrix is selected in a way that has a small
	// number of ones set in the binary representation used here.
	// A combination of precomputation and heuristics provides a
	// near-optimal matrix selection for each value of k, m.

	cauchy_init();

	// Generate Cauchy matrix
	int stride;
	u8 stack_space[CAT_CAUCHY_MATRIX_STACK_SIZE];
	bool dynamic_matrix;
	const u8 *matrix = cauchy_matrix(k, m, stride, stack_space, dynamic_matrix);

	// From the Cauchy matrix, each byte value can be expanded into
	// an 8x8 submatrix containing a minimal number of ones.
	// The rows that made it through from the original data provide
	// some of the column values for the matrix, so those can be
	// eliminated immediately.  This is useful because it conceptually
	// zeroes out those eliminated matrix elements.  And so when it
	// comes time to laborously generate the bitmatrix and solve it
	// with Gaussian elimination, that bitmatrix can be smaller since
	// it does not need to include these rows and columns.

	// Eliminate original data from recovery rows
	const int subbytes = block_bytes / 8;
	eliminate_original(original, original_count, recovery, recovery_count, matrix, stride, subbytes);

	// Now that the columns that are missing have been identified,
	// it is time to generate a bitmatrix to represent the original
	// rows that have been XOR'd together to produce the recovery data.
	// This matrix is guaranteed to be inverible as it was selected
	// from the rows/columns of a Cauchy matrix.

	// Generate square bitmatrix for erased columns from recovery rows
	int bitstride;
	u64 *bitmatrix = generate_bitmatrix(k, recovery, recovery_count, matrix,
										stride, erasures, bitstride);

	DLOG(print_matrix(bitmatrix, bitstride, recovery_count * 8);)

	// Finally, solving the matrix.
	// The most efficient approach is Gaussian elimination: An alternative
	// would be to recursively solve submatrices.  However, since the initial
	// matrix is sparse it is undesirable to add matrix rows together.
	// By working to put the matrix in upper-triangular form, the number of
	// row additions is reduced by about half.  And then a solution can be
	// immediately found without performing more row additions.

	// Gaussian elimination to put matrix in upper triangular form
	gaussian_elimination(recovery_count, recovery, bitmatrix, bitstride, subbytes);

	DLOG(print_matrix(bitmatrix, bitstride, recovery_count * 8);)

	// The matrix is now in an upper-triangular form, and can be worked from
	// right to left to conceptually produce an identity matrix.  The matrix
	// itself is not adjusted since the important result is the output values.

	// Use back-substitution to solve value for each column
	back_substitution(recovery_count, recovery, bitmatrix, bitstride, subbytes);

	// Free temporary space
	delete []bitmatrix;

	if (dynamic_matrix) {
		delete []matrix;
	}

	return true;
}

int main() {
	m_clock.OnInitialize();

	cauchy_init();

	m_clock.usec();

	cout << "Cauchy matrix solver" << endl;

	int block_bytes = 8 * 162; // a multiple of 8
	int block_count = 180;
	int recovery_block_count = 72;

	u8 *data = new u8[block_bytes * block_count];
	u8 *recovery_blocks = new u8[block_bytes * recovery_block_count];

	Abyssinian prng;
	prng.Initialize(0);
	for (int ii = 0; ii < block_bytes * block_count; ++ii) {
		data[ii] = (u8)prng.Next();
	}

	double t0 = m_clock.usec();

	assert(cauchy_encode(block_count, recovery_block_count, data, recovery_blocks, block_bytes));

	double t1 = m_clock.usec();

	cout << "Cauchy encode in " << (t1 - t0) << " usec" << endl;

	Block *blocks = new Block[block_count];

	DLOG(cout << "Original data:" << endl;)
	for (int ii = 0; ii < block_count; ++ii) {
		blocks[ii].data = data + ii * block_bytes;
		blocks[ii].row = ii;
		DLOG(print(blocks[ii].data, block_bytes);)
	}

	// Erase first block
	const int erasures_count = 7;
	int original_remaining = block_count - erasures_count;
	int erasures[256];
	for (int ii = 0; ii < erasures_count; ++ii) {
		blocks[ii].data = recovery_blocks + ii * block_bytes;
		blocks[ii].row = block_count + ii;
		erasures[ii] = ii;
	}

	t0 = m_clock.usec();

	assert(cauchy_decode(block_count, recovery_block_count, blocks, block_bytes));

	t1 = m_clock.usec();

	for (int ii = 0; ii < erasures_count; ++ii) {
		int erasure_index = erasures[ii];
		cout << "Data erasure " << ii << " and row=" << (int)blocks[erasure_index].row << endl;
		DLOG(print(blocks[erasure_index].data, block_bytes);)
	}
	for (int ii = 0; ii < erasures_count; ++ii) {
		int erasure_index = erasures[ii];
		cout << "At row " << (int)blocks[erasure_index].row << endl;
		assert(!memcmp(blocks[erasure_index].data, data + blocks[erasure_index].row * block_bytes, block_bytes));
	}

	cout << "Cauchy decode in " << (t1 - t0) << " usec" << endl;

	m_clock.OnFinalize();

	delete []data;
	delete []recovery_blocks;

	return 0;
}

