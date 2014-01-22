#include "Platform.hpp"
#include "Galois256.hpp"
using namespace cat;

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
 * For packet error correction, smaller values of w are interesting for higher speed.
 * The only good option is w = 4, which allows for up to 16 blocks.  Not bad.
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
 * Jerasure is also designed to be generic, so it has best matrices for m = 2 for all of the
 * values of w that may be of interest.  But it does not attempt to optimize for m > 2, which
 * is a huge optimization that is helpful for packet error correction use cases.
 *
 * Jerasure also only tries one generator polynomial for GF(256) instead of exploring all 16
 * of the possible generators to find minimal Cauchy matrices.  I went through the extra
 * work of evaluating all the possible generators to improve on the state of the art.
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

// TODO: Optimize for w = 4 as well
// TODO: GF(16) math

// TODO: Generate best Cauchy matrices for GF(256) generator polynomials

// Cauchy codes with lowest Hamming weight for various m

static const u8 CAUCHY_BEST_2[255] = {
	1, 2, 142, 4, 71, 8, 70, 173, 3, 35, 143, 16, 17, 67, 134, 140, 172, 6, 34, 69, 201, 216,
	5, 33, 86, 12, 65, 138, 158, 159, 175, 10, 32, 43, 66, 108, 130, 193, 234, 9, 24, 25, 50,
	68, 79, 100, 132, 174, 200, 217, 20, 21, 42, 48, 87, 169, 41, 54, 64, 84, 96, 117, 154,
	155, 165, 226, 77, 82, 135, 136, 141, 168, 192, 218, 238, 7, 18, 19, 39, 40, 78, 113, 116,
	128, 164, 180, 195, 205, 220, 232, 14, 26, 27, 58, 109, 156, 157, 203, 235, 13, 28, 29
	38, 51, 56, 75, 85, 90, 101, 110, 112, 139, 171, 11, 37, 49, 52, 76, 83, 102, 119, 131,
	150, 151, 167, 182, 184, 188, 197, 219, 224, 45, 55, 80, 94, 97, 133, 170, 194, 204, 221,
	227, 236, 36, 47, 73, 92, 98, 104, 118, 152, 153, 166, 202, 207, 239, 251, 22, 23, 44, 74,
	91, 148, 149, 161, 181, 190, 233, 46, 59, 88, 137, 146, 147, 163, 196, 208, 212, 222, 250,
	57, 81, 95, 106, 111, 129, 160, 176, 199, 243, 249, 15, 53, 72, 93, 103, 115, 125, 162,
	183, 185, 189, 206, 225, 255, 186, 210, 230, 237, 242, 248, 30, 31, 62, 89, 99, 105, 114,
	121, 124, 178, 209, 213, 223, 228, 241, 254, 60, 191, 198, 247, 120, 240, 107, 127, 144,
	145, 177, 211, 214, 246, 245, 123, 126, 187, 231, 253, 63, 179, 229, 244, 61, 122, 215, 252
};

// TODO: Generate best CRS codes for m = 3, 4, 5 once the solver is working
static const u8 CAUCHY_BEST_3[255 * 2] = {
};

// TODO: Generate best CRS codes for m = 3, 4, 5 once the solver is working
static const u8 CAUCHY_BEST_4[255 * 3] = {
};

// TODO: Generate best CRS codes for m = 3, 4, 5 once the solver is working
static const u8 CAUCHY_BEST_5[255 * 4] = {
};



void cauchy_init() {
	GF256Init();
}

/*
 * Number of 1s in Cauchy 8x8 submatrix representation
 *
 * w = 8 so the Cauchy representation is an 8x8 submatrix
 * in place of the GF(256) values of the matrix.
 *
 * To generate the 8x8 submatrix, the first column is the
 * original value in binary.  And the remaining columns
 * are the column to the left times 2 in GF(256).
 */
static const u8 CAUCHY_ONES[256] = {
	0, 8, 13, 21, 18, 22, 23, 27, 20, 28, 25, 33, 26, 30, 27, 31,
	22, 26, 29, 33, 28, 28, 35, 35, 28, 32, 31, 35, 30, 30, 29, 29,
	24, 22, 29, 27, 30, 32, 35, 37, 32, 30, 29, 27, 34, 36, 35, 37,
	30, 32, 33, 35, 32, 38, 35, 41, 32, 34, 31, 33, 30, 36, 25, 31,
	27, 31, 22, 26, 33, 33, 28, 28, 31, 35, 34, 38, 37, 37, 36, 36,
	31, 31, 32, 32, 33, 29, 26, 22, 33, 33, 38, 38, 35, 31, 36, 32,
	33, 31, 32, 30, 35, 37, 34, 36, 33, 31, 40, 38, 35, 37, 38, 40,
	33, 35, 34, 36, 31, 37, 32, 38, 31, 33, 36, 38, 21, 27, 30, 36,
	30, 30, 33, 33, 22, 26, 29, 33, 36, 36, 35, 35, 28, 32, 27, 31,
	32, 28, 37, 33, 36, 36, 37, 37, 40, 36, 37, 33, 36, 36, 33, 33,
	30, 28, 33, 31, 34, 28, 33, 27, 32, 30, 31, 29, 28, 22, 19, 13,
	32, 34, 33, 35, 40, 38, 37, 35, 36, 38, 29, 31, 36, 34, 29, 27,
	35, 35, 32, 32, 35, 39, 28, 32, 37, 37, 38, 38, 33, 37, 34, 38,
	35, 31, 30, 26, 39, 39, 38, 38, 35, 31, 38, 34, 35, 35, 38, 38,
	33, 35, 34, 36, 37, 35, 34, 32, 31, 33, 36, 38, 31, 29, 36, 34,
	29, 35, 32, 38, 37, 39, 36, 38, 17, 23, 28, 34, 29, 31, 32, 34
};

static void cauchy_improve_matrix(int k, int m, u8 *matrix) {
}

static void cauchy_prebuilt_matrix(int k, int m, u8 *matrix) {
	// Attempt to use prebuilt matrices for best performance
	switch (m) {
	case 2:
		for (int x = 0; x < k; ++x) {
			matrix[x] = 1;
			matrix[x + k] = CAUCHY_BEST_2[x];
		}
		return true;
	case 3:
		for (int x = 0; x < k; ++x) {
			matrix[x] = 1;
			matrix[x + k] = CAUCHY_BEST_3[x];
			matrix[x + 2*k] = CAUCHY_BEST_3[x + 255];
		}
		return true;
	case 4:
		for (int x = 0; x < k; ++x) {
			matrix[x] = 1;
			matrix[x + k] = CAUCHY_BEST_4[x];
			matrix[x + 2*k] = CAUCHY_BEST_4[x + 255];
			matrix[x + 3*k] = CAUCHY_BEST_4[x + 255*2];
		}
		return true;
	case 5:
		for (int x = 0; x < k; ++x) {
			matrix[x] = 1;
			matrix[x + k] = CAUCHY_BEST_5[x];
			matrix[x + 2*k] = CAUCHY_BEST_5[x + 255];
			matrix[x + 3*k] = CAUCHY_BEST_5[x + 255*2];
			matrix[x + 4*k] = CAUCHY_BEST_5[x + 255*3];
		}
		return true;
	default:
		break;
	}

	// No prebuilt matrix is available
	return false;
}

bool cauchy_matrix(int k, int m, u8 *matrix) {
	// If input is invalid,
	if (k < 1 || m < 1 || !matrix ||
		k + m > 256) {
		return false;
	}

	// If 4-bit codes are possible,
	if (k + m <= 16) {
		// TODO
	}

	// Attempt to use prebuilt matrices for best performance
	if (cauchy_prebuilt_matrix(k, m, matrix)) {
		return true;
	}

	// Matrix element x, y
	// 	= 1 / (y ^ (m + x)) in GF(256)

	// For each x, y
	for (int y = 0; y < m; ++y) {
		for (int x = 0; x < k; ++x) {
			// Compute element
			matrix[c] = GF256_INV_TABLE[y ^ (m + x)];
		}

		matrix += k;
	}

	return true;
}

