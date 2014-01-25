#include <iostream>
#include <cassert>
using namespace std;

#include "BitMath.hpp"
#include "MemXOR.hpp"
#include "AbyssinianPRNG.hpp"
#include "Clock.hpp"
using namespace cat;

static Clock m_clock;

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

// Optimal Cauchy matrices for some small values of m:

static const u8 CAUCHY_MATRIX_2[1 * 254] = {
1,195,2,4,162,81,8,194,3,97,6,163,5,10,12,20,80,40,235,16,
146,193,24,73,48,243,9,96,160,18,36,199,182,192,72,231,186,89,178,32,
176,17,166,83,234,227,69,138,7,147,161,203,21,88,65,225,13,197,11,41,
34,93,85,91,201,25,242,44,198,22,14,82,144,64,233,117,170,239,179,49,
99,164,28,183,68,167,113,121,67,226,42,84,154,207,77,168,215,230,134,152,
19,211,26,33,98,202,219,101,180,241,187,56,37,251,209,92,237,90,139,50,
130,76,75,174,229,177,112,249,46,171,145,38,150,238,100,128,116,115,232,196,
87,158,224,184,45,66,52,58,190,105,23,15,30,60,74,120,200,155,29,247,
165,181,210,136,255,240,79,142,71,214,250,206,131,43,103,169,35,104,228,148,
213,109,27,151,188,107,95,70,86,57,135,114,218,153,205,175,191,245,119,208,
51,140,132,185,236,248,39,159,143,129,125,246,172,54,78,137,94,217,53,102,
156,223,118,204,124,254,47,106,212,123,108,61,149,59,133,253,31,221,189,173,
122,62,216,127,141,157,244,222,252,111,55,126,110,220};

static const u8 CAUCHY_MATRIX_3[2 * 253] = {
4,16,81,6,5,80,162,83,8,235,48,163,9,178,195,33,1,103,40,161,
3,177,2,156,69,231,171,12,166,11,199,233,229,182,226,96,146,74,205,28,
138,97,73,58,160,46,174,18,186,84,35,113,24,17,164,21,185,10,179,201,
194,145,72,106,90,147,168,44,101,93,181,75,192,211,13,198,234,150,79,32,
124,20,219,239,215,167,142,89,228,176,183,193,36,119,65,255,130,253,202,14,
170,29,105,251,61,77,52,249,227,191,225,107,98,243,190,42,68,45,85,148,
88,236,100,53,30,49,66,165,99,137,7,245,151,188,187,196,136,76,139,125,
64,242,203,56,197,27,120,238,200,115,70,54,184,67,91,128,232,82,246,60,
153,131,26,41,92,143,38,240,254,117,118,25,204,109,34,144,39,207,108,112,
209,213,50,248,155,214,172,132,51,230,180,237,210,221,104,22,169,15,57,71,
134,47,121,19,62,149,23,152,159,31,122,127,241,154,95,59,208,140,86,158,
223,216,114,250,37,189,94,87,218,212,43,217,252,224,55,135,222,173,206,116,
244,247,102,63,78,111,141,175,126,129,110,220,157,
// For row 2:
16,4,6,81,80,5,83,162,178,48,235,9,163,8,33,195,103,1,161,40,
177,3,156,2,231,69,12,171,11,166,233,199,96,226,182,229,74,146,97,138,
28,205,58,73,46,160,18,174,84,186,24,17,35,113,21,164,10,185,201,179,
106,72,145,194,147,90,44,168,93,101,192,13,181,198,75,211,150,234,32,79,
20,124,239,219,167,215,176,228,89,142,130,253,119,36,255,65,183,193,251,29,
105,14,170,202,243,52,77,98,191,227,107,225,249,61,42,190,45,68,148,85,
236,88,196,7,139,188,76,187,151,203,53,242,99,49,165,100,56,66,30,197,
27,245,137,136,125,64,128,60,115,200,67,91,232,70,54,120,184,246,82,238,
26,92,153,254,131,117,240,38,41,143,144,204,25,112,108,118,207,39,34,109,
213,209,132,180,214,155,230,50,237,172,248,51,169,22,71,221,210,57,15,104,
47,134,62,149,121,19,159,31,23,152,241,154,122,127,208,158,95,86,140,59,
250,37,94,223,216,87,114,189,217,43,212,218,224,252,116,173,206,135,222,55,
102,63,244,247,111,78,133,126,175,110,129,157,220};

static const u8 CAUCHY_MATRIX_4[3 * 252] = {
195,2,1,65,149,99,34,81,16,163,186,72,224,243,86,148,242,246,38,25,
41,191,24,182,194,215,162,12,73,234,69,245,5,75,84,20,141,218,36,187,
96,192,6,10,205,166,139,152,37,7,198,44,17,30,174,105,201,74,144,168,
142,235,97,26,91,179,13,164,115,117,33,22,167,67,214,138,4,107,29,40,
211,66,32,171,236,82,134,160,46,76,231,239,28,60,8,209,172,189,176,146,
58,184,165,153,177,123,87,68,130,83,51,193,180,3,203,247,109,178,47,131,
19,89,132,85,219,101,95,230,248,48,54,233,199,77,112,232,79,170,225,135,
183,137,92,197,57,244,161,125,237,129,227,64,147,45,106,23,108,114,35,56,
188,252,212,185,143,49,229,61,181,70,226,196,157,93,104,120,71,202,9,18,
80,238,208,145,158,116,27,255,249,217,100,11,15,228,204,63,113,121,14,59,
150,156,254,240,200,241,220,52,39,210,88,53,206,175,253,50,154,155,124,98,
119,216,173,222,31,190,94,207,250,102,42,128,90,251,21,78,122,111,223,136,
55,213,221,126,140,118,127,110,169,159,43,133,
// For row 2:
167,6,22,235,1,96,163,138,198,241,8,160,251,97,80,162,186,72,3,85,
32,4,82,195,31,10,206,42,68,89,227,41,19,178,90,144,5,81,12,73,
74,114,106,117,69,77,255,158,16,53,109,193,15,116,2,30,93,20,65,171,
9,176,37,83,197,177,147,203,17,192,145,170,226,91,161,199,248,234,231,27,
224,153,223,48,230,215,183,148,36,44,7,188,50,23,189,243,24,194,239,156,
155,76,14,18,46,52,28,211,168,242,40,70,218,108,238,146,84,94,213,64,
164,233,112,174,205,182,33,121,99,204,184,151,217,38,254,225,240,247,61,79,
26,88,104,103,107,172,71,21,169,115,51,54,139,180,92,232,219,58,137,135,
214,179,159,166,105,45,249,49,87,86,196,236,209,185,128,35,134,102,228,124,
57,133,100,34,142,39,140,190,245,187,60,126,66,132,119,98,149,123,216,25,
122,165,136,101,200,152,130,250,120,127,220,191,202,13,207,252,208,210,229,222,
253,237,56,201,154,95,29,47,212,157,111,62,63,43,246,143,131,11,244,75,
67,150,181,173,129,55,78,113,59,125,141,110,
// For row 3:
81,154,192,4,2,162,146,34,73,163,230,186,1,46,194,182,96,195,83,235,
178,234,121,221,40,49,67,98,28,179,25,8,211,231,193,45,80,187,102,101,
85,18,89,35,243,164,97,48,145,5,24,151,138,10,62,6,68,120,19,36,
183,149,245,152,38,225,103,128,167,104,32,90,226,251,44,156,208,203,215,74,
227,16,147,86,9,158,168,87,52,209,106,17,26,199,155,94,58,115,51,116,
176,14,202,15,99,12,50,232,249,123,236,175,22,119,214,212,64,76,20,130,
153,173,65,188,233,31,170,95,207,131,72,139,105,71,166,78,144,79,42,11,
140,66,237,135,161,3,54,60,75,242,206,112,59,210,197,150,91,239,7,171,
198,69,160,191,117,129,218,241,23,93,39,174,13,57,200,184,132,100,252,53,
110,113,30,244,205,190,219,213,181,136,159,37,125,238,82,41,142,70,250,185,
177,196,229,133,223,111,88,143,228,77,224,180,253,216,114,134,254,124,107,247,
84,240,27,172,246,217,109,108,165,201,148,169,33,127,55,255,29,47,56,220,
137,126,118,21,204,92,43,61,189,222,122,248};

static const u8 CAUCHY_MATRIX_5[4 * 251] = {
81,227,178,171,4,24,46,101,67,243,83,10,195,96,194,162,43,228,32,99,
229,26,70,12,134,125,213,235,25,242,14,1,76,97,85,190,174,36,94,37,
88,112,208,48,5,8,197,209,182,84,6,80,128,79,65,53,152,234,92,11,
35,187,161,3,169,16,40,86,51,52,148,241,44,9,20,191,163,160,131,193,
93,186,90,144,55,56,72,226,71,33,202,100,206,89,41,217,2,58,95,222,
179,214,196,175,68,192,199,249,236,168,78,34,210,30,238,136,137,27,211,176,
172,215,124,147,13,50,60,28,114,138,180,248,166,183,143,164,145,207,212,75,
73,159,19,126,150,140,139,122,54,77,231,49,253,252,237,225,239,117,17,91,
198,64,69,129,155,181,7,146,121,39,15,233,221,188,156,113,29,251,245,98,
57,224,38,107,62,232,158,185,59,118,115,200,135,21,109,250,18,247,255,66,
177,230,31,45,102,151,203,105,170,104,189,204,149,167,87,216,74,219,165,184,
244,127,106,130,116,201,173,223,103,108,120,111,154,218,205,254,119,22,153,47,
157,23,133,61,42,63,240,142,132,246,110,
// For row 2:
229,4,32,162,227,213,10,102,48,44,12,46,11,6,14,171,3,49,178,112,
81,234,235,83,86,80,24,70,111,20,194,76,1,241,163,207,193,239,30,68,
167,99,176,67,197,79,5,35,34,202,96,125,91,8,98,237,50,26,233,195,
209,144,248,43,211,191,180,134,7,192,9,97,243,148,242,16,85,185,199,174,
244,251,56,187,146,90,60,108,154,166,84,117,156,188,179,224,222,138,142,2,
41,231,143,105,37,52,131,255,215,61,204,182,198,94,145,147,22,158,169,208,
23,236,245,136,223,152,72,73,139,58,40,161,33,121,196,225,238,190,177,250,
28,230,17,170,153,47,114,205,120,181,214,228,18,218,53,164,36,100,19,128,
210,59,155,232,69,77,51,55,183,42,115,92,203,89,206,149,106,186,124,65,
200,217,135,21,130,129,27,160,64,165,15,57,38,107,151,75,253,189,249,103,
212,159,132,216,101,109,221,175,126,127,247,78,113,88,63,45,116,119,118,110,
93,104,29,62,74,254,240,13,66,226,54,25,71,252,122,201,219,137,150,140,
133,172,157,168,39,87,173,95,31,123,184,
// For row 3:
48,64,8,81,227,162,88,80,166,20,134,97,139,37,198,140,33,4,65,3,
45,68,195,96,1,113,77,209,72,29,224,188,60,109,231,5,25,50,2,154,
9,28,193,24,168,32,19,192,83,10,70,243,15,249,187,112,165,235,41,172,
12,242,186,213,6,16,202,22,161,55,155,248,92,190,244,158,39,215,201,251,
197,160,93,233,179,106,62,171,163,94,13,240,255,191,133,147,79,71,18,98,
87,210,135,73,237,90,176,128,85,199,49,52,229,152,100,14,75,196,86,170,
17,184,101,102,226,91,143,51,206,53,204,111,159,167,182,253,183,219,11,34,
127,44,115,138,67,203,132,40,234,74,217,59,84,146,180,105,222,108,150,153,
47,117,252,82,211,136,151,142,241,36,103,95,46,27,178,26,218,250,194,54,
76,116,114,126,205,225,130,122,121,131,221,137,23,174,58,38,78,69,118,145,
35,149,144,7,239,212,223,43,177,148,207,238,57,129,42,104,208,156,124,21,
119,232,169,236,254,247,230,173,228,181,200,89,110,56,99,220,31,185,216,66,
107,61,246,245,125,214,30,175,120,189,123,
// For row 4:
81,8,16,233,112,25,40,1,18,101,72,98,209,183,168,12,235,138,14,21,
36,6,236,39,219,20,24,48,2,179,87,45,205,186,159,34,26,69,242,193,
128,160,22,31,188,218,214,9,35,105,95,52,146,56,227,195,194,119,197,78,
144,198,225,88,171,53,153,96,91,162,13,33,27,82,117,28,15,154,196,114,
243,122,134,116,41,5,73,163,129,85,206,65,4,178,203,83,136,67,113,217,
135,226,3,150,50,131,254,147,202,211,97,204,145,192,99,143,170,32,199,240,
77,11,80,201,177,191,237,253,231,215,94,10,165,212,71,38,79,103,152,255,
132,155,59,164,130,176,89,232,158,238,76,7,185,42,84,221,228,180,244,169,
207,58,44,200,189,190,184,100,127,61,229,86,37,109,29,74,182,216,47,106,
142,92,187,17,161,172,151,43,51,93,239,121,148,157,139,248,107,57,224,175,
174,49,90,70,108,234,241,30,46,166,23,68,210,110,75,167,125,181,64,250,
247,251,66,249,140,133,213,104,156,118,102,126,120,60,54,19,245,222,208,62,
115,124,230,55,63,246,111,220,252,137,123};

static const u8 CAUCHY_MATRIX_6[5 * 250] = {
120,3,193,22,16,87,2,233,6,239,10,101,20,65,195,179,145,175,232,38,
99,182,100,91,40,49,171,69,1,81,83,8,48,139,64,247,14,166,183,186,
19,76,25,79,39,237,93,157,188,189,184,197,177,90,227,4,77,165,89,163,
167,58,243,97,209,56,52,133,86,246,88,190,162,68,82,98,221,18,178,107,
95,240,159,147,250,12,192,146,134,80,172,201,181,13,199,29,21,155,74,169,
226,7,41,219,28,200,141,211,36,224,118,109,245,137,47,92,54,15,158,73,
115,238,198,17,127,35,185,24,26,150,70,63,206,116,229,176,71,104,75,128,
42,60,234,5,138,53,112,225,114,205,156,30,94,241,108,154,110,34,72,113,
117,187,27,67,164,228,235,170,96,23,11,216,230,161,129,106,168,135,119,144,
61,33,160,66,103,123,255,244,142,84,210,231,31,102,121,236,191,46,32,194,
130,208,105,203,140,222,85,152,51,151,149,202,122,59,153,215,214,173,37,204,
148,9,126,44,217,180,43,196,45,242,174,57,125,252,50,223,111,143,207,251,
55,212,253,249,78,254,131,62,136,248,
// For row 2:
81,162,168,178,97,194,8,4,35,152,93,24,25,10,139,243,96,161,51,40,
227,134,228,3,2,99,101,208,49,13,187,144,9,43,84,54,29,88,32,15,
56,36,193,37,125,80,75,209,171,235,72,85,112,190,172,27,195,83,117,203,
48,6,159,250,78,140,46,175,26,206,91,53,200,202,74,217,92,207,57,16,
77,181,226,28,160,68,224,240,14,64,12,86,18,174,215,67,114,20,130,179,
44,237,138,90,177,65,198,199,186,89,196,184,146,1,70,94,155,95,147,113,
124,231,79,107,176,156,45,116,197,7,158,234,41,52,210,38,249,242,252,17,
238,212,60,205,214,145,30,167,19,213,131,69,229,211,58,137,5,39,149,251,
239,109,164,241,129,182,23,62,105,143,191,136,151,55,189,128,122,135,21,108,
163,232,165,157,71,204,11,169,183,118,66,22,104,170,76,150,148,115,102,153,
233,180,173,254,120,103,98,188,87,236,50,216,126,31,247,106,255,121,132,218,
245,142,192,59,219,111,100,154,225,222,73,221,230,127,119,33,166,133,63,244,
248,34,47,253,185,42,110,201,61,123,
// For row 3:
195,48,66,64,97,41,102,71,201,6,194,100,8,15,227,19,85,138,235,237,
32,98,161,231,169,84,4,88,28,77,22,87,207,162,35,241,81,137,56,154,
20,147,79,10,80,208,242,36,170,192,163,245,9,136,179,128,74,167,89,202,
124,251,229,238,3,73,160,209,24,16,45,2,99,54,115,203,26,27,177,212,
182,113,5,250,225,200,247,59,219,232,156,7,43,104,44,72,146,83,234,149,
108,215,82,145,49,130,50,246,111,65,12,230,23,180,190,198,125,67,42,175,
199,183,224,68,206,116,159,185,86,105,213,91,34,1,205,101,197,69,168,222,
11,233,186,216,134,187,14,63,236,96,30,176,144,18,139,33,121,244,70,110,
158,165,152,132,46,140,119,13,120,155,112,181,166,164,92,193,174,217,253,95,
60,57,135,75,141,151,191,37,94,243,53,118,17,184,58,129,52,178,55,78,
223,255,90,76,31,249,47,51,171,38,117,21,25,239,228,143,142,123,254,240,
188,126,133,204,93,127,131,103,122,106,153,148,196,226,218,211,29,150,214,210,
40,62,248,109,172,252,107,157,220,61,
// For row 4:
48,140,199,8,109,198,32,227,12,165,197,162,72,97,132,20,37,186,161,202,
64,1,4,213,79,3,80,193,59,226,242,233,190,33,10,234,218,9,65,103,
106,50,251,154,113,243,34,192,81,209,62,231,28,5,17,196,139,134,108,223,
24,70,44,38,49,203,88,73,68,255,153,112,137,13,208,147,41,219,76,16,
74,136,171,51,215,237,116,30,224,117,96,22,78,25,184,166,206,244,236,87,
92,180,53,93,35,187,47,176,160,191,135,21,142,188,195,2,211,18,102,26,
101,217,249,126,170,178,7,254,19,151,130,235,133,55,229,114,128,29,146,150,
100,11,143,99,210,183,152,129,115,77,201,252,45,86,71,75,168,36,57,250,
222,58,253,248,82,83,61,205,43,182,158,14,212,179,207,15,40,23,174,181,
39,225,124,107,163,238,172,6,167,131,145,185,148,177,60,67,155,221,239,216,
95,204,230,220,200,228,54,27,42,85,91,104,138,144,69,169,118,241,120,56,
194,175,90,121,156,89,240,110,105,98,127,46,149,232,31,94,159,246,214,119,
111,52,66,84,122,125,123,247,245,141,
// For row 5:
209,69,1,80,14,194,187,235,18,2,23,242,86,85,12,231,4,3,81,16,
34,189,170,101,95,87,61,225,19,188,83,150,74,183,82,195,25,168,56,10,
160,215,99,49,5,136,67,243,75,65,186,91,50,21,11,98,106,88,233,92,
218,201,64,226,249,52,255,97,116,73,176,79,205,89,146,28,24,238,22,228,
144,40,60,17,217,123,33,178,100,252,103,132,36,120,124,47,153,181,93,6,
13,129,114,190,237,155,162,191,108,131,113,42,112,244,148,220,8,240,241,118,
152,229,180,20,72,166,182,239,165,211,9,119,198,125,234,210,66,57,142,197,
102,251,53,71,169,121,156,193,147,159,41,232,177,222,199,105,248,37,161,48,
76,77,44,117,137,216,58,130,221,175,246,203,200,54,227,59,111,96,230,164,
109,122,206,192,107,163,196,133,51,141,84,149,38,185,62,128,219,39,174,151,
43,179,184,32,138,7,250,139,212,172,157,70,207,171,254,154,167,202,26,158,
63,45,90,140,223,224,94,145,253,27,126,30,208,68,214,127,204,15,35,245,
213,173,115,55,134,135,46,110,236,104};

// TODO: Expand these precomputed optimal tables by reconstructing them given
// the optimal column selections etc during generation

void cauchy_init() {
	GFC256Init();
}

static const u8 *cauchy_matrix(int k, int m, int &stride) {
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

	// TODO: Fill this in

	return 0;
}

/*
 * BitMatrix has 8*m rows and 8*k columns
 *
 * Each 8x8 submatrix is a transposed bit matrix as in Jerasure's implementation.
 * First row is the GF(256) symbol, and each following symbol is the previous
 * one multiplied by 2.
 */

void cauchy_expand(int k, int m, const u8 *matrix, int stride, u64 *bitmatrix,
				   int bitstride, u16 *row_ones) {
	// While there are more columns to write,
	u64 *bitrow = bitmatrix;
	int x = k;
	while (x > 0) {
		int limit = x;
		if (limit > 8) {
			limit = 8;
		}
		x -= limit;

		u64 w = 0x0101010101010101ULL;

		if (limit < 8) {
			w >>= (8 - limit) * 8;
		}

		// Write 64-bit column of bitmatrix
		bitrow[0] = w;
		bitrow[bitstride] = w << 1;
		bitrow[bitstride*2] = w << 2;
		bitrow[bitstride*3] = w << 3;
		bitrow[bitstride*4] = w << 4;
		bitrow[bitstride*5] = w << 5;
		bitrow[bitstride*6] = w << 6;
		bitrow[bitstride*7] = w << 7;
		++bitrow;
	}

	// For each row of input (excludes the first row of all ones),
	const u8 *row = matrix;
	for (int y = 1; y < m; ++y, row += stride) {
		bitrow += bitstride * 7;

		// Initialize count of ones for each of the 8 rows
		row_ones += 8;
		for (int ii = 0; ii < 8; ++ii) {
			row_ones[ii] = 0;
		}

		// While there are more columns to write,
		x = k;
		while (x > 0) {
			int limit = x;
			if (limit > 8) {
				limit = 8;
			}
			x -= limit;

			// Generate low 8 bits of the word
			int shift = limit * 8;
			u8 slice = *row++;
			u64 w[8];
			w[0] = (u64)slice << shift;
			slice = GFC256Multiply(slice, 2);
			w[1] = (u64)slice << shift;
			slice = GFC256Multiply(slice, 2);
			w[2] = (u64)slice << shift;
			slice = GFC256Multiply(slice, 2);
			w[3] = (u64)slice << shift;
			slice = GFC256Multiply(slice, 2);
			w[4] = (u64)slice << shift;
			slice = GFC256Multiply(slice, 2);
			w[5] = (u64)slice << shift;
			slice = GFC256Multiply(slice, 2);
			w[6] = (u64)slice << shift;
			slice = GFC256Multiply(slice, 2);
			w[7] = (u64)slice << shift;

			// For each remaining 8 bit chunk,
			while (--limit > 0) {
				u8 slice = *row++;
				w[0] = (w[0] >> 8) | ((u64)slice << shift);
				slice = GFC256Multiply(slice, 2);
				w[1] = (w[1] >> 8) | ((u64)slice << shift);
				slice = GFC256Multiply(slice, 2);
				w[2] = (w[2] >> 8) | ((u64)slice << shift);
				slice = GFC256Multiply(slice, 2);
				w[3] = (w[3] >> 8) | ((u64)slice << shift);
				slice = GFC256Multiply(slice, 2);
				w[4] = (w[4] >> 8) | ((u64)slice << shift);
				slice = GFC256Multiply(slice, 2);
				w[5] = (w[5] >> 8) | ((u64)slice << shift);
				slice = GFC256Multiply(slice, 2);
				w[6] = (w[6] >> 8) | ((u64)slice << shift);
				slice = GFC256Multiply(slice, 2);
				w[7] = (w[7] >> 8) | ((u64)slice << shift);
			}

			// Write 64-bit column of bitmatrix
			bitrow[0] = w[0];
			bitrow[bitstride] = w[1];
			bitrow[bitstride*2] = w[2];
			bitrow[bitstride*3] = w[3];
			bitrow[bitstride*4] = w[4];
			bitrow[bitstride*5] = w[5];
			bitrow[bitstride*6] = w[6];
			bitrow[bitstride*7] = w[7];
			++bitrow;

			// Update the number of bits set in the row
			row_ones[0] += BitCount(w[0]);
			row_ones[1] += BitCount(w[1]);
			row_ones[2] += BitCount(w[2]);
			row_ones[3] += BitCount(w[3]);
			row_ones[4] += BitCount(w[4]);
			row_ones[5] += BitCount(w[5]);
			row_ones[6] += BitCount(w[6]);
			row_ones[7] += BitCount(w[7]);
		}
	}
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

	// Temporarily allocate more space if needed
	int bitstride = (k*8 + 63) / 64;
	u16 row_ones[256];
	u64 static_bitmatrix[128];
	u64 *bitmatrix;
	if (bitstride * m * 8 > 128) {
		bitmatrix = new u64[bitstride * m * 8];
	} else {
		bitmatrix = static_bitmatrix;
	}

	// Generate Cauchy matrix
	int stride;
	const u8 *matrix = cauchy_matrix(k, m, stride);

	// Expand Cauchy matrix into bit matrix
	cauchy_expand(k, m, matrix, stride, bitmatrix, bitstride, row_ones);

	// The first 8 rows of the bitmatrix are always the same, 8x8 identity
	// matrices all the way across.  So we don't even bother generating those
	// with a bitmatrix.  In fact the initial XOR for m=1 case has already
	// taken care of these bitmatrix rows.

	// Start on the second recovery block
	u8 *out = recovery_blocks + block_bytes;
	u64 *bitmatrix_row = bitmatrix + bitstride * 8;
	int subbytes = block_bytes / 8;

	// For each remaining row to generate,
	for (int y = 8; y < m*8; ++y, bitmatrix_row += bitstride, out += subbytes) {
		// Find smaller XOR count from previous rows
		u64 *prev_row = bitmatrix;
		int lowest_src_row, lowest_xors = row_ones[y];
		for (int prev = 0; prev < y; ++prev, prev_row += bitstride) {
			// Calculate the number of XORs required if we started from this row
			int xors = BitCount(prev_row[0] ^ bitmatrix_row[0]);
			for (int x = 1; x < bitstride; ++x) {
				xors += BitCount(prev_row[x] ^ bitmatrix_row[x]);
			}

			// Use it if it's lower
			if (lowest_xors > xors) {
				lowest_xors = xors;
				lowest_src_row = prev;
			}
		}

		// Number of bytes per sub-block
		const u8 *in = 0;

		// If it's better to start with an existing recovery block,
		if (lowest_xors < row_ones[y]) {
			prev_row = bitmatrix + lowest_src_row * bitstride;

			// Generate the combined matrix row
			for (int x = 0; x < bitstride; ++x) {
				bitmatrix_row[x] ^= prev_row[x];
			}

			// Copy the existing recovery block
			in = recovery_blocks + subbytes * lowest_src_row;
		}

		// XOR sub-blocks together
		const u8 *src = data;
		bool out_is_set = false;

		// For each bit in the bitmatrix row,
		for (int x = 0; x < bitstride; ++x) {
			u64 word = bitmatrix_row[x];
			for (int bit = 0; bit < 64; ++bit, src += subbytes) {
				// If bit is set,
				if (word & (1 << bit)) {
					// Set up operation
					if (!in) {
						in = src;
					} else {
						// Perform dual XOR operation
						if (!out_is_set) {
							memxor_set(out, in, src, subbytes);
							out_is_set = true;
						} else {
							memxor_add(out, in, src, subbytes);
						}
						in = 0;
					}
				}
			}
		}

		// Write whatever is left over
		if (in) {
			// If out is already written,
			if (out_is_set) {
				// XOR the remaining stuff in
				memxor(out, in, subbytes);
			} else {
				// Copy the input (there should always be at least one to copy)
				memcpy(out, in, subbytes);
			}
		}
	}

	// Free temporary space
	if (bitmatrix != static_bitmatrix) {
		delete []bitmatrix;
	}

	return true;
}

/*
 * Received blocks are passed in.  There will always be k of them.
 */

struct ReceivedBlock {
	u8 *data;
	u8 row;
};

static void cauchy_decode_m1(int k, ReceivedBlock *original_blocks, ReceivedBlock *recovery_blocks, int block_bytes) {
	// XOR all other blocks into the recovery block
	u8 *out = recovery_blocks[0].data;
	const u8 *in = 0;
	int original_block_count = k - 1;

	// For each block,
	for (int ii = 0; ii < original_block_count; ++ii) {
		const u8 *src = original_blocks[ii].data;

		if (!in) {
			in = src;
		} else {
			memxor_add(out, in, src, block_bytes);
			in = 0;
		}
	}

	// Complete XORs
	if (in) {
		memxor(out, in, block_bytes);
	}
}

/*
 * Interested in recovering some rows
 *
 * For each of the recovery rows, perform Gaussian elimination.  Other rows can
 * be skipped.
 */

bool cauchy_decode(int k, int m, ReceivedBlock *original_blocks, ReceivedBlock *recovery_blocks, int recovery_block_count, int block_bytes) {
	// If nothing is erased,
	if (recovery_block_count < 0) {
		return true;
	}

	// For the special case of one erasure,
	if (m == 1) {
		cauchy_decode_m1(k, original_blocks, recovery_blocks, block_bytes);
		return true;
	}

	int original_block_count = k - recovery_block_count;

	// TODO: Generate bitmatrix.  We only need rows that are part of the recovery
	// set, so that should cut down on the work of decoding.

	u8 pivots[8*256];

	// For solution, we can quickly set pivots for the rows that correspond to
	// original data.  This sets regions of 8 pivots quickly.  We want to avoid
	// actually generating bitmatrix rows for these.

	// Now search across from left to right looking for pivots.  If the original
	// data row is available, then use that data to eliminate columns of the
	// recovery rows.

	// If a pivot is not found, then search the recovery rows.  One should be
	// found if the matrix is invertible (100% likelihood unless i messed up
	// the matrices).  Once a pivot is found, eliminate it from all the other
	// rows immediately.

	// During forward pivot search, the source of XORs is moving forward from
	// the front of the data in a one-to-many pattern.  This is the best you
	// can do, so there's no reason to schedule it.

	// At the end of this process, the matrix is in upper-triangular form.
	// Note that none of the original bits need to be modified during this
	// process.

	// Finally, working from right to left:

	// For recovery rows, eliminate any bits that are set unless that is the
	// pivot row.

	for (int ii 
	// TODO

	return true;
}

void print(const u8 *data, int bytes) {
	for (int ii = 0; ii < bytes; ++ii) {
		cout << (int)data[ii] << " ";
	}
	cout << endl;
}

int main() {
	m_clock.OnInitialize();

	cauchy_init();

	m_clock.usec();

	cout << "Cauchy matrix solver" << endl;

	int block_bytes = 8 * 162; // a multiple of 8
	int block_count = 32;
	int recovery_block_count = 1;

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

	ReceivedBlock *original_array = new ReceivedBlock[block_count];
	ReceivedBlock *recovery_array = new ReceivedBlock[recovery_block_count];

	for (int ii = 0; ii < block_count; ++ii) {
		original_array[ii].data = data + ii * block_bytes;
		original_array[ii].row = ii;
	}

	int erased = block_count - 1;

	recovery_array[0].data = recovery_blocks;
	recovery_array[0].row = 0;

	t0 = m_clock.usec();

	assert(cauchy_decode(block_count, recovery_block_count, original_array, recovery_array, 1, block_bytes));

	t1 = m_clock.usec();

	assert(!memcmp(original_array[erased].data, recovery_array[0].data, block_bytes));

	cout << "Cauchy decode in " << (t1 - t0) << " usec" << endl;

	m_clock.OnFinalize();

	return 0;
}

