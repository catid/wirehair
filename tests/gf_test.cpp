#include <iostream>
#include <cassert>
using namespace std;

#include "Platform.hpp"
#include "MemXOR.hpp"
#include "AbyssinianPRNG.hpp"
#include "Clock.hpp"
using namespace cat;

//// Utility: GF(2^16) Math functions

#define GF_BITS 16
#define GF_SIZE ((u32)1 << GF_BITS)
#define GF_ORDER (u16)(GF_SIZE - 1)

static u16 *GF_LOG = 0;
static u16 *GF_EXP = 0;
static const u32 GF_POLY = 0x1100B;

// Based on GF-Complete 1.02 by James S. Plank (University of Tennessee)
// This optimized version was written from scratch

static void gf_init() {
	// If already initialized,
	if (GF_LOG) {
		return;
	}

	// Allocate space for tables
	GF_LOG = new u16[3 * GF_SIZE];
	GF_EXP = GF_LOG + GF_SIZE;

	// Define log(0) and log(Ord(GF)) as 0
	GF_LOG[0] = 0;
	GF_LOG[GF_ORDER] = 0;

	// Fill tables
	u32 bits = 1;
	for (int ii = 0; ii < GF_ORDER; ++ii) {
		const u32 n = bits;

		// Insert table entry
		GF_LOG[n] = (u16)ii;
		GF_EXP[ii] = n;
		GF_EXP[ii + GF_ORDER] = n;

		// LFSR iteration
		bits <<= 1;
		bits ^= GF_POLY & -(s32)(bits >> 16);
	}
}

static u16 gf_mul_ref(u16 x, u16 y) {
	u32 p = 0;

	for (int i = 0; i < GF_BITS; i++) { 
		if (x & (1 << i)) {
			p ^= (y << i);
		}
	}

	for (int i = GF_BITS*2-2; i >= GF_BITS; i--) {
		if (p & (1 << i)) {
			p ^= GF_POLY << (i - GF_BITS);
		}
	}

	return p;
}

// x * y
static CAT_INLINE u16 gf_mul(u16 x, u16 y) {
	if (x == 0 || y == 0) {
		return 0;
	}
	return GF_EXP[GF_LOG[x] + GF_LOG[y]];
}

// x / y
static CAT_INLINE u16 gf_div(u16 x, u16 y) {
	if (x == 0 || y == 0) {
		return 0;
	}
	return GF_EXP[GF_LOG[x] + GF_ORDER - GF_LOG[y]];
}

static void gf_muladd_mem_ref(u16 * CAT_RESTRICT dest, u16 n,
					  const u16 * CAT_RESTRICT src, int words)
{
	// If degenerate case of multiplying by 0,
	if (n == 0) {
		return;
	}

	// Multiply remaining words
	while (words > 0) {
		--words;

		u16 a = *src++;

		*dest++ ^= gf_mul(a, n);
	}
}

// dest += src * n
static void gf_muladd_mem(u16 * CAT_RESTRICT dest, u16 n,
					  const u16 * CAT_RESTRICT src, int words)
{
	// If degenerate case of multiplying by 0,
	if (n == 0) {
		return;
	}

	// If degenerate case of multiplying by 1,
	if (n == 1) {
		memxor(dest, src, words*2);
		return;
	}

	// Construct multiplication table
	u16 log_n = GF_LOG[n];
	u16 T[4][16];
	for (int i = 0; i < 4; ++i) {
		T[i][0] = 0;
		for (int j = 1; j < 16; ++j) {
			T[i][j] = GF_EXP[GF_LOG[j << i*4] + log_n];
		}
	}

#ifdef CAT_ENDIAN_LITTLE
	// Multiply bulk of data
	while (words >= 4) {
		words -= 4;

		u64 x = *(const u64 *)src;
		src += 4;

		u16 pa = T[0][x & 15];
		u16 pb = T[0][(x >> 16) & 15];
		u16 pc = T[0][(x >> 32) & 15];
		u16 pd = T[0][(x >> 48) & 15];

		pa ^= T[1][(x >> 4) & 15];
		pb ^= T[1][(x >> 20) & 15];
		pc ^= T[1][(x >> 36) & 15];
		pd ^= T[1][(x >> 52) & 15];

		pa ^= T[2][(x >> 8) & 15];
		pb ^= T[2][(x >> 24) & 15];
		pc ^= T[2][(x >> 40) & 15];
		pd ^= T[2][(x >> 56) & 15];

		pa ^= T[3][(x >> 12) & 15];
		pb ^= T[3][(x >> 28) & 15];
		pc ^= T[3][(x >> 44) & 15];
		pd ^= T[3][x >> 60];

		u64 r = pa | ((u32)pb << 16) | ((u64)pc << 32) | ((u64)pd << 48);

		*(u64 *)dest ^= r;
		dest += 4;
	}
#endif // CAT_ENDIAN_LITTLE

	// Multiply remaining words
	while (words > 0) {
		--words;

		u16 a = *src++;

		u16 prod = T[0][a & 15];
		prod ^= T[1][(a >> 4) & 15];
		prod ^= T[2][(a >> 8) & 15];
		prod ^= T[3][a >> 12];

		*dest++ ^= prod;
	}
}

static void gf_div_mem_ref(u16 * CAT_RESTRICT data, u16 n,
					 		int words)
{
	// If degenerate case of multiplying by 0,
	if (n == 0) {
		return;
	}

	// Multiply remaining words
	while (words > 0) {
		--words;

		u16 a = *data;

		*data++ = gf_div(a, n);
	}
}

// dest /= n
static void gf_div_mem(u16 * CAT_RESTRICT data, u16 n, int words)
{
	// If degenerate cases of dividing by 0 or 1,
	if (n <= 1) {
		return;
	}

	// Construct multiplication table
	u16 log_n = GF_ORDER - GF_LOG[n];
	u16 T[4][16];
	for (int i = 0; i < 4; i++) {
		T[i][0] = 0;
		for (int j = 1; j < 16; j++) {
			T[i][j] = GF_EXP[GF_LOG[j << i*4] + log_n];
		}
	}

#ifdef CAT_ENDIAN_LITTLE
	// Multiply bulk of data
	while (words >= 4) {
		words -= 4;

		u64 x = *(const u64 *)data;

		u16 pa = T[0][x & 15];
		u16 pb = T[0][(x >> 16) & 15];
		u16 pc = T[0][(x >> 32) & 15];
		u16 pd = T[0][(x >> 48) & 15];

		pa ^= T[1][(x >> 4) & 15];
		pb ^= T[1][(x >> 20) & 15];
		pc ^= T[1][(x >> 36) & 15];
		pd ^= T[1][(x >> 52) & 15];

		pa ^= T[2][(x >> 8) & 15];
		pb ^= T[2][(x >> 24) & 15];
		pc ^= T[2][(x >> 40) & 15];
		pd ^= T[2][(x >> 56) & 15];

		pa ^= T[3][(x >> 12) & 15];
		pb ^= T[3][(x >> 28) & 15];
		pc ^= T[3][(x >> 44) & 15];
		pd ^= T[3][x >> 60];

		u64 r = pa | ((u32)pb << 16) | ((u64)pc << 32) | ((u64)pd << 48);

		*(u64 *)data = r;
		data += 4;
	}
#endif // CAT_ENDIAN_LITTLE

	// Multiply remaining words
	while (words > 0) {
		--words;

		u16 a = *data;

		u16 p = T[0][a & 15];
		p ^= T[1][(a >> 4) & 15];
		p ^= T[2][(a >> 8) & 15];
		p ^= T[3][a >> 12];

		*data++ = p;
	}
}

int main() {
	gf_init();

	Clock m_clock;
	m_clock.OnInitialize();
	Abyssinian prng;

	prng.Initialize(0);

	static const int N = 4096;
	u16 a[N], b[N], c[N], d[N];

	for (int kk = 0; kk < N; ++kk) {
		a[kk] = (u16)prng.Next();
	}

	double t0 = m_clock.usec();
	for (int ii = 0; ii < GF_SIZE; ++ii) {
		gf_muladd_mem(c, ii, a, N);
	}
	double t1 = m_clock.usec();

	double avg = (t1 - t0) / GF_SIZE;
	cout << N * 2 / (avg) << " MB/s gf_muladd_mem" << endl;

	u16 sum = 0;
	for (int ii = 0; ii < N; ++ii) {
		sum ^= c[ii];
	}
	cout << sum << endl;

	cout << "mul ref mem test" << endl;

	for (int ii = 0; ii < GF_SIZE; ++ii) {
		for (int kk = 0; kk < N; ++kk) {
			a[kk] = (u16)prng.Next();
		}
		memset(c, 0, sizeof(c));
		gf_muladd_mem(c, ii, a, N);
		memset(d, 0, sizeof(c));
		gf_muladd_mem_ref(d, ii, a, N);
		assert(!memcmp(c, d, sizeof(a)));
	}

	cout << "div ref mem test" << endl;

	for (int ii = 0; ii < GF_SIZE; ++ii) {
		for (int kk = 0; kk < N; ++kk) {
			a[kk] = (u16)prng.Next();
		}
		memcpy(c, a, sizeof(c));
		gf_div_mem(c, ii, N);
		memcpy(d, a, sizeof(c));
		gf_div_mem_ref(d, ii, N);
		assert(!memcmp(c, d, sizeof(a)));
	}

	cout << "mul-div mem test" << endl;

	for (int ii = 1; ii < GF_SIZE; ++ii) {
		for (int kk = 0; kk < N; ++kk) {
			a[kk] = (u16)prng.Next();
		}
		memset(c, 0, sizeof(c));
		gf_muladd_mem(c, ii, a, N);
		gf_div_mem(c, ii, N);
		assert(!memcmp(c, a, sizeof(a)));
	}

	cout << "div-mul mem test" << endl;

	for (int ii = 1; ii < GF_SIZE; ++ii) {
		for (int kk = 0; kk < N; ++kk) {
			a[kk] = (u16)prng.Next();
		}
		memcpy(d, a, sizeof(d));
		gf_div_mem(d, ii, N);
		memset(c, 0, sizeof(c));
		gf_muladd_mem(c, ii, d, N);
		assert(!memcmp(c, a, sizeof(a)));
	}

	cout << "exhaustive (i/j)*(i*j) == (i*i) test" << endl;

	for (int ii = 0; ii < GF_SIZE; ++ii) {
		for (int jj = 0; jj < GF_SIZE; ++jj) {
			// x = ii / jj
			u16 x = gf_div(ii, jj);
			// y = ii * jj
			u16 y = gf_mul(ii, jj);

			if (gf_mul_ref(ii, jj) != y) {
				cout << "FAIL " << ii << "*" << jj << endl;
				return 2;
			}

			// z = ii * ii
			u16 z = gf_mul(ii, ii);
			// w = x * y
			u16 w = gf_mul(x, y);

			if (jj != 0 && w != z) {
				cout << "FAILURE " << ii << "*" << jj << endl;
				return 1;
			}
		}
	}

	cout << "Utter success." << endl;

	m_clock.OnFinalize();

	return 0;
}

