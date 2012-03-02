#include "GF2Matrix.hpp"
#include "GF256Matrix.hpp"
#include "../codec_source/SmallPRNG.hpp"
#include "../Clock.hpp"
using namespace cat;

static Clock m_clock;

#include <iostream>
#include <fstream>
using namespace std;	

#if defined(CAT_DEBUG)
#define CAT_IF_DUMP(x) x
#else
#define CAT_IF_DUMP(x)
#endif

//#define CAT_SHUFFLE_HALF

void ShuffleDeck16(CatsChoice &prng, u16 *deck, u32 count)
{
	deck[0] = 0;

	// If we can unroll 4 times,
	if (count <= 256)
	{
		for (u32 ii = 1;;)
		{
			u32 jj, rv = prng.Next();

			// 8-bit unroll
			switch (count - ii)
			{
			default:
				jj = (u8)rv % ii;
				deck[ii] = deck[jj];
				deck[jj] = ii;
				++ii;
				jj = (u8)(rv >> 8) % ii;
				deck[ii] = deck[jj];
				deck[jj] = ii;
				++ii;
				jj = (u8)(rv >> 16) % ii;
				deck[ii] = deck[jj];
				deck[jj] = ii;
				++ii;
				jj = (u8)(rv >> 24) % ii;
				deck[ii] = deck[jj];
				deck[jj] = ii;
				++ii;
				break;

			case 3:
				jj = (u8)rv % ii;
				deck[ii] = deck[jj];
				deck[jj] = ii;
				++ii;
			case 2:
				jj = (u8)(rv >> 8) % ii;
				deck[ii] = deck[jj];
				deck[jj] = ii;
				++ii;
			case 1:
				jj = (u8)(rv >> 16) % ii;
				deck[ii] = deck[jj];
				deck[jj] = ii;
			case 0:
				return;
			}
		}
	}
	else
	{
		// For each deck entry,
		for (u32 ii = 1;;)
		{
			u32 jj, rv = prng.Next();

			// 16-bit unroll
			switch (count - ii)
			{
			default:
				jj = (u16)rv % ii;
				deck[ii] = deck[jj];
				deck[jj] = ii;
				++ii;
				jj = (u16)(rv >> 16) % ii;
				deck[ii] = deck[jj];
				deck[jj] = ii;
				++ii;
				break;

			case 1:
				jj = (u16)rv % ii;
				deck[ii] = deck[jj];
				deck[jj] = ii;
			case 0:
				return;
			}
		}
	}
}

void FillMatrixColumnShuffleCode(wirehair::GF2Matrix &m, u32 seed)
{
	int check_count = m.Size();
	int pitch = m.GetPitch();
	u64 *matrix = m.GetFront();
	u32 temp_col1[100];

	CatsChoice prng;
	prng.Initialize(seed);

	u16 cols1[1000], bits1[1000];

	// Shuffle col and bit order
	ShuffleDeck16(prng, cols1, check_count);
	ShuffleDeck16(prng, bits1, check_count);

	// Initialize counters
	u16 set_count = (check_count + 1) >> 1;
	u16 *set_bits1 = bits1;
	u16 *clr_bits1 = set_bits1 + set_count;

	// Generate first col
	CAT_OBJCLR(temp_col1);
	for (int ii = 0; ii < set_count; ++ii)
	{
		// If bit is peeled,
		int bit_i1 = set_bits1[ii];
		temp_col1[bit_i1] ^= 1;
	}

	// Set up generator
	const u16 *col1 = cols1;

	// Store first row
	CAT_IF_DUMP(for (int ii = 0; ii < check_count; ++ii) cout << temp_col1[ii]; cout << " <- going to column " << *col1 << endl;)
	u16 column_i = *col1++;
	u64 *ge_dest_row = matrix + (column_i >> 6);
	for (int jj = 0; jj < check_count; ++jj) ge_dest_row[jj] ^= temp_col1[jj] << (column_i & 63);

	// Generate first half of rows
	const int loop_count = (check_count >> 1);
	for (int ii = 0; ii < loop_count; ++ii)
	{
		int bit01 = set_bits1[ii], bit11 = clr_bits1[ii];

		// Add in peeled columns
		temp_col1[bit01] ^= 1;
		temp_col1[bit11] ^= 1;

		// Store in row
		CAT_IF_DUMP(for (int ii = 0; ii < check_count; ++ii) cout << temp_col1[ii]; cout << " <- going to column " << *col1 << endl;)
		u16 column_i = *col1++;
		u64 *ge_dest_row = matrix + (column_i >> 6);
		for (int jj = 0; jj < check_count; ++jj) ge_dest_row[jj] ^= temp_col1[jj] << (column_i & 63);
	}

	// If check count is odd,
	if (check_count & 1)
	{
		// Generate middle row
		int bit01 = set_bits1[loop_count];
		temp_col1[bit01] ^= 1;

		// Store in row
		CAT_IF_DUMP(for (int ii = 0; ii < check_count; ++ii) cout << temp_col1[ii]; cout << " <- going to column " << *col1 << endl;)
		u16 column_i = *col1++;
		u64 *ge_dest_row = matrix + (column_i >> 6);
		for (int jj = 0; jj < check_count; ++jj) ge_dest_row[jj] ^= temp_col1[jj] << (column_i & 63);
	}

	const int second_loop_count = loop_count - 1;

	// Generate second half of rows
	for (int ii = 0; ii < second_loop_count; ++ii)
	{
		int bit01 = set_bits1[ii], bit11 = clr_bits1[ii];
		temp_col1[bit01] ^= 1;
		temp_col1[bit11] ^= 1;

		// Store in row
		CAT_IF_DUMP(for (int ii = 0; ii < check_count; ++ii) cout << temp_col1[ii]; cout << " <- going to column " << *col1 << endl;)
		u16 column_i = *col1++;
		u64 *ge_dest_row = matrix + (column_i >> 6);
		for (int jj = 0; jj < check_count; ++jj) ge_dest_row[jj] ^= temp_col1[jj] << (column_i & 63);
	}
}

void FillMatrixShuffleCode(wirehair::GF2Matrix &m, u32 seed)
{
	int check_count = m.Size();
	int pitch = m.GetPitch();
	u64 *matrix = m.GetFront();
	u64 temp_row1[100];

	CatsChoice prng;
	prng.Initialize(seed);

	u16 rows1[1000], bits1[1000];

	// Shuffle row and bit order
	ShuffleDeck16(prng, rows1, check_count);
	ShuffleDeck16(prng, bits1, check_count);

	// Initialize counters
	u16 set_count = (check_count + 1) >> 1;
	u16 *set_bits1 = bits1;
	u16 *clr_bits1 = set_bits1 + set_count;

	// Generate first row
	memset(temp_row1, 0, pitch * sizeof(u64));
	for (int ii = 0; ii < set_count; ++ii)
	{
		// If bit is peeled,
		int bit_i1 = set_bits1[ii];
		temp_row1[bit_i1 >> 6] ^= (u64)1 << (bit_i1 & 63);
	}

	// Set up generator
	const u16 *row1 = rows1;

	// Store first row
	CAT_IF_DUMP(for (int ii = 0; ii < check_count; ++ii) cout << ((temp_row1[ii >> 6] & ((u64)1 << (ii & 63))) ? '1' : '0'); cout << " <- going to row " << *row1 << endl;)
	u64 *ge_dest_row = matrix + pitch * *row1++;
	for (int jj = 0; jj < pitch; ++jj) ge_dest_row[jj] ^= temp_row1[jj];

	// Generate first half of rows
	const int loop_count = (check_count >> 1);
	for (int ii = 0; ii < loop_count; ++ii)
	{
		int bit01 = set_bits1[ii], bit11 = clr_bits1[ii];

		// Add in peeled columns
		temp_row1[bit01 >> 6] ^= (u64)1 << (bit01 & 63);
		temp_row1[bit11 >> 6] ^= (u64)1 << (bit11 & 63);

		// Store in row
		CAT_IF_DUMP(for (int ii = 0; ii < check_count; ++ii) cout << ((temp_row1[ii >> 6] & ((u64)1 << (ii & 63))) ? '1' : '0'); cout << " <- going to row " << *row1 << endl;)
		u64 *ge_dest_row = matrix + pitch * *row1++;
		for (int jj = 0; jj < pitch; ++jj) ge_dest_row[jj] ^= temp_row1[jj];
	}

	// If check count is odd,
	if (check_count & 1)
	{
		// Generate middle row
		int bit01 = set_bits1[loop_count];
		temp_row1[bit01 >> 6] ^= (u64)1 << (bit01 & 63);

		// Store in row
		CAT_IF_DUMP(for (int ii = 0; ii < check_count; ++ii) cout << ((temp_row1[ii >> 6] & ((u64)1 << (ii & 63))) ? '1' : '0'); cout << " <- going to row " << *row1 << endl;)
		u64 *ge_dest_row = matrix + pitch * *row1++;
		for (int jj = 0; jj < pitch; ++jj) ge_dest_row[jj] ^= temp_row1[jj];
	}

	const int second_loop_count = loop_count - 1;

	// Generate second half of rows
	for (int ii = 0; ii < second_loop_count; ++ii)
	{
		int bit01 = set_bits1[ii], bit11 = clr_bits1[ii];
		temp_row1[bit01 >> 6] ^= (u64)1 << (bit01 & 63);
		temp_row1[bit11 >> 6] ^= (u64)1 << (bit11 & 63);

		// Store in row
		CAT_IF_DUMP(for (int ii = 0; ii < check_count; ++ii) cout << ((temp_row1[ii >> 6] & ((u64)1 << (ii & 63))) ? '1' : '0'); cout << " <- going to row " << *row1 << endl;)
		u64 *ge_dest_row = matrix + pitch * *row1++;
		for (int jj = 0; jj < pitch; ++jj) ge_dest_row[jj] ^= temp_row1[jj];
	}
}

void FillMatrixShuffleCodeRand(wirehair::GF2Matrix &m, u32 seed)
{
	int check_count = m.Size();
	int pitch = m.GetPitch();
	u64 *matrix = m.GetFront();
	u64 temp_row1[100];

	CatsChoice prng;
	prng.Initialize(seed);

	u16 rows1[1000], bits1[1000];

	// Shuffle row and bit order
	ShuffleDeck16(prng, rows1, check_count);
	ShuffleDeck16(prng, bits1, check_count);

	// Initialize counters
	u16 set_count = (check_count + 1) >> 1;
	u16 *set_bits1 = bits1;
	u16 *clr_bits1 = set_bits1 + set_count;

	// Generate first row
	memset(temp_row1, 0, pitch * sizeof(u64));
	for (int ii = 0; ii < set_count; ++ii)
	{
		// If bit is peeled,
		int bit_i1 = set_bits1[ii];
		temp_row1[bit_i1 >> 6] ^= (u64)1 << (bit_i1 & 63);
	}

	// Set up generator
	const u16 *row1 = rows1;

	// Store first row
	CAT_IF_DUMP(for (int ii = 0; ii < check_count; ++ii) cout << ((temp_row1[ii >> 6] & ((u64)1 << (ii & 63))) ? '1' : '0'); cout << " <- going to row " << *row1 << endl;)
	u64 *ge_dest_row = matrix + pitch * *row1++;
	for (int jj = 0; jj < pitch; ++jj) ge_dest_row[jj] ^= temp_row1[jj];

	ShuffleDeck16(prng, bits1, check_count);

	// Generate first half of rows
	for (int ii = 1; ii < check_count; ++ii)
	{
		int bit01 = set_bits1[ii], bit11 = clr_bits1[ii];

		// Add in peeled columns
		temp_row1[bit01 >> 6] ^= (u64)1 << (bit01 & 63);
		temp_row1[bit11 >> 6] ^= (u64)1 << (bit11 & 63);

		// Store in row
		CAT_IF_DUMP(for (int ii = 0; ii < check_count; ++ii) cout << ((temp_row1[ii >> 6] & ((u64)1 << (ii & 63))) ? '1' : '0'); cout << " <- going to row " << *row1 << endl;)
		u64 *ge_dest_row = matrix + pitch * *row1++;
		for (int jj = 0; jj < pitch; ++jj) ge_dest_row[jj] ^= temp_row1[jj];
	}
}

void FindGF256GeneratorPolynomials()
{
	cout << "static const u8 GEN_POLY[] = {" << endl;
	int seen = 0;
	for (u32 taps = 0; taps < 256; ++taps)
	{
		u32 lfsr = 1;

		int count = 0;
		for (int ii = 0; ii < 255; ++ii)
		{
			u32 lsb = lfsr & 1;
			lfsr >>= 1;
			if (lsb) lfsr ^= taps;
			if (lfsr == 1) ++count;
		}

		if (lfsr == 1 && count == 1)
		{
			cout << "0x" << hex << (int)taps << dec << ", ";
			if ((++seen & 7) == 0) cout << endl;
		}
	}
	cout << "};" << endl;
}

static const u8 GEN_POLY[16] = {
	0x8e, 0x95, 0x96, 0xa6, 0xaf, 0xb1, 0xb2, 0xb4,
	0xb8, 0xc3, 0xc6, 0xd4, 0xe1, 0xe7, 0xf3, 0xfa,
};

void GenerateExpLogTables()
{
	for (int ii = 0; ii < 16; ++ii)
	{
		u32 poly = (GEN_POLY[ii] << 1) | 1;
		u16 LogTable[256], ALogTable[512*2+1];

		LogTable[0] = 512;
		ALogTable[0] = 1;
		for (u32 jj = 1; jj < 255; ++jj)
		{
			u32 next = (u32)ALogTable[jj - 1] * 2;
			if (next >= 256) next ^= poly;

			ALogTable[jj] = next;
			LogTable[ALogTable[jj]] = jj;
		}

		ALogTable[255] = ALogTable[0];
		LogTable[ALogTable[255]] = 255;

		for (u32 jj = 256; jj < 2 * 255; ++jj)
		{
			ALogTable[jj] = ALogTable[jj % 255];
		}

		ALogTable[2 * 255] = 1;

		for (u32 jj = 2 * 255 + 1; jj < 4 * 255; ++jj)
		{
			ALogTable[jj] = 0;
		}

		cout << "For generator polynomial " << hex << poly << dec << ":" << endl << endl;

		cout << "static const u16 LOG_TABLE[256] = {";
		for (int jj = 0; jj < 256; ++jj)
		{
			if ((jj & 15) == 0) cout << endl;
			cout << (int)LogTable[jj] << ", ";
		}
		cout << endl << "};" << endl << endl;

		cout << "static const u8 ALOG_TABLE[512*2+1] = {";
		for (int jj = 0; jj < 255*2+2; ++jj)
		{
			if ((jj & 15) == 0) cout << endl;
			cout << (int)ALogTable[jj] << ", ";
		}
		cout << endl << "};" << endl << endl;
	}
}

/*
	Branchless multiply and divide construction from
	"Fast Software Implementations of Finite Field Operations (Extended Abstract)"
	by Cheng Huang, Lihao Xu

	Small corrections made to paper (Q = 255):
		+ The ALOG_TABLE needs to have 512*2+1 elements to handle 0*0 = 0 case.
		+ Element 255*2 should be set to 1.

	After these corrections it works properly and reduces the execution time
	to 58% of the usual version that uses branches to handle zero input.

	This table was generated using polynomial 0x11D.
*/
static const u16 LOG_TABLE[256] = {
	512, 255, 1, 25, 2, 50, 26, 198, 3, 223, 51, 238, 27, 104, 199, 75,
	4, 100, 224, 14, 52, 141, 239, 129, 28, 193, 105, 248, 200, 8, 76, 113,
	5, 138, 101, 47, 225, 36, 15, 33, 53, 147, 142, 218, 240, 18, 130, 69,
	29, 181, 194, 125, 106, 39, 249, 185, 201, 154, 9, 120, 77, 228, 114, 166,
	6, 191, 139, 98, 102, 221, 48, 253, 226, 152, 37, 179, 16, 145, 34, 136,
	54, 208, 148, 206, 143, 150, 219, 189, 241, 210, 19, 92, 131, 56, 70, 64,
	30, 66, 182, 163, 195, 72, 126, 110, 107, 58, 40, 84, 250, 133, 186, 61,
	202, 94, 155, 159, 10, 21, 121, 43, 78, 212, 229, 172, 115, 243, 167, 87,
	7, 112, 192, 247, 140, 128, 99, 13, 103, 74, 222, 237, 49, 197, 254, 24,
	227, 165, 153, 119, 38, 184, 180, 124, 17, 68, 146, 217, 35, 32, 137, 46,
	55, 63, 209, 91, 149, 188, 207, 205, 144, 135, 151, 178, 220, 252, 190, 97,
	242, 86, 211, 171, 20, 42, 93, 158, 132, 60, 57, 83, 71, 109, 65, 162,
	31, 45, 67, 216, 183, 123, 164, 118, 196, 23, 73, 236, 127, 12, 111, 246,
	108, 161, 59, 82, 41, 157, 85, 170, 251, 96, 134, 177, 187, 204, 62, 90,
	203, 89, 95, 176, 156, 169, 160, 81, 11, 245, 22, 235, 122, 117, 44, 215,
	79, 174, 213, 233, 230, 231, 173, 232, 116, 214, 244, 234, 168, 80, 88, 175,
};

static const u8 ALOG_TABLE[512*2+1] = {
	1, 2, 4, 8, 16, 32, 64, 128, 29, 58, 116, 232, 205, 135, 19, 38,
	76, 152, 45, 90, 180, 117, 234, 201, 143, 3, 6, 12, 24, 48, 96, 192,
	157, 39, 78, 156, 37, 74, 148, 53, 106, 212, 181, 119, 238, 193, 159, 35,
	70, 140, 5, 10, 20, 40, 80, 160, 93, 186, 105, 210, 185, 111, 222, 161,
	95, 190, 97, 194, 153, 47, 94, 188, 101, 202, 137, 15, 30, 60, 120, 240,
	253, 231, 211, 187, 107, 214, 177, 127, 254, 225, 223, 163, 91, 182, 113, 226,
	217, 175, 67, 134, 17, 34, 68, 136, 13, 26, 52, 104, 208, 189, 103, 206,
	129, 31, 62, 124, 248, 237, 199, 147, 59, 118, 236, 197, 151, 51, 102, 204,
	133, 23, 46, 92, 184, 109, 218, 169, 79, 158, 33, 66, 132, 21, 42, 84,
	168, 77, 154, 41, 82, 164, 85, 170, 73, 146, 57, 114, 228, 213, 183, 115,
	230, 209, 191, 99, 198, 145, 63, 126, 252, 229, 215, 179, 123, 246, 241, 255,
	227, 219, 171, 75, 150, 49, 98, 196, 149, 55, 110, 220, 165, 87, 174, 65,
	130, 25, 50, 100, 200, 141, 7, 14, 28, 56, 112, 224, 221, 167, 83, 166,
	81, 162, 89, 178, 121, 242, 249, 239, 195, 155, 43, 86, 172, 69, 138, 9,
	18, 36, 72, 144, 61, 122, 244, 245, 247, 243, 251, 235, 203, 139, 11, 22,
	44, 88, 176, 125, 250, 233, 207, 131, 27, 54, 108, 216, 173, 71, 142, 1,
	2, 4, 8, 16, 32, 64, 128, 29, 58, 116, 232, 205, 135, 19, 38, 76,
	152, 45, 90, 180, 117, 234, 201, 143, 3, 6, 12, 24, 48, 96, 192, 157,
	39, 78, 156, 37, 74, 148, 53, 106, 212, 181, 119, 238, 193, 159, 35, 70,
	140, 5, 10, 20, 40, 80, 160, 93, 186, 105, 210, 185, 111, 222, 161, 95,
	190, 97, 194, 153, 47, 94, 188, 101, 202, 137, 15, 30, 60, 120, 240, 253,
	231, 211, 187, 107, 214, 177, 127, 254, 225, 223, 163, 91, 182, 113, 226, 217,
	175, 67, 134, 17, 34, 68, 136, 13, 26, 52, 104, 208, 189, 103, 206, 129,
	31, 62, 124, 248, 237, 199, 147, 59, 118, 236, 197, 151, 51, 102, 204, 133,
	23, 46, 92, 184, 109, 218, 169, 79, 158, 33, 66, 132, 21, 42, 84, 168,
	77, 154, 41, 82, 164, 85, 170, 73, 146, 57, 114, 228, 213, 183, 115, 230,
	209, 191, 99, 198, 145, 63, 126, 252, 229, 215, 179, 123, 246, 241, 255, 227,
	219, 171, 75, 150, 49, 98, 196, 149, 55, 110, 220, 165, 87, 174, 65, 130,
	25, 50, 100, 200, 141, 7, 14, 28, 56, 112, 224, 221, 167, 83, 166, 81,
	162, 89, 178, 121, 242, 249, 239, 195, 155, 43, 86, 172, 69, 138, 9, 18,
	36, 72, 144, 61, 122, 244, 245, 247, 243, 251, 235, 203, 139, 11, 22, 44,
	88, 176, 125, 250, 233, 207, 131, 27, 54, 108, 216, 173, 71, 142, 1, 0,
};

CAT_INLINE u8 Multiply(u8 a, u8 b)
{
	return ALOG_TABLE[LOG_TABLE[a] + LOG_TABLE[b]];
}

CAT_INLINE u8 Divide(u8 a, u8 b)
{
	// Precondition: b != 0
	return ALOG_TABLE[LOG_TABLE[a] + 255 - LOG_TABLE[b]];
}

u32 rf = 0;

void TestMultDiv()
{
	for (u32 a = 0; a < 256; ++a)
	{
		for (u32 b = 0; b < 256; ++b)
		{
			u32 r = Multiply(a, b);

			if (a == 0 && r != 0)
			{
				rf ^= r;
				cout << "FAIL for " << a << " * " << b << " = " << r << endl;
			}

			if (b == 0 && r != 0)
			{
				rf ^= r;
				cout << "FAIL for " << a << " * " << b << " = " << r << endl;
			}

			if (a != 0 && Divide(r, a) != b)
			{
				rf ^= r;
				cout << "FAIL for " << a << " * " << b << " = " << r << endl;
			}

			if (b != 0 && Divide(r, b) != a)
			{
				rf ^= r;
				cout << "FAIL for " << a << " * " << b << " = " << r << endl;
			}

			if (b != 0)
			{
				u32 d = Divide(a, b);
				rf ^= d;

				r = Multiply(d, b);
				if (r != a)
				{
					cout << "FAIL for " << a << " * " << b << " = " << r << endl;
					rf ^= r;
				}
			}
		}
	}
}

void TestInvertibleRate()
{
	wirehair::GF256Matrix m;

	int check_count = 8;

	m.Initialize(check_count);
	u64 worked = 0;

	for (int seed = 0; seed < 1000000; ++seed)
	{
		m.Zero();

		//FillMatrixShuffleCodeRand(m, seed);

		m.SetSeed(seed);m.Fill();

		//m.Print();

		cout << worked / (double)seed << endl;

		if (m.Triangle())
		{
			++worked;
			//cout << "Invertible!" << endl;
			//cin.get();
		}
		else
		{
			//cout << "Not invertible!" << endl;
		}
	}
}

int main()
{
	m_clock.OnInitialize();

	TestInvertibleRate();
	//FindGF256GeneratorPolynomials();
	GenerateExpLogTables();
	TestMultDiv();

	wirehair::GF256Matrix m;
	m.Initialize(16);
	m.Fill();
	m.Print();
	m.Triangle();
	m.Print();

	cin.get();

	m_clock.OnFinalize();
	return 0;
}
