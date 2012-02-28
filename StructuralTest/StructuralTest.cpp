#include "GF2Matrix.hpp"
#include "../codec_source/SmallPRNG.hpp"
using namespace cat;

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

int main()
{
	wirehair::GF2Matrix m;

	int check_count = 8;

	m.Initialize(check_count);
	u64 worked = 0;

	for (int seed = 0; seed < 1000000; ++seed)
	{
		m.Zero();

		FillMatrixShuffleCodeRand(m, seed);

		//m.SetSeed(seed);m.Fill();

		m.Print();

		if (m.Triangle())
		{
			cout << "Invertible!" << endl;
			cin.get();
		}
		else
		{
			cout << "Not invertible!" << endl;
		}
	}

	return 0;
}
