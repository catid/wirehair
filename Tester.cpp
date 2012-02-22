#include "Wirehair.hpp"
#include "SmallPRNG.hpp"
#include "Clock.hpp"
#include "GF2Matrix.hpp"
using namespace cat;

#include <iostream>
#include <fstream>
using namespace std;	

static Clock m_clock;




#include "WirehairUtil.hpp"
void TestInc()
{
	for (u16 b = 4; b < 256; ++b)
	{
		u16 p = cat::wirehair::NextPrime16(b);
		//cout << "p = " << p << endl;
		u16 smallest = b;

		for (u16 x = 0; x < b; ++x)
		{
			for (u16 a = 1; a < b; ++a)
			{
				u16 y = x;
				u16 y1 = y;
				cat::wirehair::IterateNextColumn(y, b, p, a);
				u16 y2 = y;

				for (u16 am = 1; am < b - 1; ++am)
				{
					u16 ap = (am * a) % p;
					if (ap >= b) ap = (((u32)a << 16) + ap - p) % a;

					u16 y3 = y;
					cat::wirehair::IterateNextColumn(y3, b, p, ap);

					if (y1 == y3 || y1 == y2 || y2 == y3)
					{
						if (am < smallest)
							smallest = am;
						//cout << "FAIL for x = " << x << ", a = " << a << ", am = " << am << endl;
					}
				}
			}
		}

		cout << b << " -> " << smallest << " for " << b * (b - 1) * (smallest - 1) << endl;
	}
}


void Print64(u64 n, int w)
{
	for (int ii = 0; ii < w; ++ii)
	{
		if (n & ((u64)1 << ii))
			cout << '1';
		else
			cout << '0';
	}
	cout << endl;
}

void GenerateMatrix(int w, int seed)
{
	u64 matrix[256];

	CatsChoice prng;
	prng.Initialize(seed);

	// Shuffle row order
	u8 rows[256];
	rows[0] = 0;
	for (int ii = 1;;)
	{
		u32 jj, rv = prng.Next();

		// Unroll to use fewer calls to prng.Next()
		jj = (u8)rv % ii;
		rows[ii] = rows[jj];
		rows[jj] = ii;
		if (++ii >= w) break;

		jj = (u8)(rv >> 8) % ii;
		rows[ii] = rows[jj];
		rows[jj] = ii;
		if (++ii >= w) break;

		jj = (u8)(rv >> 16) % ii;
		rows[ii] = rows[jj];
		rows[jj] = ii;
		if (++ii >= w) break;

		jj = (u8)(rv >> 24) % ii;
		rows[ii] = rows[jj];
		rows[jj] = ii;
		if (++ii >= w) break;
	}

	// Shuffle bit order
	u8 bits[256];
	bits[0] = 0;
	for (int ii = 1;;)
	{
		u32 jj, rv = prng.Next();

		// Unroll to use fewer calls to prng.Next()
		jj = (u8)rv % ii;
		bits[ii] = bits[jj];
		bits[jj] = ii;
		if (++ii >= w) break;

		jj = (u8)(rv >> 8) % ii;
		bits[ii] = bits[jj];
		bits[jj] = ii;
		if (++ii >= w) break;

		jj = (u8)(rv >> 16) % ii;
		bits[ii] = bits[jj];
		bits[jj] = ii;
		if (++ii >= w) break;

		jj = (u8)(rv >> 24) % ii;
		bits[ii] = bits[jj];
		bits[jj] = ii;
		if (++ii >= w) break;
	}

	// Initialize counters
	u8 set_count = (w + 1) >> 1;
	u8 *set_bits = bits;
	u8 *clr_bits = set_bits + set_count;

	// Generate first row
	u64 x0 = 0;
	for (int ii = 0; ii < set_count; ++ii)
	{
		x0 ^= (u64)1 << set_bits[ii];
	}

	// Set up generator
	const u8 *row = rows;
	int loop_count = w >> 1;

	// Store first row
	matrix[*row++] = x0;

	// Generate first half of rows
	for (int ii = 0; ii < loop_count; ++ii)
	{
		x0 ^= (u64)1 << set_bits[ii];
		x0 ^= (u64)1 << clr_bits[ii];

		matrix[*row++] = x0;
	}

	// Handle odd window sizes
	if (w & 1)
	{
		x0 ^= (u64)1 << set_bits[loop_count];

		matrix[*row++] = x0;
	}

	// Generate second half of rows
	for (int ii = 0; ii < loop_count - 1; ++ii)
	{
		x0 ^= (u64)1 << set_bits[ii];
		x0 ^= (u64)1 << clr_bits[ii];

		matrix[*row++] = x0;
	}

	for (int ii = 0; ii < w; ++ii)
	{
		//Print64(matrix[ii], w);
	}

	for (int ii = 0; ii < w; ++ii)
	{
		u64 mask = (u64)1 << ii;
		int count = 0;

		for (int jj = 0; jj < w; ++jj)
		{
			if (matrix[jj] & mask)
			{
				++count;
			}
		}

		if (count != w/2)
		{
			cout << "FAIL" << endl;
		}
	}
}

void TestDense()
{
	int seed = 0;
	for (;;)
		GenerateMatrix(20, seed++);
}



int main()
{
	m_clock.OnInitialize();

	//TestInc();
	//TestDense();

	int block_count = 4096;
	int block_bytes = 1024 + 512 + 1;
	int message_bytes = block_bytes * block_count;
	u8 *message = new u8[message_bytes];

	for (int ii = 0; ii < message_bytes; ++ii)
	{
		message[ii] = ii;
	}

	wirehair::Encoder encoder;

	g_seed = 5;
	//for (;;) encoder.Initialize(message, message_bytes, block_bytes);

	g_seed = 0;

	for (;;)
	{
		g_seed++;

		double start = m_clock.usec();
		u32 clocks = m_clock.cycles();
		bool success = encoder.Initialize(message, message_bytes, block_bytes);
		clocks = m_clock.cycles() - clocks;
		double end = m_clock.usec();

		if (success)
		{
			cout << "main: encoder.Initialize in " << clocks << " clocks and " << end - start << " usec with seed " << encoder.GetSeed() << endl;

			u8 *block = new u8[message_bytes];

			bool success = true;
			for (int ii = 0; ii < block_count; ++ii)
			{
				encoder.Generate(ii, block);

				if (block[0] != (u8)ii)
				{
					success = false;
					cout << "Block " << ii << " doesn't match: " << (int)block[0] << endl;
				}
			}

			cin.get();
		}
		else
		{
			cout << "main: encoder.Initialize FAIL in " << clocks << " clocks with seed " << encoder.GetSeed() << endl;
		}
	}

	m_clock.OnFinalize();

	return 0;
}
