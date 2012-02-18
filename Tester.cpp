#include "Wirehair.hpp"
#include "SmallPRNG.hpp"
#include "Clock.hpp"
#include "GF2Matrix.hpp"
using namespace cat;

#include <iostream>
#include <fstream>
using namespace std;	

static Clock m_clock;

static cat::wirehair::Encoder encoder;
static bool success;
static u8 *message;
static int message_bytes;
static int block_bytes;

static const u8 SQQ_TABLE[] = {
	0,  16,  22,  27,  32,  35,  39,  42,  45,  48,  50,  53,  55,  57,
	59,  61,  64,  65,  67,  69,  71,  73,  75,  76,  78,  80,  81,  83,
	84,  86,  87,  89,  90,  91,  93,  94,  96,  97,  98,  99, 101, 102,
	103, 104, 106, 107, 108, 109, 110, 112, 113, 114, 115, 116, 117, 118,
	119, 120, 121, 122, 123, 124, 125, 126, 128, 128, 129, 130, 131, 132,
	133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 144, 145,
	146, 147, 148, 149, 150, 150, 151, 152, 153, 154, 155, 155, 156, 157,
	158, 159, 160, 160, 161, 162, 163, 163, 164, 165, 166, 167, 167, 168,
	169, 170, 170, 171, 172, 173, 173, 174, 175, 176, 176, 177, 178, 178,
	179, 180, 181, 181, 182, 183, 183, 184, 185, 185, 186, 187, 187, 188,
	189, 189, 190, 191, 192, 192, 193, 193, 194, 195, 195, 196, 197, 197,
	198, 199, 199, 200, 201, 201, 202, 203, 203, 204, 204, 205, 206, 206,
	207, 208, 208, 209, 209, 210, 211, 211, 212, 212, 213, 214, 214, 215,
	215, 216, 217, 217, 218, 218, 219, 219, 220, 221, 221, 222, 222, 223,
	224, 224, 225, 225, 226, 226, 227, 227, 228, 229, 229, 230, 230, 231,
	231, 232, 232, 233, 234, 234, 235, 235, 236, 236, 237, 237, 238, 238,
	239, 240, 240, 241, 241, 242, 242, 243, 243, 244, 244, 245, 245, 246,
	246, 247, 247, 248, 248, 249, 249, 250, 250, 251, 251, 252, 252, 253,
	253, 254, 254, 255
};

static u16 cat_fred_sqrt16(u16 x)
{
	u16 r;

	if (x >= 0x100)
	{
		if (x >= 0x1000)
		{
			if (x >= 0x4000)
				r = SQQ_TABLE[x >> 8] + 1;
			else
				r = (SQQ_TABLE[x >> 6] >> 1) + 1;
		}
		else
		{
			if (x >= 0x400)
				r = (SQQ_TABLE[x >> 4] >> 2) + 1;
			else
				r = (SQQ_TABLE[x >> 2] >> 3) + 1;
		}
	}
	else
	{
		return SQQ_TABLE[x] >> 4;
	}

	// Correct rounding if necessary (compiler optimizes this form better)
	r -= (r * r > x);

	return r;
}


/*
	Hybrid Sieve of Eratosthenes Next Prime function

	It uses trial division up to the square root of the number to test.
	Uses a truncated sieve table to pick the next number to try, which
	avoids small factors 2, 3, 5, 7.  This can be considered a more
	involved version of incrementing by 2 instead of 1.

	Because of the tabular increment this is a hybrid approach.  The sieve
	would just use a very large table, but I wanted to limit the size of
	the table to something fairly small and reasonable.  210 bytes for the
	sieve table.  102 bytes for the primes list.

	It also calculates the integer square root faster than cmath sqrt()
	and uses multiplication to update the square root instead for speed.
*/
const int SIEVE_TABLE_SIZE = 2*3*5*7;
static const u8 SIEVE_TABLE[SIEVE_TABLE_SIZE] = {
	1, 0, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, 1, 0, 3, 2, 1, 0, 1, 0, 3, 2, 1, 0, 5, 4, 3, 2, 1, 0,
	1, 0, 5, 4, 3, 2, 1, 0, 3, 2, 1, 0, 1, 0, 3, 2, 1, 0, 5, 4, 3, 2, 1, 0, 5, 4, 3, 2, 1, 0,
	1, 0, 5, 4, 3, 2, 1, 0, 3, 2, 1, 0, 1, 0, 5, 4, 3, 2, 1, 0, 3, 2, 1, 0, 5, 4, 3, 2, 1, 0,
	7, 6, 5, 4, 3, 2, 1, 0, 3, 2, 1, 0, 1, 0, 3, 2, 1, 0, 1, 0, 3, 2, 1, 0, 7, 6, 5, 4, 3, 2,
	1, 0, 5, 4, 3, 2, 1, 0, 3, 2, 1, 0, 5, 4, 3, 2, 1, 0, 1, 0, 3, 2, 1, 0, 5, 4, 3, 2, 1, 0,
	1, 0, 5, 4, 3, 2, 1, 0, 5, 4, 3, 2, 1, 0, 3, 2, 1, 0, 1, 0, 3, 2, 1, 0, 5, 4, 3, 2, 1, 0,
	1, 0, 5, 4, 3, 2, 1, 0, 3, 2, 1, 0, 1, 0, 3, 2, 1, 0, 1, 0, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0,
};

static const u16 PRIMES_UNDER_256[] = {
	3, 5, 7, 
	11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61,
	67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127,
	131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191,
	193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 0x7fff
};

static u16 NextPrime16(u16 n)
{
	// Handle small n
	switch (n)
	{
	case 0:
	case 1:	return 1;
	case 2:	return 2;
	case 3:	return 3;
	case 4:
	case 5:	return 5;
	case 6:
	case 7:	return 7;
	}

	// Choose first n from table
	int offset = n % SIEVE_TABLE_SIZE;
	u32 next = SIEVE_TABLE[offset];
	offset += next + 1;
	n += next;

	// Initialize p_max to sqrt(n)
	int p_max = cat_fred_sqrt16(n);

	// For each number to try,
	for (;;)
	{
		// For each prime to test up to square root,
		const u16 *prime = PRIMES_UNDER_256 + 3;
		for (;;)
		{
			// If the next prime is above p_max we are done!
			int p = *prime;
			if (p > p_max)
				return n;

			// If composite, try next n
			if (n % p == 0)
				break;

			// Try next prime
			++prime;
		}

		// Use table to choose next trial number
		if (offset >= SIEVE_TABLE_SIZE) offset -= SIEVE_TABLE_SIZE;
		u32 next = SIEVE_TABLE[offset];
		offset += next + 1;
		n += next + 1;

		// Derivative square root iteration of p_max
		if (p_max * p_max < n)
			++p_max;
	}

	return n;
}

static u16 NextHighestPrime2(u16 n)
{
	if (n == 2) return 2;
	n |= 1;

	int p_max = cat_fred_sqrt16(n);

	for (;;)
	{
		const u16 *prime = PRIMES_UNDER_256;

		for (;;)
		{
			int p = *prime;

			if (p > p_max)
				return n;

			if (n % p == 0)
				break;

			++prime;
		}

		n += 2;

		if (p_max * p_max < n)
			++p_max;
	}

	return n;
}

void TestLoop()
{
	success = encoder.Initialize(message, message_bytes, 3 + block_bytes);
}

void TestSqrt()
{
	u32 b = 0;

	for (u32 x = 0; x < 0x10000; ++x)
	{
		u16 root1 = NextPrime16(x);
		u16 root2 = NextHighestPrime2(x);

		if (root1 != root2)
		{
			cout << "FAIL at " << x << " : " << root1 << " != " << root2 << endl;
		}

		b ^= root1;
		b += root2;
	}
}

static u32 b;

void TestSqrt1()
{
	for (u32 x = 0; x < 0x10000; ++x)
	{
		u16 root1 = NextPrime16(x);

		b ^= root1;
	}
}

void TestSqrt2()
{
	for (u32 x = 0; x < 0x10000; ++x)
	{
		u16 root2 = NextHighestPrime2(x);

		b += root2;
	}
}


static ofstream file;

double BenchmarkOverhead(int trials)
{
	u64 success_counter = 0;
	u64 trials_counter = 0;
	double clock_sum = 0;
	u32 clock_count = 0;

	for (;;)
	{
		++trials_counter;

		double start = m_clock.usec();
		TestLoop();
		double end = m_clock.usec();

		if (success)
		{
			++success_counter;
			//cout << "main: encoder.Initialize in " << clocks << " clocks and " << end - start << " usec with seed " << encoder.GetSeed() << endl;

			clock_sum += end - start;
			++clock_count;
		}

		if (trials_counter % trials == 0)
		{
			double success_rate = success_counter / (double)trials_counter;
			double speed = clock_sum / (double)clock_count;
			cout << "Check Blocks = " << g_check_block_count << ".  Success rate at 0 overhead: " << success_rate << " average " << speed << " usec" << endl;
			file << g_check_block_count << "\t" << success_rate << "\t" << speed << endl;
			clock_count = 0;
			clock_sum = 0;
			return success_rate;
		}
	}
}


static const u8 INVERTIBLE_MATRIX_SEEDS[512] = {
	0,1,2,1,5,2,3,0,0,0,7,3,3,2,0,1,1,1,2,0,0,1,1,2,5,5,4,7,4,4,6,0,
	0,0,0,1,2,4,4,4,7,0,0,3,8,1,1,0,0,0,0,0,0,6,9,2,4,4,5,4,4,8,5,1,
	3,1,0,3,1,2,0,7,0,2,5,5,0,4,1,4,4,0,0,2,0,0,8,2,4,4,1,1,2,0,4,7,
	1,1,1,2,10,12,7,3,5,3,3,0,1,1,1,4,4,3,3,2,3,3,3,0,0,0,0,4,2,11,1,0,
	1,5,2,7,3,1,1,3,1,0,1,0,0,0,1,1,1,0,2,0,0,0,0,0,3,3,2,2,1,1,0,2,
	2,1,0,0,7,4,2,2,1,0,0,0,0,0,0,0,0,2,0,0,2,6,1,1,1,3,2,1,0,0,6,0,
	0,5,3,5,0,5,3,2,1,4,1,1,1,0,3,4,4,2,0,7,7,7,2,0,0,8,1,1,6,2,0,0,
	0,3,6,3,3,3,2,7,1,1,0,2,2,1,6,0,0,0,0,2,0,0,4,2,5,0,2,3,0,2,2,5,
	3,0,1,1,1,1,1,1,1,2,0,2,2,1,1,4,9,5,2,4,4,3,0,0,4,5,2,2,5,7,1,7,
	0,4,1,2,1,3,1,10,0,0,0,1,1,0,0,0,0,2,3,1,16,3,2,0,1,2,3,0,2,3,1,7,
	7,3,7,2,2,1,1,1,1,1,0,1,2,2,0,0,4,2,6,4,8,1,2,2,7,4,1,1,2,8,0,0,
	7,4,5,4,4,1,0,0,2,1,2,2,2,1,1,1,0,1,3,1,4,4,3,5,2,2,3,3,9,0,6,11,
	0,1,4,1,1,2,2,1,1,0,2,2,6,8,5,1,1,1,4,2,2,1,6,9,0,0,7,6,6,2,2,2,
	3,4,1,3,3,3,5,7,0,4,4,1,1,1,2,3,3,3,6,2,0,0,0,3,2,6,1,2,1,1,1,4,
	1,6,1,3,4,4,4,0,0,4,8,1,4,3,3,3,1,0,0,0,1,1,1,2,2,2,2,1,6,7,2,0,
	2,3,0,0,0,0,3,1,4,0,0,5,0,0,0,0,1,12,7,0,0,0,0,0,0,0,0,4,0,0,0,0,
};

void TestGF2()
{
	wirehair::GF2Matrix gf2m;

	u8 seeds[512];

	for (int n = 1; n < 512; ++n)
	{
		if (!gf2m.Initialize(n))
		{
			cout << "FAIL INIT" << endl;
			cin.get();
		}

		u64 success_count = 0;
		u64 trial_count = 0;

		cout << "Testing with n = " << n << endl;

		for (;;)
		{
			++trial_count;
			gf2m.Fill();
			//gf2m.Print();
			if (!gf2m.Triangle())
			{
				//cout << "--> Fail... next!" << endl;
				//continue;
				cout << "Not invertible with seed " << gf2m.GetSeed() << endl;
				gf2m.NextSeed();
			}
			else
			{
				//gf2m.Print();
				++success_count;
				u32 seed = gf2m.GetSeed();
				cout << "--> Invertible with seed " << gf2m.GetSeed() << endl;
				seeds[n] = seed;
				if (seed > 255)
				{
					cout << "FAIL" << endl;
					cin.get();
				}
				break;
			}
			u32 seed = gf2m.GetSeed();
			//gf2m.Print();
			//gf2m.Diagonal();
			//gf2m.Print();
			//cout << "--> Success!" << endl;
			//cin.get();

			if (trial_count % 10000 == 0)
			{
				double success_rate = success_count / (double)trial_count;
				cout << "Success rate = " << success_rate << endl;
			}
		}
	}

	cout << "static const u8 INVERTIBLE_MATRIX_SEEDS[512] = {";
	for (int ii = 0; ii < 512; ++ii)
	{
		if (ii % 32 == 0) cout << endl;
		cout << (int)seeds[ii] << ",";
	}
	cout << endl << "};" << endl;
}

int main()
{
	m_clock.OnInitialize();

	//TestGF2();

/*
	TestSqrt();
	for (;;)
	{
		u32 clocks1 = m_clock.MeasureClocks(100, &TestSqrt1);
		u32 clocks2 = m_clock.MeasureClocks(100, &TestSqrt2);

		cout << "PRIME1 = " << clocks1 << " clocks" << endl;
		cout << "PRIME2 = " << clocks2 << " clocks" << endl;
	}
	cout << b << endl;
*/
	int block_count = 2048;
	block_bytes = 1024 + 512 + 1;
	message_bytes = block_bytes * block_count;
	message = new u8[message_bytes];

	for (int ii = 0; ii < message_bytes; ++ii)
	{
		message[ii] = ii;
	}

	file.open("tuning.txt");
	if (!file) cout << "FAIL";
	for (g_check_block_count = 6; g_check_block_count < 64; ++g_check_block_count)
	{
		for (int ii = 0; ii < 10; ++ii)
		{
			BenchmarkOverhead(100);
		}
	}

	file.close();
	cout << "DONE!" << endl;
	cin.get();

	for (;;)
	{
		//u32 clocks = m_clock.MeasureClocks(1, &TestLoop);
		double start = m_clock.usec();
		u32 clocks = m_clock.cycles();
		TestLoop();
		clocks = m_clock.cycles() - clocks;
		double end = m_clock.usec();

		if (success)
		{
			cout << "main: encoder.Initialize in " << clocks << " clocks and " << end - start << " usec with seed " << encoder.GetSeed() << endl;

			u8 *block = new u8[3 + message_bytes];

			bool success = true;
			for (int ii = 0; ii < block_count; ++ii)
			{
				encoder.Generate(block);

				if (block[3] != (u8)ii)
				{
					success = false;
					cout << "Block " << ii << " doesn't match: " << (int)block[3] << endl;
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
