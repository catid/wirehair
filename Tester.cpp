#include "Wirehair.hpp"
#include "SmallPRNG.hpp"
#include "Clock.hpp"
using namespace cat;

#include <iostream>
using namespace std;	

static Clock m_clock;

static const u16 PRIMES[] = {
	3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53,
	59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109,
	113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
	179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233,
	239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293,
	307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367,
	373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433,
	439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499,
	503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577,
	587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643,
	647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719,
	727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797,
	809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863,
	877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947,
	953, 967, 971, 977, 983, 991, 997, 1009, 1013, 1019,
	1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069,
	1087, 1091, 1093, 1097,
};

static int NextHighestPrime(int n)
{
	n |= 1;

	int p_max = (int)sqrt((float)n);

	for (;;)
	{
		const u16 *prime = PRIMES;

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


CAT_INLINE u32 BSR16(u16 x)
{
#if defined(CAT_COMPILER_MSVC) && !defined(CAT_DEBUG)

	u32 index;
	_BitScanReverse((unsigned long*)&index, x);
	return index;

#elif defined(CAT_ASM_INTEL) && defined(CAT_ISA_X86)

	u32 y = x;
	CAT_ASM_BEGIN
		BSR eax, [y]
	CAT_ASM_END

#elif defined(CAT_ASM_ATT) && defined(CAT_ISA_X86)

	register u32 y = x;
	u32 retval;

	CAT_ASM_BEGIN
		"BSRl %1, %%eax"
		: "=a" (retval)
		: "r" (y)
		: "cc"
	CAT_ASM_END

	return retval;

#else

#define CAT_NO_INTRINSIC_BSR32

	// Adapted from the Stanford Bit Twiddling Hacks collection
	u32 shift, r;

	r = (x > 0xFF) << 3; x >>= r;
	shift = (x > 0xF) << 2; x >>= shift; r |= shift;
	shift = (x > 0x3) << 1; x >>= shift; r |= shift;
	r |= (x >> 1);
	return r;

#endif
}

int IntLog3(int x)
{
	// Adapted from the Stanford Bit Twiddling Hacks collection
	u32 shift, r;

	r = (x > 0xFF) << 3; x >>= r;
	shift = (x > 0xF) << 2; x >>= shift; r |= shift;
	shift = (x > 0x3) << 1; x >>= shift; r |= shift;
	r |= (x >> 1);
	return r;
}

int IntLog2(int p)
{
	int x = 0;

	while (p)
	{
		p >>= 1;
		++x;
	}

	return x;
}

void TestLog2()
{
	u32 x = 0;

	double t0 = m_clock.usec();
	for (int y = 0; y < 1000; ++y)
	{
		for (u32 x = 1; x < 65536; ++x)
		{
			int test1 = IntLog2(x) - 1;
			x ^= test1;
		}
	}
	double t1 = m_clock.usec();
	for (int y = 0; y < 1000; ++y)
	{
		for (u32 x = 1; x < 65536; ++x)
		{
			int test2 = BSR16(x);
			x ^= test2;
		}
	}
	double t2 = m_clock.usec();
	for (int y = 0; y < 1000; ++y)
	{
		for (u32 x = 1; x < 65536; ++x)
		{
			int test3 = IntLog3(x);
			x ^= test3;
		}
	}
	double t3 = m_clock.usec();

	cout << (t1 - t0)/1000. << endl;
	cout << (t2 - t1)/1000. << endl;
	cout << (t3 - t2)/1000. << endl;
	cout << x << endl;
}

void TestOptimization()
{
	for (int b = 2; b < 100; ++b)
	{
		int p = NextHighestPrime(b);

		for (int a = 1; a < b; ++a)
		{
			for (int x = 0; x < b; ++x)
			{
				int y = (x + a) % p;
				int y0 = y;

				int z = y;
				while (z >= b)
					z = (z + a) % p;

				if (y >= b)
					y = ((a << 16) - p + y) % a;

				if (z != y)
				{
					cout << "FAIL with b=" << b << " p=" << p << " a=" << a << " x=" << x << " from y0=" << y0 << " -> " << y << " (wrong) vs " << z << " (correct)" << endl;
				}
			}
		}
	}

	u32 xx = 0;

	double start = m_clock.usec();

	for (int b = 2; b < 1024; ++b)
	{
		int p = NextHighestPrime(b);

		for (int a = 1; a < b; ++a)
		{
			for (int x = 0; x < b; ++x)
			{
				int y = (x + a) % p;
				while (y >= b) y = (y + a) % p;
				xx ^= y;
			}
		}
	}

	double end = m_clock.usec();

	for (int b = 2; b < 1024; ++b)
	{
		int p = NextHighestPrime(b);

		for (int a = 1; a < b; ++a)
		{
			for (int x = 0; x < b; ++x)
			{
				int y = (x + a) % p;
				if (y >= b)
					y = ((a << 16) - p + y) % a;
				xx ^= y;
			}
		}
	}

	double end1 = m_clock.usec();

	cout << end - start << endl;

	cout << end1 - end << endl;

	cout << xx << endl;
}

int main()
{
	m_clock.OnInitialize();

	//TestLog2();
	//TestOptimization();
	//TestOptimization();
	//cin.get();

	cat::wirehair::Encoder encoder;

	int block_count = 8;
	int block_bytes = 1;
	int message_bytes = block_bytes * block_count;
	u8 *message = new u8[message_bytes];

	if (!encoder.Initialize(message, message_bytes, 3 + block_bytes))
	{
		cout << "main: Failure on encoder initialization" << endl;
	}

	cin.get();

	m_clock.OnFinalize();

	return 0;
}
