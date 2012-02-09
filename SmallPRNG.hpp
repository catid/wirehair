/*
	Copyright (c) 2012 Christopher A. Taylor.  All rights reserved.

	Redistribution and use in source and binary forms, with or without
	modification, are permitted provided that the following conditions are :

	* Redistributions of source code must retain the above copyright notice,
	  this list of conditions and the following disclaimer.
	* Redistributions in binary form must reproduce the above copyright notice,
	  this list of conditions and the following disclaimer in the documentation
	  and/or other materials provided with the distribution.
	* Neither the name of libperFECt nor the names of its contributors may be used
	  to endorse or promote products derived from this software without
	  specific prior written permission.

	THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
	AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
	IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
	ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
	LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
	CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
	SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
	INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
	CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
	ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
	POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef CAT_SMALL_PRNG_HPP
#define CAT_SMALL_PRNG_HPP

#include "Platform.hpp"

namespace cat {


/*
	Notes on combining generators:

	All LCG, MWC, and XORS generators are safe to combine with simple addition
	since the periods of all of the generators here are relatively prime.
	In this case the overall period will be the sum of the periods.

	If you need to achieve a period of 2^^X, then the period of the generators
	should be at least 2^^(3X).  So, combine MWC with XORS or LCG to make a
	generator that would be good for 2^^32 output numbers.
*/

/*
	Performance measurements:

	Generator CatidL32_1 operates at 153 million numbers / second
	Generator Catid32_1a operates at 314 million numbers / second
	Generator Catid32_1b operates at 311 million numbers / second
	Generator Catid32_1c operates at 306 million numbers / second
	Generator Catid32_1d operates at 310 million numbers / second
	Generator Catid32_2 operates at 336 million numbers / second
	Generator Catid32_2b operates at 346 million numbers / second
	Generator Catid32_2c operates at 348 million numbers / second
	Generator JKISS32_nomult operates at 331 million numbers / second
	Generator Catid32S_1a operates at 329 million numbers / second
	Generator Catid32S_1d operates at 334 million numbers / second
	Generator Catid32S_4 operates at 426 million numbers / second

	CatsChoice is an implementation of Catid32S_5, based on Catid32S_4.
*/

/*
	Linear Congruential Generator (LCG) with power-of-two modulus

	Guidelines:
		M = 2^^b
		A = Chosen so that A - 1 is a multiple of 4, since M is a multiple of 4
			For other M, A - 1 should be divisible by all prime factors of M
		C = Relatively prime to M
			So it should be odd
			And I think it should be close to M in magnitude

	Output: b bits
	Period: 2^^b
	Issues:
		Lower bits have lower period, and lowest bit alternates
*/
template <u32 A, u32 C>
class LCG32
{
	u32 _x;

public:
	CAT_INLINE void Initialize(u32 seed)
	{
		_x = seed;
	}

	CAT_INLINE void MixSeed(u32 seed)
	{
		Next();
		_x ^= seed;
	}

	CAT_INLINE u32 Next()
	{
		return (_x = A * _x + C);
	}
};

// 64-bit version:
template <u64 A, u64 C>
class LCG64
{
	u64 _x;

public:
	CAT_INLINE void Initialize(u64 seed)
	{
		_x = seed;
	}

	CAT_INLINE void MixSeed(u64 seed)
	{
		Next();
		_x ^= seed;
	}

	CAT_INLINE u64 Next()
	{
		return (_x = A * _x + C);
	}
};


/*
	from "TABLES OF LINEAR CONGRUENTIAL GENERATORS OF DIFFERENT SIZES AND GOOD LATTICE STRUCTURE" (1999)
	by Pierre L'ecuyer
*/
typedef LCG32<2891336453, 1234567897> LecuyerLCG32_1;
typedef LCG32<29943829, 1234567897> LecuyerLCG32_2;
typedef LCG32<32310901, 1234567897> LecuyerLCG32_3;
typedef LCG64<2862933555777941757ULL, 7891234567891234567ULL> LecuyerLCG64_1;
typedef LCG64<3202034522624059733ULL, 7891234567891234567ULL> LecuyerLCG64_2;
typedef LCG64<3935559000370003845ULL, 7891234567891234567ULL> LecuyerLCG64_3;


/*
	Multiply With Carry (MWC) PRNG
	by George Marsaglia

	Guidelines:
		B = 2^^32 (base)
		A = Chosen such that A*B-1 and A*B/2 -1 are prime

	Output: 32 bits
	Period: [2^^32 * A] / 2 - 1
	Issues:
		Will get stuck if both M and C are zero
		High bits tend to be less random
*/
template <u32 A, u32 M0, u32 C0>
class MWC
{
	u32 _m, _c;

public:
	CAT_INLINE void Initialize(u32 seed)
	{
		_m = M0 ^ seed;
		_c = C0;
	}

	CAT_INLINE void MixSeed(u32 seed)
	{
		Next();
		_m ^= seed;

		if (_m == 0 && _c == 0)
			Initialize(seed);
	}

	CAT_INLINE u32 Next()
	{
		u64 t = (u64)A * _m + _c;
		_m = (u32)t;
		_c = (u32)(t >> 32);
		return _m;
	}
};

/*
	Maximal safe prime version, and the maximum period version
	from the wikipedia article

	MaxSafeMWC period = 9223371654602686463 (prime)
	MaximalMWC period = 9223371873646018559 = 773 * 1621 * 7360837163623
*/
typedef MWC<4294967118, 21987643, 1732654> MaxSafeMWC;
typedef MWC<4294967220, 21987643, 1732654> MaximalMWC;
/*
	from "Good Practice in (Pseudo) Random Number Generation for Bioinformatics Applications" (2010)
	by David Jones

	DJonesMWC1 period = 9222549758923505663 (prime)
	DJonesMWC2 period = 9119241012177272831 (prime)
*/
typedef MWC<4294584393, 43219876, 6543217> DJonesMWC1;
typedef MWC<4246477509, 21987643, 1732654> DJonesMWC2;
/*
	These are parameters that I generated:

	I want to find the largest few values of A that satisfy
	this requirements.  Furthermore, I want A to have a minimal
	number of prime factors; if this is the case, then the
	output will look random in addition to having a large period.

	So I set out to devise a fast 64-bit primality tester.
	The best approach I know of for these small sizes is the
	probabilistic Rabin-Miller primality test.  After a few
	iterations I can conclude the number is probably prime
	and then look at the prime factorization of A to select
	candidates that may go well together.  Wolfram Alpha was
	used to check that the values are really prime.

	I came up with these results with one big factor:

	-- Candidate 0xffffbe17.  Factors = 3, 1431650141
	-- Candidate 0xffff4b9f.  Factors = 3, 1431640373
	-- Candidate 0xffff0207.  Factors = 3, 1431634093
	-- Candidate 0xfffe1495.  Factors = 3, 1431613831
	-- Candidate 0xfffd8b79.  Factors = 3, 1431602131
	-- Candidate 0xfffd6389.  Factors = 3, 1431598723
	-- Candidate 0xfffd21a7.  Factors = 3, 1431593101
	-- Candidate 0xfffd1361.  Factors = 3, 1431591883
	...

	So I guess that 3 is always a factor of A...
	I then produced CatsChoice-like generators using pairs
	of these A values and tested them with BigCrush:

		A1 = 0xffffbe17, A2 = 0xffff4b9f <- Failed test 48
		A1 = 0xffff0207, A2 = 0xfffe1495 <- Passed all tests
		A1 = 0xfffd8b79, A2 = 0xfffd6389 <- Passed all tests
		A1 = 0xfffd21a7, A2 = 0xfffd1361 <- Passed all tests (chosen)
*/
typedef MWC<0xfffd21a7, 43219876, 6543217> CatMWC1;
typedef MWC<0xfffd1361, 21987643, 1732654> CatMWC2;


/*
	Type-I XOR Shift Linear Feedback Shift Register (LFSR) PRNG
	from "Xorshift RNGs" (2003)
	by George Marsaglia

	Guidelines:
		Choose shifts A,B,C from Marsaglia's comprehensive list

	Output: b bits
	Period: 2^^b - 1
		32-bit period factors = 3󬊁7�257�65537
		64-bit period factors = 3󬊁7�257�641�65537�6700417
	Issues:
		Halts on zero
		Linear relationship between blocks of b + 1 consecutive bits
*/
template <int A, int B, int C, u32 X0>
class XORShift32
{
	u32 _x;

public:
	CAT_INLINE void Initialize(u32 seed)
	{
		_x = X0 ^ seed;

		if (_x == 0)
			_x = ~(u32)0;
	}

	CAT_INLINE void MixSeed(u32 seed)
	{
		Next();
		_x += seed;

		if (_x == 0)
			Initialize(seed);
	}

	CAT_INLINE u32 Next()
	{
		register u32 x = _x;
		x ^= x << A;
		x ^= x >> B;
		x ^= x << C;
		return (_x = x);
	}
};

// 64-bit version:
template <int A, int B, int C, u64 X0>
class XORShift64
{
	u64 _x;

public:
	CAT_INLINE void Initialize(u64 seed)
	{
		_x = X0 ^ seed;

		if (_x == 0)
			_x = ~(u64)0;
	}

	CAT_INLINE void MixSeed(u64 seed)
	{
		Next();
		_x += seed;

		if (_x == 0)
			Initialize(seed);
	}

	CAT_INLINE u64 Next()
	{
		register u64 x = _x;
		x ^= x << A;
		x ^= x >> B;
		x ^= x << C;
		return (_x = x);
	}
};

/*
	Chose these at random from the list
*/
typedef XORShift32<5, 7, 22, 0x56A53625> XORShift32_1;
typedef XORShift32<8, 7, 23, 0x56A53625> XORShift32_2;
typedef XORShift32<3, 13, 7, 0x56A53625> XORShift32_3;
typedef XORShift32<5, 7, 22, 234567891> XORShift32_4;	// Used in JKISS32
typedef XORShift64<21, 17, 30, 0x4A3CE93555573AABULL> XORShift64_1;	// Used in JLKISS64
typedef XORShift64<17, 23, 29, 0x4A3CE93555573AABULL> XORShift64_2;
typedef XORShift64<16, 21, 35, 0x4A3CE93555573AABULL> XORShift64_3;


/*
	Weyl Generator PRNG
	from "Some long-period random number generators using shifts and xor" (2007)
	by Richard. P. Brent

	Guidelines:
		A = Odd, close to 2^^(b-1) * (sqrt(5) - 1)
		For b=32, close to 2654435769
		For b=64, close to 11400714819323198485
		Weak generator for combining with other generators.

	Output: b bits
	Period: 2^^b
	Issues:
		Horrible in general
*/
template <u32 A, u32 X0>
class WeylGenerator32
{
	u32 _x;

public:
	CAT_INLINE void Initialize(u32 seed)
	{
		_x = X0 ^ seed;
	}

	CAT_INLINE void MixSeed(u32 seed)
	{
		Next();
		_x ^= seed;
	}

	CAT_INLINE u32 Next()
	{
		return (_x += A);
	}
};

// 64-bit version:
template <u64 A, u64 X0>
class WeylGenerator64
{
	u64 _x;

public:
	CAT_INLINE void Initialize(u64 seed)
	{
		_x = X0 ^ seed;
	}

	CAT_INLINE void MixSeed(u64 seed)
	{
		Next();
		_x ^= seed;
	}

	CAT_INLINE u64 Next()
	{
		return (_x += A);
	}
};

// Close to choice criterion from Brent
typedef WeylGenerator32<2654435769, 1223235218> Weyl32_1;
typedef WeylGenerator64<11400714819323198485ULL, 0xFEE9095D248AB2ABULL> Weyl64_1;

/*
	from "Good Practice in (Pseudo) Random Number Generation for Bioinformatics Applications" (2010)
	by David Jones
*/
typedef WeylGenerator32<1411392427, 123456789> Weyl32_2;


/*
	Add With Carry (AWC) PRNG
	by George Marsaglia

	Weak generator for combining with other generators.

	Output: 32 bits
	Period: <2^^28 with random seeding, ~2^^31 with chosen values
	Issues:
		Cannot be seeded without seriously affecting the period
		Horrible in general
*/
template <u32 Z0, u32 W0>
class AWC
{
	u32 _z, _w, _c;

public:
	CAT_INLINE void Initialize(u32 seed)
	{
		_z = Z0;
		_w = W0;
		_c = 0;
	}

	CAT_INLINE void MixSeed(u32 seed)
	{
	}

	CAT_INLINE u32 Next()
	{
		u32 t = _z + _w + _c;
		_z = _w;
		_c = t >> 31;
		_w = t & 0x7fffffff;
		return _w;
	}
};

/*
	Factors 3󬊁7�257�641�65537�6700417 cannot be combined with XORShift

	After a short random search I came up with these values:
		(2686646964, 3741327162) period=4202554829 combo period=4202554829
		(2026632552, 1483949311) period=4150427771 combo period=4150427771
		(3631468667, 1476107563) period=3635438413 combo period=3635438413
*/
typedef AWC<2686646964, 3741327162> AWC32_1;
typedef AWC<2026632552, 1483949311> AWC32_2;
typedef AWC<3631468667, 1476107563> AWC32_3;
typedef AWC<345678912, 456789123> AWC32_4;		// from JKISS32


/*
	Single-bit LFSR PRNG

	Guidelines:
		Choose taps wisely

	Output: b bits
	Period: 2^^b - 1
	Issues:
		Halts on zero
		Horrible in general
*/
template <u32 TAP_MASK>
class SingleBitLFSR32
{
	u32 _x;

public:
	CAT_INLINE void Initialize(u32 seed)
	{
		_x = seed;

		if (_x == 0)
			_x = ~(u32)0;
	}

	CAT_INLINE void MixSeed(u32 seed)
	{
		Next();
		_x += seed;

		if (_x == 0)
			Initialize(seed);
	}

	CAT_INLINE bool Next()
	{
		_x = (_x >> 1) ^ (-(s32)(_x & 1) & TAP_MASK); 
		return (_x & 1) != 0;
	}
};

// 64-bit version:
template <u64 TAP_MASK>
class SingleBitLFSR64
{
	u64 _x;

public:
	CAT_INLINE void Initialize(u64 seed)
	{
		_x = seed;

		if (_x == 0)
			_x = ~(u64)0;
	}

	CAT_INLINE void MixSeed(u64 seed)
	{
		Next();
		_x += seed;

		if (_x == 0)
			Initialize(seed);
	}

	CAT_INLINE bool Next()
	{
		_x = (_x >> 1) ^ (-(_x & 1) & TAP_MASK); 
		return (_x & 1) != 0;
	}
};

/*
	From an LFSR taps table floating around the net
	32-bit characteristic polynomial: x^32 + x^22 + x + 1
	64-bit characteristic polynomial: x^64 + x^63 + x^61 + x^60
*/
typedef SingleBitLFSR32<0x80200003> SingleBitLFSR32_1;
typedef SingleBitLFSR64<0xD800000000000000ULL> SingleBitLFSR64_1;

/*
	From Wikipedia
	Characteristic polynomial: x^32 + x^31 + x^29 + x + 1
*/
typedef SingleBitLFSR32<0xD0000001> SingleBitLFSR32_2;


/*
	Catid's KISS with LFSR

	Period of combined generators should be about twice as long.

	Always adds in generator 1 result.
	Uses an LFSR to gate generators 2 and 3.
*/
template<typename T, class LFSR, class G1, class G2, class G3>
class CAT_EXPORT CKISSL
{
	LFSR _lfsr;
	G1 _g1;
	G2 _g2;
	G3 _g3;

public:
	void Initialize(T seed)
	{
		_lfsr.Initialize(seed);
		_g1.Initialize(seed);
		_g2.Initialize(seed);
		_g3.Initialize(seed);
	}

	void MixSeed(T seed)
	{
		_lfsr.MixSeed(seed);
		_g1.MixSeed(seed);
		_g2.MixSeed(seed);
		_g3.MixSeed(seed);
	}

	T Next()
	{
		T result = _g1.Next();

		if (_lfsr.Next())
			result += _g2.Next();
		else
			result += _g3.Next();

		return result;
	}
};

/*
	Period of ~2^^128

	Good for making the generator harder to analyze from its output.

	Passes all BigCrush tests.
*/
typedef CKISSL<u32, SingleBitLFSR32_2, MaxSafeMWC, XORShift32_1, LecuyerLCG32_1> CatidL32_1;


/*
	Catid's KISS

	Read notes on combining generators for proper usage.

	Mixes results from all generators.
*/
template<typename T, class G1, class G2, class G3>
class CAT_EXPORT CKISS
{
	G1 _g1;
	G2 _g2;
	G3 _g3;

public:
	void Initialize(T seed)
	{
		_g1.Initialize(seed);
		_g2.Initialize(seed);
		_g3.Initialize(seed);
	}

	void MixSeed(T seed)
	{
		_g1.MixSeed(seed);
		_g2.MixSeed(seed);
		_g3.MixSeed(seed);
	}

	T Next()
	{
		return _g1.Next() + _g2.Next() + _g3.Next();
	}
};

/*
	Period of ~2^^127

	Fails BigCrush tests:
		23  ClosePairs mNP2S, t = 5         0.9994
*/
typedef CKISS<u32, MaxSafeMWC, XORShift32_1, LecuyerLCG32_1> Catid32_1;
/*
	Period of ~2^^127

	Passes all BigCrush tests.
*/
typedef CKISS<u32, MaximalMWC, XORShift32_1, LecuyerLCG32_1> Catid32_1a;
/*
	Period of ~2^^127

	Passes all BigCrush tests.
*/
typedef CKISS<u32, MaxSafeMWC, XORShift32_2, LecuyerLCG32_1> Catid32_1b;
/*
	Period of ~2^^127

	Passes all BigCrush tests.
*/
typedef CKISS<u32, MaxSafeMWC, XORShift32_1, LecuyerLCG32_2> Catid32_1c;
/*
	Period of ~2^^127

	Passes all BigCrush tests.
*/
typedef CKISS<u32, MaximalMWC, XORShift32_2, LecuyerLCG32_2> Catid32_1d;
/*
	Period of ~2^^96

	Passes all BigCrush tests.
*/
typedef CKISS<u32, XORShift32_1, AWC32_1, Weyl32_1> Catid32_2;
/*
	Period of ~2^^96

	Fails BigCrush tests:
		50  SampleProd, t = 8               5.2e-4
*/
typedef CKISS<u32, XORShift32_1, AWC32_2, Weyl32_1> Catid32_2a;
/*
	Period of ~2^^96

	Passes all BigCrush tests.
*/
typedef CKISS<u32, XORShift32_2, AWC32_1, Weyl32_1> Catid32_2b;
/*
	Period of ~2^^96

	Passes all BigCrush tests.
*/
typedef CKISS<u32, XORShift32_1, AWC32_1, Weyl32_2> Catid32_2c;
/*
	Period of ~2^^96

	Fails BigCrush tests:
		38  Run, r = 0                      6.5e-7
*/
typedef CKISS<u32, XORShift32_2, AWC32_2, Weyl32_2> Catid32_2d;
/*
	Equivalent to the JKISS32 generator with no multiplies

	Passes all BigCrush tests.
*/
typedef CKISS<u32, XORShift32_4, AWC32_4, Weyl32_2> JKISS32_nomult;


/*
	Catid's Smootch

	Read notes on combining generators for proper usage.

	Mixes just two generators.
*/
template<typename T, class G1, class G2>
class CAT_EXPORT CSmootch
{
	G1 _g1;
	G2 _g2;

public:
	void Initialize(T seed)
	{
		_g1.Initialize(seed);
		_g2.Initialize(seed);
	}

	void MixSeed(T seed)
	{
		_g1.MixSeed(seed);
		_g2.MixSeed(seed);
	}

	T Next()
	{
		return _g1.Next() + _g2.Next();
	}
};

/*
	Period of ~2^^95

	Fails BigCrush tests:
		77  RandomWalk1 R (L=1000, r=20)    3.4e-4
*/
typedef CSmootch<u32, XORShift32_1, MaxSafeMWC> Catid32S_1;
/*
	Period of ~2^^95

	Passes all BigCrush tests.
*/
typedef CSmootch<u32, XORShift32_2, MaxSafeMWC> Catid32S_1a;
/*
	Period of ~2^^95

	Fails BigCrush tests:
		11  CollisionOver, t = 21          6.0e-04
*/
typedef CSmootch<u32, XORShift32_3, MaxSafeMWC> Catid32S_1b;
/*
	Period of ~2^^95

	Fails BigCrush tests:
		14  BirthdaySpacings, t = 3        3.4e-04
*/
typedef CSmootch<u32, XORShift32_1, MaximalMWC> Catid32S_1c;
/*
	Period of ~2^^95

	Passes all BigCrush tests.
*/
typedef CSmootch<u32, XORShift32_2, MaximalMWC> Catid32S_1d;
/*
	Period of ~2^^95

	Fails BigCrush tests:
		15  BirthdaySpacings, t = 4          eps
*/
typedef CSmootch<u32, MaxSafeMWC, LecuyerLCG32_1> Catid32S_2;
/*
	Period of ~2^^95

	Fails BigCrush tests:
		15  BirthdaySpacings, t = 4        2.8e-86
*/
typedef CSmootch<u32, MaxSafeMWC, LecuyerLCG32_2> Catid32S_2a;
/*
	Period of ~2^^95

	Fails BigCrush tests:
		 15  BirthdaySpacings, t = 4          eps
		 19  BirthdaySpacings, t = 8          eps
*/
typedef CSmootch<u32, MaximalMWC, LecuyerLCG32_1> Catid32S_2b;
/*
	Period of ~2^^95

	Fails BigCrush tests:
		 15  BirthdaySpacings, t = 4          eps
		 19  BirthdaySpacings, t = 8          eps
*/
typedef CSmootch<u32, MaximalMWC, LecuyerLCG32_2> Catid32S_2c;
/*
	Period of ~2^^95

	Fails BigCrush tests:
		15  BirthdaySpacings, t = 4        1.8e-98
*/
typedef CSmootch<u32, MaxSafeMWC, LecuyerLCG32_3> Catid32S_2d;
/*
	Period of ~2^^64

	Fails BigCrush tests:
		  2  SerialOver, r = 22               eps
		 19  BirthdaySpacings, t = 8          eps
		 21  BirthdaySpacings, t = 16         eps
		 69  MatrixRank, L=1000, r=26         eps
		 70  MatrixRank, L=5000               eps
		 81  LinearComp, r = 29             1 - eps1
*/
typedef CSmootch<u32, XORShift32_1, LecuyerLCG32_1> Catid32S_3;
/*
	Period of ~2^^64

	Fails BigCrush tests:
		6  CollisionOver, t = 3            6.0e-4
		8  CollisionOver, t = 7             eps
		10  CollisionOver, t = 14          7.9e-71
		12  CollisionOver, t = 21           8.8e-8
		19  BirthdaySpacings, t = 8          eps
		21  BirthdaySpacings, t = 16         eps
		27  SimpPoker, r = 27              9.1e-13
		58  AppearanceSpacings, r = 27     1 - eps1
		69  MatrixRank, L=1000, r=26         eps
		70  MatrixRank, L=5000               eps
		81  LinearComp, r = 29             1 - eps1
		87  LongestHeadRun, r = 27          1.9e-6
		102  Run of bits, r = 27              eps
*/
typedef CSmootch<u32, XORShift32_2, LecuyerLCG32_1> Catid32S_3a;
/*
	Period of ~2^^64

	Fails BigCrush tests:
		2  SerialOver, r = 22               eps
		19  BirthdaySpacings, t = 8       5.1e-167
		21  BirthdaySpacings, t = 16         eps
		69  MatrixRank, L=1000, r=26         eps
		70  MatrixRank, L=5000               eps
		81  LinearComp, r = 29             1 - eps1
*/
typedef CSmootch<u32, XORShift32_3, LecuyerLCG32_1> Catid32S_3b;
/*
	Period of ~2^^64

	Fails BigCrush tests:
		2  SerialOver, r = 22               eps
		8  CollisionOver, t = 7             eps
		10  CollisionOver, t = 14          4.0e-68
		12  CollisionOver, t = 21           8.8e-8
		19  BirthdaySpacings, t = 8          eps
		21  BirthdaySpacings, t = 16         eps
		27  SimpPoker, r = 27               4.0e-6
		58  AppearanceSpacings, r = 27     1 - eps1
		69  MatrixRank, L=1000, r=26         eps
		70  MatrixRank, L=5000               eps
		81  LinearComp, r = 29             1 - eps1
		87  LongestHeadRun, r = 27         1.3e-13
		102  Run of bits, r = 27              eps
*/
typedef CSmootch<u32, XORShift32_2, LecuyerLCG32_2> Catid32S_3c;
/*
	Period of ~2^^64

	Fails BigCrush tests:
		2  SerialOver, r = 22              9.5e-5
		8  CollisionOver, t = 7             eps
		10  CollisionOver, t = 14          4.0e-68
		12  CollisionOver, t = 21           8.8e-8
		19  BirthdaySpacings, t = 8          eps
		21  BirthdaySpacings, t = 16         eps
		27  SimpPoker, r = 27               4.0e-6
		58  AppearanceSpacings, r = 27     1 - eps1
		69  MatrixRank, L=1000, r=26         eps
		70  MatrixRank, L=5000               eps
		81  LinearComp, r = 29             1 - eps1
		87  LongestHeadRun, r = 27         6.7e-12
		102  Run of bits, r = 27              eps
*/
typedef CSmootch<u32, XORShift32_2, LecuyerLCG32_3> Catid32S_3d;
/*
	Period of ~2^^126

	Passes all BigCrush tests.
*/
typedef CSmootch<u32, MaxSafeMWC, DJonesMWC1> Catid32S_4;
/*
	Period of ~2^^126

	Fails BigCrush tests:
		76  RandomWalk1 C (L=1000, r=0)     0.9995
*/
typedef CSmootch<u32, MaxSafeMWC, MaximalMWC> Catid32S_4a;
/*
	Period of ~2^^126

	Fails BigCrush tests:
		11  CollisionOver, t = 21          6.7e-05
*/
typedef CSmootch<u32, MaxSafeMWC, DJonesMWC2> Catid32S_4b;
/*
	Period of ~2^^126

	Passes all BigCrush tests.
*/
typedef CSmootch<u32, CatMWC1, CatMWC2> Catid32S_5;


/*
	This is a unified implementation of my favorite generator
	that is designed to generate up to 2^^32 numbers per seed.

	Its period is about 2^^126 and passes all BigCrush tests.
	It is the fastest generator I could find that passes all tests.

	Furthermore, the input seeds are hashed to avoid linear
	relationships between the input seeds and the low bits of
	the first few outputs.
*/
class CAT_EXPORT CatsChoice
{
	u64 _x, _y;

public:
	CAT_INLINE void Initialize(u32 x, u32 y)
	{
		// Based on the mixing functions of MurmurHash3
		static const u64 C1 = 0xff51afd7ed558ccdULL;
		static const u64 C2 = 0xc4ceb9fe1a85ec53ULL;

		x += y;
		y += x;

		u64 seed_x = 0x9368e53c2f6af274ULL ^ x;
		u64 seed_y = 0x586dcd208f7cd3fdULL ^ y;

		seed_x *= C1;
		seed_x ^= seed_x >> 33;
		seed_x *= C2;
		seed_x ^= seed_x >> 33;

		seed_y *= C1;
		seed_y ^= seed_y >> 33;
		seed_y *= C2;
		seed_y ^= seed_y >> 33;

		_x = seed_x;
		_y = seed_y;

		Next();
	}

	CAT_INLINE void Initialize(u32 seed)
	{
		Initialize(seed, seed);
	}

	CAT_INLINE u32 Next()
	{
		_x = (u64)0xfffd21a7 * (u32)_x + (u32)(_x >> 32);
		_y = (u64)0xfffd1361 * (u32)_y + (u32)(_y >> 32);
		return (u32)_x + (u32)_y;
	}
};


} // namespace cat

#endif // CAT_SMALL_PRNG_HPP
