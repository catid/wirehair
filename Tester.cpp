#include "Wirehair.hpp"
#include "SmallPRNG.hpp"
#include "Clock.hpp"
#include "GF2Matrix.hpp"
#include "WirehairUtil.hpp"
using namespace cat;

#include <iostream>
#include <fstream>
using namespace std;	

static Clock m_clock;

#define TRIALS 100

u32 GenerateGoodCheckSeed(int block_count)
{
	const int BLOCK_BYTES = 1;
	const int message_bytes = BLOCK_BYTES * block_count;
	u8 *message = new u8[message_bytes];

	for (int ii = 0; ii < message_bytes; ++ii)
		message[ii] = ii;

	wirehair::Encoder encoder;

	cout << "Finding best check seed for " << block_count << "..." << endl;

	u32 best_work = 0;
	u32 best_seed = 0;
	for (u32 seed = 0; seed < 10; ++seed)
	{
		u32 worked = 0;
		//g_p_seed = m_clock.cycles();
		g_p_seed = 0;
		for (int ii = 0; ii < TRIALS; ++ii)
		{
			++g_p_seed;
			if (encoder.Initialize(message, message_bytes, BLOCK_BYTES))
			{
				++worked;
			}
		}

		if (worked > best_work)
		{
			best_work = worked;
			best_seed = seed;
			cout << "New best: " << best_work << " of " << TRIALS << " with seed " << seed << endl;
		}
	}

	return best_seed;
}



int main()
{
	m_clock.OnInitialize();

	//TestInc();
	//TestDense();

	int block_count = 256;
	int block_bytes = 1024 + 512 + 1;
	int message_bytes = block_bytes * block_count;
	u8 *message = new u8[message_bytes];

	for (int ii = 0; ii < message_bytes; ++ii)
	{
		message[ii] = ii;
	}

	wirehair::Encoder encoder;

#if 1

	g_c_seed = GenerateGoodCheckSeed(block_count);
	cout << "Using CSeed : " << g_c_seed << endl;

	g_p_seed = 1;
	//for (;;) encoder.Initialize(message, message_bytes, block_bytes);

	g_p_seed = m_clock.cycles();
	u64 trials = 0, worked = 0;
	while (++trials < TRIALS)
	{
		++g_p_seed;
		if (encoder.Initialize(message, message_bytes, block_bytes))
		{
			++worked;
		}
	}

	cout << "Success rate = " << worked / (double)trials << endl;

#endif

	//g_c_seed = 0;
	g_p_seed = 0;

	for (;;)
	{
		double start = m_clock.usec();
		u32 clocks = m_clock.cycles();
		bool success = encoder.Initialize(message, message_bytes, block_bytes);
		clocks = m_clock.cycles() - clocks;
		double end = m_clock.usec();

		if (success)
		{
			cout << ">> OKAY! encoder.Initialize in " << clocks << " clocks and " << end - start << " usec with PSeed " << encoder.GetPSeed() << " and CSeed " << encoder.GetCSeed() << endl;

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
			cout << "-- FAIL: encoder.Initialize in " << clocks << " clocks and " << end - start << " usec with PSeed " << encoder.GetPSeed() << " and CSeed " << encoder.GetCSeed() << endl;
		}

		g_p_seed++;
	}

	m_clock.OnFinalize();

	return 0;
}
