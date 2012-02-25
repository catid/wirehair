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

#define TRIALS 1000

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
/*
	cat::wirehair::GF2Matrix m;
	m.Initialize(4096);
	u32 seed = 0;
	for (;;)
	{
		m.SetSeed(seed++);
		m.Fill();
		if (m.Triangle())
			cout << "Invertible!" << endl;
		else
			cout << "Not invertible!" << endl;
	}
*/
	//TestInc();
	//TestDense();

	int block_count = 1024;
	int block_bytes = 1;
	int message_bytes = block_bytes * block_count;
	u8 *message = new u8[message_bytes];
	u8 *message_out = new u8[message_bytes];
	u8 *block = new u8[block_bytes];

	for (int ii = 0; ii < message_bytes; ++ii)
	{
		message[ii] = ii;
	}

	wirehair::Encoder encoder;

	g_c_seed = 0;

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

	g_p_seed = 22;

#if 1

	for (;;)
	{
		if (!encoder.Initialize(message, message_bytes, block_bytes))
		{
			++g_p_seed;
		}
		else
		{
			break;
		}
	}

	cout << "Measuring average execution time of encoder initialization with seed " << g_p_seed << " over " << TRIALS << " trials..." << endl;

	double start = m_clock.usec();
	for (int ii = 0; ii < TRIALS; ++ii)
	{
		encoder.Initialize(message, message_bytes, block_bytes);
	}
	double end = m_clock.usec();
	double avg_time = (end - start) / (double)TRIALS;
	cout << "Average time = " << avg_time << " usec, " << message_bytes / avg_time << " MB/s" << endl;

#endif

	start = m_clock.usec();
	encoder.Initialize(message, message_bytes, block_bytes);
	end = m_clock.usec();
	cout << ">> OKAY! encoder.Initialize in " << end - start << " usec, " << message_bytes / (end - start) << " MB/s with PSeed " << encoder.GetPSeed() << " and CSeed " << encoder.GetCSeed() << endl;

	CatsChoice prng;
	prng.Initialize(0);

	for (;;)
	{
		int blocks_needed = 0;

		cat::wirehair::Decoder decoder;
		if (!decoder.Initialize(message_out, message_bytes, block_bytes))
		{
			cout << "Decoder initialization failed!" << endl;
			cin.get();
		}

		for (u32 row_id = 0;; ++row_id)
		{
			if (prng.Next() & 1) continue;
			encoder.Generate(row_id, block);

			++blocks_needed;
			start = m_clock.usec();
			if (decoder.Decode(row_id, block))
			{
				end = m_clock.usec();
				//cout << ">> OKAY! decoder.Decode in " << end - start << " usec, " << message_bytes / (end - start) << " MB/s with PSeed " << decoder.GetPSeed() << " and CSeed " << decoder.GetCSeed() << endl;
				break;
			}
		}

		cout << "Overhead = " << blocks_needed - decoder.GetBlockCount() << endl;

		if (!memcmp(message, message_out, message_bytes))
		{
			//cout << "Match!" << endl;
		}
		else
		{
			cout << "FAAAAAIL!" << endl;

			for (int ii = 0; ii < message_bytes; ++ii)
			{
				cout << (int)message_out[ii] << endl;
			}

			cin.get();
		}
	}

	m_clock.OnFinalize();

	return 0;
}
