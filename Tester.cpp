#include "codec_source/Wirehair.hpp"
#include "Clock.hpp"
using namespace cat;

#include <iostream>
#include <fstream>
using namespace std;	

static Clock m_clock;

#define TRIALS 10

void GenerateWeightTable()
{
	static const int N = 64;
	double table[N] = {0};

	double p1 = 0;
	table[1] = p1;

	for (int k = 2; k < N; ++k)
	{
		double p = 1. / (k * (k - 1));

		table[k] = table[k - 1] + p;
	}

	cout << "static const u32 WEIGHT_DIST[] = {" << endl;
	for (int k = 1; k < N; ++k)
	{
		double p = table[k];

		cout << "0x" << hex;
		cout << (u32)((u64)0x100000000 * p) << dec << ", ";

		if ((k & 7) == 0) cout << endl;
	}
	cout << "0xffffffff" << endl << "};" << endl;
}

int main()
{
	m_clock.OnInitialize();

	//GenerateWeightTable();
	//cin.get();

	int block_count = 64000;
	int block_bytes = 1500;
	int message_bytes = block_bytes * block_count;
	u8 *message = new u8[message_bytes];
	u8 *message_out = new u8[message_bytes];
	u8 *block = new u8[block_bytes];

	for (int ii = 0; ii < message_bytes; ++ii)
	{
		message[ii] = ii;
	}

	g_c_seed = 7;
	g_p_seed = 5;

	wirehair::Encoder encoder;
	u64 successes = 0, trials = 0;
	//double time_sum = 0;

	for (;; ++g_p_seed)
	{
		++trials;

		double start = m_clock.usec();
		wirehair::Result r = encoder.BeginEncode(message, message_bytes, block_bytes);
		double end = m_clock.usec();

		if (r)
		{
			//if (trials % 10000 == 0)
			cout << "-- FAIL! encoder.BeginEncode error " << wirehair::GetResultString(r) << ".  Success rate = " << successes / (double)trials << endl;
			//cin.get();
			//return 1;
		}
		else
		{
			//time_sum += end - start;
			//cout << time_sum / successes << endl;

			++successes;
			//if (trials % 10000 == 0)
			cout << ">> OKAY! encoder.BeginEncode in " << end - start << " usec, " << message_bytes / (end - start) << " MB/s with seeds " << g_c_seed << " and " << g_p_seed << ".  Success rate = " << successes / (double)trials << endl;
			//cin.get();
			break;
		}
	}

	CatsChoice prng;

	u32 drop_seed = 354;
	for (;;)
	{
		int blocks_needed = 0;

		cat::wirehair::Decoder decoder;
		wirehair::Result s = decoder.BeginDecode(message_out, message_bytes, block_bytes);
		if (s)
		{
			cout << "-- FAIL! decoder.BeginEncode error " << wirehair::GetResultString(s) << endl;
			cin.get();
			return 1;
		}

		prng.Initialize(drop_seed);
		for (u32 id = 0;; ++id)
		{
			if (prng.Next() & 1) continue;
			encoder.Encode(id, block);

			++blocks_needed;
			double start = m_clock.usec();
			wirehair::Result r = decoder.Decode(id, block);
			double end = m_clock.usec();

			if (r != wirehair::R_MORE_BLOCKS)
			{
				if (r == wirehair::R_WIN)
				{
					cout << "Seed = " << drop_seed << endl;
					cout << ">> OKAY! decoder.Decode in " << end - start << " usec, " << message_bytes / (end - start) << " MB/s after " << blocks_needed - decoder.BlockCount() << " extra blocks" << endl;

					if (!memcmp(message, message_out, message_bytes))
					{
						cout << "Match!" << endl;
					}
					else
					{
						cout << "FAAAAAIL!" << endl;

						//for (int ii = 0; ii < message_bytes; ++ii)
//							cout << (int)message_out[ii] << endl;

						cin.get();
					}
				}
				else
				{
					cout << "-- FAIL! decoder.Decode error " << wirehair::GetResultString(r) << endl;

					cin.get();
				}

				//cin.get();
				break;
			}
		}

		++drop_seed;
	}

	m_clock.OnFinalize();

	return 0;
}
