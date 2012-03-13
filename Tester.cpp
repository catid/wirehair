#include "codec_source/Wirehair.hpp"
#include "Clock.hpp"
using namespace cat;

#include <iostream>
#include <fstream>
using namespace std;	

static Clock m_clock;

int main()
{
	m_clock.OnInitialize();

	for (int ii = 64000; ii <= 64000; ++ii)
	{
		int block_count = ii;
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

		double start = m_clock.usec();
		wirehair::Result r = encoder.BeginEncode(message, message_bytes, block_bytes);
		double end = m_clock.usec();

		if (r)
		{
			cout << "-- FAIL! N=" << encoder.BlockCount() << " encoder.BeginEncode error " << wirehair::GetResultString(r) << endl;
			cin.get();
			continue;
		}
		else
		{
			cout << ">> OKAY! N=" << encoder.BlockCount() << " encoder.BeginEncode in " << end - start << " usec, " << message_bytes / (end - start) << " MB/s" << endl;
			//cin.get();
		}

		CatsChoice prng;
		cat::wirehair::Decoder decoder;

		u32 overhead_sum = 0, overhead_trials = 0;
		u32 drop_seed = 44519;
		double time_sum = 0;
		const int trials = 100000;
		for (int jj = 0; jj < trials; ++jj)
		{
			int blocks_needed = 0;

			wirehair::Result s = decoder.BeginDecode(message_out, message_bytes, block_bytes);
			if (s)
			{
				cout << "-- FAIL! N=" << decoder.BlockCount() << " decoder.BeginEncode error " << wirehair::GetResultString(s) << endl;
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
						u32 overhead = blocks_needed - decoder.BlockCount();
						overhead_sum += overhead;
						++overhead_trials;

						//cout << ">> OKAY! N=" << decoder.BlockCount() << " decoder.Decode in " << end - start << " usec, " << message_bytes / (end - start) << " MB/s after " << overhead << " extra blocks.  Average extra = " << overhead_sum / (double)overhead_trials << ". Seed = " << drop_seed << endl;
						time_sum += end - start;

						if (!memcmp(message, message_out, message_bytes))
						{
							//cout << "Match!" << endl;
						}
						else
						{
							cout << "FAAAAAIL! Seed = " << drop_seed << endl;

							//for (int ii = 0; ii < message_bytes; ++ii) cout << (int)message_out[ii] << endl;

							cin.get();
						}
					}
					else
					{
						cout << "-- FAIL!  N=" << decoder.BlockCount() << " decoder.Decode error " << wirehair::GetResultString(r) << " from drop seed " << drop_seed << endl;

						overhead_sum += 1;
						++overhead_trials;

						//cin.get();
					}

					//cin.get();
					break;
				}
			}

			++drop_seed;
		}

		double avg_time = time_sum/trials;
		cout << "N=" << decoder.BlockCount() << " decoder.Decode in " << avg_time << " usec, " << message_bytes / avg_time << " MB/s.  Average overhead = " << overhead_sum / (double)overhead_trials << endl;
	}

	m_clock.OnFinalize();

	return 0;
}
