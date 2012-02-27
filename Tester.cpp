#include "codec_source/Wirehair.hpp"
#include "Clock.hpp"
using namespace cat;

#include <iostream>
#include <fstream>
using namespace std;	

static Clock m_clock;

#define TRIALS 10

int main()
{
	m_clock.OnInitialize();

	int block_count = 1024;
	int block_bytes = 1500;
	int message_bytes = block_bytes * block_count;
	u8 *message = new u8[message_bytes];
	u8 *message_out = new u8[message_bytes];
	u8 *block = new u8[block_bytes];

	for (int ii = 0; ii < message_bytes; ++ii)
	{
		message[ii] = ii;
	}

	g_c_seed = 0;
	g_p_seed = 0;

	wirehair::Encoder encoder;

	double start = m_clock.usec();
	wirehair::Result r = encoder.BeginEncode(message, message_bytes, block_bytes);
	double end = m_clock.usec();
	if (r)
	{
		cout << "-- FAIL! encoder.BeginEncode error " << wirehair::GetResultString(r) << endl;
		cin.get();
		return 1;
	}

	cout << ">> OKAY! encoder.BeginEncode in " << end - start << " usec, " << message_bytes / (end - start) << " MB/s" << endl;

	CatsChoice prng;

	u32 drop_seed = 0;
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
			start = m_clock.usec();
			wirehair::Result r = decoder.Decode(id, block);
			end = m_clock.usec();

			if (r != wirehair::R_MORE_BLOCKS)
			{
				if (r == wirehair::R_WIN)
				{
					cout << ">> OKAY! decoder.Decode in " << end - start << " usec, " << message_bytes / (end - start) << " MB/s after " << blocks_needed - decoder.BlockCount() << " extra blocks" << endl;

					if (!memcmp(message, message_out, message_bytes))
					{
						cout << "Match!" << endl;
					}
					else
					{
						cout << "Seed = " << drop_seed << endl;
						cout << "FAAAAAIL!" << endl;

						for (int ii = 0; ii < message_bytes; ++ii)
							cout << (int)message_out[ii] << endl;

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
