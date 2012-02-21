#include "Wirehair.hpp"
#include "SmallPRNG.hpp"
#include "Clock.hpp"
#include "GF2Matrix.hpp"
using namespace cat;

#include <iostream>
#include <fstream>
using namespace std;	

static Clock m_clock;

int main()
{
	m_clock.OnInitialize();

	int block_count = 4096;
	int block_bytes = 1024 + 512 + 1;
	int message_bytes = block_bytes * block_count;
	u8 *message = new u8[message_bytes];

	for (int ii = 0; ii < message_bytes; ++ii)
	{
		message[ii] = ii;
	}

	wirehair::Encoder encoder;

	//for (;;) encoder.Initialize(message, message_bytes, block_bytes);

	for (;;)
	{
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
