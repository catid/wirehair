#include "Wirehair.hpp"
#include "SmallPRNG.hpp"
#include "Clock.hpp"
using namespace cat;

#include <iostream>
using namespace std;	

static Clock m_clock;

static cat::wirehair::Encoder encoder;
static bool success;
static u8 *message;
static int message_bytes;
static int block_bytes;

void TestLoop()
{
	success = encoder.Initialize(message, message_bytes, 3 + block_bytes);
}

int main()
{
	m_clock.OnInitialize();

	int block_count = 256;
	block_bytes = 1024 + 512 + 1;
	message_bytes = block_bytes * block_count;
	message = new u8[message_bytes];

	for (int ii = 0; ii < message_bytes; ++ii)
	{
		message[ii] = ii;
	}

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
					//cout << "Block " << ii << " doesn't match: " << (int)block[3] << endl;
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
