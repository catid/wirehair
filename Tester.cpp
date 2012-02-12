#include "Wirehair.hpp"
#include "SmallPRNG.hpp"
#include "Clock.hpp"
using namespace cat;

#include <iostream>
using namespace std;	

static Clock m_clock;

int main()
{
	m_clock.OnInitialize();

	cat::wirehair::Encoder encoder;

	int block_count = 1024;
	int block_bytes = (1500 & ~255) + 1;
	int message_bytes = block_bytes * block_count;
	u8 *message = new u8[message_bytes];

	for (int ii = 0; ii < message_bytes; ++ii)
	{
		message[ii] = ii;
	}

	for (;;)
	{
		double start = m_clock.usec();
		bool success = encoder.Initialize(message, message_bytes, 3 + block_bytes);
		double end = m_clock.usec();

		cout << "main: encoder.Initialize in " << end - start << " usec" << endl;

		if (success)
		{
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
	}

	m_clock.OnFinalize();

	return 0;
}
