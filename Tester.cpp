#include "Wirehair.hpp"
#include "SmallPRNG.hpp"
#include "Clock.hpp"
#include "GF2Matrix.hpp"
using namespace cat;

#include <iostream>
#include <fstream>
using namespace std;	

static Clock m_clock;




#include "WirehairUtil.hpp"
void TestInc()
{
	for (u16 b = 4; b < 256; ++b)
	{
		u16 p = cat::wirehair::NextPrime16(b);
		//cout << "p = " << p << endl;
		u16 smallest = b;

		for (u16 x = 0; x < b; ++x)
		{
			for (u16 a = 1; a < b; ++a)
			{
				u16 y = x;
				u16 y1 = y;
				cat::wirehair::IterateNextColumn(y, b, p, a);
				u16 y2 = y;

				for (u16 am = 1; am < b - 1; ++am)
				{
					u16 ap = (am * a) % p;
					if (ap >= b) ap = (((u32)a << 16) + ap - p) % a;

					u16 y3 = y;
					cat::wirehair::IterateNextColumn(y3, b, p, ap);

					if (y1 == y3 || y1 == y2 || y2 == y3)
					{
						if (am < smallest)
							smallest = am;
						//cout << "FAIL for x = " << x << ", a = " << a << ", am = " << am << endl;
					}
				}
			}
		}

		cout << b << " -> " << smallest << " for " << b * (b - 1) * (smallest - 1) << endl;
	}
}



int main()
{
	m_clock.OnInitialize();

	//TestInc();

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

	g_seed = 0;

	for (;;)
	{
		g_seed++;

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
