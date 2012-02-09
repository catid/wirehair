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

	int block_count = 8;
	int block_bytes = 1;
	int message_bytes = block_bytes * block_count;
	u8 *message = new u8[message_bytes];

	if (!encoder.Initialize(message, message_bytes, 3 + block_bytes))
	{
		cout << "main: Failure on encoder initialization" << endl;
	}

	cin.get();

	m_clock.OnFinalize();

	return 0;
}
