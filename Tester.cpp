#include "Wirehair.hpp"
#include "SmallPRNG.hpp"
using namespace cat;

#include <iostream>
using namespace std;	

int main()
{
	cat::wirehair::Encoder encoder;

	int block_count = 16;
	int block_bytes = 1500;
	int message_bytes = block_bytes * block_count;
	u8 *message = new u8[message_bytes];

	if (!encoder.Initialize(message, message_bytes, 3 + block_bytes))
	{
		return -1;
	}

	cin.get();

	return 0;
}
