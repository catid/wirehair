#include "Wirehair.hpp"
#include "SmallPRNG.hpp"
using namespace cat;

#include <iostream>
using namespace std;	

int main()
{
	cat::wirehair::Encoder encoder;

<<<<<<< HEAD
	int block_count = 8;
	int block_bytes = 1;
=======
	int block_count = 64;
	int block_bytes = 1500;
>>>>>>> 79d3c673d612dd6ed9648aea23d6f1ef26976ad9
	int message_bytes = block_bytes * block_count;
	u8 *message = new u8[message_bytes];

	if (!encoder.Initialize(message, message_bytes, 3 + block_bytes))
	{
		cout << "main: Failure on encoder initialization" << endl;

		cin.get();

		return -1;
	}

	cin.get();

	return 0;
}
