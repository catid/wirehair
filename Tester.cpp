#include "codec_source/Wirehair.hpp"
#include "Clock.hpp"
using namespace cat;

#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;	

static Clock m_clock;

/*
	In order to support a skip list for seeds, either need to support zero additions
	in the codec or need to have a better exception list for seeds.  Having a better
	exception list would be more efficient, however the list may be larger than expected.
	The advantage of a skip list is that seeds can be more easily blacklisted.
*/

/*
	Process for selecting better heavy seeds:

	12x18 binary rows above -- too many bits to try all of them (216 bits)
	6x18 heavy rows below

	we want to choose a single set of heavy row values to maximize the
	likelihood that the matrix is invertible for all binary rows.

	If missing more than 6 columns it is not possible, so ignore that case.

	For now assume that binary rows are full rank:

	Heavy rows are eliminated left to right by XORing the left-most byte across.

	We want to reduce or eliminate the chance that doing this selectively is able
	to cause zero columns.

	This means that any XOR combination of bytes to the left cannot be equal to
	a byte to the right.  This is totally doable.  2^^18 or so

	Finally, for the last 6 columns it is important that adding the heavy rows
	together cannot cause zeroes.  This can be measured by simulation.

	Conclusion:

		There is no practical selection of matrix values that are
	better than any other selection of matrix values.  It seems as
	though the best option is to use a good peeling seed.

		Handling weak seeds therefore is just a matter of thorough
	testing.
*/

// Returns the count of bits set in the input for types up to 128 bits
template<typename T> CAT_INLINE T BitCount(T v)
{
	// From Stanford Bit Twiddling Hacks collection
	// http://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetParallel
	v = v - ((v >> 1) & (T)~(T)0/3);
	v = (v & (T)~(T)0/15*3) + ((v >> 2) & (T)~(T)0/15*3);
	v = (v + (v >> 4)) & (T)~(T)0/255*15;
	return (T)(v * ((T)~(T)0/255)) >> ((sizeof(v) - 1) * 8);
}

static const u16 LOG_TABLE[256] = {
	512, 255, 1, 122, 2, 244, 123, 181, 3, 48, 245, 224, 124, 84, 182, 111,
	4, 233, 49, 19, 246, 107, 225, 206, 125, 56, 85, 170, 183, 91, 112, 250,
	5, 117, 234, 10, 50, 156, 20, 213, 247, 203, 108, 178, 226, 37, 207, 210,
	126, 150, 57, 100, 86, 141, 171, 40, 184, 73, 92, 164, 113, 146, 251, 229,
	6, 96, 118, 15, 235, 193, 11, 13, 51, 68, 157, 195, 21, 31, 214, 237,
	248, 168, 204, 17, 109, 222, 179, 120, 227, 162, 38, 98, 208, 176, 211, 8,
	127, 188, 151, 239, 58, 132, 101, 216, 87, 80, 142, 33, 172, 27, 41, 23,
	185, 77, 74, 197, 93, 65, 165, 159, 114, 200, 147, 70, 252, 45, 230, 53,
	7, 175, 97, 161, 119, 221, 16, 167, 236, 30, 194, 67, 12, 192, 14, 95,
	52, 44, 69, 199, 158, 64, 196, 76, 22, 26, 32, 79, 215, 131, 238, 187,
	249, 90, 169, 55, 205, 106, 18, 232, 110, 83, 223, 47, 180, 243, 121, 254,
	228, 145, 163, 72, 39, 140, 99, 149, 209, 36, 177, 202, 212, 155, 9, 116,
	128, 61, 189, 218, 152, 137, 240, 103, 59, 135, 133, 134, 102, 136, 217, 60,
	88, 104, 81, 241, 143, 138, 34, 153, 173, 219, 28, 190, 42, 62, 24, 129,
	186, 130, 78, 25, 75, 63, 198, 43, 94, 191, 66, 29, 166, 220, 160, 174,
	115, 154, 201, 35, 148, 139, 71, 144, 253, 242, 46, 82, 231, 105, 54, 89
};

static const u8 EXP_TABLE[512*2+1] = {
	1, 2, 4, 8, 16, 32, 64, 128, 95, 190, 35, 70, 140, 71, 142, 67,
	134, 83, 166, 19, 38, 76, 152, 111, 222, 227, 153, 109, 218, 235, 137, 77,
	154, 107, 214, 243, 185, 45, 90, 180, 55, 110, 220, 231, 145, 125, 250, 171,
	9, 18, 36, 72, 144, 127, 254, 163, 25, 50, 100, 200, 207, 193, 221, 229,
	149, 117, 234, 139, 73, 146, 123, 246, 179, 57, 114, 228, 151, 113, 226, 155,
	105, 210, 251, 169, 13, 26, 52, 104, 208, 255, 161, 29, 58, 116, 232, 143,
	65, 130, 91, 182, 51, 102, 204, 199, 209, 253, 165, 21, 42, 84, 168, 15,
	30, 60, 120, 240, 191, 33, 66, 132, 87, 174, 3, 6, 12, 24, 48, 96,
	192, 223, 225, 157, 101, 202, 203, 201, 205, 197, 213, 245, 181, 53, 106, 212,
	247, 177, 61, 122, 244, 183, 49, 98, 196, 215, 241, 189, 37, 74, 148, 119,
	238, 131, 89, 178, 59, 118, 236, 135, 81, 162, 27, 54, 108, 216, 239, 129,
	93, 186, 43, 86, 172, 7, 14, 28, 56, 112, 224, 159, 97, 194, 219, 233,
	141, 69, 138, 75, 150, 115, 230, 147, 121, 242, 187, 41, 82, 164, 23, 46,
	92, 184, 47, 94, 188, 39, 78, 156, 103, 206, 195, 217, 237, 133, 85, 170,
	11, 22, 44, 88, 176, 63, 126, 252, 167, 17, 34, 68, 136, 79, 158, 99,
	198, 211, 249, 173, 5, 10, 20, 40, 80, 160, 31, 62, 124, 248, 175, 1,
	2, 4, 8, 16, 32, 64, 128, 95, 190, 35, 70, 140, 71, 142, 67, 134,
	83, 166, 19, 38, 76, 152, 111, 222, 227, 153, 109, 218, 235, 137, 77, 154,
	107, 214, 243, 185, 45, 90, 180, 55, 110, 220, 231, 145, 125, 250, 171, 9,
	18, 36, 72, 144, 127, 254, 163, 25, 50, 100, 200, 207, 193, 221, 229, 149,
	117, 234, 139, 73, 146, 123, 246, 179, 57, 114, 228, 151, 113, 226, 155, 105,
	210, 251, 169, 13, 26, 52, 104, 208, 255, 161, 29, 58, 116, 232, 143, 65,
	130, 91, 182, 51, 102, 204, 199, 209, 253, 165, 21, 42, 84, 168, 15, 30,
	60, 120, 240, 191, 33, 66, 132, 87, 174, 3, 6, 12, 24, 48, 96, 192,
	223, 225, 157, 101, 202, 203, 201, 205, 197, 213, 245, 181, 53, 106, 212, 247,
	177, 61, 122, 244, 183, 49, 98, 196, 215, 241, 189, 37, 74, 148, 119, 238,
	131, 89, 178, 59, 118, 236, 135, 81, 162, 27, 54, 108, 216, 239, 129, 93,
	186, 43, 86, 172, 7, 14, 28, 56, 112, 224, 159, 97, 194, 219, 233, 141,
	69, 138, 75, 150, 115, 230, 147, 121, 242, 187, 41, 82, 164, 23, 46, 92,
	184, 47, 94, 188, 39, 78, 156, 103, 206, 195, 217, 237, 133, 85, 170, 11,
	22, 44, 88, 176, 63, 126, 252, 167, 17, 34, 68, 136, 79, 158, 99, 198,
	211, 249, 173, 5, 10, 20, 40, 80, 160, 31, 62, 124, 248, 175, 1, 0,
};

static const u8 INV_TABLE[256] = {
	0, 1, 175, 202, 248, 70, 101, 114, 124, 46, 35, 77, 157, 54, 57, 247,
	62, 152, 23, 136, 190, 244, 137, 18, 225, 147, 27, 26, 179, 59, 212, 32,
	31, 213, 76, 10, 164, 182, 68, 220, 95, 144, 122, 113, 235, 195, 9, 125,
	223, 253, 230, 189, 162, 120, 13, 156, 246, 14, 178, 29, 106, 84, 16, 153,
	160, 119, 197, 198, 38, 221, 5, 249, 82, 159, 91, 207, 34, 11, 110, 166,
	128, 104, 72, 158, 61, 107, 151, 201, 218, 116, 206, 74, 171, 155, 145, 40,
	192, 139, 209, 134, 115, 6, 241, 180, 81, 129, 60, 85, 169, 176, 78, 167,
	123, 43, 7, 100, 89, 219, 161, 65, 53, 163, 42, 112, 8, 47, 227, 187,
	80, 105, 148, 232, 205, 214, 99, 208, 19, 22, 193, 97, 173, 229, 211, 238,
	41, 94, 224, 25, 130, 233, 200, 86, 17, 63, 170, 93, 55, 12, 83, 73,
	64, 118, 52, 121, 36, 183, 79, 111, 177, 108, 154, 92, 228, 140, 203, 2,
	109, 168, 58, 28, 103, 240, 37, 165, 250, 217, 226, 127, 231, 51, 20, 245,
	96, 138, 234, 45, 199, 66, 67, 196, 150, 87, 3, 174, 215, 132, 90, 75,
	135, 98, 239, 142, 30, 33, 133, 204, 251, 185, 88, 117, 39, 69, 252, 48,
	146, 24, 186, 126, 172, 141, 50, 188, 131, 149, 194, 44, 255, 243, 143, 210,
	181, 102, 254, 237, 21, 191, 56, 15, 4, 71, 184, 216, 222, 49, 242, 236
};

static CAT_INLINE u8 GF256Multiply(u8 x, u8 y)
{
	return EXP_TABLE[LOG_TABLE[x] + LOG_TABLE[y]];
}

static CAT_INLINE u8 GF256Divide(u8 x, u8 y)
{
	// Precondition: y != 0
	return EXP_TABLE[LOG_TABLE[x] + 255 - LOG_TABLE[y]];
}


void FindGoodHeavySeeds()
{
	CatsChoice prng, prng2;

	u8 heavy_matrix[6*18];

	for (u32 seed = 0xdeadbeef;;++seed)
	{
		prng.Initialize(seed);
		prng2.Initialize(seed+1);

		u8 temp_row[18];

		int fail_count = 0;

		// For each row,
		for (int ii = 0; ii < 6; ++ii)
		{
			// Fill identity columns
			for (int jj = 12; jj < 18; ++jj)
			{
				temp_row[jj] = (ii == (jj - 12)) ? 1 : 0;
			}

			bool fail;
			do 
			{
				fail = false;

				// Fill columns
				for (int jj = 0; jj < 12; ++jj)
				{
					temp_row[jj] = (u8)(prng.Next() ^ prng2.Next());
				}

				for (int jj = 0; jj < 18; ++jj)
				{
					u8 target = temp_row[jj];

					if (target == 0)
					{
						if (jj < 12)
						{
							//cout << "--Zero failure" << endl;
							fail = true;
							break;
						}

						continue;
					}

					// Generate remaining symbols
					u32 win_lim = (1 << jj);
					u32 prev_gray = 0;
					u8 state = 0;
					for (u32 kk = 1; kk < win_lim; ++kk)
					{
						u32 gray = kk ^ (kk >> 1);
						u32 diff = gray ^ prev_gray;

						u32 index;
						_BitScanForward((unsigned long*)&index, diff);

						state ^= temp_row[index];

						if (state == target)
						{
							int bci = BitCount(gray);

							if (bci < 4)
							{
								//cout << "--Derivative failure at " << jj << " hw " << bci << endl;
								jj = 18;
								fail = true;
								break;
							}
						}

						prev_gray = gray;
					}
				}

				// If row value matches an existing one,
				for (int ll = 0; ll < ii; ++ll)
				{
					for (int mm = 0; mm < 12; ++mm)
					{
						if (heavy_matrix[ll*18 + mm] == temp_row[mm])
						{
							fail = true;
							break;
						}
					}
				}

			} while (fail);

			memcpy(&heavy_matrix[ii*18], temp_row, 18);
			/*
			cout << "Found row " << ii << endl;

			for (int ll = 0; ll < 18; ++ll)
			{
				cout <<hex << setw(2) << setfill('0') << (int)temp_row[ll] << dec << " ";
			}
			cout << endl;
			*/
			// If row value matches an existing one,
			for (int ll = 0; ll < ii; ++ll)
			{
				for (int mm = 0; mm < 12; ++mm)
				{
					u8 mul_val = GF256Divide(temp_row[mm], heavy_matrix[ll*18 + mm]);

					for (int nn = mm + 1; nn < 18; ++nn)
					{
						if (temp_row[nn] == GF256Multiply(heavy_matrix[ll*18 + nn], mul_val))
						{
							++fail_count;
						}
					}
				}
			}
		}

		cout << "Completed matrix!  Had " << fail_count << " failures" << endl;
	}
}

void FindBadDenseSeeds()
{
	int block_bytes = 1;
	int max_blocks = 64000;
	int max_message_bytes = block_bytes * max_blocks;
	u8 *message = new u8[max_message_bytes];
	u8 *message_out = new u8[max_message_bytes];
	u8 *block = new u8[block_bytes];

	for (int ii = 0; ii < max_message_bytes; ++ii)
	{
		message[ii] = ii;
	}

	ofstream file("exception_list.txt");
	if (!file) return;

	file << "struct SeedException { u16 n, seed; };" << endl;
	file << "static const SeedException EXCEPTIONS[] = {" << endl;

	int fails = 0;
	int seen = 0;
	for (int ii = 2; ii <= 64000; ++ii)
	{
		++seen;
		int block_count = ii;
		int message_bytes = block_bytes * block_count;

		wirehair::Encoder encoder;

		wirehair::Result r = encoder.BeginEncode(message, message_bytes, block_bytes);

		if (r == wirehair::R_BAD_PEEL_SEED)
		{
			++fails;
			cout << "-- FAIL! N=" << encoder.BlockCount() << " encoder.BeginEncode error " << wirehair::GetResultString(r) << " at " << fails/(double)seen << endl;

			file << " {" << ii << "," << 0 << "},";

			if (fails % 16 == 0) file << endl;

			//cin.get();
		}

		if (ii % 1000 == 0) cout << ii << endl;
	}

	file << "};" << endl;

	file << "static const int EXCEPTION_COUNT = " << fails << ";" << endl;

	cin.get();
}

int main()
{
	m_clock.OnInitialize();

	//FindGoodHeavySeeds();
	//FindBadDenseSeeds();

	for (int ii = 64; ii <= 64000; ii += 1000)
	{
		int block_count = ii;
		int block_bytes = 1500;
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
			double mbytes = message_bytes / 1000000.;

			cout << ">> OKAY! N=" << encoder.BlockCount() << "(" << mbytes << " MB) encoder.BeginEncode in " << end - start << " usec, " << message_bytes / (end - start) << " MB/s" << endl;
			//cin.get();
		}

		CatsChoice prng;
		cat::wirehair::Decoder decoder;

		u32 overhead_sum = 0, overhead_trials = 0;
		u32 drop_seed = 50000;
		double time_sum = 0;
		const int trials = 1000;
		for (int jj = 0; jj < trials; ++jj)
		{
			int blocks_needed = 0;

			wirehair::Result s = decoder.BeginDecode(message_out, message_bytes, block_bytes);
			if (s)
			{
				cout << "-- FAIL! N=" << decoder.BlockCount() << " decoder.BeginDecode error " << wirehair::GetResultString(s) << endl;
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

		double avg_time = time_sum / trials;
		double avg_overhead = overhead_sum / (double)overhead_trials;
		double avg_bytes = message_bytes * (decoder.BlockCount() + avg_overhead) / (double)decoder.BlockCount() - message_bytes;
		cout << "N=" << decoder.BlockCount() << " decoder.Decode in " << avg_time << " usec, " << message_bytes / avg_time << " MB/s.  Average overhead = " << avg_overhead << " (" << avg_bytes << " bytes)" << endl;
	}

	m_clock.OnFinalize();

	return 0;
}
