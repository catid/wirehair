#include <iostream>
#include <iomanip>
using namespace std;

#include "Platform.hpp"
using namespace cat;

int main() {
	cout << "4-bit window LUT generator" << endl;

	/*
	 * The purpose of this tool is to generate a lookup table that converts
	 * 4 bits from a row to eliminate to a sequence of pivot row XOR operations.
	 *
	 * This is useful because the bitmatrix is sparse and we don't want to
	 * eliminate the upper triangle fully, which would dirty the pivot rows.
	 *
	 * For example:
	 *
	 * 1010...
	 * 0100...
	 * 0011...
	 * 0001...
	 * 1010... <- row to eliminate
	 *
	 * So take the decimal number "10" from the row to eliminate and convert it
	 * into a bit sequence that represents the XORs to perform:
	 *
	 * 1010 -> LUT -> 1011, since the third pivot row has the last bit in the
	 * window set.
	 */

	u8 sequence[64][16];

	cout << "static const u8 WIN4LUT[64][16] = {" << endl;

	for (int ii = 0; ii < 64; ++ii) {
		/*
		 * The upper triangular bits:
		 *
		 * 1xxx...
		 * 01yy...
		 * 001z...
		 * 0001...
		 *
		 * xxxyyz <- has 64 states
		 */

		cout << "{ ";

		for (int jj = 0; jj < 16; ++jj) {
			int b = jj;

			if (b & 8) {
				b ^= ii >> 3;
			}
			if (b & 4) {
				b ^= (ii >> 1) & 3;
			}
			if (b & 2) {
				b ^= ii & 1;
			}

			sequence[ii][jj] = (u8)b;

			cout << b;
			if (jj < 15) {
				cout << ",";
			}
		}

		cout << " }," << endl;
	}

	cout << "};" << endl;

	return 0;
}

