#include <iostream>
using namespace std;

#include "Platform.hpp"
#include "BitMath.hpp"
#include "Galois256.hpp"
using namespace cat;

// Calculate number of 1s in Cauchy 8x8 submatrix representation
static int cauchy_ones(u8 n) {
	/*
	 * w = 8 so the Cauchy representation is an 8x8 submatrix
	 * in place of the GF(256) values of the matrix.
	 *
	 * To generate the 8x8 submatrix, the first column is the
	 * original value in binary.  And the remaining columns
	 * are the column to the left times 2 in GF(256).
	 */

	// Count number of ones in the first column
	int ones = BIT_COUNT_TABLE[n];

	// For each remaining column,
	for (int w = 1; w < 8; ++w) {
		// Generate the column
		// NOTE: Optimizes to a single 256-byte table lookup
		n = GF256Multiply(n, 2);

		// Count the bits in it
		ones += BIT_COUNT_TABLE[n];
	}

	return ones;
}

int main() {
	GF256Init();

	cout << "static const u8 CAUCHY_ONES[256] = {" << endl;
	for (int x = 0; x < 256; ++x) {
		if ((x & 15) == 0) {
			cout << endl;
		}

		int ones = cauchy_ones(x);

		cout << ones << ", ";
	}
	cout << endl << "};" << endl;

	return 0;
}

