#include <iostream>
#include <iomanip>
using namespace std;

#include "Platform.hpp"
#include "BitMath.hpp"
#include "AbyssinianPRNG.hpp"
#include "Clock.hpp"
using namespace cat;

static Clock m_clock;

static const u8 GEN_POLY[] = {
	0x8e, 0x95, 0x96, 0xa6, 0xaf, 0xb1, 0xb2, 0xb4,
	0xb8, 0xc3, 0xc6, 0xd4, 0xe1, 0xe7, 0xf3, 0xfa,
};

// This one has average weight of 16.5 ones for the first 20 min-weight elements of GF(256)
// The one used by Jerasure has 17.55 average weight.
static const int FAVORITE_POLY = 9; // 0xC3 => 111000011

static u16 GF256_LOG_TABLE[256];
static u8 GF256_EXP_TABLE[512*2+1];
static u8 GF256_INV_TABLE[256];

u8 * CAT_RESTRICT GF256_MUL_TABLE = 0;
u8 * CAT_RESTRICT GF256_DIV_TABLE = 0;

static u8 CAUCHY_ONES[256];

static u8 MINWEIGHT_TABLE[256];

// return x * y in GF(256)
// For repeated multiplication by a constant, it is faster to put the constant in y.
static CAT_INLINE u8 GF256Multiply(u8 x, u8 y) {
	return GF256_MUL_TABLE[((u32)y << 8) + x];
}

// return x / y in GF(256)
// Memory-access optimized for constant divisors in y.
static CAT_INLINE u8 GF256Divide(u8 x, u8 y) {
	return GF256_DIV_TABLE[((u32)y << 8) + x];
}

void InitInvTable() {
	for (int x = 0; x < 256; ++x) {
		GF256_INV_TABLE[x] = GF256Divide(1, x);
	}
}

// Unpack 256x256 multiplication tables
void InitMulDivTables() {
	// If not initialized already,
	if (!GF256_MUL_TABLE) {
		GF256_MUL_TABLE = (u8 *)malloc(256 * 256 * 2);
		GF256_DIV_TABLE = GF256_MUL_TABLE + 256 * 256;
	}

	// Allocate table memory 65KB x 2
	u8 *m = GF256_MUL_TABLE, *d = GF256_DIV_TABLE;

	// Unroll y = 0 subtable
	for (int x = 0; x < 256; ++x) {
		m[x] = d[x] = 0;
	}

	// For each other y value,
	for (int y = 1; y < 256; ++y) {
		// Calculate log(y) for mult and 255 - log(y) for div
		const u8 log_y = GF256_LOG_TABLE[y];
		const u8 log_yn = 255 - log_y;

		// Next subtable
		m += 256;
		d += 256;

		// Unroll x = 0
		m[0] = 0;
		d[0] = 0;

		// Calculate x * y, x / y
		for (int x = 1; x < 256; ++x) {
			int log_x = GF256_LOG_TABLE[x];

			m[x] = GF256_EXP_TABLE[log_x + log_y];
			d[x] = GF256_EXP_TABLE[log_x + log_yn];
		}
	}
}

void GenerateExpLogTables(int ii, u16 LogTable[256], u8 ALogTable[512*2+1])
{
	u32 poly = (GEN_POLY[ii] << 1) | 1;

	LogTable[0] = 512;
	ALogTable[0] = 1;
	for (u32 jj = 1; jj < 255; ++jj)
	{
		u32 next = (u32)ALogTable[jj - 1] * 2;
		if (next >= 256) next ^= poly;

		ALogTable[jj] = next;
		LogTable[ALogTable[jj]] = jj;
	}

	ALogTable[255] = ALogTable[0];
	LogTable[ALogTable[255]] = 255;

	for (u32 jj = 256; jj < 2 * 255; ++jj)
	{
		ALogTable[jj] = ALogTable[jj % 255];
	}

	ALogTable[2 * 255] = 1;

	for (u32 jj = 2 * 255 + 1; jj < 4 * 255; ++jj)
	{
		ALogTable[jj] = 0;
	}
}

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

void GenerateCauchyOnes() {
	for (int x = 0; x < 256; ++x) {
		CAUCHY_ONES[x] = cauchy_ones(x);
	}
}

void print(int k, int m, u8 *matrix, bool cstyle = false) {
	if (cstyle) {
		cout << "static const u8 CAUCHY_MATRIX_" << m << "[" << (m - 1) << " * " << k << "] = {" << endl;
		for (int y = m > 1 ? 1 : 0; y < m; ++y) {
			if (y > 1) {
				cout << "// For row " << y << ":" << endl;
			}
			for (int x = 0; x < k; ++x) {
				cout << dec << (int)matrix[y * k + x];
				if (x > 0 && (x % 20) == 19) {
					cout << "," << endl;
				} else if (x < k - 1) {
					cout << ",";
				}
			}
			if (y < m - 1) {
				cout << "," << endl;
			}
		}
		cout << dec << "};" << endl;
	} else {
		cout << "[" << endl;
		for (int y = 0; y < m; ++y) {
			for (int x = 0; x < k; ++x) {
				cout << hex << setw(2) << setfill('0') << (int)matrix[y * k + x] << " ";
			}
			cout << endl;
		}
		cout << dec << "]" << endl;
	}
}

void SortColumns(int k, int m, u8 *matrix) {
	int *counts = new int[k];

	for (int x = 0; x < k; ++x) {
		int ones = 0;
		for (int y = 0; y < m; ++y) {
			ones += CAUCHY_ONES[matrix[y*k + x]];
		}
		counts[x] = ones;
	}

	for (int x = 0; x < k; ++x) {
		int smallest = counts[x], best_x = x;

		for (int z = x + 1; z < k; ++z) {
			int ones = counts[z];

			if (counts[z] < smallest) {
				smallest = ones;
				best_x = z;
			}
		}

		// swap counts
		counts[best_x] = counts[x];
		counts[x] = smallest;

		// swap columns
		for (int y = 1; y < m; ++y) {
			u8 temp = matrix[y*k + x];
			matrix[y*k + x] = matrix[y*k + best_x];
			matrix[y*k + best_x] = temp;
		}
	}

	delete []counts;
}

/*
 * Cauchy matrices are defined by two vectors X, Y s.t. X, Y share no elements
 * in common from the set GF(256).
 *
 * Each element i,j of the Cauchy matrix is 1/(Xi + Yj).
 *
 * Another useful property is that you can multiply each row or column of a
 * Cauchy matrix and it will still be invertible.
 *
 * Since the number of Cauchy ones is far better for 1 than other elements of
 * GF(256), all of the best options will have a 1 in at least one row of each
 * column.  And one column can be all 1s.
 */

int GenerateCauchyMatrix(int k, int m, u8 *matrix, u8 *X, u8 *Y) {
	int ones = CAUCHY_ONES[1] * k;

	u8 *row = matrix;

	for (int y = 0; y < k; ++y) {
		row[y] = 1;
	}

	for (int y = 1; y < m; ++y) {
		u8 YC = Y[y];

		row += k;

		for (int x = 0; x < k; ++x) {
			u8 XC = X[x];

			u8 D = Y[0] ^ XC;
			u8 C = GF256Multiply(GF256_INV_TABLE[XC ^ YC], D);

			row[x] = C;
		}
	}

	row = matrix;
	for (int y = 1; y < m; ++y) {
		row += k;

		int best = 0x7fffffff, best_x = 0;

		for (int x = 0; x < k; ++x) {
			u8 XC = row[x];

			int count = 0;

			for (int z = 0; z < k; ++z) {
				u8 C = GF256Divide(row[z], XC);

				count += CAUCHY_ONES[C];
			}

			if (count < best) {
				best = count;
				best_x = x;
			}
		}

		u8 XC = row[best_x];

		for (int z = 0; z < k; ++z) {
			u8 C = GF256Divide(row[z], XC);

			row[z] = C;

			ones += CAUCHY_ONES[C];
		}
	}

	return ones;
}

int ImproveMatrixRows(int k, int subk, int m, u8 *matrix) {
	for (int y = 1; y < m; ++y) {
		int best = 0x7fffffff, best_A;
		for (int x = 0; x < k; ++x) {
			u8 A = matrix[y*k + x];
			u8 IA = GF256_INV_TABLE[A];

			int ones = 0;
			for (int z = 0; z < k; ++z) {
				u8 B = matrix[y*k + z];
				u8 M = GF256Multiply(B, IA);
				ones += CAUCHY_ONES[M];
			}
			if (ones < best) {
				best = ones;
				best_A = IA;
			}
		}
		for (int z = 0; z < k; ++z) {
			u8 B = matrix[y*k + z];
			u8 M = GF256Multiply(B, best_A);
			matrix[y*k + z] = M;
		}
	}

	SortColumns(k, m, matrix);

	int total = 0;
	for (int y = 0; y < m; ++y) {
		for (int x = 0; x < subk; ++x) {
			total += CAUCHY_ONES[matrix[y*k + x]];
		}
	}
	return total;
}

void SolveBestMatrix(int m, int subk) {
	//   A B C D E
	// F 1 1 1 1 1
	// G a b c d e
	// H f g h i j

	const int k = 256 - m;
	u8 *matrix = new u8[k * m];
	u8 *best_matrix = new u8[k * m];
	int best_matrix_ones = 0x7fffffff;
	u8 X[256], Y[256];

	double t0 = m_clock.usec();

	// First row is always all ones
	for (int x = 0; x < k; ++x) {
		matrix[x] = 1;
	}

	// Choose a seed of A,F and solve the rest with a greedy algorithm
	//for (int F = 0; F < 256; ++F) {
		//for (int A = 0; A < 256; ++A) {
			//if (A != F) {
	int F = 0, A = 1; {{{
				u8 seen[256];
				for (int ii = 0; ii < 256; ++ii) {
					seen[ii] = 0;
				}
				seen[A] = 1;
				seen[F] = 1;
				X[0] = A;
				Y[0] = F;
				u8 AF = A ^ F;

				// In the order of greedy solution solve the Y values first.

				int trial_ones = 8 * subk;

				// Unroll the first column
				for (int y = 1; y < m; ++y) {
					// Pick the next worst element of GF(256) in weight
					for (int ii = 1; ii < 256; ++ii) {
						u8 a = MINWEIGHT_TABLE[ii];

						// a * (A + G) = A + F
						// aA + aG = A + F
						// G = (A + F + aA) / a
						u8 G = GF256Divide(AF ^ GF256Multiply(a, A), a);

						if (seen[G]) {
							continue;
						}

						seen[G] = 1;
						Y[y] = G;
						matrix[y * k] = a;
						trial_ones += CAUCHY_ONES[a];
						break;
					}
				}

				// Now solve the X values.
				// TODO: Maybe iterate through this inside Y value solutions

				// For each remaining column,
				for (int x = 1; x < k; ++x) {
					int best_ones = 0x7fffffff, best_B;

					// Verify that a solution is possible for all column values
					for (int B = 0; B < 256; ++B) {
						if (seen[B]) {
							continue;
						}

						int ones = 0;
						for (int y = 1; y < m; ++y) {
							u8 b = GF256Divide(B ^ F, B ^ Y[y]);
							ones += CAUCHY_ONES[b];
						}
						if (ones < best_ones) {
							best_ones = ones;
							best_B = B;
						}
					}

					int B = best_B;
					X[x] = B;
					seen[B] = 1;

					for (int y = 1; y < m; ++y) {
						u8 b = GF256Divide(B ^ F, B ^ Y[y]);
						matrix[y * k + x] = b;
					}
					if (x < subk) {
						trial_ones += best_ones;
					}
				}
				//print(k, m, matrix);

				double t1 = m_clock.usec();

				int improved_ones = ImproveMatrixRows(k, subk, m, matrix);

				double t2 = m_clock.usec();

				if (improved_ones < best_matrix_ones) {
					best_matrix_ones = improved_ones;
					memcpy(best_matrix, matrix, k * m);

					cout << "Pre-improved ones = " << trial_ones << " in " << (t1 - t0) << " usec" << endl;
					cout << "Best ones for first " << subk << " columns = " << best_matrix_ones << " in " << (t2 - t0) << " usec" << endl;
					print(k, m, best_matrix, true);

					return;
				}
			}
		}
	}
}

void SortMinWeightElements(u8 *elements) {
	int *counts = new int[256];

	for (int x = 0; x < 256; ++x) {
		counts[x] = CAUCHY_ONES[x];
	}

	for (int x = 0; x < 256; ++x) {
		int smallest = counts[x], best_x = x;

		for (int z = x + 1; z < 256; ++z) {
			int ones = counts[z];

			if (counts[z] < smallest) {
				smallest = ones;
				best_x = z;
			}
		}

		// swap counts
		counts[best_x] = counts[x];
		counts[x] = smallest;

		// swap columns
		u8 temp = elements[best_x];
		elements[best_x] = elements[x];
		elements[x] = temp;
	}

	delete []counts;
}

static void ShuffleDeck8(Abyssinian &prng, u8 * CAT_RESTRICT deck)
{
	deck[0] = 0;
	const int count = 256;

	for (u32 ii = 1;;)
	{
		u32 jj, rv = prng.Next();

		// 8-bit unroll
		switch (count - ii)
		{
			default:
				jj = (u8)rv % ii;
				deck[ii] = deck[jj];
				deck[jj] = ii;
				++ii;
				jj = (u8)(rv >> 8) % ii;
				deck[ii] = deck[jj];
				deck[jj] = ii;
				++ii;
				jj = (u8)(rv >> 16) % ii;
				deck[ii] = deck[jj];
				deck[jj] = ii;
				++ii;
				jj = (u8)(rv >> 24) % ii;
				deck[ii] = deck[jj];
				deck[jj] = ii;
				++ii;
				break;

			case 3:
				jj = (u8)rv % ii;
				deck[ii] = deck[jj];
				deck[jj] = ii;
				++ii;
			case 2:
				jj = (u8)(rv >> 8) % ii;
				deck[ii] = deck[jj];
				deck[jj] = ii;
				++ii;
			case 1:
				jj = (u8)(rv >> 16) % ii;
				deck[ii] = deck[jj];
				deck[jj] = ii;
			case 0:
				return;
		}
	}
}

void PrintMinWeights() {
	const int k = 256;
	const int m = 2;
	u8 *matrix = new u8[k * m];

	for (int ii = 0; ii < 16; ++ii) {
		cout << "*** For generator " << ii << ":" << endl;

		GenerateExpLogTables(ii, GF256_LOG_TABLE, GF256_EXP_TABLE);

		InitMulDivTables();

		InitInvTable();

		GenerateCauchyOnes();

		for (int x = 0; x < 256; ++x) {
			matrix[x] = 1;
			matrix[x + 256] = x;
		}

		cout << "Symbols in order:" << endl;

		SortColumns(k, m, matrix);

		print(k, m, matrix);

		for (int x = 1; x <= 32; ++x) {
			int ones = 0;
			for (int z = 1; z <= x; ++z) {
				ones += CAUCHY_ONES[matrix[k + z]];
			}
			cout << x << " columns = " << ones << " ones" << endl;
		}
	}
}

void Explore(int k, int m) {
	u8 *matrix = new u8[k * m];
	u8 *best = new u8[k * m];
	u8 *XY = new u8[256];
	u8 *X = XY;
	u8 *Y = XY + k;

	Abyssinian prng;

	prng.Initialize(1);

	int least = 0x7fffffff;

	for (int ii = 0; ii < 1000000000; ++ii) {
		ShuffleDeck8(prng, XY);

		int ones = GenerateCauchyMatrix(k, m, matrix, X, Y);

		if (ones < least) {
			least = ones;
			memcpy(best, matrix, k * m);

			cout << "Found a better matrix with ones = " << ones << ":" << endl;
			print(k, m, best);
		}
	}

	SortColumns(k, m, best);

	cout << "Sorted matrix:" << endl;

	print(k, m, best);

	delete []matrix;
	delete []best;
	delete []XY;
}

void PrintTables() {
	print(512*2+1, 1, GF256_EXP_TABLE, true);
	print(256, 1, GF256_INV_TABLE, true);

	cout << "static const u8 GFC256_LOG[256] = {" << endl;
	for (int x = 0; x < 256; ++x) {
		cout << dec << (int)GF256_LOG_TABLE[x];
		if (x > 0 && (x % 20) == 19) {
			cout << "," << endl;
		} else {
			cout << ",";
		}
	}
	cout << dec << "};" << endl;
}

int main() {
	m_clock.OnInitialize();

	cout << "Exploring options..." << endl;

	GenerateExpLogTables(FAVORITE_POLY, GF256_LOG_TABLE, GF256_EXP_TABLE);

	InitMulDivTables();

	InitInvTable();

	GenerateCauchyOnes();

	for (int ii = 0; ii < 256; ++ii) {
		MINWEIGHT_TABLE[ii] = ii;
	}
	SortMinWeightElements(MINWEIGHT_TABLE);

	print(256, 1, MINWEIGHT_TABLE);

	SolveBestMatrix(6, 29);
	//PrintMinWeights();
	//Explore(29, 3);

	PrintTables();

	m_clock.OnFinalize();

	return 0;
}

