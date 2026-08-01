#include "V2TinyDenseOracle.h"

#include "../gf256.h"

#include <algorithm>
#include <cstring>
#include <limits>
#include <new>
#include <stdexcept>

namespace wirehair_v2 {
namespace test {
namespace {

bool CheckedAllocation(size_t count, size_t width, size_t& bytes)
{
    if (width != 0u && count >
        kTinyOracleMaxAllocationBytes / width)
    {
        return false;
    }
    bytes = count * width;
    return bytes <= kTinyOracleMaxAllocationBytes;
}

bool IsZero(const uint8_t* data, size_t bytes)
{
    for (size_t i = 0; i < bytes; ++i) {
        if (data[i] != 0u) {
            return false;
        }
    }
    return true;
}

} // namespace

WirehairResult SolvePrecodeSystemTinyDenseOracle(
    const PrecodeSystem& system,
    const PacketRowConfig& config,
    const std::vector<SolvePacket>& packets,
    uint32_t block_bytes,
    std::vector<uint8_t>& intermediate_blocks_out)
{
    const uint32_t K = system.Params.BlockCount;
    const uint32_t S = system.Params.Staircase;
    const uint32_t D2 = system.Params.DenseRows;
    const uint32_t H = system.Params.HeavyRows;
    const uint64_t P_wide = (uint64_t)S + D2 + H;
    if (K < 2u || K > kTinyOracleMaxSourceBlocks ||
        block_bytes == 0u || block_bytes > kTinyOracleMaxBlockBytes ||
        P_wide > UINT32_MAX ||
        !IsPacketRowDomainValid(K, (uint32_t)P_wide, config.MixCount) ||
        !ValidatePrecodeSystem(system))
    {
        return Wirehair_InvalidInput;
    }
    for (const SolvePacket& packet : packets) {
        if (!packet.Data) {
            return Wirehair_InvalidInput;
        }
    }
    if (packets.size() < K) {
        return Wirehair_NeedMore;
    }
    if (gf256_init() != 0) {
        return Wirehair_UnsupportedPlatform;
    }

    const uint32_t P = (uint32_t)P_wide;
    const uint64_t L_wide = (uint64_t)K + P;
    const uint64_t row_count_wide = P_wide + packets.size();
    // These structural caps keep the cubic oracle useful and bounded even if
    // a caller supplies an otherwise-valid but fuzz-hostile custom precode.
    if (L_wide > 512u || row_count_wide > 1024u) {
        return Wirehair_InvalidInput;
    }
    const uint32_t L = (uint32_t)L_wide;
    const uint32_t row_count = (uint32_t)row_count_wide;

    size_t matrix_bytes = 0u, rhs_bytes = 0u, output_bytes = 0u;
    if (!CheckedAllocation(row_count, L, matrix_bytes) ||
        !CheckedAllocation(row_count, block_bytes, rhs_bytes) ||
        !CheckedAllocation(L, block_bytes, output_bytes))
    {
        return Wirehair_InvalidInput;
    }

    try
    {
        std::vector<uint8_t> matrix(matrix_bytes, 0u);
        std::vector<uint8_t> rhs(rhs_bytes, 0u);
        uint32_t row = 0u;
        const auto add_binary = [&](const std::vector<uint32_t>& columns,
                                    const uint8_t* data) -> bool {
            if (row >= row_count) {
                return false;
            }
            uint8_t* coefficients = matrix.data() + (size_t)row * L;
            for (uint32_t column : columns) {
                if (column >= L) {
                    return false;
                }
                coefficients[column] ^= 1u;
            }
            if (data) {
                std::memcpy(
                    rhs.data() + (size_t)row * block_bytes,
                    data,
                    block_bytes);
            }
            ++row;
            return true;
        };

        for (const std::vector<uint32_t>& columns : system.StaircaseRows) {
            if (!add_binary(columns, nullptr)) {
                return Wirehair_InvalidInput;
            }
        }
        for (const std::vector<uint32_t>& columns : system.DenseRowColumns) {
            if (!add_binary(columns, nullptr)) {
                return Wirehair_InvalidInput;
            }
        }
        for (uint32_t heavy = 0; heavy < H; ++heavy)
        {
            uint8_t* coefficients = matrix.data() + (size_t)row * L;
            for (uint32_t column = 0; column < L; ++column) {
                coefficients[column] =
                    HeavyCoefficientForParams(
                        system.Params, heavy, column);
            }
            ++row;
        }
        for (const SolvePacket& packet : packets)
        {
            const std::vector<uint32_t> columns = GeneratePacketMatrixRow(
                K, P, packet.BlockId, config);
            if (columns.empty() || !add_binary(columns, packet.Data)) {
                return Wirehair_InvalidInput;
            }
        }
        if (row != row_count) {
            return Wirehair_Error;
        }

        std::vector<uint32_t> pivot_columns;
        pivot_columns.reserve(L);
        uint32_t rank = 0u;
        for (uint32_t column = 0; column < L && rank < row_count; ++column)
        {
            uint32_t pivot = rank;
            while (pivot < row_count &&
                   matrix[(size_t)pivot * L + column] == 0u)
            {
                ++pivot;
            }
            if (pivot == row_count) {
                continue;
            }
            if (pivot != rank)
            {
                for (uint32_t j = 0; j < L; ++j) {
                    std::swap(
                        matrix[(size_t)rank * L + j],
                        matrix[(size_t)pivot * L + j]);
                }
                for (uint32_t b = 0; b < block_bytes; ++b) {
                    std::swap(
                        rhs[(size_t)rank * block_bytes + b],
                        rhs[(size_t)pivot * block_bytes + b]);
                }
            }

            const uint8_t pivot_value =
                matrix[(size_t)rank * L + column];
            if (pivot_value != 1u)
            {
                const uint8_t inverse = gf256_inv(pivot_value);
                for (uint32_t j = column; j < L; ++j) {
                    matrix[(size_t)rank * L + j] = gf256_mul(
                        matrix[(size_t)rank * L + j], inverse);
                }
                gf256_div_mem(
                    rhs.data() + (size_t)rank * block_bytes,
                    rhs.data() + (size_t)rank * block_bytes,
                    pivot_value,
                    (int)block_bytes);
            }

            for (uint32_t target = 0; target < row_count; ++target)
            {
                if (target == rank) {
                    continue;
                }
                const uint8_t scale =
                    matrix[(size_t)target * L + column];
                if (scale == 0u) {
                    continue;
                }
                for (uint32_t j = column; j < L; ++j) {
                    matrix[(size_t)target * L + j] ^= gf256_mul(
                        scale, matrix[(size_t)rank * L + j]);
                }
                gf256_muladd_mem(
                    rhs.data() + (size_t)target * block_bytes,
                    scale,
                    rhs.data() + (size_t)rank * block_bytes,
                    (int)block_bytes);
            }
            pivot_columns.push_back(column);
            ++rank;
        }

        for (uint32_t r = rank; r < row_count; ++r)
        {
            if (IsZero(matrix.data() + (size_t)r * L, L) &&
                !IsZero(rhs.data() + (size_t)r * block_bytes, block_bytes))
            {
                return Wirehair_Error;
            }
        }
        if (rank < L) {
            return Wirehair_NeedMore;
        }

        std::vector<uint8_t> output(output_bytes, 0u);
        for (uint32_t r = 0; r < rank; ++r) {
            std::memcpy(
                output.data() + (size_t)pivot_columns[r] * block_bytes,
                rhs.data() + (size_t)r * block_bytes,
                block_bytes);
        }
        intermediate_blocks_out.swap(output);
        return Wirehair_Success;
    }
    catch (const std::bad_alloc&) {
        return Wirehair_OOM;
    }
    catch (const std::length_error&) {
        return Wirehair_OOM;
    }
}

} // namespace test
} // namespace wirehair_v2
