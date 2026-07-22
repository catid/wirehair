#pragma once

#include "../gf256.h"
#include "WirehairV2Precode.h"

#include <algorithm>
#include <cstring>
#include <limits>
#include <memory>
#include <vector>

namespace wirehair_v2 {

// Initializes a block from its first source without a separate memset/read
// pass.  Later sources are XORed in batches.
class MixedBatchedBlockXorInitializer
{
public:
    MixedBatchedBlockXorInitializer(
        uint8_t* destination,
        uint32_t block_bytes,
        const uint8_t* first_source,
        bool destination_initially_zero = false)
        : Destination(destination)
        , BlockBytes(block_bytes)
        , DestinationInitiallyZero(destination_initially_zero)
        , BatchCapacity(kFusedBatchSize)
    {
        if (first_source) {
            PendingSources[PendingCount++] = first_source;
        }
    }

    void Add(const uint8_t* source)
    {
        PendingSources[PendingCount++] = source;
        if (PendingCount == BatchCapacity) {
            Flush();
        }
    }

    void Flush()
    {
        if (!Initialized)
        {
            if (PendingCount == 0u) {
                if (!DestinationInitiallyZero) {
                    std::memset(Destination, 0, BlockBytes);
                }
            }
            else {
                gf256_addset_multi_mem(
                    Destination, PendingSources, (int)PendingCount,
                    (int)BlockBytes);
            }
            Initialized = true;
            BatchCapacity = kRegularBatchSize;
        }
        else if (PendingCount != 0u)
        {
            gf256_add_multi_mem(
                Destination, PendingSources, (int)PendingCount,
                (int)BlockBytes);
        }
        PendingCount = 0u;
    }

private:
    static const uint32_t kRegularBatchSize = 8u;
    static const uint32_t kFusedBatchSize = 16u;
    uint8_t* Destination;
    uint32_t BlockBytes;
    bool DestinationInitiallyZero;
    const void* PendingSources[kFusedBatchSize];
    uint32_t PendingCount = 0u;
    uint32_t BatchCapacity;
    bool Initialized = false;
};

struct MixedJointResidueBuckets
{
    std::unique_ptr<uint8_t[]> Subfield;
    std::unique_ptr<uint8_t[]> Extension;
    uint64_t SourceXors = 0u;
    uint64_t MarginalXors = 0u;
    uint64_t MarginalCopies = 0u;
    // Bytes in the three P*block_bytes data planes.  This intentionally does
    // not claim allocator-exact accounting for the small scheduling vectors.
    uint64_t ScratchBytes = 0u;
    uint32_t ActiveDeltas = 0u;
};

// Accumulate the two independently rotated P-bucket marginals through one
// temporary P-bucket plane per active secondary-minus-A block-shift delta.
// The secondary schedule is B for the independent-extension path and C for
// grouped GF(256) rows.  `source_at`
// returns the value of one column and `is_active` excludes decoder-inactive
// columns.  Complete and partial final P-column blocks share the same path.
template<class ShiftAAt, class ShiftBAt, class SourceAt, class IsActive>
bool AccumulateMixedJointResidueBucketsWithShifts(
    uint32_t column_count,
    uint32_t coefficient_period,
    uint32_t block_bytes,
    ShiftAAt shift_a_at,
    ShiftBAt shift_b_at,
    SourceAt source_at,
    IsActive is_active,
    bool batch_sources,
    MixedJointResidueBuckets& output)
{
    output = MixedJointResidueBuckets{};
    if (coefficient_period == 0u ||
        coefficient_period > kMixedCoefficientPeriod ||
        block_bytes == 0u || block_bytes > 0x7fffffffu)
    {
        return false;
    }
    const uint64_t plane_bytes_wide =
        (uint64_t)coefficient_period * block_bytes;
    if (plane_bytes_wide > std::numeric_limits<size_t>::max() ||
        plane_bytes_wide > UINT64_MAX / 3u)
    {
        return false;
    }
    const size_t plane_bytes = (size_t)plane_bytes_wide;
    const uint32_t block_count = column_count / coefficient_period +
        (column_count % coefficient_period != 0u ? 1u : 0u);

    std::vector<uint32_t> delta_offsets(coefficient_period + 1u, 0u);
    std::vector<uint32_t> block_a_shifts(block_count);
    std::vector<uint32_t> block_deltas(block_count);
    for (uint32_t block = 0u; block < block_count; ++block)
    {
        const uint32_t a_shift = shift_a_at(block);
        const uint32_t b_shift = shift_b_at(block);
        if (a_shift >= coefficient_period || b_shift >= coefficient_period) {
            return false;
        }
        const uint32_t delta = b_shift >= a_shift ?
            b_shift - a_shift : b_shift + coefficient_period - a_shift;
        block_a_shifts[block] = a_shift;
        block_deltas[block] = delta;
        ++delta_offsets[delta + 1u];
    }
    for (uint32_t delta = 0u; delta < coefficient_period; ++delta) {
        delta_offsets[delta + 1u] += delta_offsets[delta];
    }
    std::vector<uint32_t> next = delta_offsets;
    std::vector<uint32_t> blocks_by_delta(block_count);
    for (uint32_t block = 0u; block < block_count; ++block)
    {
        const uint32_t delta = block_deltas[block];
        blocks_by_delta[next[delta]++] = block;
    }

    output.Subfield.reset(new uint8_t[plane_bytes]);
    output.Extension.reset(new uint8_t[plane_bytes]);
    std::unique_ptr<uint8_t[]> temporary(new uint8_t[plane_bytes]);
    output.ScratchBytes = plane_bytes_wide * 3u;
    if (block_count == 0u)
    {
        std::memset(output.Subfield.get(), 0, plane_bytes);
        std::memset(output.Extension.get(), 0, plane_bytes);
        return true;
    }
    bool marginals_initialized = false;
    // The active-delta loop reconstructs the initializer objects because
    // their pending counts and initialization state are delta-local.  Keep
    // the vector storage itself across iterations: P is small but a full
    // cycle can activate many deltas, and repeatedly allocating the same P
    // objects showed up in grouped-schedule solve profiles.
    std::vector<MixedBatchedBlockXorInitializer> accumulators;
    std::vector<uint8_t> initialized;
    if (batch_sources) {
        accumulators.reserve(coefficient_period);
    }
    else {
        initialized.reserve(coefficient_period);
    }

    for (uint32_t delta = 0u; delta < coefficient_period; ++delta)
    {
        const uint32_t begin = delta_offsets[delta];
        const uint32_t end = delta_offsets[delta + 1u];
        if (begin == end) continue;
        ++output.ActiveDeltas;
        if (batch_sources)
        {
            accumulators.clear();
            for (uint32_t a = 0u; a < coefficient_period; ++a)
            {
                accumulators.emplace_back(
                    temporary.get() + (size_t)a * block_bytes,
                    block_bytes, nullptr);
            }
        }
        else {
            initialized.assign(coefficient_period, uint8_t{0});
        }

        // Discover source columns in sequential P-column block order.  Each
        // destination batches the discovered pointers so the SIMD multi-XOR
        // kernel initializes and extends its temporary bucket in one pass.
        for (uint32_t position = begin; position < end; ++position)
        {
            const uint32_t block = blocks_by_delta[position];
            const uint32_t block_base = block * coefficient_period;
            const uint32_t block_columns = std::min(
                coefficient_period, column_count - block_base);
            const uint32_t a_shift = block_a_shifts[block];
            for (uint32_t offset = 0u; offset < block_columns; ++offset)
            {
                const uint32_t column = block_base + offset;
                if (!is_active(column)) continue;
                uint32_t a = offset + a_shift;
                if (a >= coefficient_period) a -= coefficient_period;
                const uint8_t* source = source_at(column);
                if (!source) return false;
                if (batch_sources) {
                    accumulators[a].Add(source);
                }
                else
                {
                    uint8_t* destination =
                        temporary.get() + (size_t)a * block_bytes;
                    if (!initialized[a])
                    {
                        MixedBatchedBlockXorInitializer initialize(
                            destination, block_bytes, source);
                        initialize.Flush();
                        initialized[a] = 1u;
                    }
                    else {
                        gf256_add_mem(
                            destination, source, (int)block_bytes);
                    }
                }
                ++output.SourceXors;
            }
        }
        if (batch_sources)
        {
            for (uint32_t a = 0u; a < coefficient_period; ++a) {
                accumulators[a].Flush();
            }
        }
        else
        {
            for (uint32_t a = 0u; a < coefficient_period; ++a) {
                if (!initialized[a]) {
                    std::memset(
                        temporary.get() + (size_t)a * block_bytes,
                        0, block_bytes);
                }
            }
        }

        if (!marginals_initialized)
        {
            // The first marginal is a write-only initialization: A retains
            // temp order and B is the same plane rotated forward by delta.
            // Coalesce those P bucket copies into at most three memcpy calls.
            std::memcpy(
                output.Subfield.get(), temporary.get(), plane_bytes);
            const size_t extension_prefix =
                (size_t)delta * block_bytes;
            const size_t extension_suffix = plane_bytes - extension_prefix;
            std::memcpy(
                output.Extension.get() + extension_prefix,
                temporary.get(), extension_suffix);
            if (extension_prefix != 0u)
            {
                std::memcpy(
                    output.Extension.get(),
                    temporary.get() + extension_suffix,
                    extension_prefix);
            }
            output.MarginalCopies = (uint64_t)2u * coefficient_period;
        }
        else
        {
            for (uint32_t a = 0u; a < coefficient_period; ++a)
            {
                uint8_t* temporary_bucket =
                    temporary.get() + (size_t)a * block_bytes;
                gf256_add_mem(
                    output.Subfield.get() + (size_t)a * block_bytes,
                    temporary_bucket, (int)block_bytes);
                uint32_t b = a + delta;
                if (b >= coefficient_period) b -= coefficient_period;
                gf256_add_mem(
                    output.Extension.get() + (size_t)b * block_bytes,
                    temporary_bucket, (int)block_bytes);
                output.MarginalXors += 2u;
            }
        }
        marginals_initialized = true;
    }
    return true;
}

// Production A/B wrapper.  Keeping the original interface centralizes the
// experiment-only choice of a different second schedule in its callers and
// leaves every existing call site and no-hooks instantiation unchanged.
template<class SourceAt, class IsActive>
bool AccumulateMixedJointResidueBuckets(
    uint32_t column_count,
    uint32_t coefficient_period,
    uint32_t block_bytes,
    SourceAt source_at,
    IsActive is_active,
    bool batch_sources,
    MixedJointResidueBuckets& output)
{
    return AccumulateMixedJointResidueBucketsWithShifts(
        column_count,
        coefficient_period,
        block_bytes,
        [](uint32_t block) {
            return ActiveMixedResidueBlockShift(block);
        },
        [](uint32_t block) {
            return ActiveMixedExtensionResidueBlockShift(block);
        },
        source_at,
        is_active,
        batch_sources,
        output);
}

} // namespace wirehair_v2
