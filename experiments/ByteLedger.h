#pragma once

#include <cstddef>
#include <cstdint>

namespace wirehair_experiments {

// LogicalBytes counts destination bytes processed by one logical operation.
// ReadBytes/WriteBytes estimate explicit algorithmic operand traffic.  They do
// not include write-allocate, cache-line rounding, prefetching, or cache hits,
// so they must not be presented as measured DRAM traffic.
struct ByteLedger
{
    uint64_t LogicalBytes = 0;
    uint64_t ReadBytes = 0;
    uint64_t WriteBytes = 0;

    ByteLedger& operator+=(const ByteLedger& other)
    {
        LogicalBytes += other.LogicalBytes;
        ReadBytes += other.ReadBytes;
        WriteBytes += other.WriteBytes;
        return *this;
    }

    uint64_t TrafficBytes() const { return ReadBytes + WriteBytes; }
};

inline bool TryAccumulateByteLedger(
    ByteLedger& total,
    const ByteLedger& increment)
{
    if (increment.LogicalBytes > UINT64_MAX - total.LogicalBytes ||
        increment.ReadBytes > UINT64_MAX - total.ReadBytes ||
        increment.WriteBytes > UINT64_MAX - total.WriteBytes)
    {
        return false;
    }
    const uint64_t next_read = total.ReadBytes + increment.ReadBytes;
    const uint64_t next_write = total.WriteBytes + increment.WriteBytes;
    if (next_read > UINT64_MAX - next_write) {
        return false;
    }
    total += increment;
    return true;
}

enum class ByteOperation
{
    Zero,
    Memcpy,
    Xor,
    Addset,
    Add2,
    Mul,
    Div,
    Muladd,
    FaninFused,
    FaninChained
};

inline ByteLedger LedgerFor(
    ByteOperation operation,
    uint64_t bytes,
    uint64_t fanin = 0)
{
    ByteLedger ledger;
    ledger.LogicalBytes = bytes;
    switch (operation)
    {
    case ByteOperation::Zero:
        ledger.WriteBytes = bytes;
        break;
    case ByteOperation::Memcpy:
    case ByteOperation::Mul:
    case ByteOperation::Div:
        ledger.ReadBytes = bytes;
        ledger.WriteBytes = bytes;
        break;
    case ByteOperation::Xor:
    case ByteOperation::Muladd:
        ledger.ReadBytes = 2 * bytes;
        ledger.WriteBytes = bytes;
        break;
    case ByteOperation::Addset:
        ledger.ReadBytes = 2 * bytes;
        ledger.WriteBytes = bytes;
        break;
    case ByteOperation::Add2:
        ledger.ReadBytes = 3 * bytes;
        ledger.WriteBytes = bytes;
        break;
    case ByteOperation::FaninFused:
        ledger.ReadBytes = (fanin + 1) * bytes;
        ledger.WriteBytes = bytes;
        break;
    case ByteOperation::FaninChained:
        ledger.ReadBytes = 2 * fanin * bytes;
        ledger.WriteBytes = fanin * bytes;
        break;
    }
    return ledger;
}

inline bool VerifyByteLedger()
{
    const uint64_t bytes = 10;
    const ByteLedger zero = LedgerFor(ByteOperation::Zero, bytes);
    const ByteLedger copy = LedgerFor(ByteOperation::Memcpy, bytes);
    const ByteLedger xorr = LedgerFor(ByteOperation::Xor, bytes);
    const ByteLedger addset = LedgerFor(ByteOperation::Addset, bytes);
    const ByteLedger add2 = LedgerFor(ByteOperation::Add2, bytes);
    const ByteLedger mul = LedgerFor(ByteOperation::Mul, bytes);
    const ByteLedger div = LedgerFor(ByteOperation::Div, bytes);
    const ByteLedger muladd = LedgerFor(ByteOperation::Muladd, bytes);
    const ByteLedger fused = LedgerFor(ByteOperation::FaninFused, bytes, 3);
    const ByteLedger chained = LedgerFor(ByteOperation::FaninChained, bytes, 3);
    if (zero.LogicalBytes != 10 || zero.ReadBytes != 0 || zero.WriteBytes != 10 ||
        copy.ReadBytes != 10 || copy.WriteBytes != 10 ||
        xorr.ReadBytes != 20 || xorr.WriteBytes != 10 ||
        addset.ReadBytes != 20 || addset.WriteBytes != 10 ||
        add2.ReadBytes != 30 || add2.WriteBytes != 10 ||
        mul.ReadBytes != 10 || mul.WriteBytes != 10 ||
        div.ReadBytes != 10 || div.WriteBytes != 10 ||
        muladd.ReadBytes != 20 || muladd.WriteBytes != 10 ||
        fused.ReadBytes != 40 || fused.WriteBytes != 10 ||
        chained.ReadBytes != 60 || chained.WriteBytes != 30)
    {
        return false;
    }
    ByteLedger mixed;
    mixed += zero;
    mixed += copy;
    mixed += xorr;
    mixed += addset;
    mixed += add2;
    mixed += mul;
    mixed += div;
    mixed += muladd;
    if (mixed.LogicalBytes != 80 || mixed.ReadBytes != 120 ||
        mixed.WriteBytes != 80 || mixed.TrafficBytes() != 200)
    {
        return false;
    }
    ByteLedger maximum;
    maximum.LogicalBytes = UINT64_MAX;
    maximum.ReadBytes = UINT64_MAX;
    maximum.WriteBytes = UINT64_MAX;
    const ByteLedger before = maximum;
    if (TryAccumulateByteLedger(maximum, zero) ||
        maximum.LogicalBytes != before.LogicalBytes ||
        maximum.ReadBytes != before.ReadBytes ||
        maximum.WriteBytes != before.WriteBytes)
    {
        return false;
    }
    ByteLedger traffic_limit;
    traffic_limit.ReadBytes = UINT64_MAX - 5;
    traffic_limit.WriteBytes = 5;
    const ByteLedger traffic_before = traffic_limit;
    ByteLedger one_read;
    one_read.ReadBytes = 1;
    return !TryAccumulateByteLedger(traffic_limit, one_read) &&
        traffic_limit.LogicalBytes == traffic_before.LogicalBytes &&
        traffic_limit.ReadBytes == traffic_before.ReadBytes &&
        traffic_limit.WriteBytes == traffic_before.WriteBytes;
}

} // namespace wirehair_experiments
