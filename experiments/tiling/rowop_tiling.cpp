// Standalone benchmark for replaying block-data row operations by byte tiles.

#include <algorithm>
#include <cerrno>
#include <chrono>
#include <climits>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <new>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "gf256.h"
#include "../ByteLedger.h"

namespace {

using Clock = std::chrono::steady_clock;

static const size_t kDefaultTraceMemoryMiB = 256u;
static const size_t kDefaultMaxTraceOperations = 1000000u;
static const size_t kMaxTraceLineBytes = 4096u;
static const size_t kTraceFixedOverheadBytes = 1024u * 1024u;
static const char kTraceCsvHeader[] =
    "stage,op_type,dst_kind,dst_block,dst_offset,src0_kind,src0_block,"
    "src0_offset,src1_kind,src1_block,src1_offset,scalar,bytes";

struct Rng
{
    uint64_t State;

    explicit Rng(uint64_t seed) : State(seed) {}

    uint64_t next()
    {
        uint64_t z = (State += UINT64_C(0x9e3779b97f4a7c15));
        z = (z ^ (z >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
        z = (z ^ (z >> 27)) * UINT64_C(0x94d049bb133111eb);
        return z ^ (z >> 31);
    }

    uint32_t uniform(uint32_t n)
    {
        return n <= 1 ? 0 : (uint32_t)(next() % n);
    }
};

struct Op
{
    uint32_t Dst;
    uint32_t Src;
};

enum class BlockKind
{
    None,
    Recovery,
    Input
};

enum class TraceOpKind
{
    Zero,
    Memcpy,
    Xor,
    Addset,
    Add2,
    Muladd,
    Div
};

struct TraceRef
{
    BlockKind Kind = BlockKind::None;
    uint32_t Block = 0;
    uint32_t Offset = 0;
};

struct TraceOp
{
    TraceOpKind Kind = TraceOpKind::Zero;
    TraceRef Dst;
    TraceRef Src0;
    TraceRef Src1;
    uint8_t Scalar = 0;
    uint32_t Bytes = 0;
};

struct TraceSchedule
{
    std::string Name;
    size_t BlockBytes = 0;
    uint32_t RecoveryBlocks = 0;
    uint32_t InputBlocks = 0;
    wirehair_experiments::ByteLedger Ledger;
    std::vector<TraceOp> Ops;
};

static bool checked_size_add(size_t a, size_t b, size_t* result)
{
    if (a > (size_t)-1 - b) {
        return false;
    }
    *result = a + b;
    return true;
}

static bool checked_size_multiply(size_t a, size_t b, size_t* result)
{
    if (a != 0u && b > (size_t)-1 / a) {
        return false;
    }
    *result = a * b;
    return true;
}

static bool trace_persistent_storage_bytes(
    const TraceSchedule& schedule,
    unsigned repeats,
    size_t* result)
{
    size_t operation_bytes = 0;
    size_t sample_bytes = 0;
    size_t name_bytes = 0;
    size_t total = kTraceFixedOverheadBytes;
    if (!checked_size_multiply(
            schedule.Ops.capacity(), sizeof(TraceOp), &operation_bytes) ||
        !checked_size_multiply(
            (size_t)repeats, sizeof(double), &sample_bytes) ||
        !checked_size_add(schedule.Name.capacity(), 1u, &name_bytes) ||
        !checked_size_add(total, operation_bytes, &total) ||
        !checked_size_add(total, name_bytes, &total) ||
        !checked_size_add(total, sample_bytes, &total))
    {
        return false;
    }
    *result = total;
    return true;
}

static bool reserve_for_trace_op(
    TraceSchedule* schedule,
    size_t max_operations,
    size_t max_memory_bytes)
{
    if (schedule->Ops.size() < schedule->Ops.capacity()) {
        return true;
    }

    const size_t old_capacity = schedule->Ops.capacity();
    size_t new_capacity = old_capacity == 0u ? 1024u : old_capacity;
    if (new_capacity > max_operations) {
        new_capacity = max_operations;
    }
    if (new_capacity < max_operations)
    {
        if (new_capacity > max_operations - new_capacity) {
            new_capacity = max_operations;
        }
        else {
            new_capacity += new_capacity;
        }
    }
    if (new_capacity <= old_capacity ||
        new_capacity > schedule->Ops.max_size())
    {
        std::fprintf(stderr,
            "schedule operation storage exceeds container limits\n");
        return false;
    }

    size_t old_bytes = 0;
    size_t new_bytes = 0;
    size_t name_bytes = 0;
    size_t transient_bytes = kTraceFixedOverheadBytes;
    if (!checked_size_multiply(old_capacity, sizeof(TraceOp), &old_bytes) ||
        !checked_size_multiply(new_capacity, sizeof(TraceOp), &new_bytes) ||
        !checked_size_add(schedule->Name.capacity(), 1u, &name_bytes) ||
        !checked_size_add(transient_bytes, old_bytes, &transient_bytes) ||
        !checked_size_add(transient_bytes, new_bytes, &transient_bytes) ||
        !checked_size_add(transient_bytes, name_bytes, &transient_bytes) ||
        transient_bytes > max_memory_bytes)
    {
        std::fprintf(stderr,
            "schedule operation storage exceeds --max-memory-mib policy\n");
        return false;
    }

    schedule->Ops.reserve(new_capacity);
    return true;
}

static wirehair_experiments::ByteOperation ledger_operation(TraceOpKind kind)
{
    using wirehair_experiments::ByteOperation;
    switch (kind)
    {
    case TraceOpKind::Zero: return ByteOperation::Zero;
    case TraceOpKind::Memcpy: return ByteOperation::Memcpy;
    case TraceOpKind::Xor: return ByteOperation::Xor;
    case TraceOpKind::Addset: return ByteOperation::Addset;
    case TraceOpKind::Add2: return ByteOperation::Add2;
    case TraceOpKind::Muladd: return ByteOperation::Muladd;
    case TraceOpKind::Div: return ByteOperation::Div;
    }
    return ByteOperation::Zero;
}

struct Case
{
    const char* Name;
    size_t BlockBytes;
    uint32_t Blocks;
    uint32_t OpsPerSource;
    uint32_t SourcePasses;
};

static double now_sec()
{
    return std::chrono::duration<double>(Clock::now().time_since_epoch()).count();
}

static void xor_words(uint64_t* dst, const uint64_t* src, size_t words)
{
    for (size_t i = 0; i < words; ++i) {
        dst[i] ^= src[i];
    }
}

static uint64_t checksum_words(const std::vector<uint64_t>& data)
{
    uint64_t h = UINT64_C(1469598103934665603);
    for (uint64_t x : data)
    {
        h ^= x;
        h *= UINT64_C(1099511628211);
    }
    return h;
}

static void fill_data(std::vector<uint64_t>& data, uint64_t seed)
{
    Rng rng(seed);
    for (uint64_t& x : data) {
        x = rng.next();
    }
}

static void fill_bytes(std::vector<uint8_t>& data, uint64_t seed)
{
    Rng rng(seed);
    uint64_t word = 0;
    for (size_t i = 0; i < data.size(); ++i)
    {
        if ((i & 7) == 0) {
            word = rng.next();
        }
        data[i] = (uint8_t)(word >> (8 * (i & 7)));
    }
}

static uint64_t checksum_bytes(const std::vector<uint8_t>& data)
{
    uint64_t h = UINT64_C(1469598103934665603);
    for (uint8_t x : data)
    {
        h ^= x;
        h *= UINT64_C(1099511628211);
    }
    return h;
}

static std::vector<Op> make_schedule(
    uint32_t blocks,
    uint32_t ops_per_source,
    uint32_t source_passes,
    uint64_t seed)
{
    Rng rng(seed);
    std::vector<Op> ops;
    ops.reserve((size_t)blocks * source_passes * ops_per_source);

    for (uint32_t pass = 0; pass < source_passes; ++pass)
    {
        for (uint32_t src = 0; src < blocks; ++src)
        {
            for (uint32_t j = 0; j < ops_per_source; ++j)
            {
                uint32_t dst = rng.uniform(blocks);
                if (dst == src) {
                    dst = (dst + 1) % blocks;
                }
                ops.push_back(Op{dst, src});
            }
        }
    }

    return ops;
}

static void replay_untiled(
    std::vector<uint64_t>& data,
    const std::vector<Op>& ops,
    size_t words_per_block)
{
    for (const Op& op : ops)
    {
        uint64_t* dst = &data[(size_t)op.Dst * words_per_block];
        const uint64_t* src = &data[(size_t)op.Src * words_per_block];
        xor_words(dst, src, words_per_block);
    }
}

static void replay_tiled(
    std::vector<uint64_t>& data,
    const std::vector<Op>& ops,
    size_t words_per_block,
    size_t tile_bytes)
{
    size_t tile_words = tile_bytes / sizeof(uint64_t);
    if (tile_words == 0) {
        tile_words = 1;
    }

    for (size_t offset = 0; offset < words_per_block; offset += tile_words)
    {
        const size_t words = std::min(tile_words, words_per_block - offset);
        for (const Op& op : ops)
        {
            uint64_t* dst = &data[(size_t)op.Dst * words_per_block + offset];
            const uint64_t* src = &data[(size_t)op.Src * words_per_block + offset];
            xor_words(dst, src, words);
        }
    }
}

static std::vector<std::string> split_csv(const std::string& line)
{
    std::vector<std::string> fields;
    size_t start = 0;
    for (;;)
    {
        const size_t comma = line.find(',', start);
        if (comma == std::string::npos)
        {
            fields.push_back(line.substr(start));
            break;
        }
        fields.push_back(line.substr(start, comma - start));
        start = comma + 1;
    }
    return fields;
}

enum class TraceLineResult
{
    Line,
    End,
    TooLong,
    ReadError
};

static TraceLineResult read_trace_line(
    std::istream& in,
    std::string* line)
{
    line->clear();
    for (;;)
    {
        const int next = in.get();
        if (next == std::char_traits<char>::eof())
        {
            if (in.bad()) {
                return TraceLineResult::ReadError;
            }
            return line->empty() ? TraceLineResult::End : TraceLineResult::Line;
        }
        if (next == '\n')
        {
            if (!line->empty() && line->back() == '\r') {
                line->pop_back();
            }
            return TraceLineResult::Line;
        }
        if (line->size() >= kMaxTraceLineBytes) {
            return TraceLineResult::TooLong;
        }
        line->push_back((char)next);
    }
}

static bool parse_u32(const std::string& s, uint32_t* out)
{
    if (s.empty() || s[0] < '0' || s[0] > '9') {
        return false;
    }
    errno = 0;
    char* end = nullptr;
    const unsigned long v = std::strtoul(s.c_str(), &end, 0);
    if (!end || *end || errno == ERANGE || v > UINT32_MAX) {
        return false;
    }
    *out = (uint32_t)v;
    return true;
}

static bool parse_int_field(const std::string& s, int* out)
{
    if (s.empty()) {
        return false;
    }
    errno = 0;
    char* end = nullptr;
    const long v = std::strtol(s.c_str(), &end, 0);
    if (!end || *end || errno == ERANGE || v < INT_MIN || v > INT_MAX) {
        return false;
    }
    *out = (int)v;
    return true;
}

static bool parse_size_field(const std::string& s, size_t* out)
{
    if (s.empty() || s[0] < '0' || s[0] > '9') {
        return false;
    }
    errno = 0;
    char* end = nullptr;
    const unsigned long long v = std::strtoull(s.c_str(), &end, 0);
    if (!end || *end || errno == ERANGE ||
        v > (unsigned long long)((size_t)-1))
    {
        return false;
    }
    *out = (size_t)v;
    return true;
}

static bool parse_size_list(const std::string& s, std::vector<size_t>* out)
{
    out->clear();
    if (s.empty() || s[s.size() - 1] == ',') {
        return false;
    }
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, ','))
    {
        size_t value = 0;
        if (!parse_size_field(item, &value) || value == 0) {
            return false;
        }
        out->push_back(value);
    }
    return !out->empty();
}

static bool parse_block_kind(const std::string& s, BlockKind* out)
{
    if (s == "none") {
        *out = BlockKind::None;
    }
    else if (s == "recovery") {
        *out = BlockKind::Recovery;
    }
    else if (s == "input") {
        *out = BlockKind::Input;
    }
    else {
        return false;
    }
    return true;
}

static bool parse_trace_op_kind(const std::string& s, TraceOpKind* out)
{
    if (s == "zero") {
        *out = TraceOpKind::Zero;
    }
    else if (s == "memcpy") {
        *out = TraceOpKind::Memcpy;
    }
    else if (s == "xor") {
        *out = TraceOpKind::Xor;
    }
    else if (s == "addset") {
        *out = TraceOpKind::Addset;
    }
    else if (s == "add2") {
        *out = TraceOpKind::Add2;
    }
    else if (s == "muladd") {
        *out = TraceOpKind::Muladd;
    }
    else if (s == "div") {
        *out = TraceOpKind::Div;
    }
    else {
        return false;
    }
    return true;
}

static bool parse_ref(
    const std::vector<std::string>& fields,
    size_t kind_i,
    size_t block_i,
    size_t offset_i,
    TraceRef* ref)
{
    int block = -1;
    if (!parse_block_kind(fields[kind_i], &ref->Kind) ||
        !parse_int_field(fields[block_i], &block) ||
        !parse_u32(fields[offset_i], &ref->Offset))
    {
        return false;
    }

    if (ref->Kind == BlockKind::None)
    {
        if (block != -1) {
            return false;
        }
        ref->Block = 0;
    }
    else
    {
        if (block < 0) {
            return false;
        }
        ref->Block = (uint32_t)block;
    }
    return true;
}

static bool validate_ref_bounds(
    const TraceRef& ref,
    size_t bytes,
    size_t block_bytes,
    uint32_t recovery_blocks,
    uint32_t input_blocks)
{
    if (ref.Kind == BlockKind::None) {
        return true;
    }
    if ((uint64_t)ref.Offset + bytes > block_bytes) {
        return false;
    }
    if (ref.Kind == BlockKind::Recovery) {
        return ref.Block < recovery_blocks;
    }
    return ref.Block < input_blocks;
}

static bool trace_ref_ranges_overlap(
    const TraceRef& a,
    const TraceRef& b,
    uint32_t bytes)
{
    if (bytes == 0u || a.Kind == BlockKind::None ||
        a.Kind != b.Kind || a.Block != b.Block)
    {
        return false;
    }
    const uint64_t a_begin = a.Offset;
    const uint64_t b_begin = b.Offset;
    return a_begin < b_begin + bytes && b_begin < a_begin + bytes;
}

static bool trace_op_has_unsupported_overlap(const TraceOp& op)
{
    const bool dst_src0_overlap =
        trace_ref_ranges_overlap(op.Dst, op.Src0, op.Bytes);
    const bool dst_src1_overlap =
        trace_ref_ranges_overlap(op.Dst, op.Src1, op.Bytes);
    const bool src0_src1_overlap =
        trace_ref_ranges_overlap(op.Src0, op.Src1, op.Bytes);

    // gf256_div_mem explicitly supports exact in-place operation.  The other
    // replay kernels use memcpy or restrict-qualified operands.
    if (op.Kind == TraceOpKind::Div)
    {
        const bool exact_in_place =
            op.Dst.Kind == op.Src0.Kind &&
            op.Dst.Block == op.Src0.Block &&
            op.Dst.Offset == op.Src0.Offset;
        return dst_src0_overlap && !exact_in_place;
    }
    if (op.Kind == TraceOpKind::Addset || op.Kind == TraceOpKind::Add2) {
        return dst_src0_overlap || dst_src1_overlap || src0_src1_overlap;
    }
    return dst_src0_overlap;
}

static bool load_trace_schedule(
    const std::string& path,
    size_t requested_block_bytes,
    size_t max_operations,
    size_t max_memory_bytes,
    TraceSchedule* schedule)
{
    std::ifstream in(path.c_str());
    if (!in) {
        std::fprintf(stderr, "failed to open schedule: %s\n", path.c_str());
        return false;
    }

    schedule->Name = path;
    schedule->BlockBytes = requested_block_bytes;
    schedule->RecoveryBlocks = 0;
    schedule->InputBlocks = 0;
    schedule->Ledger = wirehair_experiments::ByteLedger();
    schedule->Ops.clear();

    std::string line;
    const TraceLineResult header_result = read_trace_line(in, &line);
    if (header_result == TraceLineResult::End)
    {
        std::fprintf(stderr, "empty schedule: %s\n", path.c_str());
        return false;
    }
    if (header_result == TraceLineResult::TooLong)
    {
        std::fprintf(stderr,
            "%s:1: schedule line exceeds %zu-byte limit\n",
            path.c_str(), kMaxTraceLineBytes);
        return false;
    }
    if (header_result == TraceLineResult::ReadError)
    {
        std::fprintf(stderr, "failed to read schedule: %s\n", path.c_str());
        return false;
    }
    if (line != kTraceCsvHeader)
    {
        std::fprintf(stderr, "%s:1: unexpected trace CSV header\n", path.c_str());
        return false;
    }

    size_t max_bytes = 0;
    unsigned line_no = 1;
    for (;;)
    {
        const TraceLineResult line_result = read_trace_line(in, &line);
        if (line_result == TraceLineResult::End) {
            break;
        }
        ++line_no;
        if (line_result == TraceLineResult::TooLong)
        {
            std::fprintf(stderr,
                "%s:%u: schedule line exceeds %zu-byte limit\n",
                path.c_str(), line_no, kMaxTraceLineBytes);
            return false;
        }
        if (line_result == TraceLineResult::ReadError)
        {
            std::fprintf(stderr,
                "%s:%u: failed to read schedule\n", path.c_str(), line_no);
            return false;
        }
        if (line.empty()) {
            continue;
        }

        const std::vector<std::string> fields = split_csv(line);
        if (fields.size() != 13)
        {
            std::fprintf(stderr, "%s:%u: expected 13 CSV fields, got %zu\n",
                path.c_str(), line_no, fields.size());
            return false;
        }

        TraceOp op;
        uint32_t scalar = 0;
        if (!parse_trace_op_kind(fields[1], &op.Kind) ||
            !parse_ref(fields, 2, 3, 4, &op.Dst) ||
            !parse_ref(fields, 5, 6, 7, &op.Src0) ||
            !parse_ref(fields, 8, 9, 10, &op.Src1) ||
            !parse_u32(fields[11], &scalar) ||
            !parse_u32(fields[12], &op.Bytes) ||
            scalar > 255)
        {
            std::fprintf(stderr, "%s:%u: invalid trace row\n", path.c_str(), line_no);
            return false;
        }
        op.Scalar = (uint8_t)scalar;

        if (op.Dst.Kind != BlockKind::Recovery)
        {
            std::fprintf(stderr, "%s:%u: destination must be a recovery block\n",
                path.c_str(), line_no);
            return false;
        }
        if (op.Bytes > (uint32_t)INT_MAX)
        {
            std::fprintf(stderr, "%s:%u: operation byte count is too large\n",
                path.c_str(), line_no);
            return false;
        }

        const bool has_src0 = op.Src0.Kind != BlockKind::None;
        const bool has_src1 = op.Src1.Kind != BlockKind::None;
        bool valid_sources = false;
        switch (op.Kind)
        {
        case TraceOpKind::Zero:
            valid_sources = !has_src0 && !has_src1;
            break;
        case TraceOpKind::Memcpy:
        case TraceOpKind::Xor:
        case TraceOpKind::Muladd:
            valid_sources = has_src0 && !has_src1;
            break;
        case TraceOpKind::Div:
            valid_sources = has_src0 && !has_src1 && op.Scalar != 0;
            break;
        case TraceOpKind::Addset:
        case TraceOpKind::Add2:
            valid_sources = has_src0 && has_src1;
            break;
        }
        if (!valid_sources)
        {
            std::fprintf(stderr, "%s:%u: invalid source fields for op_type=%s\n",
                path.c_str(), line_no, fields[1].c_str());
            return false;
        }
        if (trace_op_has_unsupported_overlap(op))
        {
            std::fprintf(stderr,
                "%s:%u: overlapping operands are unsupported for op_type=%s\n",
                path.c_str(), line_no, fields[1].c_str());
            return false;
        }

        const TraceRef refs[] = {op.Dst, op.Src0, op.Src1};
        for (const TraceRef& ref : refs)
        {
            if (ref.Kind == BlockKind::Recovery) {
                schedule->RecoveryBlocks = std::max(schedule->RecoveryBlocks, ref.Block + 1);
            }
            else if (ref.Kind == BlockKind::Input) {
                schedule->InputBlocks = std::max(schedule->InputBlocks, ref.Block + 1);
            }
        }

        max_bytes = std::max(max_bytes, (size_t)op.Bytes);
        if (schedule->Ops.size() >= max_operations)
        {
            std::fprintf(stderr,
                "schedule exceeds the %zu-operation replay limit\n",
                max_operations);
            return false;
        }
        const wirehair_experiments::ByteLedger operation_ledger =
            wirehair_experiments::LedgerFor(
                ledger_operation(op.Kind), op.Bytes);
        if (!wirehair_experiments::TryAccumulateByteLedger(
                schedule->Ledger, operation_ledger))
        {
            std::fprintf(stderr, "schedule byte ledger overflows uint64_t\n");
            return false;
        }
        if (!reserve_for_trace_op(
                schedule, max_operations, max_memory_bytes)) {
            return false;
        }
        schedule->Ops.push_back(op);
    }

    if (schedule->Ops.empty())
    {
        std::fprintf(stderr, "schedule has no operations: %s\n", path.c_str());
        return false;
    }

    if (schedule->BlockBytes == 0) {
        schedule->BlockBytes = max_bytes;
    }
    if (schedule->BlockBytes == 0)
    {
        std::fprintf(stderr, "could not infer --block-bytes for schedule: %s\n", path.c_str());
        return false;
    }

    for (const TraceOp& op : schedule->Ops)
    {
        if (!validate_ref_bounds(op.Dst, op.Bytes, schedule->BlockBytes, schedule->RecoveryBlocks, schedule->InputBlocks) ||
            !validate_ref_bounds(op.Src0, op.Bytes, schedule->BlockBytes, schedule->RecoveryBlocks, schedule->InputBlocks) ||
            !validate_ref_bounds(op.Src1, op.Bytes, schedule->BlockBytes, schedule->RecoveryBlocks, schedule->InputBlocks))
        {
            std::fprintf(stderr, "schedule op out of bounds for block_bytes=%zu\n",
                schedule->BlockBytes);
            return false;
        }
    }

    return true;
}

static uint8_t* mutable_ref_ptr(
    std::vector<uint8_t>& recovery,
    const TraceRef& ref,
    size_t block_bytes,
    size_t delta)
{
    return &recovery[(size_t)ref.Block * block_bytes + ref.Offset + delta];
}

static const uint8_t* const_ref_ptr(
    const std::vector<uint8_t>& recovery,
    const std::vector<uint8_t>& input,
    const TraceRef& ref,
    size_t block_bytes,
    size_t delta)
{
    if (ref.Kind == BlockKind::Recovery) {
        return &recovery[(size_t)ref.Block * block_bytes + ref.Offset + delta];
    }
    if (ref.Kind == BlockKind::Input) {
        return &input[(size_t)ref.Block * block_bytes + ref.Offset + delta];
    }
    return nullptr;
}

static void apply_trace_op(
    std::vector<uint8_t>& recovery,
    const std::vector<uint8_t>& input,
    const TraceOp& op,
    size_t block_bytes,
    size_t delta,
    size_t bytes)
{
    if (bytes == 0) {
        return;
    }
    uint8_t* dst = mutable_ref_ptr(recovery, op.Dst, block_bytes, delta);
    const uint8_t* src0 = const_ref_ptr(recovery, input, op.Src0, block_bytes, delta);
    const uint8_t* src1 = const_ref_ptr(recovery, input, op.Src1, block_bytes, delta);
    const int n = (int)bytes;

    switch (op.Kind)
    {
    case TraceOpKind::Zero:
        std::memset(dst, 0, bytes);
        break;
    case TraceOpKind::Memcpy:
        std::memcpy(dst, src0, bytes);
        break;
    case TraceOpKind::Xor:
        gf256_add_mem(dst, src0, n);
        break;
    case TraceOpKind::Addset:
        gf256_addset_mem(dst, src0, src1, n);
        break;
    case TraceOpKind::Add2:
        gf256_add2_mem(dst, src0, src1, n);
        break;
    case TraceOpKind::Muladd:
        gf256_muladd_mem(dst, op.Scalar, src0, n);
        break;
    case TraceOpKind::Div:
        gf256_div_mem(dst, src0, op.Scalar, n);
        break;
    }
}

static void replay_trace_untiled(
    std::vector<uint8_t>& recovery,
    const std::vector<uint8_t>& input,
    const TraceSchedule& schedule)
{
    for (const TraceOp& op : schedule.Ops) {
        apply_trace_op(recovery, input, op, schedule.BlockBytes, 0, op.Bytes);
    }
}

static void replay_trace_tiled(
    std::vector<uint8_t>& recovery,
    const std::vector<uint8_t>& input,
    const TraceSchedule& schedule,
    size_t tile_bytes)
{
    for (size_t tile = 0; tile < schedule.BlockBytes; tile += tile_bytes)
    {
        const size_t tile_end = std::min(tile + tile_bytes, schedule.BlockBytes);
        for (const TraceOp& op : schedule.Ops)
        {
            const size_t op_begin = op.Dst.Offset;
            const size_t op_end = op_begin + op.Bytes;
            const size_t begin = std::max(tile, op_begin);
            const size_t end = std::min(tile_end, op_end);
            if (begin < end) {
                apply_trace_op(recovery, input, op, schedule.BlockBytes, begin - op_begin, end - begin);
            }
        }
    }
}

static double median(std::vector<double>& values)
{
    std::sort(values.begin(), values.end());
    return values[values.size() / 2];
}

static bool verify_case(size_t block_bytes, size_t tile_bytes)
{
    const uint32_t blocks = 64;
    const size_t words_per_block = block_bytes / sizeof(uint64_t);
    std::vector<uint64_t> base((size_t)blocks * words_per_block);
    fill_data(base, UINT64_C(0x766572696679));
    const std::vector<Op> ops = make_schedule(
        blocks, 8, 2, UINT64_C(0x7363686564756c65));

    std::vector<uint64_t> untiled = base;
    std::vector<uint64_t> tiled = base;
    replay_untiled(untiled, ops, words_per_block);
    replay_tiled(tiled, ops, words_per_block, tile_bytes);
    return untiled == tiled;
}

static double time_untiled(
    const std::vector<uint64_t>& base,
    const std::vector<Op>& ops,
    size_t words_per_block,
    unsigned repeats,
    uint64_t* checksum)
{
    std::vector<double> samples;
    samples.reserve(repeats);
    std::vector<uint64_t> data;
    for (unsigned repeat = 0; repeat < repeats; ++repeat)
    {
        data = base;
        const double start = now_sec();
        replay_untiled(data, ops, words_per_block);
        samples.push_back(now_sec() - start);
    }
    *checksum = checksum_words(data);
    return median(samples);
}

static double time_tiled(
    const std::vector<uint64_t>& base,
    const std::vector<Op>& ops,
    size_t words_per_block,
    size_t tile_bytes,
    unsigned repeats,
    uint64_t* checksum)
{
    std::vector<double> samples;
    samples.reserve(repeats);
    std::vector<uint64_t> data;
    for (unsigned repeat = 0; repeat < repeats; ++repeat)
    {
        data = base;
        const double start = now_sec();
        replay_tiled(data, ops, words_per_block, tile_bytes);
        samples.push_back(now_sec() - start);
    }
    *checksum = checksum_words(data);
    return median(samples);
}

static bool make_trace_base(
    const TraceSchedule& schedule,
    size_t max_memory_bytes,
    unsigned repeats,
    std::vector<uint8_t>* recovery,
    std::vector<uint8_t>* input)
{
    size_t recovery_bytes = 0;
    size_t input_bytes = 0;
    if ((schedule.RecoveryBlocks != 0u &&
         schedule.BlockBytes > (size_t)-1 / schedule.RecoveryBlocks) ||
        (schedule.InputBlocks != 0u &&
         schedule.BlockBytes > (size_t)-1 / schedule.InputBlocks))
    {
        std::fprintf(stderr, "schedule storage size overflows size_t\n");
        return false;
    }
    recovery_bytes = (size_t)schedule.RecoveryBlocks * schedule.BlockBytes;
    input_bytes = (size_t)schedule.InputBlocks * schedule.BlockBytes;
    if (recovery_bytes > (size_t)-1 - input_bytes)
    {
        std::fprintf(stderr, "schedule aggregate storage size overflows size_t\n");
        return false;
    }
    if (recovery_bytes > ((size_t)-1 - input_bytes) / 3u)
    {
        std::fprintf(stderr, "schedule peak replay storage size overflows size_t\n");
        return false;
    }
    const size_t peak_bytes = 3u * recovery_bytes + input_bytes;
    size_t persistent_bytes = 0;
    size_t aggregate_peak_bytes = 0;
    if (!trace_persistent_storage_bytes(
            schedule, repeats, &persistent_bytes) ||
        !checked_size_add(
            persistent_bytes, peak_bytes, &aggregate_peak_bytes))
    {
        std::fprintf(stderr,
            "schedule aggregate replay storage size overflows size_t\n");
        return false;
    }
    if (aggregate_peak_bytes > max_memory_bytes)
    {
        std::fprintf(stderr,
            "schedule aggregate storage exceeds --max-memory-mib policy\n");
        return false;
    }
    if (recovery_bytes > recovery->max_size() || input_bytes > input->max_size())
    {
        std::fprintf(stderr, "schedule storage size exceeds container limits\n");
        return false;
    }
    try
    {
        recovery->resize(recovery_bytes);
        input->resize(input_bytes);
    }
    catch (const std::bad_alloc&)
    {
        std::fprintf(stderr, "schedule storage allocation failed\n");
        return false;
    }
    catch (const std::length_error&)
    {
        std::fprintf(stderr, "schedule storage size exceeds container limits\n");
        return false;
    }
    fill_bytes(*recovery, UINT64_C(0x7265636f76657279) ^ schedule.BlockBytes);
    fill_bytes(*input, UINT64_C(0x696e7075745f5f5f) ^ schedule.BlockBytes);
    return true;
}

static bool verify_trace_schedule(
    const TraceSchedule& schedule,
    const std::vector<size_t>& tiles,
    size_t max_memory_bytes,
    unsigned repeats)
{
    try
    {
        std::vector<uint8_t> recovery_base;
        std::vector<uint8_t> input;
        if (!make_trace_base(
                schedule, max_memory_bytes, repeats,
                &recovery_base, &input)) {
            return false;
        }

        std::vector<uint8_t> untiled = recovery_base;
        replay_trace_untiled(untiled, input, schedule);

        bool checked_tile = false;
        for (size_t tile : tiles)
        {
            if (tile >= schedule.BlockBytes) {
                continue;
            }
            checked_tile = true;
            std::vector<uint8_t> tiled = recovery_base;
            replay_trace_tiled(tiled, input, schedule, tile);
            if (tiled != untiled)
            {
                std::fprintf(stderr,
                    "schedule verification failed for tile=%zu\n", tile);
                return false;
            }
        }

        if (!checked_tile)
        {
            std::fprintf(stderr,
                "no tile size below block_bytes=%zu was provided\n",
                schedule.BlockBytes);
            return false;
        }

        return true;
    }
    catch (const std::bad_alloc&)
    {
        std::fprintf(stderr, "schedule verification allocation failed\n");
        return false;
    }
    catch (const std::length_error&)
    {
        std::fprintf(stderr,
            "schedule verification storage exceeds container limits\n");
        return false;
    }
}

static double time_trace_untiled(
    const std::vector<uint8_t>& recovery_base,
    const std::vector<uint8_t>& input,
    const TraceSchedule& schedule,
    unsigned repeats,
    std::vector<uint8_t>& recovery,
    std::vector<double>& samples,
    uint64_t* checksum)
{
    samples.clear();
    for (unsigned repeat = 0; repeat < repeats; ++repeat)
    {
        std::copy(recovery_base.begin(), recovery_base.end(), recovery.begin());
        const double start = now_sec();
        replay_trace_untiled(recovery, input, schedule);
        samples.push_back(now_sec() - start);
    }
    *checksum = checksum_bytes(recovery);
    return median(samples);
}

static double time_trace_tiled(
    const std::vector<uint8_t>& recovery_base,
    const std::vector<uint8_t>& input,
    const TraceSchedule& schedule,
    size_t tile_bytes,
    unsigned repeats,
    std::vector<uint8_t>& recovery,
    std::vector<double>& samples,
    uint64_t* checksum)
{
    samples.clear();
    for (unsigned repeat = 0; repeat < repeats; ++repeat)
    {
        std::copy(recovery_base.begin(), recovery_base.end(), recovery.begin());
        const double start = now_sec();
        replay_trace_tiled(recovery, input, schedule, tile_bytes);
        samples.push_back(now_sec() - start);
    }
    *checksum = checksum_bytes(recovery);
    return median(samples);
}

static void print_result(
    const Case& c,
    const char* mode,
    size_t tile_bytes,
    size_t op_count,
    double seconds,
    uint64_t checksum)
{
    wirehair_experiments::ByteLedger ledger =
        wirehair_experiments::LedgerFor(
            wirehair_experiments::ByteOperation::Xor, c.BlockBytes);
    ledger.LogicalBytes *= op_count;
    ledger.ReadBytes *= op_count;
    ledger.WriteBytes *= op_count;
    const double logical_gib = (double)ledger.LogicalBytes / 1073741824.0;
    const double memory_gib_3stream = 3.0 * logical_gib;
    const double gib_per_s = memory_gib_3stream / seconds;
    const double read_gib = (double)ledger.ReadBytes / 1073741824.0;
    const double write_gib = (double)ledger.WriteBytes / 1073741824.0;
    const double traffic_gib = read_gib + write_gib;
    std::printf(
        "%s,%zu,%u,%zu,%s,%zu,%.6f,%.6f,%.6f,%.3f,0x%016llx,"
        "2,%.9f,%.9f,%.9f,%.9f,%.3f,%.3f\n",
        c.Name,
        c.BlockBytes,
        c.Blocks,
        op_count,
        mode,
        tile_bytes,
        seconds * 1000.0,
        logical_gib,
        memory_gib_3stream,
        gib_per_s,
        (unsigned long long)checksum,
        logical_gib,
        read_gib,
        write_gib,
        traffic_gib,
        logical_gib / seconds,
        traffic_gib / seconds);
}

static void print_trace_result(
    const TraceSchedule& schedule,
    const char* mode,
    size_t tile_bytes,
    double seconds,
    uint64_t checksum)
{
    const double logical_gib =
        (double)schedule.Ledger.LogicalBytes / 1073741824.0;
    const double memory_gib_3stream = 3.0 * logical_gib;
    const double gib_per_s = seconds > 0.0 ? memory_gib_3stream / seconds : 0.0;
    const double read_gib =
        (double)schedule.Ledger.ReadBytes / 1073741824.0;
    const double write_gib =
        (double)schedule.Ledger.WriteBytes / 1073741824.0;
    const double traffic_gib = read_gib + write_gib;
    std::printf(
        "%s,%zu,%u,%u,%zu,%s,%zu,%.6f,%.6f,%.6f,%.3f,0x%016llx,"
        "2,%.9f,%.9f,%.9f,%.9f,%.3f,%.3f\n",
        schedule.Name.c_str(),
        schedule.BlockBytes,
        schedule.RecoveryBlocks,
        schedule.InputBlocks,
        schedule.Ops.size(),
        mode,
        tile_bytes,
        seconds * 1000.0,
        logical_gib,
        memory_gib_3stream,
        gib_per_s,
        (unsigned long long)checksum,
        logical_gib,
        read_gib,
        write_gib,
        traffic_gib,
        seconds > 0.0 ? logical_gib / seconds : 0.0,
        seconds > 0.0 ? traffic_gib / seconds : 0.0);
}

} // namespace

int main(int argc, char** argv)
{
    if (gf256_init() != 0)
    {
        std::fprintf(stderr, "gf256_init failed\n");
        return 2;
    }
    if (!wirehair_experiments::VerifyByteLedger())
    {
        std::fprintf(stderr, "byte-ledger self-check failed\n");
        return 2;
    }
    std::fprintf(stderr, "byte-ledger self-check passed\n");

    const Case cases[] = {
        {"mtu1280", 1280, 32768, 8, 1},
        {"kib100", 102400, 2048, 8, 1},
        {"mib1", 1048576, 256, 8, 1},
    };
    std::vector<size_t> tiles = {
        256,        // sub-block tiles so the small-block cases measure
        512,        // a real tiled replay instead of degenerating to untiled
        16 * 1024,
        32 * 1024,
        64 * 1024,
        128 * 1024,
        256 * 1024,
        512 * 1024,
    };
    unsigned repeats = 3;
    bool verify_only = false;
    std::string schedule_path;
    size_t schedule_block_bytes = 0;
    size_t max_memory_mib = kDefaultTraceMemoryMiB;
    size_t max_operations = kDefaultMaxTraceOperations;

    for (int i = 1; i < argc; ++i)
    {
        if (!std::strcmp(argv[i], "--verify-only")) {
            verify_only = true;
        }
        else if (!std::strcmp(argv[i], "--schedule") && i + 1 < argc) {
            schedule_path = argv[++i];
        }
        else if (!std::strcmp(argv[i], "--block-bytes") && i + 1 < argc)
        {
            if (!parse_size_field(argv[++i], &schedule_block_bytes) ||
                schedule_block_bytes == 0)
            {
                std::fprintf(stderr, "--block-bytes requires a positive integer\n");
                return 1;
            }
        }
        else if (!std::strcmp(argv[i], "--tiles") && i + 1 < argc)
        {
            if (!parse_size_list(argv[++i], &tiles))
            {
                std::fprintf(stderr, "--tiles requires a non-empty CSV of positive integers\n");
                return 1;
            }
        }
        else if (!std::strcmp(argv[i], "--repeats") && i + 1 < argc)
        {
            size_t parsed = 0;
            if (!parse_size_field(argv[++i], &parsed) || parsed == 0 || parsed > UINT32_MAX)
            {
                std::fprintf(stderr, "--repeats requires a positive integer\n");
                return 1;
            }
            repeats = (unsigned)parsed;
        }
        else if (!std::strcmp(argv[i], "--max-memory-mib") && i + 1 < argc)
        {
            if (!parse_size_field(argv[++i], &max_memory_mib) ||
                max_memory_mib == 0u ||
                max_memory_mib > (size_t)-1 / (1024u * 1024u))
            {
                std::fprintf(stderr,
                    "--max-memory-mib requires a representable positive integer\n");
                return 1;
            }
        }
        else if (!std::strcmp(argv[i], "--max-operations") && i + 1 < argc)
        {
            if (!parse_size_field(argv[++i], &max_operations) ||
                max_operations == 0u ||
                max_operations > kDefaultMaxTraceOperations)
            {
                std::fprintf(stderr,
                    "--max-operations requires an integer in [1,%zu]\n",
                    kDefaultMaxTraceOperations);
                return 1;
            }
        }
        else
        {
            std::fprintf(stderr,
                "usage: rowop_tiling [--verify-only] [--schedule CSV] "
                "[--block-bytes N] [--tiles csv] [--repeats N] "
                "[--max-memory-mib N] [--max-operations N]\n");
            return 1;
        }
    }

    const size_t max_memory_bytes = max_memory_mib * 1024u * 1024u;

    if (!schedule_path.empty())
    {
        TraceSchedule schedule;
        bool loaded = false;
        try {
            loaded = load_trace_schedule(
                schedule_path, schedule_block_bytes, max_operations,
                max_memory_bytes, &schedule);
        }
        catch (const std::bad_alloc&) {
            std::fprintf(stderr, "schedule operation storage allocation failed\n");
            return 1;
        }
        catch (const std::length_error&) {
            std::fprintf(stderr,
                "schedule operation storage exceeds container limits\n");
            return 1;
        }
        if (!loaded) {
            return 1;
        }
        if (!verify_trace_schedule(
                schedule, tiles, max_memory_bytes, repeats)) {
            return 1;
        }
        if (verify_only)
        {
            std::printf("schedule verify: ok\n");
            return 0;
        }

        try
        {
            std::vector<uint8_t> recovery_base;
            std::vector<uint8_t> input;
            if (!make_trace_base(
                    schedule, max_memory_bytes, repeats,
                    &recovery_base, &input)) {
                return 1;
            }
            std::vector<uint8_t> replay_recovery(recovery_base.size());
            std::vector<double> samples;
            samples.reserve(repeats);

            std::printf(
                "case,block_bytes,recovery_blocks,input_blocks,ops,mode,tile_bytes,"
                "median_ms,logical_gib,memory_gib_3stream,gib_per_s,checksum,"
                "schema_version,logical_work_gib,estimated_read_gib,"
                "estimated_write_gib,estimated_traffic_gib,logical_work_gib_per_s,"
                "estimated_traffic_gib_per_s\n");

            uint64_t checksum = 0;
            const double untiled = time_trace_untiled(
                recovery_base, input, schedule, repeats,
                replay_recovery, samples, &checksum);
            print_trace_result(schedule, "untiled", 0, untiled, checksum);

            for (size_t tile : tiles)
            {
                if (tile >= schedule.BlockBytes) {
                    continue;
                }
                const double tiled = time_trace_tiled(
                    recovery_base, input, schedule, tile, repeats,
                    replay_recovery, samples, &checksum);
                print_trace_result(schedule, "tiled", tile, tiled, checksum);
            }
        }
        catch (const std::bad_alloc&)
        {
            std::fprintf(stderr, "schedule replay allocation failed\n");
            return 1;
        }
        catch (const std::length_error&)
        {
            std::fprintf(stderr,
                "schedule replay storage exceeds container limits\n");
            return 1;
        }

        return 0;
    }

    for (const Case& c : cases)
    {
        if (c.BlockBytes % sizeof(uint64_t) != 0)
        {
            std::fprintf(stderr, "block size must be a multiple of 8: %zu\n",
                c.BlockBytes);
            return 1;
        }
        for (size_t tile : tiles)
        {
            if (!verify_case(c.BlockBytes, tile))
            {
                std::fprintf(stderr, "verification failed for %s tile=%zu\n",
                    c.Name, tile);
                return 1;
            }
        }
    }
    if (verify_only)
    {
        std::printf("verify: ok\n");
        return 0;
    }

    std::printf(
        "case,block_bytes,blocks,ops,mode,tile_bytes,median_ms,"
        "logical_gib,memory_gib_3stream,gib_per_s,checksum,schema_version,"
        "logical_work_gib,estimated_read_gib,estimated_write_gib,"
        "estimated_traffic_gib,logical_work_gib_per_s,"
        "estimated_traffic_gib_per_s\n");

    for (const Case& c : cases)
    {
        const size_t words_per_block = c.BlockBytes / sizeof(uint64_t);
        std::vector<uint64_t> base((size_t)c.Blocks * words_per_block);
        fill_data(base, UINT64_C(0x646174615f626173) ^ (uint64_t)c.BlockBytes);
        const std::vector<Op> ops = make_schedule(
            c.Blocks, c.OpsPerSource, c.SourcePasses,
            UINT64_C(0x6f70735f73656564) ^ (uint64_t)c.BlockBytes);

        uint64_t checksum = 0;
        const double untiled = time_untiled(
            base, ops, words_per_block, repeats, &checksum);
        print_result(c, "untiled", 0, ops.size(), untiled, checksum);

        for (size_t tile : tiles)
        {
            // A tile at or above the block size replays the identical loop
            // as untiled; reporting it as "tiled" would present run-to-run
            // noise as a tiling result.
            if (tile >= c.BlockBytes) {
                continue;
            }
            const double tiled = time_tiled(
                base, ops, words_per_block, tile, repeats, &checksum);
            print_result(c, "tiled", tile, ops.size(), tiled, checksum);
        }
    }

    return 0;
}
