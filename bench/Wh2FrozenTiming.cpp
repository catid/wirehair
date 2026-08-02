#include "Wh2FrozenTiming.h"

#include <algorithm>

namespace wirehair {
namespace wh2_benchmark {
namespace {

static const uint64_t kSplitMixIncrement = UINT64_C(0x9e3779b97f4a7c15);
static const uint64_t kSplitMixMultiplier1 = UINT64_C(0xbf58476d1ce4e5b9);
static const uint64_t kSplitMixMultiplier2 = UINT64_C(0x94d049bb133111eb);
static const uint32_t kDevelopmentTimingRepetitions = 96u;
static const uint32_t kDevelopmentTimingWorkers = 8u;
static const uint32_t kDevelopmentTimingReceiveOverheadCap = 256u;
static const std::size_t kDevelopmentTimingCohorts = 264u;
static const char kDevelopmentTimingInterleavePolicy[] =
    "self-counterbalanced-repeat-major-v1";

static const uint32_t kTimingK[] = {
    8u, 32u, 100u, 128u, 512u, 1000u,
    2048u, 5000u, 8192u, 20000u, 32768u, 64000u
};
static const uint32_t kTimingWidths[] = { 64u, 1280u };
static const uint64_t kTrainingRoots[] = {
    UINT64_C(0xd1b54a32d192ed03),
    UINT64_C(0x94d049bb133111eb),
    UINT64_C(0x8538ecb5bd456ea3)
};

uint64_t SplitMix64Once(uint64_t value)
{
    value += kSplitMixIncrement;
    value = (value ^ (value >> 30)) * kSplitMixMultiplier1;
    value = (value ^ (value >> 27)) * kSplitMixMultiplier2;
    return value ^ (value >> 31);
}

std::string HexSeed(uint64_t seed)
{
    static const char hex[] = "0123456789abcdef";
    std::string text("0x0000000000000000");
    for (unsigned i = 0u; i < 16u; ++i) {
        text[17u - i] = hex[seed & 15u];
        seed >>= 4;
    }
    return text;
}

std::string BandForK(uint32_t K)
{
    if (K >= 2u && K <= 100u) {
        return "2-100";
    }
    if (K <= 1000u) {
        return "101-1000";
    }
    if (K <= 5000u) {
        return "1001-5000";
    }
    if (K <= 10000u) {
        return "5001-10000";
    }
    if (K <= 20000u) {
        return "10001-20000";
    }
    if (K <= 64000u) {
        return "20001-64000";
    }
    return std::string();
}

bool IsLowerSha256(const std::string& value)
{
    if (value.size() != 64u) {
        return false;
    }
    for (std::size_t i = 0u; i < value.size(); ++i) {
        const char ch = value[i];
        if (!((ch >= '0' && ch <= '9') || (ch >= 'a' && ch <= 'f'))) {
            return false;
        }
    }
    return true;
}

bool AppendJsonString(const std::string& value, std::string& output)
{
    static const char hex[] = "0123456789abcdef";
    output += '"';
    for (std::size_t i = 0u; i < value.size();)
    {
        const unsigned char ch = static_cast<unsigned char>(value[i]);
        if (ch >= 0x80u)
        {
            uint32_t codepoint = 0u;
            std::size_t continuation_count = 0u;
            if (ch >= 0xc2u && ch <= 0xdfu) {
                codepoint = ch & 0x1fu;
                continuation_count = 1u;
            }
            else if (ch >= 0xe0u && ch <= 0xefu) {
                codepoint = ch & 0x0fu;
                continuation_count = 2u;
            }
            else if (ch >= 0xf0u && ch <= 0xf4u) {
                codepoint = ch & 0x07u;
                continuation_count = 3u;
            }
            else {
                return false;
            }
            if (i + continuation_count >= value.size()) {
                return false;
            }
            for (std::size_t tail = 1u; tail <= continuation_count; ++tail)
            {
                const unsigned char next =
                    static_cast<unsigned char>(value[i + tail]);
                if ((next & 0xc0u) != 0x80u) {
                    return false;
                }
                codepoint = (codepoint << 6) | (next & 0x3fu);
            }
            const uint32_t minimum = continuation_count == 1u ? 0x80u :
                (continuation_count == 2u ? 0x800u : 0x10000u);
            if (codepoint < minimum || codepoint > 0x10ffffu ||
                (codepoint >= 0xd800u && codepoint <= 0xdfffu))
            {
                return false;
            }
            const auto append_escape = [&output](uint32_t unit)
            {
                output += "\\u";
                output += hex[(unit >> 12) & 15u];
                output += hex[(unit >> 8) & 15u];
                output += hex[(unit >> 4) & 15u];
                output += hex[unit & 15u];
            };
            if (codepoint <= 0xffffu) {
                append_escape(codepoint);
            }
            else
            {
                const uint32_t surrogate = codepoint - 0x10000u;
                append_escape(0xd800u + (surrogate >> 10));
                append_escape(0xdc00u + (surrogate & 0x3ffu));
            }
            i += continuation_count + 1u;
            continue;
        }
        switch (ch)
        {
        case '"': output += "\\\""; break;
        case '\\': output += "\\\\"; break;
        case '\b': output += "\\b"; break;
        case '\f': output += "\\f"; break;
        case '\n': output += "\\n"; break;
        case '\r': output += "\\r"; break;
        case '\t': output += "\\t"; break;
        default:
            if (ch < 0x20u || ch == 0x7fu) {
                output += "\\u00";
                output += hex[ch >> 4];
                output += hex[ch & 15u];
            }
            else {
                output += static_cast<char>(ch);
            }
            break;
        }
        ++i;
    }
    output += '"';
    return true;
}

bool IsCandidateArm(const std::string& arm)
{
    if (arm.empty() || arm == "wirehair1" || arm == "wirehair2_head") {
        return false;
    }
    std::string quoted;
    return AppendJsonString(arm, quoted);
}

bool ExactDevelopmentTimingBaseCell(const FrozenTimingBaseCell& cell)
{
    const uint32_t* const k_end =
        kTimingK + sizeof(kTimingK) / sizeof(kTimingK[0]);
    const uint32_t* const k_position =
        std::find(kTimingK, k_end, cell.K);
    const uint32_t* const width_end =
        kTimingWidths + sizeof(kTimingWidths) / sizeof(kTimingWidths[0]);
    const uint32_t* const width_position =
        std::find(kTimingWidths, width_end, cell.block_bytes);
    if (k_position == k_end || width_position == width_end ||
        cell.replicate >= kDevelopmentTimingRepetitions)
    {
        return false;
    }

    uint32_t expected_attempt = 0u;
    uint64_t expected_loss_seed = 0u;
    if (!DevelopmentTimingSeed(
            cell.replicate, expected_attempt, expected_loss_seed))
    {
        return false;
    }
    const uint32_t round_index =
        cell.replicate / kDevelopmentTimingWorkers;
    const uint32_t lane = cell.replicate % kDevelopmentTimingWorkers;
    const std::size_t expected_ordinal =
        static_cast<std::size_t>(round_index) * 192u +
        static_cast<std::size_t>(width_position - kTimingWidths) * 96u +
        static_cast<std::size_t>(k_position - kTimingK) * 8u + lane;
    return cell.ordinal == expected_ordinal &&
        cell.phase == "development" &&
        cell.band == BandForK(cell.K) &&
        cell.loss_ppm == 100000u &&
        cell.schedule == FrozenSchedule::Iid &&
        cell.base_seed_attempt == expected_attempt &&
        cell.base_loss_seed == expected_loss_seed &&
        cell.fixed_received_overhead == 4u &&
        cell.receive_overhead_cap ==
            kDevelopmentTimingReceiveOverheadCap &&
        cell.interleave_policy == kDevelopmentTimingInterleavePolicy &&
        cell.invocations_per_slot ==
            DevelopmentTimingInvocationsPerSlot(cell.K);
}

FrozenTimingBaseCell BaseCellFromQualified(const FrozenTimingCell& cell)
{
    FrozenTimingBaseCell base;
    base.ordinal = cell.ordinal;
    base.phase = cell.phase;
    base.band = cell.band;
    base.K = cell.K;
    base.block_bytes = cell.block_bytes;
    base.loss_ppm = cell.loss_ppm;
    base.schedule = cell.schedule;
    base.replicate = cell.replicate;
    base.base_seed_attempt = cell.base_seed_attempt;
    base.base_loss_seed = cell.base_loss_seed;
    base.fixed_received_overhead = cell.fixed_received_overhead;
    base.receive_overhead_cap = cell.receive_overhead_cap;
    base.interleave_policy = cell.interleave_policy;
    base.invocations_per_slot = cell.invocations_per_slot;
    return base;
}

bool ExactDevelopmentTimingCell(const FrozenTimingCell& cell)
{
    const FrozenTimingBaseCell base = BaseCellFromQualified(cell);
    uint64_t expected_loss_seed = 0u;
    return ExactDevelopmentTimingBaseCell(base) &&
        cell.base_cell_sha256 == TimingBaseCellSha256(base) &&
        DevelopmentTimingLossSeed(
            cell.base_loss_seed,
            cell.loss_retry_offset,
            expected_loss_seed) &&
        cell.loss_seed == expected_loss_seed;
}

bool ExactTimingPanel(const FrozenTimingPanel& panel)
{
    switch (panel.ordinal)
    {
    case 0u:
        return panel.panel_kind == FrozenPanelKind::AA &&
            panel.scope == FrozenTimingScope::DecoderSolve &&
            panel.left_arm == "wirehair2_head" &&
            panel.right_arm == "wirehair2_head";
    case 1u:
        return panel.panel_kind == FrozenPanelKind::AA &&
            panel.scope == FrozenTimingScope::ReceiveToSuccess &&
            panel.left_arm == "wirehair1" &&
            panel.right_arm == "wirehair1";
    case 2u:
        return panel.panel_kind == FrozenPanelKind::AA &&
            panel.scope == FrozenTimingScope::EncoderInitPlusFirstKSymbols &&
            panel.left_arm == "wirehair2_head" &&
            panel.right_arm == "wirehair2_head";
    case 3u:
        return panel.panel_kind == FrozenPanelKind::AA &&
            panel.scope == FrozenTimingScope::EncoderInitPlusFirstKSymbols &&
            panel.left_arm == "wirehair1" &&
            panel.right_arm == "wirehair1";
    case 4u:
        return panel.panel_kind == FrozenPanelKind::AA &&
            panel.scope == FrozenTimingScope::DecoderSolve &&
            panel.left_arm == panel.right_arm &&
            IsCandidateArm(panel.left_arm);
    case 5u:
        return panel.panel_kind == FrozenPanelKind::AA &&
            panel.scope == FrozenTimingScope::ReceiveToSuccess &&
            panel.left_arm == panel.right_arm &&
            IsCandidateArm(panel.left_arm);
    case 6u:
        return panel.panel_kind == FrozenPanelKind::AA &&
            panel.scope == FrozenTimingScope::EncoderInitPlusFirstKSymbols &&
            panel.left_arm == panel.right_arm &&
            IsCandidateArm(panel.left_arm);
    case 7u:
        return panel.panel_kind == FrozenPanelKind::AB &&
            panel.scope == FrozenTimingScope::DecoderSolve &&
            IsCandidateArm(panel.left_arm) &&
            panel.right_arm == "wirehair2_head";
    case 8u:
        return panel.panel_kind == FrozenPanelKind::AB &&
            panel.scope == FrozenTimingScope::ReceiveToSuccess &&
            IsCandidateArm(panel.left_arm) &&
            panel.right_arm == "wirehair1";
    case 9u:
        return panel.panel_kind == FrozenPanelKind::AB &&
            panel.scope == FrozenTimingScope::EncoderInitPlusFirstKSymbols &&
            IsCandidateArm(panel.left_arm) &&
            panel.right_arm == "wirehair2_head";
    case 10u:
        return panel.panel_kind == FrozenPanelKind::AB &&
            panel.scope == FrozenTimingScope::EncoderInitPlusFirstKSymbols &&
            IsCandidateArm(panel.left_arm) &&
            panel.right_arm == "wirehair1";
    default:
        return false;
    }
}

int HexNibble(char ch)
{
    if (ch >= '0' && ch <= '9') {
        return ch - '0';
    }
    if (ch >= 'a' && ch <= 'f') {
        return ch - 'a' + 10;
    }
    return -1;
}

} // namespace

FrozenTimingBaseCell::FrozenTimingBaseCell()
    : ordinal(0u)
    , K(0u)
    , block_bytes(0u)
    , loss_ppm(0u)
    , schedule(FrozenSchedule::Iid)
    , replicate(0u)
    , base_seed_attempt(0u)
    , base_loss_seed(0u)
    , fixed_received_overhead(0u)
    , receive_overhead_cap(0u)
    , interleave_policy()
    , invocations_per_slot(0u)
{
}

FrozenTimingCell::FrozenTimingCell()
    : ordinal(0u)
    , K(0u)
    , block_bytes(0u)
    , loss_ppm(0u)
    , schedule(FrozenSchedule::Iid)
    , replicate(0u)
    , base_seed_attempt(0u)
    , base_loss_seed(0u)
    , base_cell_sha256()
    , loss_retry_offset(0u)
    , loss_seed(0u)
    , fixed_received_overhead(0u)
    , receive_overhead_cap(0u)
    , interleave_policy()
    , invocations_per_slot(0u)
{
}

uint32_t DevelopmentTimingInvocationsPerSlot(uint32_t K)
{
    if (K == 0u) {
        return 0u;
    }
    const uint32_t quotient = 65536u / K;
    const uint32_t rounded_up = quotient + (65536u % K != 0u ? 1u : 0u);
    return std::max(2u, rounded_up);
}

uint32_t DevelopmentTimingWorkerCount()
{
    return kDevelopmentTimingWorkers;
}

std::size_t DevelopmentTimingCohortCount()
{
    return kDevelopmentTimingCohorts;
}

bool DevelopmentTimingWorkerSlot(
    std::size_t cohort_index,
    uint32_t replicate,
    uint32_t& worker_slot)
{
    worker_slot = 0u;
    if (cohort_index >= kDevelopmentTimingCohorts ||
        replicate >= kDevelopmentTimingRepetitions)
    {
        return false;
    }
    const uint32_t wave_index = replicate / kDevelopmentTimingWorkers;
    worker_slot = (
        replicate % kDevelopmentTimingWorkers +
        static_cast<uint32_t>(cohort_index % kDevelopmentTimingWorkers) +
        wave_index) % kDevelopmentTimingWorkers;
    return true;
}

bool DevelopmentTimingSeed(
    uint32_t replicate,
    uint32_t& base_seed_attempt,
    uint64_t& base_loss_seed)
{
    base_seed_attempt = 0u;
    base_loss_seed = 0u;
    if (replicate >= kDevelopmentTimingRepetitions) {
        return false;
    }
    const uint32_t pair_index = replicate % 3u;
    base_seed_attempt = 0u;
    const uint64_t salt =
        static_cast<uint64_t>(replicate) * kSplitMixIncrement;
    base_loss_seed = SplitMix64Once(kTrainingRoots[pair_index] ^ salt);
    return true;
}

std::vector<FrozenTimingBaseCell> EnumerateDevelopmentTimingBaseCells()
{
    std::vector<FrozenTimingBaseCell> cells;
    cells.reserve(2304u);
    for (uint32_t round_index = 0u;
        round_index < kDevelopmentTimingRepetitions /
            kDevelopmentTimingWorkers;
        ++round_index)
    {
        for (std::size_t width_index = 0u;
            width_index < sizeof(kTimingWidths) / sizeof(kTimingWidths[0]);
            ++width_index)
        {
            for (std::size_t k_index = 0u;
                k_index < sizeof(kTimingK) / sizeof(kTimingK[0]);
                ++k_index)
            {
                for (uint32_t lane = 0u;
                    lane < kDevelopmentTimingWorkers;
                    ++lane)
                {
                    const uint32_t replicate =
                        round_index * kDevelopmentTimingWorkers + lane;
                    uint32_t attempt = 0u;
                    uint64_t base_loss_seed = 0u;
                    if (!DevelopmentTimingSeed(
                            replicate, attempt, base_loss_seed))
                    {
                        return std::vector<FrozenTimingBaseCell>();
                    }
                    FrozenTimingBaseCell cell;
                    cell.ordinal = cells.size();
                    cell.phase = "development";
                    cell.band = BandForK(kTimingK[k_index]);
                    cell.K = kTimingK[k_index];
                    cell.block_bytes = kTimingWidths[width_index];
                    cell.loss_ppm = 100000u;
                    cell.schedule = FrozenSchedule::Iid;
                    cell.replicate = replicate;
                    cell.base_seed_attempt = attempt;
                    cell.base_loss_seed = base_loss_seed;
                    cell.fixed_received_overhead = 4u;
                    cell.receive_overhead_cap =
                        kDevelopmentTimingReceiveOverheadCap;
                    cell.interleave_policy =
                        kDevelopmentTimingInterleavePolicy;
                    cell.invocations_per_slot =
                        DevelopmentTimingInvocationsPerSlot(cell.K);
                    if (cell.invocations_per_slot == 0u) {
                        return std::vector<FrozenTimingBaseCell>();
                    }
                    cells.push_back(cell);
                }
            }
        }
    }
    return cells;
}

std::string CanonicalTimingBaseCellJson(const FrozenTimingBaseCell& cell)
{
    if (!ExactDevelopmentTimingBaseCell(cell)) {
        return std::string();
    }
    std::string json;
    json.reserve(320u);
    json += "{\"K\":";
    json += std::to_string(cell.K);
    json += ",\"band\":\"";
    json += cell.band;
    json += "\",\"base_loss_seed\":\"";
    json += HexSeed(cell.base_loss_seed);
    json += "\",\"base_seed_attempt\":";
    json += std::to_string(cell.base_seed_attempt);
    json += ",\"block_bytes\":";
    json += std::to_string(cell.block_bytes);
    json += ",\"fixed_received_overhead\":";
    json += std::to_string(cell.fixed_received_overhead);
    json += ",\"interleave_policy\":\"";
    json += cell.interleave_policy;
    json += "\"";
    json += ",\"invocations_per_slot\":";
    json += std::to_string(cell.invocations_per_slot);
    json += ",\"loss_ppm\":";
    json += std::to_string(cell.loss_ppm);
    json += ",\"phase\":\"development\",";
    json += "\"receive_overhead_cap\":";
    json += std::to_string(cell.receive_overhead_cap);
    json += ",\"replicate\":";
    json += std::to_string(cell.replicate);
    json += ",\"schedule\":\"iid\"}";
    return json;
}

std::string TimingBaseCellSha256(const FrozenTimingBaseCell& cell)
{
    const std::string json = CanonicalTimingBaseCellJson(cell);
    return json.empty() ? std::string() : Sha256Hex(json);
}

std::string DevelopmentTimingBaseDomainSha256()
{
    const std::vector<FrozenTimingBaseCell> cells =
        EnumerateDevelopmentTimingBaseCells();
    if (cells.size() != 2304u) {
        return std::string();
    }
    std::string stream;
    stream.reserve(700000u);
    for (std::size_t i = 0u; i < cells.size(); ++i)
    {
        const std::string json = CanonicalTimingBaseCellJson(cells[i]);
        if (json.empty()) {
            return std::string();
        }
        stream += json;
        stream += '\n';
    }
    return Sha256Hex(stream);
}

bool DevelopmentTimingLossSeed(
    uint64_t base_loss_seed,
    uint32_t loss_retry_offset,
    uint64_t& loss_seed)
{
    loss_seed = 0u;
    if (loss_retry_offset > 255u) {
        return false;
    }
    loss_seed = base_loss_seed +
        static_cast<uint64_t>(loss_retry_offset) * kSplitMixIncrement;
    return true;
}

bool QualifyDevelopmentTimingCell(
    const FrozenTimingBaseCell& base_cell,
    uint32_t loss_retry_offset,
    FrozenTimingCell& cell)
{
    cell = FrozenTimingCell();
    uint64_t loss_seed = 0u;
    if (!ExactDevelopmentTimingBaseCell(base_cell) ||
        !DevelopmentTimingLossSeed(
            base_cell.base_loss_seed, loss_retry_offset, loss_seed))
    {
        return false;
    }
    const std::string base_cell_sha256 = TimingBaseCellSha256(base_cell);
    if (!IsLowerSha256(base_cell_sha256)) {
        return false;
    }

    cell.ordinal = base_cell.ordinal;
    cell.phase = base_cell.phase;
    cell.band = base_cell.band;
    cell.K = base_cell.K;
    cell.block_bytes = base_cell.block_bytes;
    cell.loss_ppm = base_cell.loss_ppm;
    cell.schedule = base_cell.schedule;
    cell.replicate = base_cell.replicate;
    cell.base_seed_attempt = base_cell.base_seed_attempt;
    cell.base_loss_seed = base_cell.base_loss_seed;
    cell.base_cell_sha256 = base_cell_sha256;
    cell.loss_retry_offset = loss_retry_offset;
    cell.loss_seed = loss_seed;
    cell.fixed_received_overhead = base_cell.fixed_received_overhead;
    cell.receive_overhead_cap = base_cell.receive_overhead_cap;
    cell.interleave_policy = base_cell.interleave_policy;
    cell.invocations_per_slot = base_cell.invocations_per_slot;
    return true;
}

std::vector<FrozenTimingCell> EnumerateDevelopmentTimingCells(
    const std::vector<uint32_t>& retry_offsets)
{
    if (retry_offsets.size() != 2304u) {
        return std::vector<FrozenTimingCell>();
    }
    const std::vector<FrozenTimingBaseCell> base_cells =
        EnumerateDevelopmentTimingBaseCells();
    if (base_cells.size() != retry_offsets.size()) {
        return std::vector<FrozenTimingCell>();
    }
    std::vector<FrozenTimingCell> cells;
    cells.reserve(base_cells.size());
    for (std::size_t i = 0u; i < base_cells.size(); ++i)
    {
        FrozenTimingCell cell;
        if (!QualifyDevelopmentTimingCell(
                base_cells[i], retry_offsets[i], cell))
        {
            return std::vector<FrozenTimingCell>();
        }
        cells.push_back(cell);
    }
    return cells;
}

std::string CanonicalTimingCellJson(const FrozenTimingCell& cell)
{
    if (!ExactDevelopmentTimingCell(cell)) {
        return std::string();
    }
    std::string json;
    json.reserve(450u);
    json += "{\"K\":";
    json += std::to_string(cell.K);
    json += ",\"band\":\"";
    json += cell.band;
    json += "\",\"base_cell_sha256\":\"";
    json += cell.base_cell_sha256;
    json += "\",\"base_loss_seed\":\"";
    json += HexSeed(cell.base_loss_seed);
    json += "\",\"base_seed_attempt\":";
    json += std::to_string(cell.base_seed_attempt);
    json += ",\"block_bytes\":";
    json += std::to_string(cell.block_bytes);
    json += ",\"fixed_received_overhead\":";
    json += std::to_string(cell.fixed_received_overhead);
    json += ",\"interleave_policy\":\"";
    json += cell.interleave_policy;
    json += "\"";
    json += ",\"invocations_per_slot\":";
    json += std::to_string(cell.invocations_per_slot);
    json += ",\"loss_ppm\":";
    json += std::to_string(cell.loss_ppm);
    json += ",\"loss_retry_offset\":";
    json += std::to_string(cell.loss_retry_offset);
    json += ",\"loss_seed\":\"";
    json += HexSeed(cell.loss_seed);
    json += "\",\"phase\":\"development\",";
    json += "\"receive_overhead_cap\":";
    json += std::to_string(cell.receive_overhead_cap);
    json += ",\"replicate\":";
    json += std::to_string(cell.replicate);
    json += ",\"schedule\":\"iid\"}";
    return json;
}

std::string TimingCellSha256(const FrozenTimingCell& cell)
{
    const std::string json = CanonicalTimingCellJson(cell);
    return json.empty() ? std::string() : Sha256Hex(json);
}

std::string DevelopmentTimingDomainSha256(
    const std::vector<uint32_t>& retry_offsets)
{
    const std::vector<FrozenTimingCell> cells =
        EnumerateDevelopmentTimingCells(retry_offsets);
    if (cells.size() != 2304u) {
        return std::string();
    }
    std::string stream;
    stream.reserve(1000000u);
    for (std::size_t i = 0u; i < cells.size(); ++i)
    {
        const std::string json = CanonicalTimingCellJson(cells[i]);
        if (json.empty()) {
            return std::string();
        }
        stream += json;
        stream += '\n';
    }
    return Sha256Hex(stream);
}

std::string CanonicalTimingSourceJson(const FrozenTimingCell& cell)
{
    if (!ExactDevelopmentTimingCell(cell)) {
        return std::string();
    }
    return CanonicalTimingBaseCellJson(BaseCellFromQualified(cell));
}

std::string TimingSourceIdentitySha256(const FrozenTimingCell& cell)
{
    const std::string source_json = CanonicalTimingSourceJson(cell);
    if (source_json.empty()) {
        return std::string();
    }
    const std::string identity = Sha256Hex(source_json);
    return identity == cell.base_cell_sha256 ? identity : std::string();
}

FrozenTimingTraceReceipt::FrozenTimingTraceReceipt()
    : ordinal(0u)
    , K(0u)
    , attempted_candidates(0u)
    , candidate_limit(0u)
{
}

FrozenTraceStatus GenerateDevelopmentTimingTrace(
    const FrozenTimingCell& cell,
    FrozenPacketTrace& trace,
    FrozenTimingTraceReceipt& receipt)
{
    trace = FrozenPacketTrace();
    receipt = FrozenTimingTraceReceipt();
    if (!ExactDevelopmentTimingCell(cell)) {
        return FrozenTraceStatus::InvalidInput;
    }

    receipt.ordinal = cell.ordinal;
    receipt.K = cell.K;
    receipt.cell_sha256 = TimingCellSha256(cell);
    const FrozenTraceStatus status = GenerateFrozenPacketTrace(
        cell.K,
        cell.block_bytes,
        cell.loss_ppm,
        cell.loss_seed,
        FrozenSchedule::Iid,
        cell.receive_overhead_cap,
        trace);
    receipt.trace_sha256 = trace.trace_sha256;
    receipt.attempted_candidates = trace.attempted_candidates;
    receipt.candidate_limit = trace.candidate_limit;
    return status;
}

std::string CanonicalTimingTraceManifestRow(
    const FrozenTimingCell& cell,
    const FrozenTimingTraceReceipt& receipt)
{
    if (!ExactDevelopmentTimingCell(cell) ||
        receipt.ordinal != cell.ordinal)
    {
        return std::string();
    }
    FrozenPacketTrace expected_trace;
    FrozenTimingTraceReceipt expected_receipt;
    if (GenerateDevelopmentTimingTrace(
            cell, expected_trace, expected_receipt) !=
            FrozenTraceStatus::Complete ||
        receipt.K != expected_receipt.K ||
        receipt.cell_sha256 != expected_receipt.cell_sha256 ||
        receipt.trace_sha256 != expected_receipt.trace_sha256 ||
        receipt.attempted_candidates !=
            expected_receipt.attempted_candidates ||
        receipt.candidate_limit != expected_receipt.candidate_limit ||
        !IsLowerSha256(receipt.cell_sha256) ||
        !IsLowerSha256(receipt.trace_sha256))
    {
        return std::string();
    }

    std::string json;
    json.reserve(190u);
    json += "{\"cell_sha256\":\"";
    json += receipt.cell_sha256;
    json += "\",\"ordinal\":";
    json += std::to_string(receipt.ordinal);
    json += ",\"trace_sha256\":\"";
    json += receipt.trace_sha256;
    json += "\"}";
    return json;
}

const char* FrozenPanelKindName(FrozenPanelKind kind)
{
    switch (kind)
    {
    case FrozenPanelKind::AA: return "AA";
    case FrozenPanelKind::AB: return "AB";
    }
    return "";
}

const char* FrozenTimingScopeName(FrozenTimingScope scope)
{
    switch (scope)
    {
    case FrozenTimingScope::DecoderSolve:
        return "decoder_solve";
    case FrozenTimingScope::ReceiveToSuccess:
        return "receive_to_success";
    case FrozenTimingScope::EncoderInitPlusFirstKSymbols:
        return "encoder_init_plus_first_K_symbols";
    }
    return "";
}

FrozenTimingPanel::FrozenTimingPanel()
    : ordinal(0u)
    , panel_kind(FrozenPanelKind::AA)
    , scope(FrozenTimingScope::DecoderSolve)
{
}

std::vector<FrozenTimingPanel> EnumerateOneCandidateTimingPanels(
    const std::string& candidate_arm)
{
    if (!IsCandidateArm(candidate_arm)) {
        return std::vector<FrozenTimingPanel>();
    }

    std::vector<FrozenTimingPanel> panels;
    panels.reserve(11u);
    const auto add = [&panels](
        FrozenPanelKind kind,
        FrozenTimingScope scope,
        const std::string& left,
        const std::string& right)
    {
        FrozenTimingPanel panel;
        panel.ordinal = panels.size();
        panel.panel_kind = kind;
        panel.scope = scope;
        panel.left_arm = left;
        panel.right_arm = right;
        panels.push_back(panel);
    };

    add(FrozenPanelKind::AA, FrozenTimingScope::DecoderSolve,
        "wirehair2_head", "wirehair2_head");
    add(FrozenPanelKind::AA, FrozenTimingScope::ReceiveToSuccess,
        "wirehair1", "wirehair1");
    add(FrozenPanelKind::AA,
        FrozenTimingScope::EncoderInitPlusFirstKSymbols,
        "wirehair2_head", "wirehair2_head");
    add(FrozenPanelKind::AA,
        FrozenTimingScope::EncoderInitPlusFirstKSymbols,
        "wirehair1", "wirehair1");

    add(FrozenPanelKind::AA, FrozenTimingScope::DecoderSolve,
        candidate_arm, candidate_arm);
    add(FrozenPanelKind::AA, FrozenTimingScope::ReceiveToSuccess,
        candidate_arm, candidate_arm);
    add(FrozenPanelKind::AA,
        FrozenTimingScope::EncoderInitPlusFirstKSymbols,
        candidate_arm, candidate_arm);

    add(FrozenPanelKind::AB, FrozenTimingScope::DecoderSolve,
        candidate_arm, "wirehair2_head");
    add(FrozenPanelKind::AB, FrozenTimingScope::ReceiveToSuccess,
        candidate_arm, "wirehair1");
    add(FrozenPanelKind::AB,
        FrozenTimingScope::EncoderInitPlusFirstKSymbols,
        candidate_arm, "wirehair2_head");
    add(FrozenPanelKind::AB,
        FrozenTimingScope::EncoderInitPlusFirstKSymbols,
        candidate_arm, "wirehair1");
    return panels;
}

std::string CanonicalTimingPanelJson(const FrozenTimingPanel& panel)
{
    if (!ExactTimingPanel(panel)) {
        return std::string();
    }
    std::string json;
    json.reserve(180u);
    json += "{\"left_arm\":";
    if (!AppendJsonString(panel.left_arm, json)) {
        return std::string();
    }
    json += ",\"panel_kind\":\"";
    json += FrozenPanelKindName(panel.panel_kind);
    json += "\",\"right_arm\":";
    if (!AppendJsonString(panel.right_arm, json)) {
        return std::string();
    }
    json += ",\"scope\":\"";
    json += FrozenTimingScopeName(panel.scope);
    json += "\"}";
    return json;
}

std::string TimingPanelSha256(const FrozenTimingPanel& panel)
{
    const std::string json = CanonicalTimingPanelJson(panel);
    return json.empty() ? std::string() : Sha256Hex(json);
}

const char* FrozenTimingOrderName(FrozenTimingOrder order)
{
    switch (order)
    {
    case FrozenTimingOrder::ABBA: return "ABBA";
    case FrozenTimingOrder::BAAB: return "BAAB";
    case FrozenTimingOrder::Invalid: return "";
    }
    return "";
}

FrozenTimingOrder TimingPanelOrder(
    const FrozenTimingPanel& panel,
    uint32_t replicate)
{
    const std::string digest = TimingPanelSha256(panel);
    if (!IsLowerSha256(digest)) {
        return FrozenTimingOrder::Invalid;
    }
    const int low_nibble = HexNibble(digest[digest.size() - 1u]);
    if (low_nibble < 0) {
        return FrozenTimingOrder::Invalid;
    }
    const uint32_t phase_bit = static_cast<uint32_t>(low_nibble) & 1u;
    return (replicate & 1u) == phase_bit ?
        FrozenTimingOrder::ABBA : FrozenTimingOrder::BAAB;
}

} // namespace wh2_benchmark
} // namespace wirehair
