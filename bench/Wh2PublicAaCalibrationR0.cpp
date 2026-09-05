// The r1 worker is immutable. Reuse its public-API oracle/encoding helpers in
// this additive translation unit, but expose only the calibration entrypoint.
#define main Wh2PublicBorrowedR1UnusedMain
#include "Wh2PublicBorrowedCurrentVsWh1R1.cpp"
#undef main

namespace wh2_aa_calibration {

using wirehair::wh2_benchmark::Sha256Hex;
static const char kName[] = "wh2-public-aa-calibration-r0";
static const char kPrefix[] = "wirehair.wh2.public-aa-calibration-r0";
static const char kDirectory[] = "/var/tmp/wh2-public-aa-calibration-r0";
static const uint32_t kReplicates = 12u;
static const uint32_t kConditions[4][4] = {
    {0u, 1u, 3u, 2u}, {1u, 2u, 0u, 3u},
    {2u, 3u, 1u, 0u}, {3u, 0u, 2u, 1u}
};
static const CellShape kShapes[4] = {
    {8u, 64u}, {8u, 1280u}, {128u, 1280u}, {8192u, 1280u}
};
struct Tuple { uint32_t CellIndex; Scope ScopeValue; };
static const Tuple kTuples[8] = {
    {0u, Scope::PrebuiltSystematic}, {1u, Scope::PrebuiltSystematic},
    {1u, Scope::PrebuiltRepair}, {2u, Scope::PrebuiltSystematic},
    {2u, Scope::FreshRepair}, {3u, Scope::PrebuiltRepair},
    {3u, Scope::PrebuiltSystematic}, {2u, Scope::FreshSystematic}
};

uint32_t ScopeInvocations(const Tuple& tuple)
{
    return IsFresh(tuple.ScopeValue) ? 1u :
        (8192u + kShapes[tuple.CellIndex].K - 1u) / kShapes[tuple.CellIndex].K;
}

const char* Mapping(uint32_t condition) { return condition < 2u ? "direct" : "swapped"; }
const char* Order(uint32_t condition) { return condition % 2u == 0u ? "ABBA" : "BAAB"; }
const char* LogicalName(uint32_t lane) { return lane == 0u ? "left" : "right"; }
uint32_t PhysicalIndex(uint32_t logical, uint32_t condition)
{
    return logical ^ (condition / 2u);
}

std::array<uint32_t, 8> SlotOrder(uint32_t condition)
{
    const std::array<uint32_t, 8> abba = {{0u, 1u, 1u, 0u, 1u, 0u, 0u, 1u}};
    std::array<uint32_t, 8> result = abba;
    for (uint32_t& lane : result) lane ^= condition % 2u;
    return result;
}

std::string UnitKey(uint32_t replicate, uint32_t unit)
{
    return Sha256Hex("aa-cal-r0:rep=" + std::to_string(replicate) +
        ":unit=" + std::to_string(unit));
}

std::array<uint32_t, 16> UnitOrder(uint32_t replicate)
{
    std::array<std::pair<std::string, uint32_t>, 16> keyed;
    for (uint32_t unit = 0u; unit < keyed.size(); ++unit)
        keyed[unit] = std::make_pair(UnitKey(replicate, unit), unit);
    std::sort(keyed.begin(), keyed.end());
    std::array<uint32_t, 16> result;
    for (size_t i = 0u; i < result.size(); ++i) result[i] = keyed[i].second;
    return result;
}

std::string Describe()
{
    std::ostringstream out;
    out << "{\"campaign\":\"" << kName << "\",\"condition_sequences\":[";
    for (size_t row = 0u; row < 4u; ++row) {
        if (row) out << ',';
        out << '[';
        for (size_t col = 0u; col < 4u; ++col) {
            if (col) out << ',';
            out << kConditions[row][col];
        }
        out << ']';
    }
    out << "],\"conditions\":[";
    for (uint32_t condition = 0u; condition < 4u; ++condition) {
        if (condition) out << ',';
        const uint32_t first = condition % 2u;
        out << "{\"condition\":" << condition
            << ",\"mapping\":\"" << Mapping(condition)
            << "\",\"order\":\"" << Order(condition)
            << "\",\"warmup\":[\"" << LogicalName(first)
            << "\",\"" << LogicalName(first ^ 1u) << "\"]}";
    }
    out << "],\"expected_encode_calls\":47431680"
        << ",\"expected_measured_batches\":6144"
        << ",\"expected_measured_encode_calls\":37945344"
        << ",\"expected_measured_invocations\":2411520"
        << ",\"expected_panels\":768,\"expected_records\":770"
        << ",\"expected_warmup_batches\":1536"
        << ",\"expected_warmup_encode_calls\":9486336"
        << ",\"expected_warmup_invocations\":602880"
        << ",\"internal_deadline_seconds\":110,\"panel_replicates\":12"
        << ",\"roles\":[\"C\",\"L\"],\"schema\":\"" << kPrefix
        << ".describe.v1\",\"sibling_cpu\":56,\"source_seed_base\":\""
        << Hex64(kSourceSeedBase)
        << "\",\"target_cpu\":120,\"timing_granularity\":\"whole-batch\",\"tuples\":[";
    for (size_t i = 0u; i < 8u; ++i) {
        if (i) out << ',';
        const Tuple& tuple = kTuples[i];
        const CellShape& shape = kShapes[tuple.CellIndex];
        out << "{\"K\":" << shape.K << ",\"block_bytes\":" << shape.BlockBytes
            << ",\"scope\":\"" << ScopeName(tuple.ScopeValue)
            << "\",\"scope_invocations_per_batch\":" << ScopeInvocations(tuple) << '}';
    }
    out << "],\"unit_order\":\"sha256:aa-cal-r0:rep=R:unit=U\"}\n";
    return out.str();
}

class PhysicalLane
{
public:
    PhysicalLane(const ScreenCell& cell, Arm role, Scope scope, uint32_t invocations)
        : Cell(cell), Role(role), ScopeValue(scope), Invocations(invocations),
          Oracle(FindOracle(cell, role)) {}

    bool Initialize()
    {
        if (!Oracle || !PrepareScratch(Cell.Shape.K, Cell.Shape.BlockBytes, Scratch, Lengths))
            return false;
        if (IsFresh(ScopeValue)) return Invocations == 1u;
        return CreateEncoder(Role, Cell.Source, Cell.Shape.BlockBytes, Prebuilt, Artifacts) == 0 &&
            CheckProfile(Artifacts);
    }

    // Warmups and the portable selftest use the same work without any clocks.
    bool Batch(bool timed, uint64_t& elapsed)
    {
        std::fill(Scratch.begin(), Scratch.end(), kScratchCanary);
        std::fill(Lengths.begin(), Lengths.end(), kLengthCanary);
        EncoderHandle fresh;
        ConstructionArtifacts artifacts;
        EncoderHandle* encoder = &Prebuilt;
        int64_t result = 0;
        const uint32_t first_id = IsSystematic(ScopeValue) ? 0u : Cell.Shape.K;
        Clock::time_point start;
        if (timed) start = Clock::now();
        if (IsFresh(ScopeValue)) {
            result = CreateEncoder(Role, Cell.Source, Cell.Shape.BlockBytes, fresh, artifacts);
            encoder = &fresh;
        }
        for (uint32_t i = 0u; result == 0 && i < Invocations; ++i)
            result = EncodeBatch(*encoder, first_id, Cell.Shape.K, Cell.Shape.BlockBytes,
                Cell.Shape.BlockBytes, Scratch, Lengths);
        elapsed = timed ? PositiveNanoseconds(Clock::now() - start) : 0u;
        fresh.Reset(); // The constructor-inclusive estimand excludes free.
        return result == 0 && (!timed || elapsed != 0u) &&
            (!IsFresh(ScopeValue) || CheckProfile(artifacts));
    }

    bool Verify(std::string& hash) const
    {
        return VerifyScratch(Cell.Source, Cell.Shape.K, Cell.Shape.BlockBytes,
            Cell.Shape.BlockBytes, IsSystematic(ScopeValue), Scratch, Lengths, hash) &&
            hash == (IsSystematic(ScopeValue) ? Oracle->SystematicSha256 : Oracle->FullRepairSha256);
    }

    std::string Address() const
    {
        std::ostringstream out;
        out << "0x" << std::hex << reinterpret_cast<uintptr_t>(Scratch.data());
        return out.str();
    }

private:
    bool CheckProfile(const ConstructionArtifacts& artifacts) const
    {
        return Role == Arm::L ? artifacts.ProfileBytes == 0u :
            ValidateProfile(Cell.Source, Cell.Shape.BlockBytes, artifacts) &&
            Sha256Hex(artifacts.Profile.data(), artifacts.ProfileBytes) == Oracle->ProfileSha256;
    }

    const ScreenCell& Cell;
    Arm Role;
    Scope ScopeValue;
    uint32_t Invocations;
    const OracleReceipt* Oracle;
    EncoderHandle Prebuilt;
    ConstructionArtifacts Artifacts;
    std::vector<uint8_t> Scratch;
    std::vector<uint32_t> Lengths;
};

// Extract only the immutable oracle serialization. Its legacy per-replicate
// invocation field is an oracle-format receipt, never calibration scheduling.
std::string CellJson(const std::vector<std::shared_ptr<const ScreenCell> >& cells)
{
    std::ostringstream legacy;
    const std::string zero(64u, '0');
    if (!EmitConfig(legacy, cells, zero, "/synthetic/libwirehair.so.2.0.0",
            kTargetCpu, zero, "{}")) return std::string();
    const std::string text = legacy.str();
    const std::string begin_marker = "\"cells\":";
    const size_t begin = text.find(begin_marker);
    if (begin == std::string::npos) return std::string();
    const size_t first = begin + begin_marker.size();
    const size_t end = text.find(",\"comparisons\":", first);
    if (end == std::string::npos) return std::string();
    return text.substr(first, end - first);
}

bool EmitCalibrationConfig(const std::vector<std::shared_ptr<const ScreenCell> >& cells,
    const std::string& maps_hash, const std::string& library, const std::string& identity)
{
    const std::string cells_json = CellJson(cells);
    if (cells_json.empty()) return false;
    std::cout << "{\"campaign\":\"" << kName << "\",\"cells\":" << cells_json
        << ",\"compile\":{\"compiler_path\":\"" << WIREHAIR_WH2_CXX_COMPILER_PATH
        << "\",\"compiler_sha256\":\"" << WIREHAIR_WH2_CXX_COMPILER_SHA256
        << "\",\"compiler_version\":\"" << WIREHAIR_WH2_CMAKE_CXX_COMPILER_VERSION
        << "\",\"harness_git_commit\":\"" << WIREHAIR_WH2_SOURCE_GIT_COMMIT
        << "\",\"implementation_git_commit\":\"" << WIREHAIR_WH2_SOURCE_GIT_COMMIT
        << "\"},\"description_sha256\":\"" << Sha256Hex(Describe())
        << "\",\"runtime_library_maps_sha256\":\"" << maps_hash
        << "\",\"runtime_library_path\":\"" << library
        << "\",\"schema\":\"" << kPrefix << ".config.v1\",\"target_identity\":"
        << identity << "}\n";
    std::cout.flush();
    return std::cout.good();
}

bool Guard(const Clock::time_point& deadline)
{
    if (Clock::now() >= deadline || !VerifyTargetCpu(kTargetCpu)) {
        std::cerr << "AA calibration deadline or target CPU drift\n";
        return false;
    }
    return true;
}

bool RunCondition(const ScreenCell& cell, const Tuple& tuple, uint32_t replicate,
    uint32_t unit, uint32_t condition, uint32_t sequence, PhysicalLane& lane0,
    PhysicalLane& lane1, const std::string& maps_text, const Clock::time_point& deadline)
{
    PhysicalLane* lanes[2] = {&lane0, &lane1};
    uint64_t ignored = 0u;
    for (uint32_t i = 0u; i < 2u; ++i) {
        const uint32_t logical = (condition % 2u) ^ i;
        if (!Guard(deadline) || !lanes[PhysicalIndex(logical, condition)]->Batch(false, ignored) ||
            !Guard(deadline)) return false;
    }
    const std::array<uint32_t, 8> order = SlotOrder(condition);
    std::array<uint64_t, 8> elapsed = {{0u}};
    for (size_t slot = 0u; slot < order.size(); ++slot) {
        if (!Guard(deadline) ||
            !lanes[PhysicalIndex(order[slot], condition)]->Batch(true, elapsed[slot]) ||
            !Guard(deadline)) return false;
    }
    std::string hashes[2];
    if (!lane0.Verify(hashes[0]) || !lane1.Verify(hashes[1]) ||
        Sha256Hex(cell.Source.data(), cell.Source.size()) != cell.SourceSha256 ||
        RuntimeWirehairMapsText() != maps_text || !Guard(deadline)) return false;
    std::cout << "{\"K\":" << cell.Shape.K << ",\"block_bytes\":" << cell.Shape.BlockBytes
        << ",\"condition\":" << condition << ",\"description_sha256\":\"" << Sha256Hex(Describe())
        << "\",\"left_output_sha256\":\"" << hashes[PhysicalIndex(0u, condition)]
        << "\",\"mapping\":\"" << Mapping(condition) << "\",\"order\":\"" << Order(condition)
        << "\",\"physical_scratch_addresses\":[\"" << lane0.Address() << "\",\"" << lane1.Address()
        << "\"],\"replicate\":" << replicate
        << ",\"right_output_sha256\":\"" << hashes[PhysicalIndex(1u, condition)]
        << "\",\"role\":\"" << (unit % 2u == 0u ? "C" : "L")
        << "\",\"runtime_library_maps_sha256\":\"" << Sha256Hex(maps_text)
        << "\",\"schema\":\"" << kPrefix << ".panel.v1\",\"scope\":\"" << ScopeName(tuple.ScopeValue)
        << "\",\"scope_invocations_per_batch\":" << ScopeInvocations(tuple)
        << ",\"sequence\":" << sequence << ",\"slots\":[";
    for (size_t slot = 0u; slot < order.size(); ++slot) {
        if (slot) std::cout << ',';
        std::cout << "{\"elapsed_ns\":" << elapsed[slot]
            << ",\"logical_lane\":\"" << LogicalName(order[slot])
            << "\",\"physical_lane\":" << PhysicalIndex(order[slot], condition) << '}';
    }
    std::cout << "],\"source_immutable_verified\":true,\"target_cpu\":120,\"tuple_index\":" << unit / 2u
        << ",\"unit_index\":" << unit << ",\"unit_key_sha256\":\"" << UnitKey(replicate, unit) << "\"}\n";
    std::cout.flush();
    return std::cout.good();
}

bool CaptureIdentity(std::string& identity)
{
    wirehair_wh2_bench::TargetIdentityReceiptV2 receipt;
    std::string diagnostic;
    if (!wirehair_wh2_bench::CapturePublicBorrowedTargetIdentity(kTargetCpu, receipt, diagnostic)) {
        std::cerr << "AA calibration target identity: " << diagnostic << '\n';
        return false;
    }
    identity = TargetIdentityJson(receipt, diagnostic);
    return !identity.empty();
}

bool Run(const Clock::time_point& deadline)
{
    std::string diagnostic, identity, later_identity;
    if (!PinAndVerifyTargetCpu(kTargetCpu, diagnostic) || !Guard(deadline) ||
        !CaptureIdentity(identity) || wirehair_init() != Wirehair_Success) return false;
    std::string maps_text, library, binding;
    if (!RuntimeWirehairMapsReceipt(maps_text, library) ||
        !RuntimeTimedBindingPath(binding) || binding != library) return false;
    std::vector<std::shared_ptr<const ScreenCell> > cells;
    for (const CellShape& shape : kShapes) {
        std::shared_ptr<ScreenCell> cell(new ScreenCell);
        if (!Guard(deadline) || !BuildCell(shape, *cell) || !Guard(deadline)) return false;
        cells.push_back(cell);
    }
    if (!CaptureIdentity(later_identity) || identity != later_identity ||
        RuntimeWirehairMapsText() != maps_text || !Guard(deadline) ||
        !EmitCalibrationConfig(cells, Sha256Hex(maps_text), library, identity)) return false;
    uint32_t panels = 0u;
    uint64_t measured_invocations = 0u, warmup_invocations = 0u, encode_calls = 0u;
    for (uint32_t replicate = 0u; replicate < kReplicates; ++replicate) {
        for (uint32_t unit : UnitOrder(replicate)) {
            const Tuple& tuple = kTuples[unit / 2u];
            const ScreenCell& cell = *cells[tuple.CellIndex];
            const Arm role = unit % 2u == 0u ? Arm::C : Arm::L;
            const uint32_t count = ScopeInvocations(tuple);
            PhysicalLane lane0(cell, role, tuple.ScopeValue, count);
            PhysicalLane lane1(cell, role, tuple.ScopeValue, count);
            if (!Guard(deadline) || !lane0.Initialize() || !lane1.Initialize() ||
                lane0.Address() == lane1.Address()) return false;
            for (uint32_t condition : kConditions[replicate % 4u]) {
                if (!RunCondition(cell, tuple, replicate, unit, condition, panels,
                        lane0, lane1, maps_text, deadline)) return false;
                ++panels;
                measured_invocations += 8u * count;
                warmup_invocations += 2u * count;
                encode_calls += UINT64_C(10) * count * cell.Shape.K;
            }
        }
    }
    for (const std::shared_ptr<const ScreenCell>& cell : cells)
        if (Sha256Hex(cell->Source.data(), cell->Source.size()) != cell->SourceSha256) return false;
    if (panels != 768u || measured_invocations != UINT64_C(2411520) ||
        warmup_invocations != UINT64_C(602880) || encode_calls != UINT64_C(47431680) ||
        !CaptureIdentity(later_identity) || identity != later_identity ||
        RuntimeWirehairMapsText() != maps_text || !Guard(deadline)) return false;
    std::cout << "{\"campaign\":\"" << kName << "\",\"encode_call_count\":" << encode_calls
        << ",\"measured_batch_count\":6144,\"measured_invocation_count\":" << measured_invocations
        << ",\"panel_count\":" << panels << ",\"record_count\":770,\"schema\":\"" << kPrefix
        << ".terminal.v1\",\"status\":\"complete\",\"warmup_batch_count\":1536"
        << ",\"warmup_invocation_count\":" << warmup_invocations << "}\n";
    std::cout.flush();
    return std::cout.good() && Guard(deadline);
}

#if defined(__linux__)
class Fd
{
public:
    explicit Fd(int fd) : Value(fd) {}
    ~Fd() { if (Value >= 0) ::close(Value); }
    Fd(const Fd&) = delete;
    Fd& operator=(const Fd&) = delete;
    int Value;
};

bool Authorize(const char* directory, const char* claim_hash, const char* worker_hash,
    const char* description_hash, std::string& diagnostic)
{
    if (std::strcmp(directory, kDirectory) != 0 || !IsLowerSha256(claim_hash) ||
        !IsLowerSha256(worker_hash) || !IsLowerSha256(description_hash) ||
        Sha256Hex(Describe()) != description_hash) return false;
    Fd dir(::open(directory, O_RDONLY | O_DIRECTORY | O_CLOEXEC | O_NOFOLLOW));
    struct stat info;
    if (dir.Value < 0 || ::fstat(dir.Value, &info) != 0 || !S_ISDIR(info.st_mode) ||
        (info.st_mode & 07777) != 0700) return false;
    Fd claim(::openat(dir.Value, "CLAIM", O_RDONLY | O_CLOEXEC | O_NOFOLLOW));
    Fd executable(::open("/proc/self/exe", O_RDONLY | O_CLOEXEC));
    std::string bytes;
    if (claim.Value < 0 || executable.Value < 0 ||
        !VerifyRegularAuthorizationFile(claim.Value, 0400, 1024u * 1024u, bytes, diagnostic) ||
        Sha256Hex(bytes) != claim_hash ||
        !VerifyRegularAuthorizationFile(executable.Value, 0, 64u * 1024u * 1024u, bytes, diagnostic) ||
        Sha256Hex(bytes) != worker_hash) return false;
    Fd marker(::openat(dir.Value, "WORKER_STARTED",
        O_WRONLY | O_CREAT | O_EXCL | O_CLOEXEC | O_NOFOLLOW, 0400));
    if (marker.Value < 0 || ::fchmod(marker.Value, 0400) != 0) return false;
    std::ostringstream json;
    json << "{\"campaign\":\"" << kName << "\",\"claim_sha256\":\"" << claim_hash
        << "\",\"description_sha256\":\"" << description_hash << "\",\"schema\":\""
        << kPrefix << ".worker-started.v1\",\"source_commit\":\"" << WIREHAIR_WH2_SOURCE_GIT_COMMIT
        << "\",\"status\":\"started\",\"worker_sha256\":\"" << worker_hash << "\"}\n";
    return WriteFdAll(marker.Value, json.str(), diagnostic) &&
        ::fsync(marker.Value) == 0 && ::fsync(dir.Value) == 0;
}
#else
bool Authorize(const char*, const char*, const char*, const char*, std::string&) { return false; }
#endif

bool Selftest()
{
    uint64_t measured = 0u, warmup = 0u, calls = 0u;
    uint32_t panels = 0u;
    for (uint32_t replicate = 0u; replicate < kReplicates; ++replicate) {
        const std::array<uint32_t, 16> order = UnitOrder(replicate);
        if (std::set<uint32_t>(order.begin(), order.end()).size() != 16u) return false;
        std::string last;
        for (uint32_t unit : order) {
            const std::string key = UnitKey(replicate, unit);
            if (!last.empty() && key <= last) return false;
            last = key;
            const Tuple& tuple = kTuples[unit / 2u];
            for (uint32_t condition : kConditions[replicate % 4u]) {
                uint32_t same[2] = {0u, 0u}, other[2] = {0u, 0u};
                uint32_t prior = (condition % 2u) ^ 1u;
                for (uint32_t logical : SlotOrder(condition)) {
                    ++(logical == prior ? same[logical] : other[logical]);
                    prior = logical;
                }
                if (same[0] != 1u || same[1] != 1u || other[0] != 3u || other[1] != 3u ||
                    PhysicalIndex(0u, condition) == PhysicalIndex(1u, condition)) return false;
                ++panels;
                measured += UINT64_C(8) * ScopeInvocations(tuple);
                warmup += UINT64_C(2) * ScopeInvocations(tuple);
                calls += UINT64_C(10) * ScopeInvocations(tuple) * kShapes[tuple.CellIndex].K;
            }
        }
    }
    if (panels != 768u || measured != UINT64_C(2411520) || warmup != UINT64_C(602880) ||
        calls != UINT64_C(47431680) || wirehair_init() != Wirehair_Success) return false;
    ScreenCell small;
    if (!BuildCell(kShapes[0], small)) return false;
    for (Arm role : {Arm::C, Arm::L}) {
        for (Scope scope : kScopes) {
            PhysicalLane lane(small, role, scope, IsFresh(scope) ? 1u : 2u);
            uint64_t ignored = 1u;
            std::string hash;
            if (!lane.Initialize() || !lane.Batch(false, ignored) || ignored != 0u ||
                !lane.Verify(hash)) return false;
        }
    }
    return Sha256Hex(small.Source.data(), small.Source.size()) == small.SourceSha256 &&
        Describe().back() == '\n';
}

} // namespace wh2_aa_calibration

int main(int argc, char** argv)
{
    using namespace wh2_aa_calibration;
    try {
        if (argc == 2 && std::strcmp(argv[1], "--describe") == 0) {
            std::cout << Describe();
            return std::cout.good() ? 0 : 1;
        }
        if (argc == 2 && std::strcmp(argv[1], "--selftest") == 0) {
            if (!Selftest()) { std::cerr << "AA calibration selftest failed\n"; return 1; }
            std::cout << "WH2 public AA calibration r0 selftest passed\n";
            return 0;
        }
        const Clock::time_point deadline = Clock::now() + std::chrono::seconds(110);
        if (argc != 12 || std::strcmp(argv[1], "--run") != 0 ||
            std::strcmp(argv[2], "--cpu") != 0 || std::strcmp(argv[3], "120") != 0 ||
            std::strcmp(argv[4], "--evidence-dir") != 0 ||
            std::strcmp(argv[6], "--claim-sha256") != 0 ||
            std::strcmp(argv[8], "--worker-sha256") != 0 ||
            std::strcmp(argv[10], "--config-identity-sha256") != 0) {
            std::cerr << "AA calibration expects --describe, --selftest, or exact authorized --run arguments\n";
            return 1;
        }
        std::string diagnostic;
        if (!Authorize(argv[5], argv[7], argv[9], argv[11], diagnostic)) {
            std::cerr << "AA calibration launch authorization failed: " << diagnostic << '\n';
            return 1;
        }
        if (!Run(deadline)) { std::cerr << "AA calibration failed\n"; return 1; }
        return 0;
    } catch (const std::exception& exc) {
        std::cerr << "AA calibration exception: " << exc.what() << '\n';
    } catch (...) { std::cerr << "AA calibration unknown exception\n"; }
    return 1;
}
