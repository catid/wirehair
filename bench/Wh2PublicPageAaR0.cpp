// Additive page-placement diagnostic. The r1 worker and public codec stay frozen.
#define main Wh2PublicBorrowedR1UnusedMain
#include "Wh2PublicBorrowedCurrentVsWh1R1.cpp"
#undef main
#include <new>

namespace wh2_page_aa {

using wirehair::wh2_benchmark::Sha256Hex;
static const char kName[] = "wh2-public-page-aa-r0";
static const char kPrefix[] = "wirehair.wh2.public-page-aa-r0";
static const char kDirectory[] = "/var/tmp/wh2-public-page-aa-r0";
static const uint32_t kReplicates = 12u;
static const uint32_t kConditions[4][4] = {
    {0u, 1u, 3u, 2u}, {1u, 2u, 0u, 3u},
    {2u, 3u, 1u, 0u}, {3u, 0u, 2u, 1u}
};
static const CellShape kShapes[3] = {{8u, 64u}, {8u, 1280u}, {128u, 1280u}};
static const size_t kPhases[4][2] = {
    {0x920u, 0x920u}, {0x820u, 0x820u},
    {0x920u, 0x820u}, {0x920u, 0x920u}
};

std::string Address(uintptr_t value)
{
    std::ostringstream out;
    out << "0x" << std::hex << value;
    return out.str();
}
const char* RoleName(uint32_t role) { return role == 0u ? "C" : role == 1u ? "L" : "M"; }
const char* Mapping(uint32_t condition) { return condition < 2u ? "direct" : "swapped"; }
const char* Order(uint32_t condition) { return condition % 2u == 0u ? "ABBA" : "BAAB"; }
const char* LogicalName(uint32_t lane) { return lane == 0u ? "left" : "right"; }
uint32_t PhysicalIndex(uint32_t logical, uint32_t condition) { return logical ^ (condition / 2u); }
std::array<uint32_t, 8> SlotOrder(uint32_t condition)
{
    std::array<uint32_t, 8> result = {{0u, 1u, 1u, 0u, 1u, 0u, 0u, 1u}};
    for (uint32_t& lane : result) lane ^= condition % 2u;
    return result;
}
std::string UnitKey(uint32_t replicate, uint32_t unit)
{
    return Sha256Hex("page-aa-r0:rep=" + std::to_string(replicate) +
        ":unit=" + std::to_string(unit));
}
std::array<uint32_t, 9> UnitOrder(uint32_t replicate)
{
    std::array<std::pair<std::string, uint32_t>, 9> keyed;
    for (uint32_t unit = 0u; unit < keyed.size(); ++unit)
        keyed[unit] = std::make_pair(UnitKey(replicate, unit), unit);
    std::sort(keyed.begin(), keyed.end());
    std::array<uint32_t, 9> result;
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
        out << "{\"condition\":" << condition << ",\"mapping\":\"" << Mapping(condition)
            << "\",\"order\":\"" << Order(condition) << "\",\"warmup\":[\""
            << LogicalName(first) << "\",\"" << LogicalName(first ^ 1u) << "\"]}";
    }
    out << "],\"counter_phase\":1024,\"expected_encode_calls\":141557760"
        << ",\"expected_measured_batches\":13824,\"expected_measured_encode_calls\":113246208"
        << ",\"expected_measured_invocations\":9732096,\"expected_panels\":1728"
        << ",\"expected_records\":1730,\"expected_warmup_batches\":3456"
        << ",\"expected_warmup_encode_calls\":28311552,\"expected_warmup_invocations\":2433024"
        << ",\"internal_deadline_seconds\":110,\"page_bytes\":4096,\"panel_replicates\":12"
        << ",\"primary_roles\":[\"C\",\"L\"],\"primary_scenario\":0,\"roles\":[\"C\",\"L\",\"M\"]"
        << ",\"scenario_order\":\"(replicate+unit+index)%4\",\"scenarios\":[";
    for (uint32_t scenario = 0u; scenario < 4u; ++scenario) {
        if (scenario) out << ',';
        out << "{\"output_phases\":[" << kPhases[scenario][0] << ',' << kPhases[scenario][1]
            << "],\"scenario\":" << scenario << ",\"shared_handle\":"
            << (scenario == 3u ? "true" : "false") << '}';
    }
    out << "],\"schema\":\"" << kPrefix << ".describe.v1\",\"sibling_cpu\":56"
        << ",\"source_phase\":2048,\"source_seed_base\":\"" << Hex64(kSourceSeedBase)
        << "\",\"target_cpu\":120,\"timing_granularity\":\"whole-batch\",\"tuples\":[";
    for (size_t i = 0u; i < 3u; ++i) {
        if (i) out << ',';
        out << "{\"K\":" << kShapes[i].K << ",\"block_bytes\":" << kShapes[i].BlockBytes
            << ",\"scope\":\"prebuilt-systematic\",\"scope_invocations_per_batch\":"
            << 8192u / kShapes[i].K << '}';
    }
    out << "],\"unit_order\":\"sha256:page-aa-r0:rep=R:unit=U\"}\n";
    return out.str();
}

struct CopyBinding
{
    typedef void* (*FunctionType)(void*, const void*, size_t);
    FunctionType Function = nullptr;
    uintptr_t Location = 0u;
    uintptr_t ElfOffset = 0u;
    std::string Path;

    bool Resolve()
    {
#if defined(__linux__)
        (void)dlerror();
        void* symbol = dlsym(RTLD_DEFAULT, "memcpy");
        const char* error = dlerror();
        Dl_info info = {};
        if (error || !symbol || dladdr(symbol, &info) == 0 || !info.dli_fname ||
            !CanonicalExistingPath(info.dli_fname, Path)) return false;
        static_assert(sizeof(Function) == sizeof(symbol), "POSIX function pointer size");
        std::memcpy(&Function, &symbol, sizeof(Function));
        Location = reinterpret_cast<uintptr_t>(symbol);
        const uintptr_t base = reinterpret_cast<uintptr_t>(info.dli_fbase);
        if (!base || Location < base) return false;
        ElfOffset = Location - base;
#else
        // Portable untimed selftest only. The authorized run is Linux-only.
        Function = &std::memcpy;
        Path = "portable-selftest";
#endif
        return Function != nullptr;
    }
    std::string Json() const
    {
        return "{\"address\":\"" + Address(Location) + "\",\"elf_offset\":\"" +
            Address(ElfOffset) + "\",\"path\":\"" + Path + "\"}";
    }
    bool Stable() const
    {
        CopyBinding current;
        return current.Resolve() && current.Json() == Json();
    }
};

// A fixed typed object gives the arena's counters a C++11 lifetime, not just alignment.
struct CounterBlock { uint32_t Values[128]; };
class Arena
{
public:
    ~Arena() { for (CounterBlock* counter : Counters) if (counter) counter->~CounterBlock(); }
    bool Initialize(const ScreenCell& cell)
    {
        Bytes = cell.Source.size();
        if (cell.Shape.K > 128u || Bytes != static_cast<size_t>(cell.Shape.K) * cell.Shape.BlockBytes)
            return false;
        Span = ((Bytes + 4096u + 4095u) / 4096u) * 4096u;
        Storage.reset(new uint8_t[5u * Span + 4095u]);
        const uintptr_t raw = reinterpret_cast<uintptr_t>(Storage.get());
        Base = reinterpret_cast<uint8_t*>((raw + 4095u) & ~uintptr_t(4095u));
        std::fill(Base, Base + 5u * Span, kScratchCanary);
        Source = Base + 0x800u;
        std::memcpy(Source, cell.Source.data(), Bytes);
        for (size_t lane = 0u; lane < 2u; ++lane) {
            void* address = Base + (3u + lane) * Span + 0x400u;
            Counters[lane] = new (address) CounterBlock;
            std::fill(Counters[lane]->Values, Counters[lane]->Values + 128u, kLengthCanary);
        }
        return Select(0u);
    }
    bool Select(uint32_t scenario)
    {
        if (!Base || scenario >= 4u) return false;
        // Clear prior scenarios' output footprints, never the attached source.
        std::fill(Base + Span, Base + 3u * Span, kScratchCanary);
        for (size_t lane = 0u; lane < 2u; ++lane)
            Outputs[lane] = Base + (1u + lane) * Span + kPhases[scenario][lane];
        return Source + Bytes <= Outputs[0] && Outputs[0] + Bytes <= Outputs[1] &&
            Outputs[1] + Bytes <= reinterpret_cast<uint8_t*>(Counters[0]) &&
            reinterpret_cast<uint8_t*>(Counters[0]) + sizeof(CounterBlock) <=
                reinterpret_cast<uint8_t*>(Counters[1]);
    }
    bool Gap(const uint8_t* begin, const uint8_t* end) const
    {
        return begin <= end && std::all_of(begin, end, [](uint8_t value) { return value == kScratchCanary; });
    }
    bool Verify(const ScreenCell& cell, std::string hashes[2]) const
    {
        if (Sha256Hex(Source, Bytes) != cell.SourceSha256 ||
            std::memcmp(Source, cell.Source.data(), Bytes) != 0 ||
            !Gap(Base, Source) || !Gap(Source + Bytes, Outputs[0]) ||
            !Gap(Outputs[0] + Bytes, Outputs[1]) ||
            !Gap(Outputs[1] + Bytes, reinterpret_cast<uint8_t*>(Counters[0]))) return false;
        for (size_t lane = 0u; lane < 2u; ++lane) {
            if (std::memcmp(Outputs[lane], Source, Bytes) != 0) return false;
            hashes[lane] = Sha256Hex(Outputs[lane], Bytes);
            if (hashes[lane] != cell.Oracles[0].SystematicSha256 ||
                hashes[lane] != cell.Oracles[1].SystematicSha256) return false;
            for (uint32_t i = 0u; i < 128u; ++i)
                if (Counters[lane]->Values[i] != (i < cell.Shape.K ? cell.Shape.BlockBytes : kLengthCanary))
                    return false;
            const uint8_t* next = lane == 0u ? reinterpret_cast<uint8_t*>(Counters[1]) : Base + 5u * Span;
            if (!Gap(reinterpret_cast<uint8_t*>(Counters[lane]) + sizeof(CounterBlock), next)) return false;
        }
        return true;
    }
    std::unique_ptr<uint8_t[]> Storage;
    uint8_t* Base = nullptr;
    uint8_t* Source = nullptr;
    uint8_t* Outputs[2] = {nullptr, nullptr};
    CounterBlock* Counters[2] = {nullptr, nullptr};
    size_t Bytes = 0u;
    size_t Span = 0u;
};

int64_t CreateAt(uint32_t role, const uint8_t* source, size_t bytes, uint32_t block_bytes,
    EncoderHandle& encoder, ConstructionArtifacts& artifacts)
{
    if (role == 0u) {
        WirehairV2EncoderOptions options = {};
        options.struct_bytes = sizeof(options);
        options.options_version = WIREHAIR_V2_ENCODER_OPTIONS_VERSION;
        options.source_policy = WirehairV2EncoderSource_BorrowedImmutable;
        WirehairV2Codec codec = nullptr;
        const WirehairV2Result result = wirehair_v2_encoder_create_with_options(
            source, bytes, block_bytes, &options, artifacts.Profile.data(),
            static_cast<uint32_t>(artifacts.Profile.size()), &artifacts.ProfileBytes, &codec);
        if (result != WirehairV2_Success || !codec) {
            wirehair_v2_free(codec);
            return result == WirehairV2_Success ? WirehairV2_Error : result;
        }
        encoder.Adopt(codec);
        return 0;
    }
    WirehairCodec codec = nullptr;
    const WirehairResult result = wirehair_encoder_create_ex(nullptr, source, bytes, block_bytes, &codec);
    if (result != Wirehair_Success || !codec) {
        wirehair_free(codec);
        return result == Wirehair_Success ? Wirehair_Error : result;
    }
    encoder.Adopt(codec);
    return 0;
}

class PhysicalPair
{
public:
    PhysicalPair(const ScreenCell& cell, uint32_t role, const CopyBinding& copy)
        : Cell(cell), Role(role), Copy(copy) {}
    bool Initialize()
    {
        if (Role > 2u || !Memory.Initialize(Cell)) return false;
        if (Role == 2u) return true;
        for (size_t lane = 0u; lane < 2u; ++lane) {
            if (CreateAt(Role, Memory.Source, Memory.Bytes, Cell.Shape.BlockBytes,
                    Handles[lane], Artifacts[lane]) != 0) return false;
            if (Role == 0u && (!ValidateProfile(Cell.Source, Cell.Shape.BlockBytes, Artifacts[lane]) ||
                    Sha256Hex(Artifacts[lane].Profile.data(), Artifacts[lane].ProfileBytes) !=
                        Cell.Oracles[0].ProfileSha256)) return false;
        }
        return HandleAddress(0u) != HandleAddress(1u);
    }
    bool Select(uint32_t scenario) { Scenario = scenario; return Memory.Select(scenario); }
    uintptr_t HandleAddress(uint32_t lane) const
    {
        if (Role == 2u) return 0u;
        const EncoderHandle& handle = Handles[Scenario == 3u ? 0u : lane];
        return Role == 0u ? reinterpret_cast<uintptr_t>(handle.V2Codec()) :
            reinterpret_cast<uintptr_t>(handle.LegacyCodec());
    }
    bool Batch(uint32_t lane, bool timed, uint32_t invocations, uint64_t& elapsed)
    {
        if (lane > 1u || invocations == 0u) return false;
        uint8_t* const output = Memory.Outputs[lane];
        uint32_t* const lengths = Memory.Counters[lane]->Values;
        std::fill(output, output + Memory.Bytes, kScratchCanary);
        std::fill(lengths, lengths + Cell.Shape.K, kLengthCanary);
        const EncoderHandle& handle = Handles[Scenario == 3u ? 0u : lane];
        int64_t result = 0;
        Clock::time_point start;
        if (timed) start = Clock::now();
        for (uint32_t invocation = 0u; result == 0 && invocation < invocations; ++invocation) {
            for (uint32_t id = 0u; id < Cell.Shape.K; ++id) {
                const size_t offset = static_cast<size_t>(id) * Cell.Shape.BlockBytes;
                if (Role == 2u) {
                    if (Copy.Function(output + offset, Memory.Source + offset, Cell.Shape.BlockBytes) != output + offset)
                        result = WirehairV2_Error;
                    lengths[id] = Cell.Shape.BlockBytes;
                } else {
                    result = EncodeOne(handle, id, output + offset, Cell.Shape.BlockBytes, lengths[id]);
                }
                if (result != 0 || lengths[id] != Cell.Shape.BlockBytes) {
                    if (result == 0) result = WirehairV2_Error;
                    break;
                }
            }
        }
        elapsed = timed ? PositiveNanoseconds(Clock::now() - start) : 0u;
        return result == 0 && (!timed || elapsed != 0u);
    }
    std::string Addresses() const
    {
        std::ostringstream out;
        out << "{\"arena\":\"" << Address(reinterpret_cast<uintptr_t>(Memory.Base))
            << "\",\"counters\":[\"" << Address(reinterpret_cast<uintptr_t>(Memory.Counters[0]))
            << "\",\"" << Address(reinterpret_cast<uintptr_t>(Memory.Counters[1]))
            << "\"],\"handles\":[\"" << Address(HandleAddress(0u)) << "\",\"" << Address(HandleAddress(1u))
            << "\"],\"outputs\":[\"" << Address(reinterpret_cast<uintptr_t>(Memory.Outputs[0]))
            << "\",\"" << Address(reinterpret_cast<uintptr_t>(Memory.Outputs[1]))
            << "\"],\"source\":\"" << Address(reinterpret_cast<uintptr_t>(Memory.Source))
            << "\",\"span\":" << Memory.Span << '}';
        return out.str();
    }
    bool Verify(std::string hashes[2]) const { return Memory.Verify(Cell, hashes); }
    Arena Memory; // Declared before handles: borrowed source outlives both destructors.
private:
    EncoderHandle Handles[2];
    ConstructionArtifacts Artifacts[2];
    const ScreenCell& Cell;
    uint32_t Role;
    const CopyBinding& Copy;
    uint32_t Scenario = 0u;
};

std::string CellJson(const std::vector<std::shared_ptr<const ScreenCell> >& cells)
{
    std::ostringstream legacy;
    const std::string zero(64u, '0');
    if (!EmitConfig(legacy, cells, zero, "/synthetic/libwirehair.so.2.0.0", kTargetCpu, zero, "{}"))
        return std::string();
    const std::string text = legacy.str(), marker = "\"cells\":";
    const size_t begin = text.find(marker);
    if (begin == std::string::npos) return std::string();
    const size_t first = begin + marker.size(), end = text.find(",\"comparisons\":", first);
    return end == std::string::npos ? std::string() : text.substr(first, end - first);
}
bool EmitPageConfig(const std::vector<std::shared_ptr<const ScreenCell> >& cells,
    const std::string& maps_hash, const std::string& library, const std::string& identity, const CopyBinding& copy)
{
    const std::string cells_json = CellJson(cells);
    if (cells_json.empty()) return false;
    std::cout << "{\"campaign\":\"" << kName << "\",\"cells\":" << cells_json
        << ",\"compile\":{\"compiler_path\":\"" << WIREHAIR_WH2_CXX_COMPILER_PATH
        << "\",\"compiler_sha256\":\"" << WIREHAIR_WH2_CXX_COMPILER_SHA256
        << "\",\"compiler_version\":\"" << WIREHAIR_WH2_CMAKE_CXX_COMPILER_VERSION
        << "\",\"harness_git_commit\":\"" << WIREHAIR_WH2_SOURCE_GIT_COMMIT
        << "\",\"implementation_git_commit\":\"" << WIREHAIR_WH2_SOURCE_GIT_COMMIT
        << "\"},\"copy_binding\":" << copy.Json() << ",\"description_sha256\":\"" << Sha256Hex(Describe())
        << "\",\"runtime_library_maps_sha256\":\"" << maps_hash << "\",\"runtime_library_path\":\"" << library
        << "\",\"schema\":\"" << kPrefix << ".config.v1\",\"target_identity\":" << identity << "}\n";
    std::cout.flush();
    return std::cout.good();
}
bool Guard(const Clock::time_point& deadline)
{
    if (Clock::now() >= deadline || !VerifyTargetCpu(kTargetCpu)) {
        std::cerr << "Page AA deadline or target CPU drift\n";
        return false;
    }
    return true;
}
bool CaptureIdentity(std::string& identity)
{
    wirehair_wh2_bench::TargetIdentityReceiptV2 receipt;
    std::string diagnostic;
    if (!wirehair_wh2_bench::CapturePublicBorrowedTargetIdentity(kTargetCpu, receipt, diagnostic)) return false;
    identity = TargetIdentityJson(receipt, diagnostic);
    return !identity.empty();
}
bool RunCondition(const ScreenCell& cell, uint32_t replicate, uint32_t unit, uint32_t scenario,
    uint32_t condition, uint32_t sequence, PhysicalPair& pair, const std::string& maps_text,
    const CopyBinding& copy, const Clock::time_point& deadline)
{
    const uint32_t count = 8192u / cell.Shape.K;
    uint64_t ignored = 0u;
    for (uint32_t i = 0u; i < 2u; ++i)
        if (!Guard(deadline) || !pair.Batch(PhysicalIndex((condition % 2u) ^ i, condition), false, count, ignored) ||
            !Guard(deadline)) return false;
    const std::array<uint32_t, 8> order = SlotOrder(condition);
    std::array<uint64_t, 8> elapsed = {{0u}};
    for (size_t slot = 0u; slot < order.size(); ++slot)
        if (!Guard(deadline) || !pair.Batch(PhysicalIndex(order[slot], condition), true, count, elapsed[slot]) ||
            !Guard(deadline)) return false;
    std::string hashes[2];
    if (!pair.Verify(hashes) || Sha256Hex(cell.Source.data(), cell.Source.size()) != cell.SourceSha256 ||
        RuntimeWirehairMapsText() != maps_text || !copy.Stable() || !Guard(deadline)) return false;
    std::cout << "{\"K\":" << cell.Shape.K << ",\"addresses\":" << pair.Addresses()
        << ",\"block_bytes\":" << cell.Shape.BlockBytes << ",\"condition\":" << condition
        << ",\"copy_binding_sha256\":\"" << Sha256Hex(copy.Json())
        << "\",\"description_sha256\":\"" << Sha256Hex(Describe())
        << "\",\"left_output_sha256\":\"" << hashes[PhysicalIndex(0u, condition)]
        << "\",\"mapping\":\"" << Mapping(condition) << "\",\"order\":\"" << Order(condition)
        << "\",\"replicate\":" << replicate << ",\"right_output_sha256\":\"" << hashes[PhysicalIndex(1u, condition)]
        << "\",\"role\":\"" << RoleName(unit % 3u) << "\",\"runtime_library_maps_sha256\":\"" << Sha256Hex(maps_text)
        << "\",\"scenario\":" << scenario << ",\"schema\":\"" << kPrefix << ".panel.v1"
        << "\",\"scope\":\"prebuilt-systematic\",\"scope_invocations_per_batch\":" << count
        << ",\"sequence\":" << sequence << ",\"slots\":[";
    for (size_t slot = 0u; slot < order.size(); ++slot) {
        if (slot) std::cout << ',';
        std::cout << "{\"elapsed_ns\":" << elapsed[slot] << ",\"logical_lane\":\"" << LogicalName(order[slot])
            << "\",\"physical_lane\":" << PhysicalIndex(order[slot], condition) << '}';
    }
    std::cout << "],\"source_immutable_verified\":true,\"target_cpu\":120,\"tuple_index\":" << unit / 3u
        << ",\"unit_index\":" << unit << ",\"unit_key_sha256\":\"" << UnitKey(replicate, unit) << "\"}\n";
    std::cout.flush();
    return std::cout.good();
}

bool Run(const Clock::time_point& deadline)
{
    std::string diagnostic, identity, later_identity, maps_text, library, binding;
    CopyBinding copy;
    if (!PinAndVerifyTargetCpu(kTargetCpu, diagnostic) || !Guard(deadline) || !CaptureIdentity(identity) ||
        wirehair_init() != Wirehair_Success || !RuntimeWirehairMapsReceipt(maps_text, library) ||
        !RuntimeTimedBindingPath(binding) || binding != library || !copy.Resolve()) return false;
    std::vector<std::shared_ptr<const ScreenCell> > cells;
    for (const CellShape& shape : kShapes) {
        std::shared_ptr<ScreenCell> cell(new ScreenCell);
        if (!Guard(deadline) || !BuildCell(shape, *cell) || !Guard(deadline)) return false;
        cells.push_back(cell);
    }
    if (!CaptureIdentity(later_identity) || identity != later_identity || RuntimeWirehairMapsText() != maps_text ||
        !copy.Stable() || !Guard(deadline) || !EmitPageConfig(cells, Sha256Hex(maps_text), library, identity, copy)) return false;
    uint32_t panels = 0u;
    uint64_t measured = 0u, warmup = 0u, calls = 0u;
    for (uint32_t replicate = 0u; replicate < kReplicates; ++replicate)
        for (uint32_t unit : UnitOrder(replicate)) {
            const ScreenCell& cell = *cells[unit / 3u];
            PhysicalPair pair(cell, unit % 3u, copy);
            if (!Guard(deadline) || !pair.Initialize() || !Guard(deadline)) return false;
            for (uint32_t index = 0u; index < 4u; ++index) {
                const uint32_t scenario = (replicate + unit + index) % 4u;
                if (!pair.Select(scenario)) return false;
                for (uint32_t condition : kConditions[replicate % 4u]) {
                    if (!RunCondition(cell, replicate, unit, scenario, condition, panels, pair, maps_text, copy, deadline))
                        return false;
                    ++panels;
                    measured += UINT64_C(8) * (8192u / cell.Shape.K);
                    warmup += UINT64_C(2) * (8192u / cell.Shape.K);
                    calls += UINT64_C(10) * 8192u;
                }
            }
        }
    if (panels != 1728u || measured != UINT64_C(9732096) || warmup != UINT64_C(2433024) ||
        calls != UINT64_C(141557760) || !CaptureIdentity(later_identity) || identity != later_identity ||
        RuntimeWirehairMapsText() != maps_text || !copy.Stable() || !Guard(deadline)) return false;
    std::cout << "{\"campaign\":\"" << kName << "\",\"encode_call_count\":" << calls
        << ",\"measured_batch_count\":13824,\"measured_invocation_count\":" << measured
        << ",\"panel_count\":" << panels << ",\"record_count\":1730,\"schema\":\"" << kPrefix
        << ".terminal.v1\",\"status\":\"complete\",\"warmup_batch_count\":3456"
        << ",\"warmup_invocation_count\":" << warmup << "}\n";
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
    if (std::strcmp(directory, kDirectory) != 0 || !IsLowerSha256(claim_hash) || !IsLowerSha256(worker_hash) ||
        !IsLowerSha256(description_hash) || Sha256Hex(Describe()) != description_hash) return false;
    Fd dir(::open(directory, O_RDONLY | O_DIRECTORY | O_CLOEXEC | O_NOFOLLOW));
    struct stat info;
    if (dir.Value < 0 || ::fstat(dir.Value, &info) != 0 || !S_ISDIR(info.st_mode) || (info.st_mode & 07777) != 0700)
        return false;
    Fd claim(::openat(dir.Value, "CLAIM", O_RDONLY | O_CLOEXEC | O_NOFOLLOW));
    Fd executable(::open("/proc/self/exe", O_RDONLY | O_CLOEXEC));
    std::string bytes;
    if (claim.Value < 0 || executable.Value < 0 ||
        !VerifyRegularAuthorizationFile(claim.Value, 0400, 1024u * 1024u, bytes, diagnostic) || Sha256Hex(bytes) != claim_hash ||
        !VerifyRegularAuthorizationFile(executable.Value, 0, 64u * 1024u * 1024u, bytes, diagnostic) ||
        Sha256Hex(bytes) != worker_hash) return false;
    Fd marker(::openat(dir.Value, "WORKER_STARTED", O_WRONLY | O_CREAT | O_EXCL | O_CLOEXEC | O_NOFOLLOW, 0400));
    if (marker.Value < 0 || ::fchmod(marker.Value, 0400) != 0) return false;
    std::ostringstream json;
    json << "{\"campaign\":\"" << kName << "\",\"claim_sha256\":\"" << claim_hash
        << "\",\"description_sha256\":\"" << description_hash << "\",\"schema\":\"" << kPrefix
        << ".worker-started.v1\",\"source_commit\":\"" << WIREHAIR_WH2_SOURCE_GIT_COMMIT
        << "\",\"status\":\"started\",\"worker_sha256\":\"" << worker_hash << "\"}\n";
    return WriteFdAll(marker.Value, json.str(), diagnostic) && ::fsync(marker.Value) == 0 && ::fsync(dir.Value) == 0;
}
#else
bool Authorize(const char*, const char*, const char*, const char*, std::string&) { return false; }
#endif

bool Selftest()
{
    uint32_t panels = 0u;
    uint64_t measured = 0u, warmup = 0u, calls = 0u;
    for (uint32_t replicate = 0u; replicate < kReplicates; ++replicate) {
        const std::array<uint32_t, 9> units = UnitOrder(replicate);
        if (std::set<uint32_t>(units.begin(), units.end()).size() != 9u) return false;
        std::string last;
        for (uint32_t unit : units) {
            const std::string key = UnitKey(replicate, unit);
            if (!last.empty() && key <= last) return false;
            last = key;
            std::set<uint32_t> scenarios;
            for (uint32_t index = 0u; index < 4u; ++index) {
                scenarios.insert((replicate + unit + index) % 4u);
                for (uint32_t condition : kConditions[replicate % 4u]) {
                    uint32_t count[2] = {0u, 0u};
                    for (uint32_t logical : SlotOrder(condition)) ++count[PhysicalIndex(logical, condition)];
                    if (count[0] != 4u || count[1] != 4u) return false;
                    ++panels;
                    measured += UINT64_C(8) * (8192u / kShapes[unit / 3u].K);
                    warmup += UINT64_C(2) * (8192u / kShapes[unit / 3u].K);
                    calls += UINT64_C(10) * 8192u;
                }
            }
            if (scenarios.size() != 4u) return false;
        }
    }
    CopyBinding copy;
    if (panels != 1728u || measured != UINT64_C(9732096) || warmup != UINT64_C(2433024) ||
        calls != UINT64_C(141557760) || wirehair_init() != Wirehair_Success || !copy.Resolve()) return false;
    for (const CellShape& shape : kShapes) {
        ScreenCell cell;
        if (!BuildCell(shape, cell)) return false;
        for (uint32_t role = 0u; role < 3u; ++role) {
            PhysicalPair pair(cell, role, copy);
            if (!pair.Initialize()) return false;
            for (uint32_t scenario = 0u; scenario < 4u; ++scenario) {
                if (!pair.Select(scenario)) return false;
                for (uint32_t lane = 0u; lane < 2u; ++lane) {
                    uint64_t elapsed = 1u;
                    if (!pair.Batch(lane, false, 2u, elapsed) || elapsed != 0u ||
                        (reinterpret_cast<uintptr_t>(pair.Memory.Outputs[lane]) & 4095u) != kPhases[scenario][lane])
                        return false;
                }
                std::string hashes[2];
                if (!pair.Verify(hashes) || (reinterpret_cast<uintptr_t>(pair.Memory.Source) & 4095u) != 0x800u ||
                    ((pair.HandleAddress(0u) == pair.HandleAddress(1u)) != (role == 2u || scenario == 3u))) return false;
            }
        }
    }
    return copy.Stable() && Describe().back() == '\n';
}

} // namespace wh2_page_aa

int main(int argc, char** argv)
{
    using namespace wh2_page_aa;
    try {
        if (argc == 2 && std::strcmp(argv[1], "--describe") == 0) {
            std::cout << Describe();
            return std::cout.good() ? 0 : 1;
        }
        if (argc == 2 && std::strcmp(argv[1], "--selftest") == 0) {
            if (!Selftest()) { std::cerr << "Page AA selftest failed\n"; return 1; }
            std::cout << "WH2 public page AA r0 selftest passed\n";
            return 0;
        }
        const Clock::time_point deadline = Clock::now() + std::chrono::seconds(110);
        if (argc != 12 || std::strcmp(argv[1], "--run") != 0 || std::strcmp(argv[2], "--cpu") != 0 ||
            std::strcmp(argv[3], "120") != 0 || std::strcmp(argv[4], "--evidence-dir") != 0 ||
            std::strcmp(argv[6], "--claim-sha256") != 0 || std::strcmp(argv[8], "--worker-sha256") != 0 ||
            std::strcmp(argv[10], "--config-identity-sha256") != 0) {
            std::cerr << "Page AA expects --describe, --selftest, or exact authorized --run arguments\n";
            return 1;
        }
        std::string diagnostic;
        if (!Authorize(argv[5], argv[7], argv[9], argv[11], diagnostic)) {
            std::cerr << "Page AA launch authorization failed: " << diagnostic << '\n';
            return 1;
        }
        if (!Run(deadline)) { std::cerr << "Page AA failed\n"; return 1; }
        return 0;
    } catch (const std::exception& exc) {
        std::cerr << "Page AA exception: " << exc.what() << '\n';
    } catch (...) { std::cerr << "Page AA unknown exception\n"; }
    return 1;
}
