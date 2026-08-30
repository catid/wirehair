#include "Wh2NativeCodec.h"

#include "../codec/WirehairV2PrecodeDecode.h"

#include <algorithm>
#include <chrono>
#include <cstring>
#include <limits>
#include <new>
#include <stdexcept>
#include <utility>

#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS) && \
    defined(WIREHAIR_V2_ENABLE_BENCHMARK_EQUATIONS)
#error "native WH2 adapter cannot mix test hooks with benchmark equations"
#elif !defined(WIREHAIR_V2_ENABLE_TEST_HOOKS) && \
    !defined(WIREHAIR_V2_ENABLE_BENCHMARK_EQUATIONS)
#error "native WH2 adapter requires one explicit equation capability"
#endif

namespace wirehair_wh2_bench {
namespace {

static const uint32_t kMaxConstructionAttempt = 255u;
static const uint32_t kRecoveryOverheadCap = 4u;
static const char kNativeArmOwnedSource[] = "native-arm-owned-kxb-v1";
static const char kFixtureSourceAcquisition[] =
    "fixture-copy-before-clock-move-v1";

#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
thread_local int64_t NativeCounterBatchAllocationFailureCountdown = -1;

void TriggerNativeCounterBatchAllocationFailureForTesting()
{
    if (NativeCounterBatchAllocationFailureCountdown == 0) {
        throw std::bad_alloc();
    }
    if (NativeCounterBatchAllocationFailureCountdown > 0) {
        --NativeCounterBatchAllocationFailureCountdown;
    }
}
#endif

template <bool Measure>
class FixtureTimer;

template <>
class FixtureTimer<true>
{
public:
    FixtureTimer()
        : Start(std::chrono::steady_clock::now())
    {
    }

    uint64_t Stop() const
    {
        const std::chrono::steady_clock::time_point finish =
            std::chrono::steady_clock::now();
        const std::chrono::nanoseconds elapsed =
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                finish - Start);
        return elapsed.count() > 0 ? (uint64_t)elapsed.count() : 0u;
    }

private:
    std::chrono::steady_clock::time_point Start;
};

template <>
class FixtureTimer<false>
{
public:
    uint64_t Stop() const
    {
        return 0u;
    }
};

template <bool Measure>
struct FixtureElapsedPolicy;

template <>
struct FixtureElapsedPolicy<true>
{
    static bool Valid(uint64_t elapsed)
    {
        return elapsed > 0u;
    }
};

template <>
struct FixtureElapsedPolicy<false>
{
    static bool Valid(uint64_t elapsed)
    {
        return elapsed == 0u;
    }
};

uint64_t SplitMix64(uint64_t& state)
{
    state += UINT64_C(0x9e3779b97f4a7c15);
    uint64_t value = state;
    value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
    return value ^ (value >> 31);
}

bool CheckedProduct(uint64_t left, uint64_t right, size_t& product_out)
{
    product_out = 0u;
    if (left != 0u && right >
            (uint64_t)std::numeric_limits<size_t>::max() / left)
    {
        return false;
    }
    product_out = (size_t)(left * right);
    return true;
}

bool CertifiedEquationOptions(
    const wirehair_v2::MessagePrecodeEncoderOptions& options)
{
    const wirehair_v2::MessagePrecodeEncoderOptions defaults;
    return options.RecoveryMixCount == defaults.RecoveryMixCount &&
        options.DenseAnchors == defaults.DenseAnchors &&
        options.DenseIdentityCorner == defaults.DenseIdentityCorner &&
        options.PrecodeSeedSalt == defaults.PrecodeSeedSalt &&
        options.RecoveryRowSeedSalt == defaults.RecoveryRowSeedSalt;
}

bool DefaultWh2Options(
    const wirehair_v2::MessagePrecodeEncoderOptions& options)
{
    const wirehair_v2::MessagePrecodeEncoderOptions defaults;
    return CertifiedEquationOptions(options) &&
        options.CacheSystematicSource == defaults.CacheSystematicSource &&
        options.CacheReceivedSystematicPackets ==
            defaults.CacheReceivedSystematicPackets;
}

bool ValidArmSpec(const NativeArmSpec& spec)
{
    if (spec.SystematicEmission !=
            NativeSystematicEmission::EquationEvaluation &&
        spec.SystematicEmission !=
            NativeSystematicEmission::RetainedSourceDirect)
    {
        return false;
    }
    if (spec.Wh2Options.CacheSystematicSource ||
        spec.Wh2Options.CacheReceivedSystematicPackets)
    {
        // These flags affect only facade storage, while this adapter operates
        // directly on the exact equations.  Reject rather than silently ignore
        // an apparent timing axis.
        return false;
    }
    if (spec.Kind == NativeArmKind::Wirehair1)
    {
        return spec.ConstructionAttempt == 0u && !spec.Transform &&
            !spec.TransformContext && DefaultWh2Options(spec.Wh2Options) &&
            spec.BaseKind == NativeWh2BaseKind::ProductionProfile &&
            spec.SystematicEmission ==
                NativeSystematicEmission::EquationEvaluation;
    }
    if (spec.Kind == NativeArmKind::Wirehair2Certified)
    {
        return spec.ConstructionAttempt <= kMaxConstructionAttempt &&
            !spec.Transform && !spec.TransformContext &&
            CertifiedEquationOptions(spec.Wh2Options) &&
            spec.BaseKind == NativeWh2BaseKind::ProductionProfile;
    }
    if (spec.Kind == NativeArmKind::Wirehair2Experiment)
    {
        if (spec.ConstructionAttempt > kMaxConstructionAttempt ||
            (!spec.Transform && spec.TransformContext))
        {
            return false;
        }
        if (spec.BaseKind == NativeWh2BaseKind::ProductionProfile) {
            return true;
        }
        if (spec.BaseKind ==
                NativeWh2BaseKind::CanonicalCertifiedStructure)
        {
            return spec.Transform && DefaultWh2Options(spec.Wh2Options);
        }
        return false;
    }
    return false;
}

bool ValidSourceShape(
    uint32_t block_count,
    uint32_t block_bytes,
    const std::vector<uint8_t>& source)
{
    size_t expected_source_bytes = 0u;
    return block_count >= 2u && block_count <= 64000u &&
        block_bytes != 0u && block_bytes <= UINT32_C(0x7fffffff) &&
        CheckedProduct(block_count, block_bytes, expected_source_bytes) &&
        source.size() == expected_source_bytes;
}

bool ValidPrecodeParamsShape(
    const wirehair_v2::PrecodeParams& params,
    uint32_t expected_block_count)
{
    const uint64_t binary_span = (uint64_t)params.BlockCount +
        params.Staircase + params.DenseRows;
    const uint64_t total_span = binary_span + params.HeavyRows;
    const uint64_t known_span =
        (uint64_t)params.BlockCount + params.Staircase;
    bool valid_dense_anchors = params.DenseAnchors ==
        wirehair_v2::DenseAnchorLayout::Disabled;
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS) || \
    defined(WIREHAIR_V2_ENABLE_BENCHMARK_EQUATIONS)
    valid_dense_anchors = valid_dense_anchors ||
        ((params.DenseAnchors == wirehair_v2::DenseAnchorLayout::Two07 ||
          params.DenseAnchors ==
              wirehair_v2::DenseAnchorLayout::Four0369) &&
         params.DenseRows == 12u && params.HeavyRows == 12u &&
         !params.DenseIdentityCorner &&
         params.HeavyFamily ==
             wirehair_v2::HeavyCoefficientFamily::PeriodicCauchy);
#endif
    return params.BlockCount == expected_block_count &&
        params.BlockCount >= 2u && params.BlockCount <= 64000u &&
        params.Staircase != 0u &&
        params.SourceHits >= 1u && params.SourceHits <= 8u &&
        params.DenseRows <= 64u && params.HeavyRows <= 128u &&
        (params.HeavyFamily ==
             wirehair_v2::HeavyCoefficientFamily::PeriodicCauchy ||
         params.HeavyFamily ==
             wirehair_v2::HeavyCoefficientFamily::HashedNonzero) &&
        valid_dense_anchors &&
        binary_span <= UINT16_MAX && total_span <= UINT16_MAX &&
        (!params.DenseIdentityCorner ||
         known_span >= 2u * (uint64_t)(params.DenseRows >> 1));
}

bool SamePrecodeParams(
    const wirehair_v2::PrecodeParams& left,
    const wirehair_v2::PrecodeParams& right)
{
    return left.BlockCount == right.BlockCount &&
        left.Staircase == right.Staircase &&
        left.DenseRows == right.DenseRows &&
        left.HeavyRows == right.HeavyRows &&
        left.SourceHits == right.SourceHits &&
        left.DenseIdentityCorner == right.DenseIdentityCorner &&
        left.HeavyFamily == right.HeavyFamily &&
        left.DenseAnchors == right.DenseAnchors &&
        left.Seed == right.Seed;
}

bool SameNormalizedSolveStats(
    const wirehair_v2::PrecodeSolveStats& left,
    const wirehair_v2::PrecodeSolveStats& right)
{
    return left.PacketRows == right.PacketRows &&
        left.PeeledColumns == right.PeeledColumns &&
        left.InactivatedColumns == right.InactivatedColumns &&
        left.ResidualRows == right.ResidualRows &&
        left.ResidualRank == right.ResidualRank &&
        left.BinaryResidualRank == right.BinaryResidualRank &&
        left.BinaryRowReferences == right.BinaryRowReferences &&
        left.BinaryRowStorageBytes == right.BinaryRowStorageBytes &&
        left.BinaryAdjacencyStorageBytes ==
            right.BinaryAdjacencyStorageBytes &&
        left.BinaryRowStorageAllocations ==
            right.BinaryRowStorageAllocations &&
        left.BinaryAdjacencyStorageAllocations ==
            right.BinaryAdjacencyStorageAllocations &&
        left.BlockXors == right.BlockXors &&
        left.BlockMulAdds == right.BlockMulAdds &&
        left.PacketSeedAttempt == right.PacketSeedAttempt;
}

uint32_t PrecodeCount(const wirehair_v2::PrecodeSystem& system)
{
    return system.Params.Staircase + system.Params.DenseRows +
        system.Params.HeavyRows;
}

RecoveryCellResult ClassifyConstructionResult(WirehairResult result)
{
    RecoveryCellResult classified;
    classified.Result = result;
    classified.FirstOverhead = UINT32_MAX;
    switch (result)
    {
    case Wirehair_BadDenseSeed:
    case Wirehair_BadPeelSeed:
        classified.Outcome = RecoveryOutcome::ConstructFailed;
        break;
    case Wirehair_BadInput_SmallN:
    case Wirehair_BadInput_LargeN:
    case Wirehair_UnsupportedPlatform:
        classified.Outcome = RecoveryOutcome::Unsupported;
        break;
    default:
        classified.Outcome = RecoveryOutcome::Fatal;
        break;
    }
    return classified;
}

struct LegacyCodecOwner
{
    WirehairCodec Codec = nullptr;

    ~LegacyCodecOwner()
    {
        wirehair_free(Codec);
    }

    LegacyCodecOwner() = default;
    LegacyCodecOwner(const LegacyCodecOwner&) = delete;
    LegacyCodecOwner& operator=(const LegacyCodecOwner&) = delete;
};

bool UniquePacketPrefix(
    uint32_t block_count,
    uint32_t overhead,
    const std::vector<uint32_t>& packet_ids,
    bool require_exact_size)
{
    const uint64_t required_u64 = (uint64_t)block_count + overhead;
    if (required_u64 > UINT32_MAX ||
        required_u64 > (uint64_t)std::numeric_limits<size_t>::max())
    {
        return false;
    }
    const size_t required = (size_t)required_u64;
    if (packet_ids.size() < required ||
        (require_exact_size && packet_ids.size() != required))
    {
        return false;
    }
    std::vector<uint32_t> sorted(
        packet_ids.begin(), packet_ids.begin() + required);
    std::sort(sorted.begin(), sorted.end());
    return std::adjacent_find(sorted.begin(), sorted.end()) == sorted.end();
}

bool VerifyWh2DecodedBytes(
    const wirehair_v2::PrecodeSystem& system,
    const wirehair_v2::PacketRowConfig& config,
    const wirehair_v2::PacketRowRuntime& runtime,
    uint32_t block_bytes,
    const std::vector<uint8_t>& expected_source,
    const std::vector<uint8_t>& intermediate)
{
    const uint32_t K = system.Params.BlockCount;
    size_t source_bytes = 0u;
    size_t intermediate_bytes = 0u;
    if (!CheckedProduct(K, block_bytes, source_bytes) ||
        !CheckedProduct(
            (uint64_t)K + PrecodeCount(system),
            block_bytes,
            intermediate_bytes) ||
        expected_source.size() != source_bytes ||
        intermediate.size() != intermediate_bytes)
    {
        return false;
    }

    std::vector<uint8_t> block(block_bytes, 0u);
    for (uint32_t block_id = 0u; block_id < K; ++block_id)
    {
        if (!wirehair_v2::EvaluatePacketBlockForValidatedSystemWithRuntime(
                system,
                config,
                runtime,
                intermediate.data(),
                block_bytes,
                block_id,
                block.data()) ||
            std::memcmp(
                block.data(),
                expected_source.data() + (size_t)block_id * block_bytes,
                block_bytes) != 0)
        {
            return false;
        }
    }
    return true;
}

WirehairResult EncodePrefix(
    const NativeArm& arm,
    const std::vector<uint32_t>& packet_ids,
    std::vector<uint8_t>& packet_storage_out)
{
    size_t packet_bytes = 0u;
    if (!CheckedProduct(
            packet_ids.size(), arm.BlockBytes(), packet_bytes))
    {
        return Wirehair_InvalidInput;
    }
    try
    {
        std::vector<uint8_t> storage(packet_bytes, 0u);
        std::vector<uint8_t> packet;
        for (size_t i = 0u; i < packet_ids.size(); ++i)
        {
            const WirehairResult result = arm.Encode(packet_ids[i], packet);
            if (result != Wirehair_Success) {
                return result;
            }
            if (packet.size() != arm.BlockBytes()) {
                return Wirehair_Error;
            }
            std::memcpy(
                storage.data() + i * arm.BlockBytes(),
                packet.data(),
                arm.BlockBytes());
        }
        packet_storage_out.swap(storage);
        return Wirehair_Success;
    }
    catch (const std::bad_alloc&) {
        return Wirehair_OOM;
    }
    catch (const std::length_error&) {
        return Wirehair_OOM;
    }
}

} // namespace

struct NativeArm::Impl
{
    ~Impl()
    {
        wirehair_free(LegacyEncoder);
    }

    NativeArmKind Kind = NativeArmKind::Invalid;
    uint32_t K = 0u;
    uint32_t BlockBytes = 0u;
    uint32_t Attempt = 0u;
    NativeSystematicEmission SystematicEmission =
        NativeSystematicEmission::Invalid;
    std::vector<uint8_t> Source;
    WirehairCodec LegacyEncoder = nullptr;
    wirehair_v2::PrecodeSystem System = {};
    wirehair_v2::PacketRowConfig PacketConfig = {};
    wirehair_v2::PacketRowRuntime PacketRuntime = {};
    wirehair_v2::PrecodeSolveStats Stats = {};
    std::vector<uint8_t> Intermediate;
    bool Initialized = false;
};

const char* NativeSystematicEmissionName(
    NativeSystematicEmission emission)
{
    switch (emission)
    {
    case NativeSystematicEmission::EquationEvaluation:
        return "equation-eval-v1";
    case NativeSystematicEmission::RetainedSourceDirect:
        return "retained-source-direct-v1";
    default:
        return nullptr;
    }
}

NativeArmSpec MakeWirehair1Arm()
{
    NativeArmSpec spec;
    spec.Kind = NativeArmKind::Wirehair1;
    spec.SystematicEmission = NativeSystematicEmission::EquationEvaluation;
    return spec;
}

NativeArmSpec MakeCertifiedWh2Arm(uint32_t construction_attempt)
{
    NativeArmSpec spec;
    spec.Kind = NativeArmKind::Wirehair2Certified;
    spec.SystematicEmission = NativeSystematicEmission::EquationEvaluation;
    spec.ConstructionAttempt = construction_attempt;
    return spec;
}

NativeArmSpec MakeExperimentalWh2Arm(
    uint32_t construction_attempt,
    Wh2EquationTransform transform,
    void* transform_context,
    const wirehair_v2::MessagePrecodeEncoderOptions* options,
    NativeWh2BaseKind base_kind)
{
    NativeArmSpec spec;
    spec.Kind = NativeArmKind::Wirehair2Experiment;
    spec.SystematicEmission = NativeSystematicEmission::EquationEvaluation;
    spec.BaseKind = base_kind;
    spec.ConstructionAttempt = construction_attempt;
    spec.Transform = transform;
    spec.TransformContext = transform_context;
    if (options) {
        spec.Wh2Options = *options;
    }
    return spec;
}

bool ResolveNativeWh2Configuration(
    const NativeArmSpec& spec,
    uint32_t block_count,
    uint32_t block_bytes,
    ResolvedNativeWh2Configuration& result_out)
{
    if (!ValidArmSpec(spec) ||
        spec.Kind == NativeArmKind::Wirehair1 ||
        block_count < 2u || block_count > 64000u ||
        block_bytes == 0u || block_bytes > UINT32_C(0x7fffffff))
    {
        return false;
    }

    wirehair_v2::PrecodeParams base_params;
    wirehair_v2::PacketRowConfig base_packet_config;
    if (spec.BaseKind == NativeWh2BaseKind::ProductionProfile)
    {
        const wirehair_v2::SeedProfile profile =
            wirehair_v2::SelectSeedProfile(block_count, block_bytes);
        wirehair_v2::MessagePrecodeEncoderOptions options;
        if (!wirehair_v2::ResolveMessagePrecodeOptions(
                profile, &spec.Wh2Options, options) ||
            !wirehair_v2::ResolveMessagePrecodeConfiguration(
                profile, options, base_params, base_packet_config))
        {
            return false;
        }
    }
    else if (spec.BaseKind ==
                 NativeWh2BaseKind::CanonicalCertifiedStructure)
    {
        base_params = wirehair_v2::MakeCertifiedParams(block_count, 0u);
        base_packet_config.PeelSeed = 0u;
        base_packet_config.MixCount =
            wirehair_v2::kCertifiedPacketMixCount;
    }
    else {
        return false;
    }

    if (spec.Transform &&
        !spec.Transform(
            block_count,
            block_bytes,
            base_params,
            base_packet_config,
            spec.TransformContext))
    {
        return false;
    }

    ResolvedNativeWh2Configuration resolved;
    resolved.Params = wirehair_v2::PrecodeParamsForAttempt(
        base_params, spec.ConstructionAttempt);
    resolved.PacketConfig = wirehair_v2::PacketConfigForAttempt(
        base_packet_config, spec.ConstructionAttempt);
    resolved.PrecodeAttempt = spec.ConstructionAttempt;
    resolved.PacketAttempt = spec.ConstructionAttempt;
    if (!ValidPrecodeParamsShape(resolved.Params, block_count)) {
        return false;
    }
    result_out = resolved;
    return true;
}

bool MakeDeterministicSource(
    uint32_t block_count,
    uint32_t block_bytes,
    uint64_t source_seed,
    std::vector<uint8_t>& source_out)
{
    size_t bytes = 0u;
    if (block_count < 2u || block_count > 64000u ||
        block_bytes == 0u || block_bytes > UINT32_C(0x7fffffff) ||
        !CheckedProduct(block_count, block_bytes, bytes))
    {
        return false;
    }

    try
    {
        std::vector<uint8_t> source(bytes, 0u);
        uint64_t state = source_seed ^
            ((uint64_t)block_count * UINT64_C(0xd6e8feb86659fd93)) ^
            ((uint64_t)block_bytes * UINT64_C(0xa0761d6478bd642f));
        size_t offset = 0u;
        while (offset < source.size())
        {
            const uint64_t word = SplitMix64(state);
            for (uint32_t byte = 0u; byte < 8u && offset < source.size();
                 ++byte, ++offset)
            {
                source[offset] = (uint8_t)(word >> (byte * 8u));
            }
        }
        source_out.swap(source);
        return true;
    }
    catch (const std::bad_alloc&) {
        return false;
    }
    catch (const std::length_error&) {
        return false;
    }
}

NativeArm::NativeArm()
{
}

NativeArm::~NativeArm()
{
}

NativeArm::NativeArm(NativeArm&& other) noexcept = default;
NativeArm& NativeArm::operator=(NativeArm&& other) noexcept = default;

WirehairResult NativeArm::Initialize(
    const NativeArmSpec& spec,
    uint32_t block_count,
    uint32_t block_bytes,
    const std::vector<uint8_t>& source)
{
    if (!ValidSourceShape(block_count, block_bytes, source) ||
        !ValidArmSpec(spec))
    {
        return Wirehair_InvalidInput;
    }
    const WirehairResult init_result = wirehair_init();
    if (init_result != Wirehair_Success) {
        return init_result;
    }
    try {
        std::vector<uint8_t> owned_source(source);
        return InitializeOwnedSourceAfterGlobalInit(
            spec,
            block_count,
            block_bytes,
            std::move(owned_source));
    }
    catch (const std::bad_alloc&) {
        return Wirehair_OOM;
    }
    catch (const std::length_error&) {
        return Wirehair_OOM;
    }
}

WirehairResult NativeArm::InitializeOwnedSourceAfterGlobalInit(
    const NativeArmSpec& spec,
    uint32_t block_count,
    uint32_t block_bytes,
    std::vector<uint8_t>&& source)
{
    if (!ValidSourceShape(block_count, block_bytes, source) ||
        !ValidArmSpec(spec))
    {
        return Wirehair_InvalidInput;
    }

    try
    {
        std::shared_ptr<Impl> next(new Impl);
        next->Kind = spec.Kind;
        next->K = block_count;
        next->BlockBytes = block_bytes;
        next->Attempt = spec.ConstructionAttempt;
        next->SystematicEmission = spec.SystematicEmission;
        next->Source = std::move(source);

        if (spec.Kind == NativeArmKind::Wirehair1)
        {
            const WirehairResult result = wirehair_encoder_create_ex(
                nullptr,
                next->Source.data(),
                (uint64_t)next->Source.size(),
                block_bytes,
                &next->LegacyEncoder);
            if (result != Wirehair_Success) {
                return result;
            }
            next->Initialized = true;
            ImplValue.swap(next);
            return Wirehair_Success;
        }

        ResolvedNativeWh2Configuration resolved;
        if (!ResolveNativeWh2Configuration(
                spec, block_count, block_bytes, resolved))
        {
            return Wirehair_InvalidInput;
        }
        const wirehair_v2::PrecodeParams& params = resolved.Params;
        const wirehair_v2::PacketRowConfig& packet_config =
            resolved.PacketConfig;

        wirehair_v2::PrecodeSystem system;
        if (!wirehair_v2::BuildPrecodeSystem(params, system)) {
            return Wirehair_BadPeelSeed;
        }
        const uint64_t precode_count_wide = (uint64_t)params.Staircase +
            params.DenseRows + params.HeavyRows;
        if (precode_count_wide > UINT32_MAX ||
            !next->PacketRuntime.Initialize(
                block_count,
                (uint32_t)precode_count_wide,
                packet_config.MixCount))
        {
            return Wirehair_InvalidInput;
        }

        std::vector<wirehair_v2::SolvePacket> packets;
        packets.reserve(block_count);
        for (uint32_t block_id = 0u; block_id < block_count; ++block_id)
        {
            wirehair_v2::SolvePacket packet;
            packet.BlockId = block_id;
            packet.Data = next->Source.data() +
                (size_t)block_id * block_bytes;
            packets.push_back(packet);
        }

        wirehair_v2::PrecodeSolveStats stats;
        std::vector<uint8_t> intermediate;
        const WirehairResult solve_result =
            wirehair_v2::SolvePrecodeSystemForValidatedSystemWithRuntime(
                system,
                packet_config,
                next->PacketRuntime,
                packets,
                block_bytes,
                intermediate,
                &stats);
        if (solve_result != Wirehair_Success) {
            return wirehair_v2::ClassifyExactSystematicConstructionFailure(
                system,
                packet_config,
                next->PacketRuntime,
                solve_result);
        }

        stats.PacketSeedAttempt = spec.ConstructionAttempt;
        next->System = std::move(system);
        next->PacketConfig = packet_config;
        next->Stats = stats;
        next->Intermediate.swap(intermediate);
        next->Initialized = true;
        ImplValue.swap(next);
        return Wirehair_Success;
    }
    catch (const std::bad_alloc&) {
        return Wirehair_OOM;
    }
    catch (const std::length_error&) {
        return Wirehair_OOM;
    }
    catch (...) {
        return Wirehair_Error;
    }
}

bool NativeArm::IsInitialized() const
{
    return ImplValue && ImplValue->Initialized;
}

NativeArmKind NativeArm::Kind() const
{
    return IsInitialized() ? ImplValue->Kind : NativeArmKind::Invalid;
}

uint32_t NativeArm::BlockCount() const
{
    return IsInitialized() ? ImplValue->K : 0u;
}

uint32_t NativeArm::BlockBytes() const
{
    return IsInitialized() ? ImplValue->BlockBytes : 0u;
}

uint32_t NativeArm::ConstructionAttempt() const
{
    return IsInitialized() ? ImplValue->Attempt : 0u;
}

NativeSystematicEmission NativeArm::SystematicEmission() const
{
    return IsInitialized() ? ImplValue->SystematicEmission :
        NativeSystematicEmission::Invalid;
}

bool NativeArm::UsesRetainedSourceFor(uint32_t block_id) const
{
    return IsInitialized() && ImplValue->Kind != NativeArmKind::Wirehair1 &&
        UsesRetainedSourceForInitializedWh2(block_id);
}

bool NativeArm::UsesRetainedSourceForInitializedWh2(
    uint32_t block_id) const
{
    return block_id < ImplValue->K &&
        ImplValue->SystematicEmission ==
            NativeSystematicEmission::RetainedSourceDirect;
}

bool NativeArm::HasSameWh2SolvedState(const NativeArm& other) const
{
    if (!IsInitialized() || !other.IsInitialized() ||
        ImplValue->Kind == NativeArmKind::Wirehair1 ||
        other.ImplValue->Kind == NativeArmKind::Wirehair1 ||
        ImplValue->Kind != other.ImplValue->Kind ||
        ImplValue->K != other.ImplValue->K ||
        ImplValue->BlockBytes != other.ImplValue->BlockBytes ||
        ImplValue->Attempt != other.ImplValue->Attempt ||
        !SamePrecodeParams(
            ImplValue->System.Params, other.ImplValue->System.Params) ||
        ImplValue->System.StaircaseRows !=
            other.ImplValue->System.StaircaseRows ||
        ImplValue->System.DenseBasisRowColumns !=
            other.ImplValue->System.DenseBasisRowColumns ||
        ImplValue->PacketConfig.PeelSeed !=
            other.ImplValue->PacketConfig.PeelSeed ||
        ImplValue->PacketConfig.MixCount !=
            other.ImplValue->PacketConfig.MixCount ||
        ImplValue->PacketRuntime.SourcePrime() !=
            other.ImplValue->PacketRuntime.SourcePrime() ||
        ImplValue->PacketRuntime.PrecodePrime() !=
            other.ImplValue->PacketRuntime.PrecodePrime() ||
        ImplValue->Source != other.ImplValue->Source ||
        ImplValue->Intermediate != other.ImplValue->Intermediate ||
        !SameNormalizedSolveStats(ImplValue->Stats, other.ImplValue->Stats))
    {
        return false;
    }
    const uint32_t precode_count = PrecodeCount(ImplValue->System);
    return ImplValue->PacketRuntime.IsValidFor(
               ImplValue->K,
               precode_count,
               ImplValue->PacketConfig.MixCount) &&
        other.ImplValue->PacketRuntime.IsValidFor(
            other.ImplValue->K,
            precode_count,
            other.ImplValue->PacketConfig.MixCount);
}

WirehairResult NativeArm::Encode(
    uint32_t block_id,
    std::vector<uint8_t>& packet_out) const
{
    if (!IsInitialized()) {
        return Wirehair_InvalidInput;
    }
    try
    {
        std::vector<uint8_t> packet(ImplValue->BlockBytes, 0u);
        const WirehairResult result = EncodeInto(
            block_id, packet.data(), ImplValue->BlockBytes);
        if (result != Wirehair_Success) {
            return result;
        }
        packet_out.swap(packet);
        return Wirehair_Success;
    }
    catch (const std::bad_alloc&) {
        return Wirehair_OOM;
    }
    catch (const std::length_error&) {
        return Wirehair_OOM;
    }
}

WirehairResult NativeArm::EncodeInto(
    uint32_t block_id,
    uint8_t* packet_out,
    uint32_t packet_capacity) const
{
    if (!IsInitialized() || !packet_out ||
        packet_capacity < ImplValue->BlockBytes)
    {
        return Wirehair_InvalidInput;
    }
    if (ImplValue->Kind == NativeArmKind::Wirehair1)
    {
        uint32_t data_bytes = 0u;
        const WirehairResult result = wirehair_encode(
            ImplValue->LegacyEncoder,
            block_id,
            packet_out,
            packet_capacity,
            &data_bytes);
        if (result != Wirehair_Success) {
            return result;
        }
        return data_bytes == ImplValue->BlockBytes ?
            Wirehair_Success : Wirehair_Error;
    }
    if (UsesRetainedSourceForInitializedWh2(block_id))
    {
        std::memmove(
            packet_out,
            ImplValue->Source.data() +
                (size_t)block_id * ImplValue->BlockBytes,
            ImplValue->BlockBytes);
        return Wirehair_Success;
    }
    return wirehair_v2::EvaluatePacketBlockForValidatedSystemWithRuntime(
        ImplValue->System,
        ImplValue->PacketConfig,
        ImplValue->PacketRuntime,
        ImplValue->Intermediate.data(),
        ImplValue->BlockBytes,
        block_id,
        packet_out) ? Wirehair_Success : Wirehair_Error;
}

struct NativeRecoveryFixture::Impl
{
    std::shared_ptr<const NativeArm::Impl> ArmState;
    std::vector<uint32_t> PacketIds;
    std::vector<uint8_t> PacketStorage;
    bool Initialized = false;
};

NativeRecoveryFixture::NativeRecoveryFixture()
{
}

NativeRecoveryFixture::~NativeRecoveryFixture()
{
}

NativeRecoveryFixture::NativeRecoveryFixture(
    NativeRecoveryFixture&& other) noexcept = default;
NativeRecoveryFixture& NativeRecoveryFixture::operator=(
    NativeRecoveryFixture&& other) noexcept = default;

WirehairResult NativeRecoveryFixture::Initialize(
    const NativeArm& arm,
    const std::vector<uint32_t>& packet_ids)
{
    if (!arm.IsInitialized()) {
        return Wirehair_InvalidInput;
    }
    try
    {
        if (!UniquePacketPrefix(
                arm.ImplValue->K,
                kRecoveryOverheadCap,
                packet_ids,
                true))
        {
            return Wirehair_InvalidInput;
        }
        std::unique_ptr<Impl> next(new Impl);
        next->ArmState = arm.ImplValue;
        next->PacketIds = packet_ids;
        const WirehairResult encode_result =
            EncodePrefix(arm, packet_ids, next->PacketStorage);
        if (encode_result != Wirehair_Success) {
            return encode_result;
        }
        next->Initialized = true;
        ImplValue.swap(next);
        return Wirehair_Success;
    }
    catch (const std::bad_alloc&) {
        return Wirehair_OOM;
    }
    catch (const std::length_error&) {
        return Wirehair_OOM;
    }
}

bool NativeRecoveryFixture::IsInitialized() const
{
    return ImplValue && ImplValue->Initialized;
}

RecoveryCellResult NativeRecoveryFixture::RunNested() const
{
    RecoveryCellResult cell;
    if (!IsInitialized()) {
        return cell;
    }

    try
    {
        for (uint32_t overhead = 0u; overhead <= kRecoveryOverheadCap;
             ++overhead)
        {
            const NativeArm::Impl& arm = *ImplValue->ArmState;
            const uint32_t received = arm.K + overhead;
            if (arm.Kind == NativeArmKind::Wirehair1)
            {
                LegacyCodecOwner decoder;
                const WirehairResult create_result = wirehair_decoder_create_ex(
                    nullptr,
                    (uint64_t)arm.Source.size(),
                    arm.BlockBytes,
                    &decoder.Codec);
                if (create_result != Wirehair_Success) {
                    return ClassifyConstructionResult(create_result);
                }

                WirehairResult decode_result = Wirehair_NeedMore;
                for (uint32_t packet_index = 0u;
                     packet_index < received;
                     ++packet_index)
                {
                    decode_result = wirehair_decode(
                        decoder.Codec,
                        ImplValue->PacketIds[packet_index],
                        ImplValue->PacketStorage.data() +
                            (size_t)packet_index * arm.BlockBytes,
                        arm.BlockBytes);
                    if (decode_result == Wirehair_Success)
                    {
                        if (packet_index + 1u < arm.K) {
                            cell.Result = Wirehair_Error;
                            return cell;
                        }
                        std::vector<uint8_t> recovered(
                            arm.Source.size(), 0u);
                        const WirehairResult recover_result = wirehair_recover(
                            decoder.Codec,
                            recovered.data(),
                            recovered.size());
                        if (recover_result != Wirehair_Success ||
                            recovered != arm.Source)
                        {
                            cell.Result = recover_result == Wirehair_Success ?
                                Wirehair_Error : recover_result;
                            return cell;
                        }
                        cell.Outcome = RecoveryOutcome::Success;
                        cell.FirstOverhead = packet_index + 1u - arm.K;
                        cell.Result = Wirehair_Success;
                        return cell;
                    }
                    if (decode_result != Wirehair_NeedMore) {
                        cell.Result = decode_result;
                        return cell;
                    }
                }
            }
            else
            {
                std::vector<wirehair_v2::SolvePacket> packets;
                packets.reserve(received);
                for (uint32_t packet_index = 0u;
                     packet_index < received;
                     ++packet_index)
                {
                    wirehair_v2::SolvePacket packet;
                    packet.BlockId = ImplValue->PacketIds[packet_index];
                    packet.Data = ImplValue->PacketStorage.data() +
                        (size_t)packet_index * arm.BlockBytes;
                    packets.push_back(packet);
                }
                std::vector<uint8_t> intermediate;
                const WirehairResult solve_result =
                    wirehair_v2::
                        SolvePrecodeSystemForValidatedSystemWithRuntime(
                            arm.System,
                            arm.PacketConfig,
                            arm.PacketRuntime,
                            packets,
                            arm.BlockBytes,
                            intermediate);
                if (solve_result == Wirehair_Success)
                {
                    if (!VerifyWh2DecodedBytes(
                            arm.System,
                            arm.PacketConfig,
                            arm.PacketRuntime,
                            arm.BlockBytes,
                            arm.Source,
                            intermediate))
                    {
                        cell.Result = Wirehair_Error;
                        return cell;
                    }
                    cell.Outcome = RecoveryOutcome::Success;
                    cell.FirstOverhead = overhead;
                    cell.Result = Wirehair_Success;
                    return cell;
                }
                if (solve_result != Wirehair_NeedMore) {
                    cell.Result = solve_result;
                    return cell;
                }
            }
        }

        cell.Outcome = RecoveryOutcome::NeedMoreAtCap;
        cell.FirstOverhead = UINT32_MAX;
        cell.Result = Wirehair_NeedMore;
        return cell;
    }
    catch (const std::bad_alloc&) {
        cell.Result = Wirehair_OOM;
        return cell;
    }
    catch (const std::length_error&) {
        cell.Result = Wirehair_OOM;
        return cell;
    }
}

RecoveryCellResult RunRecoveryCell(
    const NativeArmSpec& spec,
    uint32_t block_count,
    uint32_t block_bytes,
    const std::vector<uint8_t>& source,
    const std::vector<uint32_t>& packet_ids)
{
    NativeArm arm;
    const WirehairResult construction =
        arm.Initialize(spec, block_count, block_bytes, source);
    if (construction != Wirehair_Success) {
        return ClassifyConstructionResult(construction);
    }
    NativeRecoveryFixture fixture;
    const WirehairResult fixture_result = fixture.Initialize(arm, packet_ids);
    if (fixture_result != Wirehair_Success)
    {
        RecoveryCellResult result;
        result.Result = fixture_result;
        return result;
    }
    return fixture.RunNested();
}

struct NativeTimingControlProbe::Impl
{
    uint32_t K = 0u;
    wirehair_v2::PrecodeSystem HeadSystem;
    wirehair_v2::PacketRowConfig HeadPacketConfig = {};
    wirehair_v2::PacketRowRuntime HeadPacketRuntime = {};
    bool Initialized = false;
};

NativeTimingControlProbe::NativeTimingControlProbe()
{
}

NativeTimingControlProbe::~NativeTimingControlProbe()
{
}

NativeTimingControlProbe::NativeTimingControlProbe(
    NativeTimingControlProbe&& other) noexcept = default;
NativeTimingControlProbe& NativeTimingControlProbe::operator=(
    NativeTimingControlProbe&& other) noexcept = default;

WirehairResult NativeTimingControlProbe::Initialize(
    const NativeArmSpec& wirehair2_head_spec,
    uint32_t block_count,
    uint32_t block_bytes)
{
    if (wirehair2_head_spec.Kind != NativeArmKind::Wirehair2Certified ||
        wirehair2_head_spec.SystematicEmission !=
            NativeSystematicEmission::EquationEvaluation ||
        !ValidArmSpec(wirehair2_head_spec) ||
        block_count < 2u || block_count > 64000u ||
        block_bytes == 0u || block_bytes > UINT32_C(0x7fffffff))
    {
        return Wirehair_InvalidInput;
    }

    const WirehairResult init_result = wirehair_init();
    if (init_result != Wirehair_Success) {
        return init_result;
    }

    try
    {
        ResolvedNativeWh2Configuration resolved;
        if (!ResolveNativeWh2Configuration(
                wirehair2_head_spec,
                block_count,
                block_bytes,
                resolved))
        {
            return Wirehair_InvalidInput;
        }

        std::unique_ptr<Impl> next(new Impl);
        next->K = block_count;
        if (!wirehair_v2::BuildPrecodeSystem(
                resolved.Params, next->HeadSystem))
        {
            return Wirehair_BadPeelSeed;
        }
        next->HeadPacketConfig = resolved.PacketConfig;
        const uint64_t precode_count_wide =
            static_cast<uint64_t>(resolved.Params.Staircase) +
            resolved.Params.DenseRows + resolved.Params.HeavyRows;
        if (precode_count_wide > UINT32_MAX ||
            !next->HeadPacketRuntime.Initialize(
                block_count,
                static_cast<uint32_t>(precode_count_wide),
                resolved.PacketConfig.MixCount))
        {
            return Wirehair_InvalidInput;
        }

        // NativeArm construction must solve this exact systematic system
        // before a timing solve fixture can exist.  Prove that once per
        // reusable (K, width, attempt) probe with a consistent one-byte RHS.
        const uint8_t zero = 0u;
        std::vector<wirehair_v2::SolvePacket> systematic(block_count);
        for (uint32_t block_id = 0u; block_id < block_count; ++block_id)
        {
            systematic[block_id].BlockId = block_id;
            systematic[block_id].Data = &zero;
        }
        std::vector<uint8_t> intermediate;
        const WirehairResult construction_result =
            wirehair_v2::SolvePrecodeSystemForValidatedSystemWithRuntime(
                next->HeadSystem,
                next->HeadPacketConfig,
                next->HeadPacketRuntime,
                systematic,
                1u,
                intermediate);
        if (construction_result != Wirehair_Success)
        {
            return construction_result == Wirehair_NeedMore ?
                Wirehair_BadPeelSeed : construction_result;
        }
        const uint64_t intermediate_blocks =
            static_cast<uint64_t>(block_count) + precode_count_wide;
        if (intermediate.size() != intermediate_blocks ||
            std::find_if(
                intermediate.begin(), intermediate.end(),
                [](uint8_t value) { return value != 0u; }) !=
                    intermediate.end())
        {
            return Wirehair_Error;
        }

        next->Initialized = true;
        ImplValue.swap(next);
        return Wirehair_Success;
    }
    catch (const std::bad_alloc&) {
        return Wirehair_OOM;
    }
    catch (const std::length_error&) {
        return Wirehair_OOM;
    }
    catch (...) {
        return Wirehair_Error;
    }
}

bool NativeTimingControlProbe::IsInitialized() const
{
    return ImplValue && ImplValue->Initialized;
}

NativeTimingControlQualificationResult NativeTimingControlProbe::Run(
    const std::vector<uint32_t>& packet_ids)
{
    NativeTimingControlQualificationResult result;
    static const uint32_t kReceiveOverheadCap = 256u;
    if (!IsInitialized())
    {
        result.Wirehair1Result = Wirehair_InvalidInput;
        result.Wirehair2HeadResult = Wirehair_InvalidInput;
        return result;
    }

    try
    {
        if (!UniquePacketPrefix(
                ImplValue->K, kReceiveOverheadCap, packet_ids, true))
        {
            result.Wirehair1Result = Wirehair_InvalidInput;
            result.Wirehair2HeadResult = Wirehair_InvalidInput;
            return result;
        }

        const uint8_t zero = 0u;
        const size_t solve_packet_count =
            static_cast<size_t>(ImplValue->K) + kRecoveryOverheadCap;
        std::vector<wirehair_v2::SolvePacket> solve_packets(
            solve_packet_count);
        for (size_t i = 0u; i < solve_packet_count; ++i)
        {
            solve_packets[i].BlockId = packet_ids[i];
            solve_packets[i].Data = &zero;
        }
        std::vector<uint8_t> intermediate;
        result.Wirehair2HeadResult =
            wirehair_v2::SolvePrecodeSystemForValidatedSystemWithRuntime(
                ImplValue->HeadSystem,
                ImplValue->HeadPacketConfig,
                ImplValue->HeadPacketRuntime,
                solve_packets,
                1u,
                intermediate);
        const uint64_t precode_count_wide =
            static_cast<uint64_t>(ImplValue->HeadSystem.Params.Staircase) +
            ImplValue->HeadSystem.Params.DenseRows +
            ImplValue->HeadSystem.Params.HeavyRows;
        const uint64_t expected_intermediate =
            static_cast<uint64_t>(ImplValue->K) + precode_count_wide;
        if (result.Wirehair2HeadResult == Wirehair_Success &&
            (intermediate.size() != expected_intermediate ||
             std::find_if(
                 intermediate.begin(), intermediate.end(),
                 [](uint8_t value) { return value != 0u; }) !=
                 intermediate.end()))
        {
            result.Wirehair2HeadResult = Wirehair_Error;
        }

        // LegacyCurrent's equation selection depends on K, not payload width.
        // K one-byte blocks therefore select the exact timing-cell matrix
        // while retaining the real public decode, identity, recovery-value,
        // and ReconstructOutput paths at minimum payload cost.  A fresh owner
        // on every call prevents a retry from retaining prior decoder state.
        LegacyCodecOwner legacy;
        const WirehairResult create_result = wirehair_decoder_create_ex(
            nullptr,
            ImplValue->K,
            1u,
            &legacy.Codec);
        result.Wirehair1Result = create_result;
        if (create_result != Wirehair_Success) {
            return result;
        }
        if (!legacy.Codec) {
            result.Wirehair1Result = Wirehair_Error;
            return result;
        }

        result.Wirehair1Result = Wirehair_NeedMore;
        for (size_t i = 0u; i < packet_ids.size(); ++i)
        {
            result.Wirehair1Result = wirehair_decode(
                legacy.Codec,
                packet_ids[i],
                &zero,
                1u);
            if (result.Wirehair1Result == Wirehair_Success)
            {
                const uint64_t received = static_cast<uint64_t>(i) + 1u;
                if (received < ImplValue->K ||
                    received >
                        static_cast<uint64_t>(ImplValue->K) +
                            kReceiveOverheadCap)
                {
                    result.Wirehair1Result = Wirehair_Error;
                }
                else {
                    result.Wirehair1DecodedOverhead =
                        static_cast<uint32_t>(received - ImplValue->K);
                    std::vector<uint8_t> recovered(ImplValue->K, 0xffu);
                    const WirehairResult recover_result = wirehair_recover(
                        legacy.Codec,
                        recovered.data(),
                        recovered.size());
                    if (recover_result != Wirehair_Success ||
                        std::find_if(
                            recovered.begin(), recovered.end(),
                            [](uint8_t value) { return value != 0u; }) !=
                            recovered.end())
                    {
                        result.Wirehair1Result =
                            recover_result == Wirehair_Success ?
                                Wirehair_Error : recover_result;
                        result.Wirehair1DecodedOverhead = UINT32_MAX;
                    }
                }
                break;
            }
            if (result.Wirehair1Result != Wirehair_NeedMore) {
                break;
            }
        }

        const bool wh1_known =
            result.Wirehair1Result == Wirehair_Success ||
            result.Wirehair1Result == Wirehair_NeedMore;
        const bool wh2_known =
            result.Wirehair2HeadResult == Wirehair_Success ||
            result.Wirehair2HeadResult == Wirehair_NeedMore;
        if (result.Wirehair1Result == Wirehair_Success &&
            result.Wirehair2HeadResult == Wirehair_Success)
        {
            result.Qualification =
                NativeTimingControlQualification::Success;
        }
        else if (wh1_known && wh2_known)
        {
            result.Qualification =
                NativeTimingControlQualification::NeedMore;
        }
        return result;
    }
    catch (const std::bad_alloc&)
    {
        result.Qualification = NativeTimingControlQualification::Fatal;
        result.Wirehair1Result = Wirehair_OOM;
        result.Wirehair1DecodedOverhead = UINT32_MAX;
        result.Wirehair2HeadResult = Wirehair_OOM;
        return result;
    }
    catch (const std::length_error&)
    {
        result.Qualification = NativeTimingControlQualification::Fatal;
        result.Wirehair1Result = Wirehair_OOM;
        result.Wirehair1DecodedOverhead = UINT32_MAX;
        result.Wirehair2HeadResult = Wirehair_OOM;
        return result;
    }
    catch (...)
    {
        result.Qualification = NativeTimingControlQualification::Fatal;
        result.Wirehair1Result = Wirehair_Error;
        result.Wirehair1DecodedOverhead = UINT32_MAX;
        result.Wirehair2HeadResult = Wirehair_Error;
        return result;
    }
}

struct NativeEncoderFixture::Impl
{
    NativeArmSpec Spec = {};
    uint32_t K = 0u;
    uint32_t BlockBytes = 0u;
    std::vector<uint8_t> Source;
    bool Initialized = false;
};

Batch8ArmResult::Batch8ArmResult()
    : ElapsedNanoseconds(0u)
    , VerificationMask(0u)
{
    Results.fill(Wirehair_Error);
    DecodedOverheads.fill(UINT32_MAX);
    DirectSystematicPackets.fill(0u);
}

Batch128CounterArmResult::Batch128CounterArmResult()
    : CounterStart(0u)
    , CounterEnd(0u)
    , CounterDelta(0u)
    , BeginReadSucceeded(false)
    , EndReadSucceeded(false)
{
    Results.fill(Wirehair_Error);
    VerificationMask.fill(0u);
    DecodedOverheads.fill(UINT32_MAX);
    DirectSystematicPackets.fill(0u);
}

#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
void SetNativeCounterBatchAllocationFailureCountdownForTesting(
    int64_t countdown)
{
    NativeCounterBatchAllocationFailureCountdown = countdown;
}
#endif

NativeEncoderFixture::NativeEncoderFixture()
{
}

NativeEncoderFixture::~NativeEncoderFixture()
{
}

NativeEncoderFixture::NativeEncoderFixture(
    NativeEncoderFixture&& other) noexcept = default;
NativeEncoderFixture& NativeEncoderFixture::operator=(
    NativeEncoderFixture&& other) noexcept = default;

WirehairResult NativeEncoderFixture::Initialize(
    const NativeArmSpec& spec,
    uint32_t block_count,
    uint32_t block_bytes,
    const std::vector<uint8_t>& source)
{
    try
    {
        if (!ValidSourceShape(block_count, block_bytes, source) ||
            !ValidArmSpec(spec))
        {
            return Wirehair_InvalidInput;
        }
        // Global CPU/GF table dispatch is process setup, not fresh encoder
        // construction.  Run it once here so it cannot enter the clock in
        // Run(); the private arm path below assumes this has succeeded.
        const WirehairResult init_result = wirehair_init();
        if (init_result != Wirehair_Success) {
            return init_result;
        }

        std::unique_ptr<Impl> next(new Impl);
        next->Spec = spec;
        next->K = block_count;
        next->BlockBytes = block_bytes;
        next->Source = source;
        next->Initialized = true;
        ImplValue.swap(next);
        return Wirehair_Success;
    }
    catch (const std::bad_alloc&) {
        return Wirehair_OOM;
    }
    catch (const std::length_error&) {
        return Wirehair_OOM;
    }
}

bool NativeEncoderFixture::IsInitialized() const
{
    return ImplValue && ImplValue->Initialized;
}

bool NativeEncoderFixture::GetStorageReceipt(
    NativeEncoderStorageReceipt& receipt_out) const
{
    receipt_out = NativeEncoderStorageReceipt();
    if (!IsInitialized() ||
        (ImplValue->Spec.Kind != NativeArmKind::Wirehair2Certified &&
         ImplValue->Spec.Kind != NativeArmKind::Wirehair2Experiment))
    {
        return false;
    }
    NativeEncoderStorageReceipt receipt;
    receipt.SourceStorage = kNativeArmOwnedSource;
    receipt.SourceAcquisition = kFixtureSourceAcquisition;
    receipt.SystematicEmission = ImplValue->Spec.SystematicEmission;
    receipt_out = receipt;
    return true;
}

template <bool Measure>
TimedArmResult NativeEncoderFixture::RunImpl() const
{
    TimedArmResult result;
    if (!IsInitialized()) {
        return result;
    }

    try
    {
        // Output storage and the empty handle exist before the measured
        // interval.  Exact construction and every requested symbol remain in
        // the interval; codec teardown occurs after the stop timestamp.
        std::vector<uint8_t> encoded(ImplValue->Source.size(), 0u);
        std::vector<uint8_t> owned_source(ImplValue->Source);
        NativeArm arm;
        WirehairResult invocation_result = Wirehair_Error;
        const FixtureTimer<Measure> timer;
        invocation_result = arm.InitializeOwnedSourceAfterGlobalInit(
            ImplValue->Spec,
            ImplValue->K,
            ImplValue->BlockBytes,
            std::move(owned_source));
        if (invocation_result == Wirehair_Success)
        {
            for (uint32_t block_id = 0u;
                 block_id < ImplValue->K;
                 ++block_id)
            {
                invocation_result = arm.EncodeInto(
                    block_id,
                    encoded.data() +
                        (size_t)block_id * ImplValue->BlockBytes,
                    ImplValue->BlockBytes);
                if (invocation_result != Wirehair_Success) {
                    break;
                }
            }
        }
        result.ElapsedNanoseconds = timer.Stop();
        result.Result = invocation_result;
        if (invocation_result == Wirehair_Success)
        {
            result.BytesVerified = encoded == ImplValue->Source;
            const uint64_t expected_direct =
                ImplValue->Spec.SystematicEmission ==
                    NativeSystematicEmission::RetainedSourceDirect ?
                    ImplValue->K : 0u;
            const bool direct_roster =
                arm.UsesRetainedSourceFor(0u) &&
                arm.UsesRetainedSourceFor(ImplValue->K - 1u) &&
                !arm.UsesRetainedSourceFor(ImplValue->K);
            const bool equation_roster =
                !arm.UsesRetainedSourceFor(0u) &&
                !arm.UsesRetainedSourceFor(ImplValue->K - 1u) &&
                !arm.UsesRetainedSourceFor(ImplValue->K);
            result.DirectSystematicPackets = direct_roster ?
                ImplValue->K : 0u;
            if (!result.BytesVerified ||
                !FixtureElapsedPolicy<Measure>::Valid(
                    result.ElapsedNanoseconds) ||
                (!direct_roster && !equation_roster) ||
                result.DirectSystematicPackets != expected_direct)
            {
                result.Result = Wirehair_Error;
            }
        }
        if (result.Result != Wirehair_Success) {
            result.ElapsedNanoseconds = 0u;
            result.BytesVerified = false;
            result.DirectSystematicPackets = 0u;
        }
        return result;
    }
    catch (const std::bad_alloc&) {
        result.Result = Wirehair_OOM;
        return result;
    }
    catch (const std::length_error&) {
        result.Result = Wirehair_OOM;
        return result;
    }
}

TimedArmResult NativeEncoderFixture::Preflight() const
{
    return RunImpl<false>();
}

TimedArmResult NativeEncoderFixture::Run() const
{
    return RunImpl<true>();
}

template <bool Measure>
Batch8ArmResult NativeEncoderFixture::RunBatch8Impl() const
{
    Batch8ArmResult result;
    if (!IsInitialized()) {
        return result;
    }

    size_t failure_invocation = 0u;
    try
    {
        // All output buffers, source copies, and empty arm handles are
        // distinct and prepared before the one aggregate interval.
        std::array<std::vector<uint8_t>, kNativeBatch8Size> encoded;
        std::array<std::vector<uint8_t>, kNativeBatch8Size> owned_sources;
        std::array<NativeArm, kNativeBatch8Size> arms;
        for (size_t invocation = 0u;
             invocation < kNativeBatch8Size;
             ++invocation)
        {
            failure_invocation = invocation;
            encoded[invocation].assign(ImplValue->Source.size(), 0u);
            owned_sources[invocation] = ImplValue->Source;
        }

        const FixtureTimer<Measure> timer;
        for (size_t invocation = 0u;
             invocation < kNativeBatch8Size;
             ++invocation)
        {
            failure_invocation = invocation;
            WirehairResult invocation_result =
                arms[invocation].InitializeOwnedSourceAfterGlobalInit(
                    ImplValue->Spec,
                    ImplValue->K,
                    ImplValue->BlockBytes,
                    std::move(owned_sources[invocation]));
            if (invocation_result == Wirehair_Success)
            {
                for (uint32_t block_id = 0u;
                     block_id < ImplValue->K;
                     ++block_id)
                {
                    invocation_result = arms[invocation].EncodeInto(
                        block_id,
                        encoded[invocation].data() +
                            (size_t)block_id * ImplValue->BlockBytes,
                        ImplValue->BlockBytes);
                    if (invocation_result != Wirehair_Success) {
                        break;
                    }
                }
            }
            result.Results[invocation] = invocation_result;
        }
        result.ElapsedNanoseconds = timer.Stop();

        bool batch_valid = FixtureElapsedPolicy<Measure>::Valid(
            result.ElapsedNanoseconds);
        const uint64_t expected_direct =
            ImplValue->Spec.SystematicEmission ==
                NativeSystematicEmission::RetainedSourceDirect ?
                ImplValue->K : 0u;
        for (size_t invocation = 0u;
             invocation < kNativeBatch8Size;
             ++invocation)
        {
            if (result.Results[invocation] != Wirehair_Success) {
                batch_valid = false;
                continue;
            }

            const NativeArm& arm = arms[invocation];
            const bool bytes_verified =
                encoded[invocation] == ImplValue->Source;
            const bool direct_roster =
                arm.UsesRetainedSourceFor(0u) &&
                arm.UsesRetainedSourceFor(ImplValue->K - 1u) &&
                !arm.UsesRetainedSourceFor(ImplValue->K);
            const bool equation_roster =
                !arm.UsesRetainedSourceFor(0u) &&
                !arm.UsesRetainedSourceFor(ImplValue->K - 1u) &&
                !arm.UsesRetainedSourceFor(ImplValue->K);
            const uint64_t direct_count = direct_roster ?
                ImplValue->K : 0u;
            if (!bytes_verified ||
                (!direct_roster && !equation_roster) ||
                direct_count != expected_direct)
            {
                batch_valid = false;
                continue;
            }
            result.VerificationMask |= UINT32_C(1) << invocation;
            result.DirectSystematicPackets[invocation] = direct_count;
        }
        if (!batch_valid) {
            result.ElapsedNanoseconds = 0u;
        }
        return result;
    }
    catch (const std::bad_alloc&) {
        result.Results[failure_invocation] = Wirehair_OOM;
        return result;
    }
    catch (const std::length_error&) {
        result.Results[failure_invocation] = Wirehair_OOM;
        return result;
    }
}

Batch8ArmResult NativeEncoderFixture::PreflightBatch8() const
{
    return RunBatch8Impl<false>();
}

Batch8ArmResult NativeEncoderFixture::RunBatch8() const
{
    return RunBatch8Impl<true>();
}

template <bool ReadCounter>
Batch128CounterArmResult
NativeEncoderFixture::RunWh1Batch128CounterImpl(
    const NativeCounterBracket* bracket) const
{
    Batch128CounterArmResult result;
    NativeCounterBracket frozen_bracket;
    if (!IsInitialized()) {
        return result;
    }
    // This counter seam is deliberately WH1-only.  Reject other identities
    // before inspecting or invoking any caller-supplied callback.
    if (ImplValue->Spec.Kind != NativeArmKind::Wirehair1)
    {
        result.Results[0] = Wirehair_InvalidInput;
        return result;
    }
    if (ReadCounter && (!bracket || !bracket->Read))
    {
        result.Results[0] = Wirehair_InvalidInput;
        return result;
    }
    if (ReadCounter) {
        frozen_bracket = *bracket;
    }

    size_t failure_invocation = 0u;
    try
    {
        // Allocate 128 independent output/source/handle states before Begin.
        // Fresh WH1 construction plus IDs [0,K) remains the existing encoder
        // timed core and is the only work between the two injected reads.
        std::array<std::vector<uint8_t>, kNativeCounterBatchSize> encoded;
        std::array<std::vector<uint8_t>, kNativeCounterBatchSize>
            owned_sources;
        std::array<NativeArm, kNativeCounterBatchSize> arms;
        for (size_t invocation = 0u;
             invocation < kNativeCounterBatchSize;
             ++invocation)
        {
            failure_invocation = invocation;
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
            TriggerNativeCounterBatchAllocationFailureForTesting();
#endif
            encoded[invocation].assign(ImplValue->Source.size(), 0u);
            owned_sources[invocation] = ImplValue->Source;
        }

        if (ReadCounter)
        {
            uint64_t counter_start = 0u;
            try
            {
                if (!frozen_bracket.Read(
                        frozen_bracket.Context,
                        frozen_bracket.Selector,
                        &counter_start))
                {
                    return result;
                }
            }
            catch (...)
            {
                return result;
            }
            result.CounterStart = counter_start;
            result.BeginReadSucceeded = true;
        }

        bool core_threw = false;
        try
        {
            for (size_t invocation = 0u;
                 invocation < kNativeCounterBatchSize;
                 ++invocation)
            {
                failure_invocation = invocation;
                WirehairResult invocation_result =
                    arms[invocation].InitializeOwnedSourceAfterGlobalInit(
                        ImplValue->Spec,
                        ImplValue->K,
                        ImplValue->BlockBytes,
                        std::move(owned_sources[invocation]));
                if (invocation_result == Wirehair_Success)
                {
                    for (uint32_t block_id = 0u;
                         block_id < ImplValue->K;
                         ++block_id)
                    {
                        invocation_result = arms[invocation].EncodeInto(
                            block_id,
                            encoded[invocation].data() +
                                (size_t)block_id * ImplValue->BlockBytes,
                            ImplValue->BlockBytes);
                        if (invocation_result != Wirehair_Success) {
                            break;
                        }
                    }
                }
                result.Results[invocation] = invocation_result;
            }
        }
        catch (const std::bad_alloc&)
        {
            result.Results[failure_invocation] = Wirehair_OOM;
            core_threw = true;
        }
        catch (const std::length_error&)
        {
            result.Results[failure_invocation] = Wirehair_OOM;
            core_threw = true;
        }
        catch (...)
        {
            result.Results[failure_invocation] = Wirehair_Error;
            core_threw = true;
        }

        if (ReadCounter)
        {
            uint64_t counter_end = 0u;
            try
            {
                if (!frozen_bracket.Read(
                        frozen_bracket.Context,
                        frozen_bracket.Selector,
                        &counter_end))
                {
                    return result;
                }
            }
            catch (...)
            {
                return result;
            }
            result.CounterEnd = counter_end;
            result.EndReadSucceeded = true;
            if (counter_end > result.CounterStart) {
                result.CounterDelta = counter_end - result.CounterStart;
            }
        }
        if (core_threw) {
            return result;
        }

        const uint64_t expected_direct = 0u;
        for (size_t invocation = 0u;
             invocation < kNativeCounterBatchSize;
             ++invocation)
        {
            if (result.Results[invocation] != Wirehair_Success) {
                continue;
            }

            const NativeArm& arm = arms[invocation];
            const bool bytes_verified =
                encoded[invocation] == ImplValue->Source;
            const bool direct_roster =
                arm.UsesRetainedSourceFor(0u) &&
                arm.UsesRetainedSourceFor(ImplValue->K - 1u) &&
                !arm.UsesRetainedSourceFor(ImplValue->K);
            const bool equation_roster =
                !arm.UsesRetainedSourceFor(0u) &&
                !arm.UsesRetainedSourceFor(ImplValue->K - 1u) &&
                !arm.UsesRetainedSourceFor(ImplValue->K);
            const uint64_t direct_count = direct_roster ?
                ImplValue->K : 0u;
            if (!bytes_verified ||
                (!direct_roster && !equation_roster) ||
                direct_count != expected_direct)
            {
                continue;
            }
            result.VerificationMask[invocation >> 6] |=
                UINT64_C(1) << (invocation & 63u);
            result.DirectSystematicPackets[invocation] = direct_count;
        }
        return result;
    }
    catch (const std::bad_alloc&)
    {
        result.Results[failure_invocation] = Wirehair_OOM;
        return result;
    }
    catch (const std::length_error&)
    {
        result.Results[failure_invocation] = Wirehair_OOM;
        return result;
    }
}

Batch128CounterArmResult NativeEncoderFixture::PreflightBatch128() const
{
    return RunWh1Batch128CounterImpl<false>(nullptr);
}

Batch128CounterArmResult NativeEncoderFixture::RunWh1Batch128Counter(
    const NativeCounterBracket& bracket) const
{
    return RunWh1Batch128CounterImpl<true>(&bracket);
}

struct NativeReceiveFixture::Impl
{
    std::shared_ptr<const NativeArm::Impl> ArmState;
    std::vector<uint32_t> PacketIds;
    std::vector<uint8_t> PacketStorage;
    uint32_t ReceiveOverheadCap = 0u;
    bool Initialized = false;
};

NativeReceiveFixture::NativeReceiveFixture()
{
}

NativeReceiveFixture::~NativeReceiveFixture()
{
}

NativeReceiveFixture::NativeReceiveFixture(
    NativeReceiveFixture&& other) noexcept = default;
NativeReceiveFixture& NativeReceiveFixture::operator=(
    NativeReceiveFixture&& other) noexcept = default;

WirehairResult NativeReceiveFixture::Initialize(
    const NativeArm& arm,
    const std::vector<uint32_t>& packet_ids,
    uint32_t receive_overhead_cap)
{
    if (!arm.IsInitialized() ||
        arm.SystematicEmission() !=
            NativeSystematicEmission::EquationEvaluation)
    {
        return Wirehair_InvalidInput;
    }
    try
    {
        if (!UniquePacketPrefix(
                arm.ImplValue->K,
                receive_overhead_cap,
                packet_ids,
                true))
        {
            return Wirehair_InvalidInput;
        }
        std::unique_ptr<Impl> next(new Impl);
        next->ArmState = arm.ImplValue;
        next->PacketIds = packet_ids;
        next->ReceiveOverheadCap = receive_overhead_cap;
        const WirehairResult encode_result =
            EncodePrefix(arm, packet_ids, next->PacketStorage);
        if (encode_result != Wirehair_Success) {
            return encode_result;
        }

        next->Initialized = true;
        ImplValue.swap(next);
        return Wirehair_Success;
    }
    catch (const std::bad_alloc&) {
        return Wirehair_OOM;
    }
    catch (const std::length_error&) {
        return Wirehair_OOM;
    }
}

bool NativeReceiveFixture::IsInitialized() const
{
    return ImplValue && ImplValue->Initialized;
}

template <bool Measure>
TimedArmResult NativeReceiveFixture::RunImpl() const
{
    TimedArmResult result;
    if (!IsInitialized()) {
        return result;
    }

    try
    {
        const NativeArm::Impl& arm = *ImplValue->ArmState;
        const size_t receive_packet_count =
            (size_t)arm.K + ImplValue->ReceiveOverheadCap;
        std::vector<uint8_t> recovered(arm.Source.size(), 0u);
        WirehairResult receive_result = Wirehair_NeedMore;
        uint32_t success_packet_count = 0u;

        if (arm.Kind == NativeArmKind::Wirehair1)
        {
            // Decoder construction is explicitly outside this scope.
            LegacyCodecOwner decoder;
            const WirehairResult create_result = wirehair_decoder_create_ex(
                nullptr,
                (uint64_t)arm.Source.size(),
                arm.BlockBytes,
                &decoder.Codec);
            if (create_result != Wirehair_Success) {
                result.Result = create_result;
                return result;
            }

            const FixtureTimer<Measure> timer;
            for (size_t packet_index = 0u;
                 packet_index < receive_packet_count;
                 ++packet_index)
            {
                receive_result = wirehair_decode(
                    decoder.Codec,
                    ImplValue->PacketIds[packet_index],
                    ImplValue->PacketStorage.data() +
                        (size_t)packet_index * arm.BlockBytes,
                    arm.BlockBytes);
                if (receive_result == Wirehair_Success)
                {
                    success_packet_count = (uint32_t)packet_index + 1u;
                    receive_result = wirehair_recover(
                        decoder.Codec,
                        recovered.data(),
                        recovered.size());
                    break;
                }
                if (receive_result != Wirehair_NeedMore) {
                    break;
                }
            }
            result.ElapsedNanoseconds = timer.Stop();
        }
        else
        {
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS) || \
    defined(WIREHAIR_V2_ENABLE_BENCHMARK_EQUATIONS)
            // Construct and initialize a fresh real decoder outside the
            // measured interval, matching Wirehair1.  The exact-system seam
            // preserves experimental systems/configurations/attempts while
            // DecodeResult and RecoverResult perform all timed receive work.
            wirehair_v2::MessagePrecodeDecoder decoder;
            const WirehairResult initialize_result =
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
                decoder.InitializeForValidatedSystemForTesting(
#else
                decoder.InitializeForValidatedSystemForBenchmark(
#endif
                    (uint64_t)arm.Source.size(),
                    arm.BlockBytes,
                    arm.System,
                    arm.PacketConfig,
                    arm.Attempt);
            if (initialize_result != Wirehair_Success) {
                result.Result = initialize_result;
                return result;
            }
            const FixtureTimer<Measure> timer;
            for (size_t packet_index = 0u;
                 packet_index < receive_packet_count;
                 ++packet_index)
            {
                receive_result = decoder.DecodeResult(
                    ImplValue->PacketIds[packet_index],
                    ImplValue->PacketStorage.data() +
                        packet_index * arm.BlockBytes,
                    arm.BlockBytes);
                if (receive_result == Wirehair_Success)
                {
                    success_packet_count = (uint32_t)packet_index + 1u;
                    receive_result = decoder.RecoverResult(
                        recovered.data(), recovered.size());
                    break;
                }
                if (receive_result != Wirehair_NeedMore) {
                    break;
                }
            }
            result.ElapsedNanoseconds = timer.Stop();
#else
            result.Result = Wirehair_UnsupportedPlatform;
            return result;
#endif
        }

        result.Result = receive_result;
        if (receive_result == Wirehair_Success)
        {
            if (success_packet_count < arm.K ||
                success_packet_count >
                    (uint64_t)arm.K + ImplValue->ReceiveOverheadCap)
            {
                result.Result = Wirehair_Error;
                result.ElapsedNanoseconds = 0u;
                return result;
            }
            result.DecodedOverhead = success_packet_count - arm.K;
            result.BytesVerified = recovered == arm.Source;
            if (!result.BytesVerified ||
                !FixtureElapsedPolicy<Measure>::Valid(
                    result.ElapsedNanoseconds))
            {
                result.Result = Wirehair_Error;
            }
        }
        if (result.Result != Wirehair_Success)
        {
            result.ElapsedNanoseconds = 0u;
            result.BytesVerified = false;
            result.DecodedOverhead = UINT32_MAX;
        }
        return result;
    }
    catch (const std::bad_alloc&) {
        result.Result = Wirehair_OOM;
        result.ElapsedNanoseconds = 0u;
        result.DecodedOverhead = UINT32_MAX;
        return result;
    }
    catch (const std::length_error&) {
        result.Result = Wirehair_OOM;
        result.ElapsedNanoseconds = 0u;
        result.DecodedOverhead = UINT32_MAX;
        return result;
    }
}

TimedArmResult NativeReceiveFixture::Preflight() const
{
    return RunImpl<false>();
}

TimedArmResult NativeReceiveFixture::Run() const
{
    return RunImpl<true>();
}

template <bool Measure>
Batch8ArmResult NativeReceiveFixture::RunBatch8Impl() const
{
    Batch8ArmResult result;
    if (!IsInitialized()) {
        return result;
    }

    size_t failure_invocation = 0u;
    try
    {
        const NativeArm::Impl& arm = *ImplValue->ArmState;
        const size_t receive_packet_count =
            (size_t)arm.K + ImplValue->ReceiveOverheadCap;

        // Scope dispatch and all fresh-state preparation happen before the
        // timer.  Each branch then runs one homogeneous core eight times.
        if (arm.Kind == NativeArmKind::Wirehair1)
        {
            std::array<LegacyCodecOwner, kNativeBatch8Size> decoders;
            std::array<std::vector<uint8_t>, kNativeBatch8Size> recovered;
            std::array<uint32_t, kNativeBatch8Size> success_packet_counts = {};
            for (size_t invocation = 0u;
                 invocation < kNativeBatch8Size;
                 ++invocation)
            {
                failure_invocation = invocation;
                recovered[invocation].assign(arm.Source.size(), 0u);
                const WirehairResult create_result = wirehair_decoder_create_ex(
                    nullptr,
                    (uint64_t)arm.Source.size(),
                    arm.BlockBytes,
                    &decoders[invocation].Codec);
                if (create_result != Wirehair_Success)
                {
                    result.Results[invocation] = create_result;
                    return result;
                }
            }

            const FixtureTimer<Measure> timer;
            for (size_t invocation = 0u;
                 invocation < kNativeBatch8Size;
                 ++invocation)
            {
                failure_invocation = invocation;
                WirehairResult receive_result = Wirehair_NeedMore;
                for (size_t packet_index = 0u;
                     packet_index < receive_packet_count;
                     ++packet_index)
                {
                    receive_result = wirehair_decode(
                        decoders[invocation].Codec,
                        ImplValue->PacketIds[packet_index],
                        ImplValue->PacketStorage.data() +
                            packet_index * arm.BlockBytes,
                        arm.BlockBytes);
                    if (receive_result == Wirehair_Success)
                    {
                        success_packet_counts[invocation] =
                            (uint32_t)packet_index + 1u;
                        receive_result = wirehair_recover(
                            decoders[invocation].Codec,
                            recovered[invocation].data(),
                            recovered[invocation].size());
                        break;
                    }
                    if (receive_result != Wirehair_NeedMore) {
                        break;
                    }
                }
                result.Results[invocation] = receive_result;
            }
            result.ElapsedNanoseconds = timer.Stop();

            bool batch_valid = FixtureElapsedPolicy<Measure>::Valid(
                result.ElapsedNanoseconds);
            for (size_t invocation = 0u;
                 invocation < kNativeBatch8Size;
                 ++invocation)
            {
                if (result.Results[invocation] != Wirehair_Success ||
                    success_packet_counts[invocation] < arm.K ||
                    success_packet_counts[invocation] >
                        (uint64_t)arm.K + ImplValue->ReceiveOverheadCap ||
                    recovered[invocation] != arm.Source)
                {
                    batch_valid = false;
                    continue;
                }
                result.VerificationMask |= UINT32_C(1) << invocation;
                result.DecodedOverheads[invocation] =
                    success_packet_counts[invocation] - arm.K;
            }
            if (!batch_valid) {
                result.ElapsedNanoseconds = 0u;
            }
            return result;
        }

#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS) || \
    defined(WIREHAIR_V2_ENABLE_BENCHMARK_EQUATIONS)
        std::array<wirehair_v2::MessagePrecodeDecoder,
            kNativeBatch8Size> decoders;
        std::array<std::vector<uint8_t>, kNativeBatch8Size> recovered;
        std::array<uint32_t, kNativeBatch8Size> success_packet_counts = {};
        for (size_t invocation = 0u;
             invocation < kNativeBatch8Size;
             ++invocation)
        {
            failure_invocation = invocation;
            recovered[invocation].assign(arm.Source.size(), 0u);
            const WirehairResult initialize_result =
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
                decoders[invocation].InitializeForValidatedSystemForTesting(
#else
                decoders[invocation].InitializeForValidatedSystemForBenchmark(
#endif
                    (uint64_t)arm.Source.size(),
                    arm.BlockBytes,
                    arm.System,
                    arm.PacketConfig,
                    arm.Attempt);
            if (initialize_result != Wirehair_Success)
            {
                result.Results[invocation] = initialize_result;
                return result;
            }
        }

        const FixtureTimer<Measure> timer;
        for (size_t invocation = 0u;
             invocation < kNativeBatch8Size;
             ++invocation)
        {
            failure_invocation = invocation;
            WirehairResult receive_result = Wirehair_NeedMore;
            for (size_t packet_index = 0u;
                 packet_index < receive_packet_count;
                 ++packet_index)
            {
                receive_result = decoders[invocation].DecodeResult(
                    ImplValue->PacketIds[packet_index],
                    ImplValue->PacketStorage.data() +
                        packet_index * arm.BlockBytes,
                    arm.BlockBytes);
                if (receive_result == Wirehair_Success)
                {
                    success_packet_counts[invocation] =
                        (uint32_t)packet_index + 1u;
                    receive_result = decoders[invocation].RecoverResult(
                        recovered[invocation].data(),
                        recovered[invocation].size());
                    break;
                }
                if (receive_result != Wirehair_NeedMore) {
                    break;
                }
            }
            result.Results[invocation] = receive_result;
        }
        result.ElapsedNanoseconds = timer.Stop();

        bool batch_valid = FixtureElapsedPolicy<Measure>::Valid(
            result.ElapsedNanoseconds);
        for (size_t invocation = 0u;
             invocation < kNativeBatch8Size;
             ++invocation)
        {
            if (result.Results[invocation] != Wirehair_Success ||
                success_packet_counts[invocation] < arm.K ||
                success_packet_counts[invocation] >
                    (uint64_t)arm.K + ImplValue->ReceiveOverheadCap ||
                recovered[invocation] != arm.Source)
            {
                batch_valid = false;
                continue;
            }
            result.VerificationMask |= UINT32_C(1) << invocation;
            result.DecodedOverheads[invocation] =
                success_packet_counts[invocation] - arm.K;
        }
        if (!batch_valid) {
            result.ElapsedNanoseconds = 0u;
        }
        return result;
#else
        result.Results[0] = Wirehair_UnsupportedPlatform;
        return result;
#endif
    }
    catch (const std::bad_alloc&) {
        result.Results[failure_invocation] = Wirehair_OOM;
        return result;
    }
    catch (const std::length_error&) {
        result.Results[failure_invocation] = Wirehair_OOM;
        return result;
    }
}

Batch8ArmResult NativeReceiveFixture::PreflightBatch8() const
{
    return RunBatch8Impl<false>();
}

Batch8ArmResult NativeReceiveFixture::RunBatch8() const
{
    return RunBatch8Impl<true>();
}

template <bool ReadCounter>
Batch128CounterArmResult
NativeReceiveFixture::RunWh1Batch128CounterImpl(
    const NativeCounterBracket* bracket) const
{
    Batch128CounterArmResult result;
    NativeCounterBracket frozen_bracket;
    if (!IsInitialized()) {
        return result;
    }

    const NativeArm::Impl& arm = *ImplValue->ArmState;
    // Reject WH2 and invalid callbacks before any caller-controlled read.
    if (arm.Kind != NativeArmKind::Wirehair1)
    {
        result.Results[0] = Wirehair_InvalidInput;
        return result;
    }
    if (ReadCounter && (!bracket || !bracket->Read))
    {
        result.Results[0] = Wirehair_InvalidInput;
        return result;
    }
    if (ReadCounter) {
        frozen_bracket = *bracket;
    }

    size_t failure_invocation = 0u;
    try
    {
        const size_t receive_packet_count =
            (size_t)arm.K + ImplValue->ReceiveOverheadCap;
        std::array<LegacyCodecOwner, kNativeCounterBatchSize> decoders;
        std::array<std::vector<uint8_t>, kNativeCounterBatchSize> recovered;
        std::array<uint32_t, kNativeCounterBatchSize>
            success_packet_counts = {};

        // All decoder construction and output allocation precede Begin.
        for (size_t invocation = 0u;
             invocation < kNativeCounterBatchSize;
             ++invocation)
        {
            failure_invocation = invocation;
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
            TriggerNativeCounterBatchAllocationFailureForTesting();
#endif
            recovered[invocation].assign(arm.Source.size(), 0u);
            const WirehairResult create_result = wirehair_decoder_create_ex(
                nullptr,
                (uint64_t)arm.Source.size(),
                arm.BlockBytes,
                &decoders[invocation].Codec);
            if (create_result != Wirehair_Success)
            {
                result.Results[invocation] = create_result;
                return result;
            }
        }

        if (ReadCounter)
        {
            uint64_t counter_start = 0u;
            try
            {
                if (!frozen_bracket.Read(
                        frozen_bracket.Context,
                        frozen_bracket.Selector,
                        &counter_start))
                {
                    return result;
                }
            }
            catch (...)
            {
                return result;
            }
            result.CounterStart = counter_start;
            result.BeginReadSucceeded = true;
        }

        bool core_threw = false;
        try
        {
            for (size_t invocation = 0u;
                 invocation < kNativeCounterBatchSize;
                 ++invocation)
            {
                failure_invocation = invocation;
                WirehairResult receive_result = Wirehair_NeedMore;
                for (size_t packet_index = 0u;
                     packet_index < receive_packet_count;
                     ++packet_index)
                {
                    receive_result = wirehair_decode(
                        decoders[invocation].Codec,
                        ImplValue->PacketIds[packet_index],
                        ImplValue->PacketStorage.data() +
                            packet_index * arm.BlockBytes,
                        arm.BlockBytes);
                    if (receive_result == Wirehair_Success)
                    {
                        success_packet_counts[invocation] =
                            (uint32_t)packet_index + 1u;
                        receive_result = wirehair_recover(
                            decoders[invocation].Codec,
                            recovered[invocation].data(),
                            recovered[invocation].size());
                        break;
                    }
                    if (receive_result != Wirehair_NeedMore) {
                        break;
                    }
                }
                result.Results[invocation] = receive_result;
            }
        }
        catch (const std::bad_alloc&)
        {
            result.Results[failure_invocation] = Wirehair_OOM;
            core_threw = true;
        }
        catch (const std::length_error&)
        {
            result.Results[failure_invocation] = Wirehair_OOM;
            core_threw = true;
        }
        catch (...)
        {
            result.Results[failure_invocation] = Wirehair_Error;
            core_threw = true;
        }

        if (ReadCounter)
        {
            uint64_t counter_end = 0u;
            try
            {
                if (!frozen_bracket.Read(
                        frozen_bracket.Context,
                        frozen_bracket.Selector,
                        &counter_end))
                {
                    return result;
                }
            }
            catch (...)
            {
                return result;
            }
            result.CounterEnd = counter_end;
            result.EndReadSucceeded = true;
            if (counter_end > result.CounterStart) {
                result.CounterDelta = counter_end - result.CounterStart;
            }
        }
        if (core_threw) {
            return result;
        }

        for (size_t invocation = 0u;
             invocation < kNativeCounterBatchSize;
             ++invocation)
        {
            if (result.Results[invocation] != Wirehair_Success ||
                success_packet_counts[invocation] < arm.K ||
                success_packet_counts[invocation] >
                    (uint64_t)arm.K + ImplValue->ReceiveOverheadCap ||
                recovered[invocation] != arm.Source)
            {
                continue;
            }
            result.VerificationMask[invocation >> 6] |=
                UINT64_C(1) << (invocation & 63u);
            result.DecodedOverheads[invocation] =
                success_packet_counts[invocation] - arm.K;
        }
        return result;
    }
    catch (const std::bad_alloc&)
    {
        result.Results[failure_invocation] = Wirehair_OOM;
        return result;
    }
    catch (const std::length_error&)
    {
        result.Results[failure_invocation] = Wirehair_OOM;
        return result;
    }
}

Batch128CounterArmResult NativeReceiveFixture::PreflightBatch128() const
{
    return RunWh1Batch128CounterImpl<false>(nullptr);
}

Batch128CounterArmResult NativeReceiveFixture::RunWh1Batch128Counter(
    const NativeCounterBracket& bracket) const
{
    return RunWh1Batch128CounterImpl<true>(&bracket);
}

struct NativeSolveFixture::Impl
{
    std::shared_ptr<const NativeArm::Impl> ArmState;
    std::vector<uint32_t> PacketIds;
    std::vector<uint8_t> PacketStorage;
    std::vector<wirehair_v2::SolvePacket> Packets;
    bool Initialized = false;
};

NativeSolveFixture::NativeSolveFixture()
{
}

NativeSolveFixture::~NativeSolveFixture()
{
}

NativeSolveFixture::NativeSolveFixture(
    NativeSolveFixture&& other) noexcept = default;
NativeSolveFixture& NativeSolveFixture::operator=(
    NativeSolveFixture&& other) noexcept = default;

WirehairResult NativeSolveFixture::Initialize(
    const NativeArm& arm,
    const std::vector<uint32_t>& packet_ids,
    uint32_t fixed_received_overhead)
{
    if (!arm.IsInitialized() ||
        arm.Kind() == NativeArmKind::Wirehair1 ||
        arm.SystematicEmission() !=
            NativeSystematicEmission::EquationEvaluation ||
        fixed_received_overhead != kRecoveryOverheadCap)
    {
        return Wirehair_InvalidInput;
    }
    try
    {
        if (!UniquePacketPrefix(
                arm.ImplValue->K,
                fixed_received_overhead,
                packet_ids,
                false))
        {
            return Wirehair_InvalidInput;
        }
        std::unique_ptr<Impl> next(new Impl);
        next->ArmState = arm.ImplValue;
        const size_t received =
            (size_t)next->ArmState->K + fixed_received_overhead;
        next->PacketIds.assign(
            packet_ids.begin(), packet_ids.begin() + received);
        const WirehairResult encode_result =
            EncodePrefix(arm, next->PacketIds, next->PacketStorage);
        if (encode_result != Wirehair_Success) {
            return encode_result;
        }

        next->Packets.reserve(received);
        for (size_t packet_index = 0u;
             packet_index < received;
             ++packet_index)
        {
            wirehair_v2::SolvePacket packet;
            packet.BlockId = next->PacketIds[packet_index];
            packet.Data = next->PacketStorage.data() +
                (size_t)packet_index * next->ArmState->BlockBytes;
            next->Packets.push_back(packet);
        }
        next->Initialized = true;
        ImplValue.swap(next);
        return Wirehair_Success;
    }
    catch (const std::bad_alloc&) {
        return Wirehair_OOM;
    }
    catch (const std::length_error&) {
        return Wirehair_OOM;
    }
}

bool NativeSolveFixture::IsInitialized() const
{
    return ImplValue && ImplValue->Initialized;
}

IsolatedSolveResult NativeSolveFixture::Run() const
{
    IsolatedSolveResult result;
    if (!IsInitialized()) {
        return result;
    }

    try
    {
        const NativeArm::Impl& arm = *ImplValue->ArmState;
        std::vector<uint8_t> intermediate;
        wirehair_v2::PrecodeSolveStats stats;
        const std::chrono::steady_clock::time_point start =
            std::chrono::steady_clock::now();
        const WirehairResult solve_result =
            wirehair_v2::SolvePrecodeSystemForValidatedSystemWithRuntime(
                arm.System,
                arm.PacketConfig,
                arm.PacketRuntime,
                ImplValue->Packets,
                arm.BlockBytes,
                intermediate,
                &stats);
        const std::chrono::steady_clock::time_point finish =
            std::chrono::steady_clock::now();
        const std::chrono::nanoseconds elapsed =
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                finish - start);
        result.ElapsedNanoseconds = elapsed.count() > 0 ?
            (uint64_t)elapsed.count() : 0u;
        result.Result = solve_result;
        result.Stats = stats;
        if (solve_result == Wirehair_Success)
        {
            result.BytesVerified = VerifyWh2DecodedBytes(
                arm.System,
                arm.PacketConfig,
                arm.PacketRuntime,
                arm.BlockBytes,
                arm.Source,
                intermediate);
            if (!result.BytesVerified || result.ElapsedNanoseconds == 0u) {
                result.Result = Wirehair_Error;
            }
        }
        if (result.Result != Wirehair_Success) {
            result.ElapsedNanoseconds = 0u;
            result.BytesVerified = false;
        }
        return result;
    }
    catch (const std::bad_alloc&) {
        result.Result = Wirehair_OOM;
        result.ElapsedNanoseconds = 0u;
        result.BytesVerified = false;
        return result;
    }
    catch (const std::length_error&) {
        result.Result = Wirehair_OOM;
        result.ElapsedNanoseconds = 0u;
        result.BytesVerified = false;
        return result;
    }
}

} // namespace wirehair_wh2_bench
