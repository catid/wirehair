#include "Wh2NativeCodec.h"

#include <algorithm>
#include <chrono>
#include <cstring>
#include <limits>
#include <new>
#include <stdexcept>
#include <utility>

namespace wirehair_wh2_bench {
namespace {

static const uint32_t kMaxConstructionAttempt = 255u;
static const uint32_t kRecoveryOverheadCap = 4u;

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
            spec.BaseKind == NativeWh2BaseKind::ProductionProfile;
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
    return params.BlockCount == expected_block_count &&
        params.BlockCount >= 2u && params.BlockCount <= 64000u &&
        params.Staircase != 0u &&
        params.SourceHits >= 1u && params.SourceHits <= 8u &&
        params.DenseRows <= 64u && params.HeavyRows <= 128u &&
        (params.HeavyFamily ==
             wirehair_v2::HeavyCoefficientFamily::PeriodicCauchy ||
         params.HeavyFamily ==
             wirehair_v2::HeavyCoefficientFamily::HashedNonzero) &&
        binary_span <= UINT16_MAX && total_span <= UINT16_MAX &&
        (!params.DenseIdentityCorner ||
         known_span >= 2u * (uint64_t)(params.DenseRows >> 1));
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
    const std::vector<uint32_t>& packet_ids)
{
    if (packet_ids.size() != (size_t)block_count + kRecoveryOverheadCap) {
        return false;
    }
    std::vector<uint32_t> sorted(packet_ids);
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

bool RecoverWh2DecodedBytesInto(
    const wirehair_v2::PrecodeSystem& system,
    const wirehair_v2::PacketRowConfig& config,
    const wirehair_v2::PacketRowRuntime& runtime,
    uint32_t block_bytes,
    const std::vector<uint8_t>& intermediate,
    uint8_t* recovered_out,
    size_t recovered_bytes)
{
    const uint32_t K = system.Params.BlockCount;
    size_t expected_source_bytes = 0u;
    size_t expected_intermediate_bytes = 0u;
    if (!recovered_out ||
        !CheckedProduct(K, block_bytes, expected_source_bytes) ||
        !CheckedProduct(
            (uint64_t)K + PrecodeCount(system),
            block_bytes,
            expected_intermediate_bytes) ||
        recovered_bytes != expected_source_bytes ||
        intermediate.size() != expected_intermediate_bytes)
    {
        return false;
    }
    for (uint32_t block_id = 0u; block_id < K; ++block_id)
    {
        if (!wirehair_v2::EvaluatePacketBlockForValidatedSystemWithRuntime(
                system,
                config,
                runtime,
                intermediate.data(),
                block_bytes,
                block_id,
                recovered_out + (size_t)block_id * block_bytes))
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
    std::vector<uint8_t> Source;
    WirehairCodec LegacyEncoder = nullptr;
    wirehair_v2::PrecodeSystem System = {};
    wirehair_v2::PacketRowConfig PacketConfig = {};
    wirehair_v2::PacketRowRuntime PacketRuntime = {};
    std::vector<uint8_t> Intermediate;
    bool Initialized = false;
};

NativeArmSpec MakeWirehair1Arm()
{
    NativeArmSpec spec;
    spec.Kind = NativeArmKind::Wirehair1;
    return spec;
}

NativeArmSpec MakeCertifiedWh2Arm(uint32_t construction_attempt)
{
    NativeArmSpec spec;
    spec.Kind = NativeArmKind::Wirehair2Certified;
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
    spec.BaseKind = base_kind;
    spec.ConstructionAttempt = construction_attempt;
    spec.Transform = transform;
    spec.TransformContext = transform_context;
    if (options) {
        spec.Wh2Options = *options;
    }
    return spec;
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
                return Wirehair_InvalidInput;
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
            return Wirehair_InvalidInput;
        }
        if (spec.Transform &&
            !spec.Transform(
                block_count,
                block_bytes,
                base_params,
                base_packet_config,
                spec.TransformContext))
        {
            return Wirehair_InvalidInput;
        }

        const wirehair_v2::PrecodeParams params =
            wirehair_v2::PrecodeParamsForAttempt(
                base_params, spec.ConstructionAttempt);
        const wirehair_v2::PacketRowConfig packet_config =
            wirehair_v2::PacketConfigForAttempt(
                base_packet_config, spec.ConstructionAttempt);
        if (!ValidPrecodeParamsShape(params, block_count)) {
            return Wirehair_InvalidInput;
        }

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
        if (!UniquePacketPrefix(arm.ImplValue->K, packet_ids)) {
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

struct NativeEncoderFixture::Impl
{
    NativeArmSpec Spec = {};
    uint32_t K = 0u;
    uint32_t BlockBytes = 0u;
    std::vector<uint8_t> Source;
    bool Initialized = false;
};

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

TimedArmResult NativeEncoderFixture::Run() const
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
        const std::chrono::steady_clock::time_point start =
            std::chrono::steady_clock::now();
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
        const std::chrono::steady_clock::time_point finish =
            std::chrono::steady_clock::now();
        const std::chrono::nanoseconds elapsed =
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                finish - start);
        result.ElapsedNanoseconds = elapsed.count() > 0 ?
            (uint64_t)elapsed.count() : 0u;
        result.Result = invocation_result;
        if (invocation_result == Wirehair_Success)
        {
            result.BytesVerified = encoded == ImplValue->Source;
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
        return result;
    }
    catch (const std::length_error&) {
        result.Result = Wirehair_OOM;
        return result;
    }
}

struct NativeReceiveFixture::Impl
{
    std::shared_ptr<const NativeArm::Impl> ArmState;
    std::vector<uint32_t> PacketIds;
    std::vector<uint8_t> PacketStorage;
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
    const std::vector<uint32_t>& packet_ids)
{
    if (!arm.IsInitialized()) {
        return Wirehair_InvalidInput;
    }
    try
    {
        if (!UniquePacketPrefix(arm.ImplValue->K, packet_ids)) {
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

bool NativeReceiveFixture::IsInitialized() const
{
    return ImplValue && ImplValue->Initialized;
}

TimedArmResult NativeReceiveFixture::Run() const
{
    TimedArmResult result;
    if (!IsInitialized()) {
        return result;
    }

    try
    {
        const NativeArm::Impl& arm = *ImplValue->ArmState;
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

            const std::chrono::steady_clock::time_point start =
                std::chrono::steady_clock::now();
            for (uint32_t packet_index = 0u;
                 packet_index < arm.K + kRecoveryOverheadCap;
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
                    success_packet_count = packet_index + 1u;
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
            const std::chrono::steady_clock::time_point finish =
                std::chrono::steady_clock::now();
            const std::chrono::nanoseconds elapsed =
                std::chrono::duration_cast<std::chrono::nanoseconds>(
                    finish - start);
            result.ElapsedNanoseconds = elapsed.count() > 0 ?
                (uint64_t)elapsed.count() : 0u;
        }
        else
        {
            // This is the benchmark-local fresh decoder state.  Every storage
            // reservation and the recovery destination are prepared before the
            // clock; solver/resume allocations remain decoder work.
            std::vector<wirehair_v2::SolvePacket> received_packets;
            received_packets.reserve((size_t)arm.K + kRecoveryOverheadCap);
            wirehair_v2::PrecodeSolveResumeState resume_state;
            wirehair_v2::PrecodeSolveResumeState cold_resume;
            std::vector<uint8_t> intermediate;

            const std::chrono::steady_clock::time_point start =
                std::chrono::steady_clock::now();
            for (uint32_t packet_index = 0u;
                 packet_index < arm.K + kRecoveryOverheadCap;
                 ++packet_index)
            {
                wirehair_v2::SolvePacket packet;
                packet.BlockId = ImplValue->PacketIds[packet_index];
                packet.Data = ImplValue->PacketStorage.data() +
                    (size_t)packet_index * arm.BlockBytes;
                received_packets.push_back(packet);
                if (received_packets.size() < arm.K) {
                    continue;
                }

                if (resume_state.Active)
                {
                    receive_result = wirehair_v2::ResumePrecodeSystem(
                        arm.System,
                        arm.PacketConfig,
                        packet.BlockId,
                        packet.Data,
                        arm.BlockBytes,
                        resume_state,
                        intermediate);
                }
                else
                {
                    cold_resume.Clear();
                    receive_result = wirehair_v2::
                        SolvePrecodeSystemForValidatedSystemWithRuntime(
                            arm.System,
                            arm.PacketConfig,
                            arm.PacketRuntime,
                            received_packets,
                            arm.BlockBytes,
                            intermediate,
                            nullptr,
                            &cold_resume);
                    if (receive_result == Wirehair_NeedMore &&
                        cold_resume.Active)
                    {
                        resume_state.Swap(cold_resume);
                    }
                }

                if (receive_result == Wirehair_Success)
                {
                    success_packet_count = packet_index + 1u;
                    if (!RecoverWh2DecodedBytesInto(
                            arm.System,
                            arm.PacketConfig,
                            arm.PacketRuntime,
                            arm.BlockBytes,
                            intermediate,
                            recovered.data(),
                            recovered.size()))
                    {
                        receive_result = Wirehair_Error;
                    }
                    break;
                }
                if (receive_result != Wirehair_NeedMore) {
                    break;
                }
            }
            const std::chrono::steady_clock::time_point finish =
                std::chrono::steady_clock::now();
            const std::chrono::nanoseconds elapsed =
                std::chrono::duration_cast<std::chrono::nanoseconds>(
                    finish - start);
            result.ElapsedNanoseconds = elapsed.count() > 0 ?
                (uint64_t)elapsed.count() : 0u;
        }

        result.Result = receive_result;
        if (receive_result == Wirehair_Success)
        {
            if (success_packet_count < arm.K ||
                success_packet_count > arm.K + kRecoveryOverheadCap)
            {
                result.Result = Wirehair_Error;
                result.ElapsedNanoseconds = 0u;
                return result;
            }
            result.DecodedOverhead = success_packet_count - arm.K;
            result.BytesVerified = recovered == arm.Source;
            if (!result.BytesVerified || result.ElapsedNanoseconds == 0u) {
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
        fixed_received_overhead > kRecoveryOverheadCap)
    {
        return Wirehair_InvalidInput;
    }
    try
    {
        if (!UniquePacketPrefix(arm.ImplValue->K, packet_ids)) {
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

        const uint32_t received = next->ArmState->K +
            fixed_received_overhead;
        next->Packets.reserve(received);
        for (uint32_t packet_index = 0u;
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
