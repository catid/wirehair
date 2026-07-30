/*
    Test-hook-only repair-v1 timing receipt.

    This file is included from WirehairV2Bench.cpp after the bandtiming
    implementation so it can reuse that command's private, independently
    authenticated loss, payload, timing, and hashing helpers without changing
    a byte of the frozen bandtiming v2 protocol.
*/

enum class RepairTimingRootDerivation
{
    Derived,
    Fixed
};

struct RepairTimingOptions
{
    uint32_t BlockCount = 0u;
    uint32_t BlockBytes = 0u;
    const wirehair_v2::RepairV1Contract* Contract = nullptr;
    const wirehair_v2::V2EquationContract* DispatchContract = nullptr;
    RepairTimingRootDerivation RootDerivation =
        RepairTimingRootDerivation::Derived;
    uint32_t ConstructionSeedBase = 0u;
    double Loss = 0.0;
    uint64_t LossSeedBase = 0u;
    PacketScheduleKind Schedule = PacketScheduleKind::Iid;
    uint32_t WarmupReplicates = 0u;
    uint32_t Replicates = 0u;
    uint32_t InnerReps = 0u;
    uint32_t MaxOverhead = 0u;
    PeelTimingCacheState CacheState = PeelTimingCacheState::Warm;
    bool SystematicCache = false;
    uint64_t EvictBytes = 0u;
    std::string ContextSha256;
    double RequiredMargin = 0.0;
    std::string DescriptorFault;
    bool TimedScopeFault = false;
    bool ForcedWitnessFault = false;
};

const char* RepairTimingRootDerivationName(
    RepairTimingRootDerivation derivation)
{
    return derivation == RepairTimingRootDerivation::Fixed ?
        "fixed_base_u32_v1" :
        "base_xor_d192ed03_times_rep_plus_1_mod2^32_v1";
}

uint32_t RepairTimingConstructionRoot(
    const RepairTimingOptions& options,
    uint32_t cell)
{
    return options.RootDerivation == RepairTimingRootDerivation::Fixed ?
        options.ConstructionSeedBase :
        BandTimingDerivedConstructionSeed(
            options.ConstructionSeedBase, cell);
}

bool ParseRepairTimingOptions(
    int argc,
    char** argv,
    RepairTimingOptions& options)
{
    bool have_N = false;
    bool have_bb = false;
    bool have_arm = false;
    bool have_policy = false;
    bool have_dispatch = false;
    bool have_derivation = false;
    bool have_construction_seed = false;
    bool have_loss = false;
    bool have_loss_seed = false;
    bool have_schedule = false;
    bool have_warmups = false;
    bool have_replicates = false;
    bool have_inner_reps = false;
    bool have_max_overhead = false;
    bool have_cache_state = false;
    bool have_systematic_cache = false;
    bool have_evict_bytes = false;
    bool have_context = false;
    bool have_margin = false;
    bool have_descriptor_fault = false;
    bool have_timed_scope_fault = false;
    bool have_forced_witness_fault = false;

    for (int i = 0; i < argc; ++i)
    {
        const char* value = nullptr;
        if (!std::strcmp(argv[i], "--N"))
        {
            if (!AcceptOptionOnce("repairtiming", "--N", have_N) ||
                !TakeArg("repairtiming", "--N", argc, argv, i, value) ||
                !ParseU32Arg("--N", value, options.BlockCount))
            {
                return false;
            }
        }
        else if (!std::strcmp(argv[i], "--bb"))
        {
            if (!AcceptOptionOnce("repairtiming", "--bb", have_bb) ||
                !TakeArg("repairtiming", "--bb", argc, argv, i, value) ||
                !ParseU32Arg("--bb", value, options.BlockBytes))
            {
                return false;
            }
        }
        else if (!std::strcmp(argv[i], "--repair-arm"))
        {
            if (!AcceptOptionOnce(
                    "repairtiming", "--repair-arm", have_arm) ||
                !TakeArg(
                    "repairtiming", "--repair-arm",
                    argc, argv, i, value))
            {
                return false;
            }
            options.Contract =
                wirehair_v2::FindRepairV1Contract(value);
            if (!options.Contract)
            {
                std::fprintf(stderr,
                    "repairtiming --repair-arm must be "
                    "pure8_s0m1_d3_repair_v1 or "
                    "pure9_s0m1_d3_repair_v1\n");
                return false;
            }
        }
        else if (!std::strcmp(argv[i], "--repair-policy"))
        {
            if (!AcceptOptionOnce(
                    "repairtiming", "--repair-policy", have_policy) ||
                !TakeArg(
                    "repairtiming", "--repair-policy",
                    argc, argv, i, value) ||
                std::strcmp(value, "repair-v1") != 0)
            {
                std::fprintf(stderr,
                    "repairtiming --repair-policy must be repair-v1\n");
                return false;
            }
        }
        else if (!std::strcmp(argv[i], "--dispatch-profile"))
        {
            if (!AcceptOptionOnce(
                    "repairtiming", "--dispatch-profile", have_dispatch) ||
                !TakeArg(
                    "repairtiming", "--dispatch-profile",
                    argc, argv, i, value))
            {
                return false;
            }
            options.DispatchContract =
                wirehair_v2::FindV2EquationContract(value);
            if (!options.DispatchContract ||
                options.DispatchContract->PublicProfile ||
                options.DispatchContract->ProfileId !=
                    wirehair_v2::kWh2DispatchV1ContractId)
            {
                std::fprintf(stderr,
                    "repairtiming --dispatch-profile must be dispatch-v1\n");
                return false;
            }
        }
        else if (!std::strcmp(
                    argv[i], "--construction-seed-derivation"))
        {
            if (!AcceptOptionOnce(
                    "repairtiming", "--construction-seed-derivation",
                    have_derivation) ||
                !TakeArg(
                    "repairtiming", "--construction-seed-derivation",
                    argc, argv, i, value))
            {
                return false;
            }
            if (!std::strcmp(value, "derived")) {
                options.RootDerivation =
                    RepairTimingRootDerivation::Derived;
            }
            else if (!std::strcmp(value, "fixed")) {
                options.RootDerivation =
                    RepairTimingRootDerivation::Fixed;
            }
            else
            {
                std::fprintf(stderr,
                    "repairtiming --construction-seed-derivation must "
                    "be derived or fixed\n");
                return false;
            }
        }
        else if (!std::strcmp(argv[i], "--construction-seed"))
        {
            if (!AcceptOptionOnce(
                    "repairtiming", "--construction-seed",
                    have_construction_seed) ||
                !TakeArg(
                    "repairtiming", "--construction-seed",
                    argc, argv, i, value) ||
                !ParseU32Arg(
                    "--construction-seed", value,
                    options.ConstructionSeedBase))
            {
                return false;
            }
        }
        else if (!std::strcmp(argv[i], "--loss"))
        {
            if (!AcceptOptionOnce("repairtiming", "--loss", have_loss) ||
                !TakeArg(
                    "repairtiming", "--loss", argc, argv, i, value) ||
                !ParseDoubleArg("--loss", value, options.Loss))
            {
                return false;
            }
        }
        else if (!std::strcmp(argv[i], "--loss-seed"))
        {
            if (!AcceptOptionOnce(
                    "repairtiming", "--loss-seed", have_loss_seed) ||
                !TakeArg(
                    "repairtiming", "--loss-seed",
                    argc, argv, i, value) ||
                !ParseU64Arg(
                    "--loss-seed", value, options.LossSeedBase))
            {
                return false;
            }
        }
        else if (!std::strcmp(argv[i], "--schedule"))
        {
            if (!AcceptOptionOnce(
                    "repairtiming", "--schedule", have_schedule) ||
                !TakeArg(
                    "repairtiming", "--schedule",
                    argc, argv, i, value) ||
                !ParsePacketSchedule(value, options.Schedule))
            {
                std::fprintf(stderr,
                    "repairtiming --schedule must be iid, burst, "
                    "permutation, systematic-first, repair-only, or "
                    "adversarial\n");
                return false;
            }
        }
        else if (!std::strcmp(argv[i], "--warmup-replicates"))
        {
            if (!AcceptOptionOnce(
                    "repairtiming", "--warmup-replicates", have_warmups) ||
                !TakeArg(
                    "repairtiming", "--warmup-replicates",
                    argc, argv, i, value) ||
                !ParseU32Arg(
                    "--warmup-replicates", value,
                    options.WarmupReplicates))
            {
                return false;
            }
        }
        else if (!std::strcmp(argv[i], "--replicates"))
        {
            if (!AcceptOptionOnce(
                    "repairtiming", "--replicates", have_replicates) ||
                !TakeArg(
                    "repairtiming", "--replicates",
                    argc, argv, i, value) ||
                !ParseU32Arg(
                    "--replicates", value, options.Replicates))
            {
                return false;
            }
        }
        else if (!std::strcmp(argv[i], "--inner-reps"))
        {
            if (!AcceptOptionOnce(
                    "repairtiming", "--inner-reps", have_inner_reps) ||
                !TakeArg(
                    "repairtiming", "--inner-reps",
                    argc, argv, i, value) ||
                !ParseU32Arg(
                    "--inner-reps", value, options.InnerReps))
            {
                return false;
            }
        }
        else if (!std::strcmp(argv[i], "--max-overhead"))
        {
            if (!AcceptOptionOnce(
                    "repairtiming", "--max-overhead",
                    have_max_overhead) ||
                !TakeArg(
                    "repairtiming", "--max-overhead",
                    argc, argv, i, value) ||
                !ParseU32Arg(
                    "--max-overhead", value, options.MaxOverhead))
            {
                return false;
            }
        }
        else if (!std::strcmp(argv[i], "--cache-state"))
        {
            if (!AcceptOptionOnce(
                    "repairtiming", "--cache-state", have_cache_state) ||
                !TakeArg(
                    "repairtiming", "--cache-state",
                    argc, argv, i, value))
            {
                return false;
            }
            if (!std::strcmp(value, "cold")) {
                options.CacheState = PeelTimingCacheState::Cold;
            }
            else if (!std::strcmp(value, "warm")) {
                options.CacheState = PeelTimingCacheState::Warm;
            }
            else
            {
                std::fprintf(stderr,
                    "repairtiming --cache-state must be cold or warm\n");
                return false;
            }
        }
        else if (!std::strcmp(argv[i], "--systematic-cache"))
        {
            if (!AcceptOptionOnce(
                    "repairtiming", "--systematic-cache",
                    have_systematic_cache) ||
                !TakeArg(
                    "repairtiming", "--systematic-cache",
                    argc, argv, i, value))
            {
                return false;
            }
            if (!std::strcmp(value, "off")) {
                options.SystematicCache = false;
            }
            else if (!std::strcmp(value, "on")) {
                options.SystematicCache = true;
            }
            else
            {
                std::fprintf(stderr,
                    "repairtiming --systematic-cache must be off or on\n");
                return false;
            }
        }
        else if (!std::strcmp(argv[i], "--evict-bytes"))
        {
            if (!AcceptOptionOnce(
                    "repairtiming", "--evict-bytes", have_evict_bytes) ||
                !TakeArg(
                    "repairtiming", "--evict-bytes",
                    argc, argv, i, value) ||
                !ParseU64Arg(
                    "--evict-bytes", value, options.EvictBytes))
            {
                return false;
            }
        }
        else if (!std::strcmp(argv[i], "--context-sha256"))
        {
            if (!AcceptOptionOnce(
                    "repairtiming", "--context-sha256", have_context) ||
                !TakeArg(
                    "repairtiming", "--context-sha256",
                    argc, argv, i, value))
            {
                return false;
            }
            options.ContextSha256 = value;
        }
        else if (!std::strcmp(argv[i], "--required-margin"))
        {
            if (!AcceptOptionOnce(
                    "repairtiming", "--required-margin", have_margin) ||
                !TakeArg(
                    "repairtiming", "--required-margin",
                    argc, argv, i, value) ||
                !ParseDoubleArg(
                    "--required-margin", value,
                    options.RequiredMargin))
            {
                return false;
            }
        }
        else if (!std::strcmp(argv[i], "--test-descriptor-fault"))
        {
            if (!AcceptOptionOnce(
                    "repairtiming", "--test-descriptor-fault",
                    have_descriptor_fault) ||
                !TakeArg(
                    "repairtiming", "--test-descriptor-fault",
                    argc, argv, i, value))
            {
                return false;
            }
            options.DescriptorFault = value;
            if (options.DescriptorFault != "bad-magic" &&
                options.DescriptorFault != "bad-id" &&
                options.DescriptorFault != "bad-attempt" &&
                options.DescriptorFault != "trailing")
            {
                std::fprintf(stderr,
                    "repairtiming --test-descriptor-fault must be "
                    "bad-magic, bad-id, bad-attempt, or trailing\n");
                return false;
            }
        }
        else if (!std::strcmp(argv[i], "--test-timed-scope-drift"))
        {
            if (!AcceptOptionOnce(
                    "repairtiming", "--test-timed-scope-drift",
                    have_timed_scope_fault) ||
                !TakeArg(
                    "repairtiming", "--test-timed-scope-drift",
                    argc, argv, i, value) ||
                std::strcmp(value, "hash") != 0)
            {
                std::fprintf(stderr,
                    "repairtiming --test-timed-scope-drift must be "
                    "hash\n");
                return false;
            }
            options.TimedScopeFault = true;
        }
        else if (!std::strcmp(
                    argv[i], "--test-forced-witness-drift"))
        {
            if (!AcceptOptionOnce(
                    "repairtiming", "--test-forced-witness-drift",
                    have_forced_witness_fault) ||
                !TakeArg(
                    "repairtiming", "--test-forced-witness-drift",
                    argc, argv, i, value) ||
                std::strcmp(value, "intermediate") != 0)
            {
                std::fprintf(stderr,
                    "repairtiming --test-forced-witness-drift must be "
                    "intermediate\n");
                return false;
            }
            options.ForcedWitnessFault = true;
        }
        else {
            return UnknownArg("repairtiming", argv[i]);
        }
    }

    if (!have_N || !have_bb || !have_arm || !have_policy ||
        !have_dispatch || !have_derivation ||
        !have_construction_seed || !have_loss || !have_loss_seed ||
        !have_schedule || !have_warmups || !have_replicates ||
        !have_inner_reps || !have_max_overhead || !have_cache_state ||
        !have_systematic_cache || !have_evict_bytes || !have_context ||
        !have_margin)
    {
        std::fprintf(stderr,
            "repairtiming requires --N, --bb, --repair-arm, "
            "--repair-policy, --dispatch-profile, "
            "--construction-seed-derivation, --construction-seed, "
            "--loss, --loss-seed, --schedule, --warmup-replicates, "
            "--replicates, --inner-reps, --max-overhead, --cache-state, "
            "--systematic-cache, --evict-bytes, --context-sha256, and "
            "--required-margin\n");
        return false;
    }

    const uint64_t cells =
        (uint64_t)options.WarmupReplicates + options.Replicates;
    const uint64_t delivered =
        (uint64_t)options.BlockCount + options.MaxOverhead;
    const long double work =
        (long double)cells * 112.0L *
        (long double)options.InnerReps * (long double)delivered;
    if (options.BlockCount < 2u || options.BlockCount > 100u ||
        options.BlockBytes < 2u || options.BlockBytes > 4096u ||
        (options.BlockBytes & 1u) != 0u ||
        !std::isfinite(options.Loss) ||
        options.Loss < 0.0 ||
        options.Loss > 0.99 ||
        (options.Loss == 0.0 && std::signbit(options.Loss)) ||
        options.Replicates < 3u ||
        cells == 0u || cells > 768u ||
        options.InnerReps == 0u || options.InnerReps > 1024u ||
        (options.CacheState == PeelTimingCacheState::Cold &&
         options.InnerReps != 1u) ||
        options.MaxOverhead != 64u ||
        delivered > UINT32_MAX ||
        options.EvictBytes < 4096u ||
        options.EvictBytes > kPeelTimingEvictBytesMax ||
        !IsPeelTimingLowerHexSha256(options.ContextSha256) ||
        !std::isfinite(options.RequiredMargin) ||
        options.RequiredMargin < 0.0 ||
        options.RequiredMargin > 1.0 ||
        (options.RequiredMargin == 0.0 &&
         std::signbit(options.RequiredMargin)) ||
        work > 1000000000000.0L)
    {
        std::fprintf(stderr, "repairtiming argument domain mismatch\n");
        return false;
    }
    return true;
}

bool ValidateRepairTimingEnvironment()
{
    static const char* const kForbidden[] = {
        "WIREHAIR_V2_PEEL_DEGREES",
        "WIREHAIR_V2_STAIRCASE_DEGREES",
        "WIREHAIR_V2_STAIRCASE_ROW_DEGREES",
        "WIREHAIR_V2_STAIRCASE_DEGREE_SCALE",
        "WIREHAIR_V2_BAND_TRACKING_X"
    };
    for (const char* name : kForbidden)
    {
        const wirehair::EnvironmentValue value(name);
        if (value.IsSet())
        {
            std::fprintf(stderr,
                "repairtiming forbids ambient %s\n", name);
            return false;
        }
    }
    return ValidateTargetExperimentEnvironment("repairtiming");
}

struct RepairTimingDescriptor
{
    const wirehair_v2::RepairV1Contract* Contract = nullptr;
    uint32_t ConstructionRoot = 0u;
    uint32_t SelectedAttempt = wirehair_v2::kRepairV1NoAttempt;
    std::vector<uint8_t> Bytes;
    std::string Hex;
    std::string Sha256;
};

void RepairTimingAppendLe64(
    std::vector<uint8_t>& bytes,
    uint64_t value)
{
    for (unsigned i = 0u; i < 8u; ++i) {
        bytes.push_back((uint8_t)(value >> (8u * i)));
    }
}

void RepairTimingAppendLe32(
    std::vector<uint8_t>& bytes,
    uint32_t value)
{
    for (unsigned i = 0u; i < 4u; ++i) {
        bytes.push_back((uint8_t)(value >> (8u * i)));
    }
}

uint64_t RepairTimingReadLe64(const uint8_t* bytes)
{
    uint64_t value = 0u;
    for (unsigned i = 0u; i < 8u; ++i) {
        value |= (uint64_t)bytes[i] << (8u * i);
    }
    return value;
}

uint32_t RepairTimingReadLe32(const uint8_t* bytes)
{
    uint32_t value = 0u;
    for (unsigned i = 0u; i < 4u; ++i) {
        value |= (uint32_t)bytes[i] << (8u * i);
    }
    return value;
}

std::string RepairTimingHex(const std::vector<uint8_t>& bytes)
{
    static const char kHex[] = "0123456789abcdef";
    std::string text;
    text.reserve(bytes.size() * 2u);
    for (uint8_t value : bytes)
    {
        text.push_back(kHex[value >> 4]);
        text.push_back(kHex[value & 15u]);
    }
    return text;
}

bool BuildRepairTimingDescriptor(
    const wirehair_v2::RepairV1Contract& contract,
    uint32_t construction_root,
    uint32_t selected_attempt,
    RepairTimingDescriptor& descriptor)
{
    if (selected_attempt >= wirehair_v2::kRepairV1AttemptCount ||
        wirehair_v2::FindRepairV1Contract(contract.ProvisionalId) !=
            &contract ||
        (contract.ProvisionalId !=
             UINT64_C(0x19cccf775ce0bf09) &&
         contract.ProvisionalId !=
             UINT64_C(0xa530f9105beaa450)))
    {
        return false;
    }
    RepairTimingDescriptor candidate;
    candidate.Bytes.reserve(17u);
    candidate.Contract = &contract;
    candidate.ConstructionRoot = construction_root;
    candidate.SelectedAttempt = selected_attempt;
    candidate.Bytes.push_back((uint8_t)'W');
    candidate.Bytes.push_back((uint8_t)'2');
    candidate.Bytes.push_back((uint8_t)'R');
    candidate.Bytes.push_back((uint8_t)'3');
    RepairTimingAppendLe64(
        candidate.Bytes, contract.ProvisionalId);
    RepairTimingAppendLe32(
        candidate.Bytes, construction_root);
    candidate.Bytes.push_back((uint8_t)selected_attempt);
    if (candidate.Bytes.size() != 17u) return false;
    candidate.Hex = RepairTimingHex(candidate.Bytes);
    candidate.Sha256 = Sha256HexBytes(
        candidate.Bytes.data(), candidate.Bytes.size());
    if (candidate.Hex.size() != 34u ||
        !IsPeelTimingLowerHexSha256(candidate.Sha256))
    {
        return false;
    }
    descriptor = std::move(candidate);
    return true;
}

bool ParseRepairTimingDescriptor(
    const std::vector<uint8_t>& bytes,
    RepairTimingDescriptor& descriptor)
{
    if (bytes.size() != 17u ||
        bytes[0] != (uint8_t)'W' ||
        bytes[1] != (uint8_t)'2' ||
        bytes[2] != (uint8_t)'R' ||
        bytes[3] != (uint8_t)'3')
    {
        return false;
    }
    const uint64_t id = RepairTimingReadLe64(bytes.data() + 4u);
    const uint32_t root = RepairTimingReadLe32(bytes.data() + 12u);
    const uint32_t attempt = bytes[16];
    const wirehair_v2::RepairV1Contract* contract =
        wirehair_v2::FindRepairV1Contract(id);
    if (!contract ||
        (id != UINT64_C(0x19cccf775ce0bf09) &&
         id != UINT64_C(0xa530f9105beaa450)) ||
        attempt >= wirehair_v2::kRepairV1AttemptCount)
    {
        return false;
    }
    RepairTimingDescriptor candidate;
    candidate.Contract = contract;
    candidate.ConstructionRoot = root;
    candidate.SelectedAttempt = attempt;
    candidate.Bytes = bytes;
    candidate.Hex = RepairTimingHex(bytes);
    candidate.Sha256 =
        Sha256HexBytes(bytes.data(), bytes.size());
    if (candidate.Hex.size() != 34u ||
        !IsPeelTimingLowerHexSha256(candidate.Sha256))
    {
        return false;
    }
    descriptor = std::move(candidate);
    return true;
}

bool ParseRepairTimingDescriptorView(
    const std::vector<uint8_t>& bytes,
    const wirehair_v2::RepairV1Contract*& contract,
    uint32_t& construction_root,
    uint32_t& selected_attempt)
{
    if (bytes.size() != 17u ||
        bytes[0] != (uint8_t)'W' ||
        bytes[1] != (uint8_t)'2' ||
        bytes[2] != (uint8_t)'R' ||
        bytes[3] != (uint8_t)'3')
    {
        return false;
    }
    const uint64_t id = RepairTimingReadLe64(bytes.data() + 4u);
    const wirehair_v2::RepairV1Contract* parsed =
        wirehair_v2::FindRepairV1Contract(id);
    const uint32_t attempt = bytes[16];
    if (!parsed ||
        (id != UINT64_C(0x19cccf775ce0bf09) &&
         id != UINT64_C(0xa530f9105beaa450)) ||
        attempt >= wirehair_v2::kRepairV1AttemptCount)
    {
        return false;
    }
    contract = parsed;
    construction_root =
        RepairTimingReadLe32(bytes.data() + 12u);
    selected_attempt = attempt;
    return true;
}

enum class RepairTimingEncoderKind
{
    RawAttemptZero,
    LazySelector,
    ForcedSelected
};

struct RepairTimingRecoveryEvidence
{
    WirehairResult EncodeResult = Wirehair_InvalidInput;
    BandTimingPayload Payload;
    BandTimingRecovery Recovery;
};

WirehairResult BuildRepairTimingPayload(
    const RepairTimingOptions& options,
    uint32_t construction_root,
    const BandTimingAlignedBytes& message,
    const std::vector<uint32_t>& ids,
    RepairTimingEncoderKind kind,
    uint32_t selected_attempt,
    BandTimingPayload& payload,
    wirehair_v2::RepairV1SelectorTelemetry* telemetry,
    wirehair_v2::PrecodeSolveStats* solve_stats,
    bool output_preallocated = false)
{
    BandTimingAlignedBytes prepared_bytes;
    std::vector<uint32_t> prepared_data_bytes;
    if (output_preallocated)
    {
        prepared_bytes = std::move(payload.Bytes);
        prepared_data_bytes = std::move(payload.DataBytes);
    }
    payload = BandTimingPayload();
    if (output_preallocated)
    {
        payload.Bytes = std::move(prepared_bytes);
        payload.DataBytes = std::move(prepared_data_bytes);
    }
    if (telemetry) {
        *telemetry = wirehair_v2::RepairV1SelectorTelemetry();
    }
    if (solve_stats) {
        *solve_stats = wirehair_v2::PrecodeSolveStats();
    }
    if (!options.Contract ||
        ids.size() >
            std::numeric_limits<size_t>::max() /
                (size_t)options.BlockBytes)
    {
        payload.ConstructResult = Wirehair_InvalidInput;
        payload.Result = Wirehair_InvalidInput;
        return payload.Result;
    }
    const size_t required_payload_bytes =
        ids.size() * (size_t)options.BlockBytes;
    if (output_preallocated)
    {
        if (payload.Bytes.Size() != required_payload_bytes ||
            payload.DataBytes.size() != ids.size())
        {
            payload.ConstructResult = Wirehair_InvalidInput;
            payload.Result = Wirehair_InvalidInput;
            return payload.Result;
        }
    }
    else if (!payload.Bytes.Resize(required_payload_bytes))
    {
        payload.ConstructResult = Wirehair_OOM;
        payload.Result = Wirehair_OOM;
        return payload.Result;
    }
    try
    {
        if (!output_preallocated) {
            payload.DataBytes.resize(ids.size());
        }
    }
    catch (const std::bad_alloc&)
    {
        payload.ConstructResult = Wirehair_OOM;
        payload.Result = Wirehair_OOM;
        return payload.Result;
    }
    catch (const std::length_error&)
    {
        payload.ConstructResult = Wirehair_OOM;
        payload.Result = Wirehair_OOM;
        return payload.Result;
    }

    wirehair_v2::MessagePrecodeEncoder encoder;
    switch (kind)
    {
    case RepairTimingEncoderKind::RawAttemptZero:
        payload.ConstructResult =
            encoder.InitializeRepairV1SelectedResultForTesting(
                message.Data(), message.Size(), options.BlockBytes,
                *options.Contract, construction_root, 0u,
                options.SystematicCache);
        break;
    case RepairTimingEncoderKind::LazySelector:
        payload.ConstructResult =
            encoder.InitializeRepairV1ResultForTesting(
                message.Data(), message.Size(), options.BlockBytes,
                *options.Contract, construction_root, telemetry,
                options.SystematicCache);
        break;
    case RepairTimingEncoderKind::ForcedSelected:
        payload.ConstructResult =
            encoder.InitializeRepairV1SelectedResultForTesting(
                message.Data(), message.Size(), options.BlockBytes,
                *options.Contract, construction_root, selected_attempt,
                options.SystematicCache);
        break;
    }
    payload.Result = payload.ConstructResult;
    if (payload.Result != Wirehair_Success) {
        return payload.Result;
    }
    if (encoder.ProvisionalRepairContractIdForTesting() !=
            options.Contract->ProvisionalId ||
        !encoder.IntermediateBlocks())
    {
        payload.Result = Wirehair_Error;
        return payload.Result;
    }
    if (solve_stats) {
        *solve_stats = encoder.SolveStats();
    }
    const wirehair_v2::PrecodeParams& params =
        encoder.BlockEncoder().System().Params;
    const uint64_t intermediate_blocks =
        (uint64_t)params.BlockCount + params.Staircase +
        params.DenseRows + params.HeavyRows;
    const uint64_t intermediate_bytes =
        intermediate_blocks * options.BlockBytes;
    if (intermediate_bytes >
        (uint64_t)std::numeric_limits<size_t>::max())
    {
        payload.Result = Wirehair_OOM;
        return payload.Result;
    }
    payload.IntermediateBytes = (size_t)intermediate_bytes;
    payload.IntermediateSha256 = Sha256HexBytes(
        encoder.IntermediateBlocks(), payload.IntermediateBytes);
    for (size_t i = 0u; i < ids.size(); ++i)
    {
        uint32_t bytes = 0u;
        payload.Result = encoder.EncodeResult(
            ids[i],
            payload.Bytes.Data() + i * (size_t)options.BlockBytes,
            options.BlockBytes, &bytes);
        const uint32_t expected =
            ids[i] + 1u == options.BlockCount ?
                1u : options.BlockBytes;
        if (payload.Result != Wirehair_Success || bytes != expected)
        {
            if (payload.Result == Wirehair_Success) {
                payload.Result = Wirehair_Error;
            }
            return payload.Result;
        }
        payload.DataBytes[i] = bytes;
    }
    payload.Sha256 =
        Sha256HexBytes(payload.Bytes.Data(), payload.Bytes.Size());
    return payload.Result;
}

BandTimingRecovery RunRepairTimingRecovery(
    const RepairTimingOptions& options,
    const RepairTimingDescriptor& descriptor,
    const BandTimingAlignedBytes& message,
    const std::vector<uint32_t>& ids,
    const BandTimingPayload& payload)
{
    BandTimingRecovery recovery;
    if (!descriptor.Contract ||
        descriptor.Contract != options.Contract ||
        payload.Result != Wirehair_Success ||
        payload.DataBytes.size() != ids.size() ||
        payload.Bytes.Size() !=
            ids.size() * (size_t)options.BlockBytes)
    {
        recovery.ConstructResult = payload.ConstructResult;
        recovery.Result = payload.Result;
        return recovery;
    }
    RepairTimingDescriptor parsed;
    if (!ParseRepairTimingDescriptor(descriptor.Bytes, parsed) ||
        parsed.Contract != options.Contract ||
        parsed.ConstructionRoot != descriptor.ConstructionRoot ||
        parsed.SelectedAttempt != descriptor.SelectedAttempt)
    {
        recovery.ConstructResult = Wirehair_InvalidInput;
        recovery.Result = Wirehair_InvalidInput;
        return recovery;
    }
    wirehair_v2::MessagePrecodeDecoder decoder;
    recovery.ConstructResult =
        decoder.InitializeRepairV1ResultForTesting(
            message.Size(), options.BlockBytes, *parsed.Contract,
            parsed.ConstructionRoot, parsed.SelectedAttempt,
            options.SystematicCache);
    recovery.Result = recovery.ConstructResult;
    if (recovery.Result != Wirehair_Success) {
        return recovery;
    }
    recovery.Result = Wirehair_NeedMore;
    for (size_t i = 0u; i < ids.size(); ++i)
    {
        recovery.Result = decoder.DecodeResult(
            ids[i],
            payload.Bytes.Data() + i * (size_t)options.BlockBytes,
            payload.DataBytes[i]);
        ++recovery.ReceivedSymbols;
        if (recovery.Result != Wirehair_NeedMore) break;
    }
    if (recovery.Result == Wirehair_Success)
    {
        recovery.Stats = decoder.SolveStats();
        recovery.StatsAvailable = true;
        std::vector<uint8_t> recovered(message.Size(), uint8_t{0});
        recovery.RecoveryResult = decoder.RecoverResult(
            recovered.data(), recovered.size());
        if (recovery.RecoveryResult == Wirehair_Success &&
            std::memcmp(
                recovered.data(), message.Data(), message.Size()) == 0)
        {
            recovery.RecoveryOk = true;
            recovery.Overhead =
                (int32_t)recovery.ReceivedSymbols -
                (int32_t)options.BlockCount;
            recovery.RecoveredSha256 =
                Sha256HexBytes(recovered.data(), recovered.size());
        }
        else if (recovery.RecoveryResult == Wirehair_Success) {
            recovery.RecoveryResult = Wirehair_Error;
        }
    }
    return recovery;
}

bool BuildRepairTimingExplicitDescriptor(
    const RepairTimingOptions& options,
    const RepairTimingDescriptor& wire_descriptor,
    BandTimingDescriptor& descriptor)
{
    descriptor = BandTimingDescriptor();
    if (!wire_descriptor.Contract ||
        wire_descriptor.Contract != options.Contract)
    {
        return false;
    }
    const wirehair_v2::RepairV1Contract& contract =
        *wire_descriptor.Contract;
    const wirehair_v2::ScopedRepairV1ContractStateForTesting scope(
        contract);
    if (!scope.IsValid() ||
        !wirehair_v2::MakeRepairV1ExplicitConfigForTesting(
            contract, options.BlockCount,
            wire_descriptor.ConstructionRoot,
            wire_descriptor.SelectedAttempt,
            options.SystematicCache,
            options.SystematicCache,
            descriptor.Explicit))
    {
        return false;
    }
    descriptor.GF256Rows = contract.GF256Rows;
    descriptor.GF16Rows = contract.GF16Rows;
    descriptor.Period = wirehair_v2::kMixedCoefficientPeriod;
    descriptor.X = BandTimingX::Frozen;
    return true;
}

bool SameRepairTimingStatsStructural(
    const wirehair_v2::PrecodeSolveStats& a,
    const wirehair_v2::PrecodeSolveStats& b)
{
    return
        a.PacketRows == b.PacketRows &&
        a.PeeledColumns == b.PeeledColumns &&
        a.InactivatedColumns == b.InactivatedColumns &&
        a.ResidualRows == b.ResidualRows &&
        a.ResidualRank == b.ResidualRank &&
        a.BinaryResidualRank == b.BinaryResidualRank &&
        a.BinaryRowReferences == b.BinaryRowReferences &&
        a.BinaryRowStorageBytes == b.BinaryRowStorageBytes &&
        a.BinaryAdjacencyStorageBytes ==
            b.BinaryAdjacencyStorageBytes &&
        a.BinaryRowStorageAllocations ==
            b.BinaryRowStorageAllocations &&
        a.BinaryAdjacencyStorageAllocations ==
            b.BinaryAdjacencyStorageAllocations &&
        a.BlockXors == b.BlockXors &&
        a.BlockMulAdds == b.BlockMulAdds &&
        a.BlockCopies == b.BlockCopies &&
        a.BlockZeroFills == b.BlockZeroFills &&
        a.BlockAddSets == b.BlockAddSets &&
        a.BlockAddSetSources == b.BlockAddSetSources &&
        a.PeelAdjacencyVisits == b.PeelAdjacencyVisits &&
        a.PeelRowScanSteps == b.PeelRowScanSteps &&
        a.PeelHeapOperations == b.PeelHeapOperations &&
        a.PeelHeapCompactionRebuildColumnProbes ==
            b.PeelHeapCompactionRebuildColumnProbes &&
        a.PeelHeapCompactionHeapifyInputKeys ==
            b.PeelHeapCompactionHeapifyInputKeys &&
        a.ProjectionWordXors == b.ProjectionWordXors &&
        a.ResidualCoeffWordXors == b.ResidualCoeffWordXors &&
        a.ResidualCoeffByteOps == b.ResidualCoeffByteOps &&
        a.PacketSeedAttempt == b.PacketSeedAttempt &&
        a.TinyMixedFastPathAcceptances ==
            b.TinyMixedFastPathAcceptances &&
        a.MixedJointSourceXors == b.MixedJointSourceXors &&
        a.MixedJointMarginalXors == b.MixedJointMarginalXors &&
        a.MixedJointMarginalCopies == b.MixedJointMarginalCopies &&
        a.MixedJointScratchBytes == b.MixedJointScratchBytes &&
        a.MixedJointActiveDeltas == b.MixedJointActiveDeltas &&
        a.MixedDualSourceColumns == b.MixedDualSourceColumns;
}

template<typename T>
bool RepairTimingAddStat(T& total, T value)
{
    if (value > std::numeric_limits<T>::max() - total) {
        return false;
    }
    total += value;
    return true;
}

bool RepairTimingAddStats(
    wirehair_v2::PrecodeSolveStats& total,
    const wirehair_v2::PrecodeSolveStats& value)
{
#define WIREHAIR_REPAIR_TIMING_ADD_STAT(field) \
    if (!RepairTimingAddStat(total.field, value.field)) return false
    WIREHAIR_REPAIR_TIMING_ADD_STAT(PacketRows);
    WIREHAIR_REPAIR_TIMING_ADD_STAT(PeeledColumns);
    WIREHAIR_REPAIR_TIMING_ADD_STAT(InactivatedColumns);
    WIREHAIR_REPAIR_TIMING_ADD_STAT(ResidualRows);
    WIREHAIR_REPAIR_TIMING_ADD_STAT(ResidualRank);
    WIREHAIR_REPAIR_TIMING_ADD_STAT(BinaryResidualRank);
    WIREHAIR_REPAIR_TIMING_ADD_STAT(BinaryRowReferences);
    WIREHAIR_REPAIR_TIMING_ADD_STAT(BinaryRowStorageBytes);
    WIREHAIR_REPAIR_TIMING_ADD_STAT(BinaryAdjacencyStorageBytes);
    WIREHAIR_REPAIR_TIMING_ADD_STAT(BinaryRowStorageAllocations);
    WIREHAIR_REPAIR_TIMING_ADD_STAT(
        BinaryAdjacencyStorageAllocations);
    WIREHAIR_REPAIR_TIMING_ADD_STAT(BlockXors);
    WIREHAIR_REPAIR_TIMING_ADD_STAT(BlockMulAdds);
    WIREHAIR_REPAIR_TIMING_ADD_STAT(BlockCopies);
    WIREHAIR_REPAIR_TIMING_ADD_STAT(BlockZeroFills);
    WIREHAIR_REPAIR_TIMING_ADD_STAT(BlockAddSets);
    WIREHAIR_REPAIR_TIMING_ADD_STAT(BlockAddSetSources);
    WIREHAIR_REPAIR_TIMING_ADD_STAT(PeelAdjacencyVisits);
    WIREHAIR_REPAIR_TIMING_ADD_STAT(PeelRowScanSteps);
    WIREHAIR_REPAIR_TIMING_ADD_STAT(PeelHeapOperations);
    // Preserve the common structural receipt for in-memory A/A equality.
    // The frozen repairtiming v3 row and digest schemas intentionally do not
    // acquire columns for these newly production-visible counters.
    WIREHAIR_REPAIR_TIMING_ADD_STAT(
        PeelHeapCompactionRebuildColumnProbes);
    WIREHAIR_REPAIR_TIMING_ADD_STAT(
        PeelHeapCompactionHeapifyInputKeys);
    WIREHAIR_REPAIR_TIMING_ADD_STAT(ProjectionWordXors);
    WIREHAIR_REPAIR_TIMING_ADD_STAT(ResidualCoeffWordXors);
    WIREHAIR_REPAIR_TIMING_ADD_STAT(ResidualCoeffByteOps);
    WIREHAIR_REPAIR_TIMING_ADD_STAT(BuildNanoseconds);
    WIREHAIR_REPAIR_TIMING_ADD_STAT(PeelNanoseconds);
    WIREHAIR_REPAIR_TIMING_ADD_STAT(ProjectNanoseconds);
    WIREHAIR_REPAIR_TIMING_ADD_STAT(ResidualNanoseconds);
    WIREHAIR_REPAIR_TIMING_ADD_STAT(BackSubNanoseconds);
    WIREHAIR_REPAIR_TIMING_ADD_STAT(PacketSeedAttempt);
    WIREHAIR_REPAIR_TIMING_ADD_STAT(
        TinyMixedFastPathAcceptances);
    WIREHAIR_REPAIR_TIMING_ADD_STAT(MixedJointSourceXors);
    WIREHAIR_REPAIR_TIMING_ADD_STAT(MixedJointMarginalXors);
    WIREHAIR_REPAIR_TIMING_ADD_STAT(MixedJointMarginalCopies);
    WIREHAIR_REPAIR_TIMING_ADD_STAT(MixedJointScratchBytes);
    WIREHAIR_REPAIR_TIMING_ADD_STAT(MixedJointActiveDeltas);
    WIREHAIR_REPAIR_TIMING_ADD_STAT(MixedDualSourceColumns);
#undef WIREHAIR_REPAIR_TIMING_ADD_STAT
    return true;
}

bool RepairTimingAggregateSelectorStats(
    const wirehair_v2::RepairV1SelectorTelemetry& telemetry,
    wirehair_v2::PrecodeSolveStats& total)
{
    total = wirehair_v2::PrecodeSolveStats();
    uint32_t calls = 0u;
    uint32_t real_calls = 0u;
    uint32_t probe_calls = 0u;
    for (uint32_t attempt = 0u;
         attempt < wirehair_v2::kRepairV1AttemptCount;
         ++attempt)
    {
        const wirehair_v2::RepairV1AttemptTelemetry& observation =
            telemetry.Attempts[attempt];
        if (observation.ProbeExecuted)
        {
            if (!observation.ProbeStatsAvailable ||
                !RepairTimingAddStats(
                    total, observation.ProbeStats))
            {
                return false;
            }
            ++calls;
            ++probe_calls;
        }
        else if (observation.ProbeStatsAvailable) {
            return false;
        }
        if (observation.RealExecuted)
        {
            if (!observation.RealStatsAvailable ||
                !RepairTimingAddStats(
                    total, observation.RealStats))
            {
                return false;
            }
            ++calls;
            ++real_calls;
        }
        else if (observation.RealStatsAvailable) {
            return false;
        }
    }
    return calls == telemetry.CallsExecuted &&
        real_calls == telemetry.RealConfigurationCalls &&
        probe_calls == telemetry.StructuralProbeCalls;
}

void AppendRepairTimingStructuralStats(
    std::ostringstream& text,
    const wirehair_v2::PrecodeSolveStats* stats,
    bool available)
{
    // Frozen selector-structural-v1 subset.  Heap-compaction C/M remain
    // internal until a future repairtiming protocol and digest-domain bump.
    const auto value = [&](uint64_t v) -> std::string {
        return available ?
            std::to_string(v) : std::string("-1");
    };
    text
        << ",probe_packet_rows="
            << value(stats ? stats->PacketRows : 0u)
        << ",probe_peeled_columns="
            << value(stats ? stats->PeeledColumns : 0u)
        << ",probe_inactivated_columns="
            << value(stats ? stats->InactivatedColumns : 0u)
        << ",probe_residual_rows="
            << value(stats ? stats->ResidualRows : 0u)
        << ",probe_residual_rank="
            << value(stats ? stats->ResidualRank : 0u)
        << ",probe_binary_residual_rank="
            << value(stats ? stats->BinaryResidualRank : 0u)
        << ",probe_binary_row_references="
            << value(stats ? stats->BinaryRowReferences : 0u)
        << ",probe_binary_row_storage_bytes="
            << value(stats ? stats->BinaryRowStorageBytes : 0u)
        << ",probe_binary_adjacency_storage_bytes="
            << value(stats ? stats->BinaryAdjacencyStorageBytes : 0u)
        << ",probe_binary_row_storage_allocations="
            << value(stats ? stats->BinaryRowStorageAllocations : 0u)
        << ",probe_binary_adjacency_storage_allocations="
            << value(
                stats ?
                    stats->BinaryAdjacencyStorageAllocations : 0u)
        << ",probe_block_xors="
            << value(stats ? stats->BlockXors : 0u)
        << ",probe_block_muladds="
            << value(stats ? stats->BlockMulAdds : 0u)
        << ",probe_block_copies="
            << value(stats ? stats->BlockCopies : 0u)
        << ",probe_block_zero_fills="
            << value(stats ? stats->BlockZeroFills : 0u)
        << ",probe_block_add_sets="
            << value(stats ? stats->BlockAddSets : 0u)
        << ",probe_block_add_set_sources="
            << value(stats ? stats->BlockAddSetSources : 0u)
        << ",probe_peel_adjacency_visits="
            << value(stats ? stats->PeelAdjacencyVisits : 0u)
        << ",probe_peel_row_scan_steps="
            << value(stats ? stats->PeelRowScanSteps : 0u)
        << ",probe_peel_heap_operations="
            << value(stats ? stats->PeelHeapOperations : 0u)
        << ",probe_projection_word_xors="
            << value(stats ? stats->ProjectionWordXors : 0u)
        << ",probe_residual_coeff_word_xors="
            << value(stats ? stats->ResidualCoeffWordXors : 0u)
        << ",probe_residual_coeff_byte_ops="
            << value(stats ? stats->ResidualCoeffByteOps : 0u)
        << ",probe_packet_seed_attempt="
            << value(stats ? stats->PacketSeedAttempt : 0u)
        << ",probe_tiny_mixed_fastpath_acceptances="
            << value(
                stats ?
                    stats->TinyMixedFastPathAcceptances : 0u)
        << ",probe_mixed_joint_source_xors="
            << value(stats ? stats->MixedJointSourceXors : 0u)
        << ",probe_mixed_joint_marginal_xors="
            << value(stats ? stats->MixedJointMarginalXors : 0u)
        << ",probe_mixed_joint_marginal_copies="
            << value(stats ? stats->MixedJointMarginalCopies : 0u)
        << ",probe_mixed_joint_scratch_bytes="
            << value(stats ? stats->MixedJointScratchBytes : 0u)
        << ",probe_mixed_joint_active_deltas="
            << value(stats ? stats->MixedJointActiveDeltas : 0u)
        << ",probe_mixed_dual_source_columns="
            << value(stats ? stats->MixedDualSourceColumns : 0u);
}

std::string RepairTimingStructuralDigest(
    const RepairTimingOptions& options,
    uint32_t construction_root,
    const wirehair_v2::RepairV1SelectorTelemetry& telemetry,
    const RepairTimingDescriptor* descriptor)
{
    std::ostringstream text;
    text
        << "wirehair-wh2-repair-v1-selector-structural-v1\n"
        << "arm_id=" << std::hex << std::setfill('0')
        << std::setw(16) << options.Contract->ProvisionalId
        << std::dec << "\n"
        << "K=" << options.BlockCount << "\n"
        << "construction_root=" << construction_root << "\n"
        << "repair_policy_sha256="
        << wirehair_v2::kRepairV1PolicySha256 << "\n"
        << "selector_result=" << (int)telemetry.Result << "\n"
        << "selected_attempt="
        << (telemetry.SelectedAttempt ==
                wirehair_v2::kRepairV1NoAttempt ?
                -1 : (int64_t)telemetry.SelectedAttempt)
        << "\n"
        << "attempts_executed=" << telemetry.AttemptsExecuted << "\n"
        << "calls_executed=" << telemetry.CallsExecuted << "\n"
        << "real_configuration_calls="
        << telemetry.RealConfigurationCalls << "\n"
        << "structural_probe_calls="
        << telemetry.StructuralProbeCalls << "\n"
        << "cap_exhausted=" << (telemetry.CapExhausted ? 1 : 0)
        << "\n"
        << "fatal_attempt_zero_mismatch="
        << (telemetry.FatalAttemptZeroMismatch ? 1 : 0) << "\n"
        << "oom=" << (telemetry.Oom ? 1 : 0) << "\n"
        << "committed=" << (telemetry.Committed ? 1 : 0) << "\n"
        << "descriptor_sha256="
        << (descriptor ? descriptor->Sha256 : "not_applicable")
        << "\n";
    for (uint32_t attempt = 0u;
         attempt < wirehair_v2::kRepairV1AttemptCount;
         ++attempt)
    {
        const wirehair_v2::RepairV1AttemptTelemetry& observation =
            telemetry.Attempts[attempt];
        text
            << "attempt=" << attempt
            << ",probe_executed="
            << (observation.ProbeExecuted ? 1 : 0)
            << ",probe_result="
            << (observation.ProbeExecuted ?
                (int)observation.ProbeResult : -1)
            << ",probe_stats_available="
            << (observation.ProbeStatsAvailable ? 1 : 0);
        AppendRepairTimingStructuralStats(
            text,
            observation.ProbeStatsAvailable ?
                &observation.ProbeStats : nullptr,
            observation.ProbeStatsAvailable);
        text << '\n';
    }
    return Sha256Hex(text.str());
}

enum class RepairTimingRole
{
    RepairSelectorEncoder,
    RepairForcedEncoder,
    Wh1Encoder,
    RepairDecoderFeed,
    Wh1DecoderFeed,
    RepairDecoderFull,
    Wh1DecoderFull,
    RepairDirect,
    DispatchDirect
};

static_assert(
    (size_t)RepairTimingRole::DispatchDirect + 1u == 9u,
    "repairtiming preflight role array requires contiguous roles");

const char* RepairTimingRoleName(RepairTimingRole role)
{
    switch (role)
    {
    case RepairTimingRole::RepairSelectorEncoder:
        return "repair_selector_encoder";
    case RepairTimingRole::RepairForcedEncoder:
        return "repair_forced_encoder";
    case RepairTimingRole::Wh1Encoder:
        return "wh1_encoder";
    case RepairTimingRole::RepairDecoderFeed:
        return "repair_decoder_feed";
    case RepairTimingRole::Wh1DecoderFeed:
        return "wh1_decoder_feed";
    case RepairTimingRole::RepairDecoderFull:
        return "repair_decoder_full";
    case RepairTimingRole::Wh1DecoderFull:
        return "wh1_decoder_full";
    case RepairTimingRole::RepairDirect:
        return "repair_direct";
    case RepairTimingRole::DispatchDirect:
        return "dispatch_direct";
    }
    return "unknown";
}

const char* RepairTimingScopeName(RepairTimingRole role)
{
    switch (role)
    {
    case RepairTimingRole::RepairSelectorEncoder:
        return "encoder_selector";
    case RepairTimingRole::RepairForcedEncoder:
        return "encoder_forced";
    case RepairTimingRole::Wh1Encoder:
        return "encoder_wh1";
    case RepairTimingRole::RepairDecoderFeed:
    case RepairTimingRole::Wh1DecoderFeed:
        return "decoder_feed";
    case RepairTimingRole::RepairDecoderFull:
    case RepairTimingRole::Wh1DecoderFull:
        return "decoder_full";
    case RepairTimingRole::RepairDirect:
    case RepairTimingRole::DispatchDirect:
        return "direct";
    }
    return "unknown";
}

struct RepairTimingPanel
{
    const char* Name;
    RepairTimingRole A;
    RepairTimingRole B;
};

static const RepairTimingPanel kRepairTimingPanels[] = {
    {
        "encoder_selector_forced",
        RepairTimingRole::RepairSelectorEncoder,
        RepairTimingRole::RepairForcedEncoder
    },
    {
        "encoder_selector_aa",
        RepairTimingRole::RepairSelectorEncoder,
        RepairTimingRole::RepairSelectorEncoder
    },
    {
        "encoder_forced_aa",
        RepairTimingRole::RepairForcedEncoder,
        RepairTimingRole::RepairForcedEncoder
    },
    {
        "encoder_selector_wh1",
        RepairTimingRole::RepairSelectorEncoder,
        RepairTimingRole::Wh1Encoder
    },
    {
        "encoder_wh1_aa",
        RepairTimingRole::Wh1Encoder,
        RepairTimingRole::Wh1Encoder
    },
    {
        "decoder_feed_wh1",
        RepairTimingRole::RepairDecoderFeed,
        RepairTimingRole::Wh1DecoderFeed
    },
    {
        "decoder_feed_aa",
        RepairTimingRole::RepairDecoderFeed,
        RepairTimingRole::RepairDecoderFeed
    },
    {
        "decoder_feed_wh1_aa",
        RepairTimingRole::Wh1DecoderFeed,
        RepairTimingRole::Wh1DecoderFeed
    },
    {
        "decoder_full_wh1",
        RepairTimingRole::RepairDecoderFull,
        RepairTimingRole::Wh1DecoderFull
    },
    {
        "decoder_full_aa",
        RepairTimingRole::RepairDecoderFull,
        RepairTimingRole::RepairDecoderFull
    },
    {
        "decoder_full_wh1_aa",
        RepairTimingRole::Wh1DecoderFull,
        RepairTimingRole::Wh1DecoderFull
    },
    {
        "direct_dispatch",
        RepairTimingRole::RepairDirect,
        RepairTimingRole::DispatchDirect
    },
    {
        "direct_aa",
        RepairTimingRole::RepairDirect,
        RepairTimingRole::RepairDirect
    },
    {
        "direct_dispatch_aa",
        RepairTimingRole::DispatchDirect,
        RepairTimingRole::DispatchDirect
    }
};

static_assert(
    sizeof(kRepairTimingPanels) / sizeof(kRepairTimingPanels[0]) == 14u,
    "repairtiming v3 freezes fourteen timing panels");

struct RepairTimingCell
{
    uint32_t Index = 0u;
    uint32_t ConstructionRoot = 0u;
    uint64_t LossSeed = 0u;
    std::string TraceSha256;
    std::string MessageSha256;
    BandTimingAlignedBytes Message;
    std::vector<uint32_t> Ids;
    std::vector<uint32_t> FirstKIds;

    wirehair_v2::RepairV1SelectorTelemetry Selector;
    RepairTimingDescriptor Descriptor;
    bool HaveDescriptor = false;
    std::string StructuralSha256;

    BandTimingPayload RawPayload;
    BandTimingPayload RepairedPayload;
    BandTimingPayload ForcedAPayload;
    BandTimingPayload ForcedBPayload;
    BandTimingPayload DispatchPayload;
    BandTimingPayload Wh1Payload;

    BandTimingRecovery RawRecovery;
    BandTimingRecovery RepairedRecovery;
    BandTimingRecovery DispatchRecovery;
    BandTimingRecovery Wh1Recovery;

    BandTimingDescriptor RepairExplicit;
    BandTimingDescriptor DispatchExplicit;
    bool HaveRepairExplicit = false;
    bool HaveDispatchExplicit = false;
    uint32_t DirectPrefixSymbols = 0u;
    bool RepairDirectExecuted = false;
    bool DispatchDirectExecuted = false;
    BandTimingDirect RepairDirect;
    BandTimingDirect DispatchDirect;
};

struct RepairTimingOperation
{
    bool Available = false;
    WirehairResult ConstructResult = Wirehair_InvalidInput;
    WirehairResult Result = Wirehair_InvalidInput;
    WirehairResult RecoveryResult = Wirehair_InvalidInput;
    bool RecoveryOk = false;
    uint32_t ReceivedSymbols = 0u;
    wirehair_v2::PrecodeSolveStats Stats;
    bool StatsAvailable = false;
    size_t IntermediateBytes = 0u;
    std::string IntermediateSha256;
    std::string SelectorStructuralSha256;
    std::string SignatureSha256;
    uint64_t ElapsedNanoseconds = 0u;
    int CpuBefore = -1;
    int CpuAfter = -1;
    int64_t MinorFaults = -1;
    int64_t MajorFaults = -1;
};

thread_local bool RepairTimingIntervalActive = false;

bool BeginRepairTimingInterval(
    const RepairTimingOptions& options,
    std::vector<uint8_t>& eviction,
    PeelTimingUsage& usage_before,
    int& cpu_before,
    uint64_t& begin_ns)
{
    if (RepairTimingIntervalActive) return false;
    if (options.CacheState == PeelTimingCacheState::Cold) {
        EvictPeelTimingCache(eviction);
    }
    usage_before = ReadPeelTimingUsage();
    cpu_before = PeelTimingCurrentCpu();
    if (cpu_before < 0 ||
        !ReadPeelTimingMonotonicNanoseconds(begin_ns))
    {
        return false;
    }
    RepairTimingIntervalActive = true;
    return true;
}

bool FinishRepairTimingInterval(
    const PeelTimingUsage& usage_before,
    int cpu_before,
    uint64_t begin_ns,
    RepairTimingOperation& operation)
{
    uint64_t end_ns = 0u;
    const bool clock_ok =
        ReadPeelTimingMonotonicNanoseconds(end_ns);
    RepairTimingIntervalActive = false;
    if (!clock_ok ||
        end_ns <= begin_ns)
    {
        return false;
    }
    operation.CpuBefore = cpu_before;
    operation.CpuAfter = PeelTimingCurrentCpu();
    const PeelTimingUsage usage_after = ReadPeelTimingUsage();
    if (operation.CpuAfter < 0 ||
        !BandTimingUsageDeltas(
            usage_before, usage_after,
            operation.MinorFaults, operation.MajorFaults))
    {
        return false;
    }
    operation.ElapsedNanoseconds = end_ns - begin_ns;
    return true;
}

bool RepairTimingWitnessOutsideInterval(
    const RepairTimingOptions& options)
{
    return !RepairTimingIntervalActive &&
        !options.TimedScopeFault;
}

std::string RepairTimingOperationSignature(
    const RepairTimingOperation& operation,
    const std::string& payload_sha256,
    const std::string& recovered_sha256,
    const std::string& descriptor_sha256)
{
    std::ostringstream text;
    text
        << "wirehair-wh2-repairtiming-operation-v1\n"
        << (int)operation.ConstructResult << ','
        << (int)operation.Result << ','
        << (int)operation.RecoveryResult << ','
        << (operation.RecoveryOk ? 1 : 0) << ','
        << operation.ReceivedSymbols << ','
        << operation.IntermediateBytes << ','
        << (operation.IntermediateSha256.empty() ?
            "not_applicable" : operation.IntermediateSha256) << ','
        << (payload_sha256.empty() ?
            "not_applicable" : payload_sha256) << ','
        << (recovered_sha256.empty() ?
            "not_applicable" : recovered_sha256) << ','
        << (descriptor_sha256.empty() ?
            "not_applicable" : descriptor_sha256) << ','
        << (operation.SelectorStructuralSha256.empty() ?
            "not_applicable" :
            operation.SelectorStructuralSha256) << '\n';
    return Sha256Hex(text.str());
}

bool RunRepairTimingEncoderOperation(
    const RepairTimingOptions& options,
    const RepairTimingCell& cell,
    RepairTimingRole role,
    std::vector<uint8_t>& eviction,
    bool measure,
    RepairTimingOperation& operation)
{
    operation = RepairTimingOperation();
    if (role != RepairTimingRole::RepairSelectorEncoder &&
        role != RepairTimingRole::RepairForcedEncoder &&
        role != RepairTimingRole::Wh1Encoder)
    {
        return false;
    }
    if (role == RepairTimingRole::RepairForcedEncoder &&
        !cell.HaveDescriptor)
    {
        return true;
    }

    BandTimingPayload payload;
    if (!payload.Bytes.Resize(
            cell.FirstKIds.size() *
                (size_t)options.BlockBytes))
    {
        return false;
    }
    try {
        payload.DataBytes.resize(cell.FirstKIds.size());
    }
    catch (const std::bad_alloc&) {
        return false;
    }
    catch (const std::length_error&) {
        return false;
    }
    wirehair_v2::RepairV1SelectorTelemetry telemetry;
    RepairTimingDescriptor descriptor;
    descriptor.Bytes.reserve(17u);
    wirehair_v2::MessagePrecodeEncoder repair_encoder;
    wirehair::Codec wh1_encoder;
    const uint32_t wh1_dense_count =
        role == RepairTimingRole::Wh1Encoder ?
            wirehair::GetDenseCount(options.BlockCount) : 0u;
    if (wh1_dense_count > UINT16_MAX) {
        return false;
    }
    PeelTimingUsage usage_before;
    int cpu_before = -1;
    uint64_t begin_ns = 0u;
    if (measure &&
        !BeginRepairTimingInterval(
            options, eviction, usage_before, cpu_before, begin_ns))
    {
        return false;
    }

    if (role == RepairTimingRole::Wh1Encoder)
    {
        wh1_encoder.SetWireProfile(
            wirehair::WireProfile::LegacyCurrent);
        wh1_encoder.OverrideSeeds(
            (uint16_t)wh1_dense_count,
            (uint16_t)(cell.ConstructionRoot & 0xffffu),
            (uint16_t)(cell.ConstructionRoot >> 16));
        operation.ConstructResult =
            wh1_encoder.InitializeEncoder(
                cell.Message.Size(), options.BlockBytes);
        operation.Result = operation.ConstructResult;
        if (operation.Result == Wirehair_Success)
        {
            operation.Result = wh1_encoder.EncodeFeed(
                cell.Message.Data(), options.SystematicCache);
        }
        for (size_t i = 0u;
             operation.Result == Wirehair_Success &&
             i < cell.FirstKIds.size(); ++i)
        {
            const uint32_t id = cell.FirstKIds[i];
            payload.DataBytes[i] = wh1_encoder.Encode(
                id,
                payload.Bytes.Data() +
                    i * (size_t)options.BlockBytes,
                options.BlockBytes);
            const uint32_t expected =
                id + 1u == options.BlockCount ?
                    1u : options.BlockBytes;
            if (payload.DataBytes[i] != expected) {
                operation.Result = Wirehair_Error;
            }
        }
    }
    else
    {
        if (role == RepairTimingRole::RepairSelectorEncoder)
        {
            operation.ConstructResult =
                repair_encoder.InitializeRepairV1ResultForTesting(
                    cell.Message.Data(), cell.Message.Size(),
                    options.BlockBytes, *options.Contract,
                    cell.ConstructionRoot, &telemetry,
                    options.SystematicCache);
        }
        else
        {
            operation.ConstructResult =
                repair_encoder.
                    InitializeRepairV1SelectedResultForTesting(
                        cell.Message.Data(), cell.Message.Size(),
                        options.BlockBytes, *options.Contract,
                        cell.ConstructionRoot,
                        cell.Descriptor.SelectedAttempt,
                        options.SystematicCache);
        }
        operation.Result = operation.ConstructResult;
        for (size_t i = 0u;
             operation.Result == Wirehair_Success &&
             i < cell.FirstKIds.size(); ++i)
        {
            const uint32_t id = cell.FirstKIds[i];
            uint32_t bytes = 0u;
            operation.Result = repair_encoder.EncodeResult(
                id,
                payload.Bytes.Data() +
                    i * (size_t)options.BlockBytes,
                options.BlockBytes, &bytes);
            const uint32_t expected =
                id + 1u == options.BlockCount ?
                    1u : options.BlockBytes;
            if (operation.Result == Wirehair_Success &&
                bytes != expected)
            {
                operation.Result = Wirehair_Error;
            }
            payload.DataBytes[i] = bytes;
        }
        if (operation.Result == Wirehair_Success)
        {
            const uint32_t attempt =
                role == RepairTimingRole::RepairSelectorEncoder ?
                    telemetry.SelectedAttempt :
                    cell.Descriptor.SelectedAttempt;
            descriptor.Contract = options.Contract;
            descriptor.ConstructionRoot = cell.ConstructionRoot;
            descriptor.SelectedAttempt = attempt;
            descriptor.Bytes.push_back((uint8_t)'W');
            descriptor.Bytes.push_back((uint8_t)'2');
            descriptor.Bytes.push_back((uint8_t)'R');
            descriptor.Bytes.push_back((uint8_t)'3');
            RepairTimingAppendLe64(
                descriptor.Bytes,
                options.Contract->ProvisionalId);
            RepairTimingAppendLe32(
                descriptor.Bytes, cell.ConstructionRoot);
            descriptor.Bytes.push_back((uint8_t)attempt);
            if (descriptor.Bytes.size() != 17u)
            {
                operation.Result = Wirehair_Error;
            }
        }
    }

    if (measure &&
        !FinishRepairTimingInterval(
            usage_before, cpu_before, begin_ns, operation))
    {
        return false;
    }
    if (!RepairTimingWitnessOutsideInterval(options)) {
        return false;
    }
    payload.ConstructResult = operation.ConstructResult;
    payload.Result = operation.Result;
    if (role != RepairTimingRole::Wh1Encoder &&
        operation.Result == Wirehair_Success)
    {
        if (role == RepairTimingRole::RepairSelectorEncoder)
        {
            if (!RepairTimingAggregateSelectorStats(
                    telemetry, operation.Stats))
            {
                return false;
            }
        }
        else {
            operation.Stats = repair_encoder.SolveStats();
        }
        operation.StatsAvailable = true;
        const wirehair_v2::PrecodeParams& params =
            repair_encoder.BlockEncoder().System().Params;
        const uint64_t intermediate_bytes =
            (
                (uint64_t)params.BlockCount +
                params.Staircase + params.DenseRows +
                params.HeavyRows
            ) * options.BlockBytes;
        if (intermediate_bytes >
                (uint64_t)std::numeric_limits<size_t>::max() ||
            !repair_encoder.IntermediateBlocks())
        {
            return false;
        }
        payload.IntermediateBytes = (size_t)intermediate_bytes;
        payload.IntermediateSha256 = Sha256HexBytes(
            repair_encoder.IntermediateBlocks(),
            payload.IntermediateBytes);
        descriptor.Hex = RepairTimingHex(descriptor.Bytes);
        descriptor.Sha256 = Sha256HexBytes(
            descriptor.Bytes.data(), descriptor.Bytes.size());
        const bool cell_witness_matches =
            cell.HaveDescriptor &&
            payload.IntermediateBytes ==
                cell.RepairedPayload.IntermediateBytes &&
            payload.IntermediateSha256 ==
                cell.RepairedPayload.IntermediateSha256 &&
            descriptor.Sha256 == cell.Descriptor.Sha256;
        if (!measure && !cell_witness_matches)
        {
            return false;
        }
        if (role == RepairTimingRole::RepairSelectorEncoder)
        {
            operation.SelectorStructuralSha256 =
                RepairTimingStructuralDigest(
                    options, cell.ConstructionRoot,
                    telemetry, &descriptor);
        }
    }
    if (operation.Result == Wirehair_Success)
    {
        payload.Sha256 =
            Sha256HexBytes(payload.Bytes.Data(), payload.Bytes.Size());
    }
    operation.IntermediateBytes = payload.IntermediateBytes;
    operation.IntermediateSha256 = payload.IntermediateSha256;
    operation.SignatureSha256 = RepairTimingOperationSignature(
        operation, payload.Sha256, std::string(),
        descriptor.Sha256);
    operation.Available =
        operation.ConstructResult == Wirehair_Success &&
        operation.Result == Wirehair_Success &&
        IsPeelTimingLowerHexSha256(payload.Sha256) &&
        (
            role == RepairTimingRole::Wh1Encoder ||
            (
                operation.StatsAvailable &&
                IsPeelTimingLowerHexSha256(
                    operation.IntermediateSha256) &&
                IsPeelTimingLowerHexSha256(descriptor.Sha256)
            )
        );
    return true;
}

bool RunRepairTimingDecoderOperation(
    const RepairTimingOptions& options,
    const RepairTimingCell& cell,
    RepairTimingRole role,
    std::vector<uint8_t>& eviction,
    bool measure,
    RepairTimingOperation& operation)
{
    operation = RepairTimingOperation();
    const bool repair =
        role == RepairTimingRole::RepairDecoderFeed ||
        role == RepairTimingRole::RepairDecoderFull;
    const bool full =
        role == RepairTimingRole::RepairDecoderFull ||
        role == RepairTimingRole::Wh1DecoderFull;
    if (!repair &&
        role != RepairTimingRole::Wh1DecoderFeed &&
        role != RepairTimingRole::Wh1DecoderFull)
    {
        return false;
    }
    if (repair &&
        (!cell.HaveDescriptor ||
         cell.RepairedPayload.Result != Wirehair_Success))
    {
        return true;
    }
    if (!repair && cell.Wh1Payload.Result != Wirehair_Success) {
        return true;
    }

    PeelTimingUsage usage_before;
    int cpu_before = -1;
    uint64_t begin_ns = 0u;
    const wirehair_v2::RepairV1Contract* parsed_contract = nullptr;
    uint32_t parsed_root = 0u;
    uint32_t parsed_attempt = wirehair_v2::kRepairV1NoAttempt;
    wirehair_v2::MessagePrecodeDecoder repair_decoder;
    wirehair::Codec wh1_decoder;
    const BandTimingPayload& payload =
        repair ? cell.RepairedPayload : cell.Wh1Payload;
    if (payload.DataBytes.size() != cell.Ids.size() ||
        payload.Bytes.Size() !=
            cell.Ids.size() * (size_t)options.BlockBytes)
    {
        return false;
    }

    std::vector<uint8_t> recovered;
    if (full)
    {
        try {
            recovered.assign(
                cell.Message.Size(), uint8_t{0});
        }
        catch (const std::bad_alloc&) {
            return false;
        }
        catch (const std::length_error&) {
            return false;
        }
    }

    if (full && measure &&
        !BeginRepairTimingInterval(
            options, eviction, usage_before, cpu_before, begin_ns))
    {
        return false;
    }

    if (repair)
    {
        if (!ParseRepairTimingDescriptorView(
                cell.Descriptor.Bytes, parsed_contract,
                parsed_root, parsed_attempt) ||
            parsed_contract != options.Contract ||
            parsed_root != cell.ConstructionRoot ||
            parsed_attempt != cell.Descriptor.SelectedAttempt)
        {
            operation.ConstructResult = Wirehair_InvalidInput;
            operation.Result = Wirehair_InvalidInput;
        }
        else
        {
            operation.ConstructResult =
                repair_decoder.InitializeRepairV1ResultForTesting(
                    cell.Message.Size(), options.BlockBytes,
                    *parsed_contract, parsed_root,
                    parsed_attempt, options.SystematicCache);
            operation.Result = operation.ConstructResult;
        }
    }
    else
    {
        wh1_decoder.SetWireProfile(
            wirehair::WireProfile::LegacyCurrent);
        const uint32_t dense_count =
            wirehair::GetDenseCount(options.BlockCount);
        if (dense_count > UINT16_MAX)
        {
            operation.ConstructResult = Wirehair_InvalidInput;
            operation.Result = Wirehair_InvalidInput;
        }
        else
        {
            wh1_decoder.OverrideSeeds(
                (uint16_t)dense_count,
                (uint16_t)(cell.ConstructionRoot & 0xffffu),
                (uint16_t)(cell.ConstructionRoot >> 16));
            operation.ConstructResult =
                wh1_decoder.InitializeDecoder(
                    cell.Message.Size(), options.BlockBytes);
            operation.Result = operation.ConstructResult;
        }
    }

    if (!full && measure &&
        !BeginRepairTimingInterval(
            options, eviction, usage_before, cpu_before, begin_ns))
    {
        return false;
    }

    uint32_t received = 0u;
    if (operation.Result == Wirehair_Success)
    {
        operation.Result = Wirehair_NeedMore;
        for (size_t i = 0u; i < cell.Ids.size(); ++i)
        {
            if (repair)
            {
                operation.Result = repair_decoder.DecodeResult(
                    cell.Ids[i],
                    payload.Bytes.Data() +
                        i * (size_t)options.BlockBytes,
                    payload.DataBytes[i]);
            }
            else
            {
                operation.Result = wh1_decoder.DecodeFeed(
                    cell.Ids[i],
                    payload.Bytes.Data() +
                        i * (size_t)options.BlockBytes,
                    payload.DataBytes[i]);
            }
            ++received;
            if (operation.Result != Wirehair_NeedMore) break;
        }
    }

    if (full && operation.Result == Wirehair_Success)
    {
        operation.RecoveryResult = repair ?
            repair_decoder.RecoverResult(
                recovered.data(), recovered.size()) :
            wh1_decoder.ReconstructOutput(
                recovered.data(), recovered.size());
    }

    if (measure &&
        !FinishRepairTimingInterval(
            usage_before, cpu_before, begin_ns, operation))
    {
        return false;
    }
    if (!RepairTimingWitnessOutsideInterval(options)) {
        return false;
    }
    if (repair && operation.Result == Wirehair_Success)
    {
        operation.Stats = repair_decoder.SolveStats();
        operation.StatsAvailable = true;
    }
    if (full && operation.RecoveryResult == Wirehair_Success)
    {
        if (std::memcmp(
                recovered.data(), cell.Message.Data(),
                cell.Message.Size()) == 0)
        {
            operation.RecoveryOk = true;
        }
        else {
            operation.RecoveryResult = Wirehair_Error;
        }
    }
    const std::string recovered_sha256 =
        operation.RecoveryOk ?
            Sha256HexBytes(recovered.data(), recovered.size()) :
            std::string();
    operation.ReceivedSymbols = received;
    operation.SignatureSha256 = RepairTimingOperationSignature(
        operation, payload.Sha256, recovered_sha256,
        repair ? cell.Descriptor.Sha256 : std::string());
    operation.Available =
        operation.ConstructResult == Wirehair_Success &&
        operation.Result == Wirehair_Success &&
        (!repair || operation.StatsAvailable) &&
        (!full ||
         (
             operation.RecoveryResult == Wirehair_Success &&
             operation.RecoveryOk
         ));
    return true;
}

bool RunRepairTimingDirectOperation(
    const RepairTimingOptions& options,
    const RepairTimingCell& cell,
    RepairTimingRole role,
    std::vector<uint8_t>& eviction,
    bool measure,
    RepairTimingOperation& operation)
{
    operation = RepairTimingOperation();
    const bool repair = role == RepairTimingRole::RepairDirect;
    if (!repair && role != RepairTimingRole::DispatchDirect) {
        return false;
    }
    const BandTimingDescriptor& descriptor =
        repair ? cell.RepairExplicit : cell.DispatchExplicit;
    const BandTimingPayload& payload =
        repair ? cell.RepairedPayload : cell.DispatchPayload;
    if ((repair ? !cell.HaveRepairExplicit :
                  !cell.HaveDispatchExplicit) ||
        payload.Result != Wirehair_Success)
    {
        return true;
    }
    if (payload.DataBytes.size() != cell.Ids.size() ||
        payload.Bytes.Size() !=
            cell.Ids.size() * (size_t)options.BlockBytes)
    {
        return false;
    }
    if (!ConfigureBandTimingWh2(
            descriptor.GF256Rows, descriptor.GF16Rows,
            descriptor.Period, descriptor.X))
    {
        return false;
    }
    wirehair_v2::PrecodeSystem system;
    if (!wirehair_v2::BuildPrecodeSystem(
            descriptor.Explicit.Params, system))
    {
        return false;
    }
    const uint64_t precode_count =
        (uint64_t)system.Params.Staircase +
        system.Params.DenseRows + system.Params.HeavyRows;
    wirehair_v2::PacketRowRuntime runtime;
    if (precode_count > UINT32_MAX ||
        !runtime.Initialize(
            options.BlockCount, (uint32_t)precode_count,
            descriptor.Explicit.Packet.MixCount))
    {
        return false;
    }
    if (cell.DirectPrefixSymbols < options.BlockCount ||
        cell.DirectPrefixSymbols >
            options.BlockCount + 64u ||
        cell.DirectPrefixSymbols > cell.Ids.size())
    {
        return false;
    }
    std::vector<wirehair_v2::SolvePacket> packets(
        cell.DirectPrefixSymbols);
    for (size_t i = 0u; i < cell.DirectPrefixSymbols; ++i)
    {
        packets[i].BlockId = cell.Ids[i];
        packets[i].Data = payload.Bytes.Data() +
            i * (size_t)options.BlockBytes;
    }
    std::vector<uint8_t> intermediate;
    PeelTimingUsage usage_before;
    int cpu_before = -1;
    uint64_t begin_ns = 0u;
    if (measure &&
        !BeginRepairTimingInterval(
            options, eviction, usage_before, cpu_before, begin_ns))
    {
        return false;
    }
    operation.ConstructResult = Wirehair_Success;
    operation.Result =
        wirehair_v2::
        SolvePrecodeSystemForValidatedSystemWithRuntime(
            system, descriptor.Explicit.Packet, runtime, packets,
            options.BlockBytes, intermediate, &operation.Stats);
    operation.StatsAvailable = true;
    if (measure &&
        !FinishRepairTimingInterval(
            usage_before, cpu_before, begin_ns, operation))
    {
        return false;
    }
    if (!RepairTimingWitnessOutsideInterval(options)) {
        return false;
    }
    operation.IntermediateBytes = intermediate.size();
    operation.IntermediateSha256 =
        Sha256HexBytes(intermediate.data(), intermediate.size());
    operation.SignatureSha256 = RepairTimingOperationSignature(
        operation, payload.Sha256, std::string(), std::string());
    operation.Available =
        operation.ConstructResult == Wirehair_Success &&
        operation.Result == Wirehair_Success &&
        operation.StatsAvailable &&
        operation.IntermediateBytes == payload.IntermediateBytes &&
        operation.IntermediateSha256 ==
            payload.IntermediateSha256;
    return true;
}

bool RunRepairTimingOperation(
    const RepairTimingOptions& options,
    const RepairTimingCell& cell,
    RepairTimingRole role,
    std::vector<uint8_t>& eviction,
    bool measure,
    RepairTimingOperation& operation)
{
    switch (role)
    {
    case RepairTimingRole::RepairSelectorEncoder:
    case RepairTimingRole::RepairForcedEncoder:
    case RepairTimingRole::Wh1Encoder:
        return RunRepairTimingEncoderOperation(
            options, cell, role, eviction, measure, operation);
    case RepairTimingRole::RepairDecoderFeed:
    case RepairTimingRole::Wh1DecoderFeed:
    case RepairTimingRole::RepairDecoderFull:
    case RepairTimingRole::Wh1DecoderFull:
        return RunRepairTimingDecoderOperation(
            options, cell, role, eviction, measure, operation);
    case RepairTimingRole::RepairDirect:
    case RepairTimingRole::DispatchDirect:
        return RunRepairTimingDirectOperation(
            options, cell, role, eviction, measure, operation);
    }
    return false;
}

bool SameRepairTimingOperation(
    const RepairTimingOperation& left,
    const RepairTimingOperation& right)
{
    return left.Available == right.Available &&
        left.ConstructResult == right.ConstructResult &&
        left.Result == right.Result &&
        left.RecoveryResult == right.RecoveryResult &&
        left.RecoveryOk == right.RecoveryOk &&
        left.ReceivedSymbols == right.ReceivedSymbols &&
        left.StatsAvailable == right.StatsAvailable &&
        (!left.StatsAvailable ||
         SameRepairTimingStatsStructural(
            left.Stats, right.Stats)) &&
        left.IntermediateBytes == right.IntermediateBytes &&
        left.IntermediateSha256 == right.IntermediateSha256 &&
        left.SelectorStructuralSha256 ==
            right.SelectorStructuralSha256 &&
        left.SignatureSha256 == right.SignatureSha256;
}

struct RepairTimingField
{
    std::string Name;
    bool Numeric = false;

    RepairTimingField(const std::string& name, bool numeric)
        : Name(name)
        , Numeric(numeric)
    {
    }
};

const std::vector<RepairTimingField>& RepairTimingSchema()
{
    static const std::vector<RepairTimingField> kFields = []() {
        std::vector<RepairTimingField> fields;
        const auto add = [&](const std::string& name, bool numeric) {
            fields.push_back({name, numeric});
        };
        add("row_kind", false);
        add("cell_index", true);
        add("measured", true);
        add("construction_root", true);
        add("loss_seed", true);
        add("trace_sha256", false);
        add("message_sha256", false);

        add("selector_result", true);
        add("selected_attempt", true);
        add("attempts_executed", true);
        add("calls_executed", true);
        add("real_configuration_calls", true);
        add("structural_probe_calls", true);
        add("cap_exhausted", true);
        add("fatal_attempt_zero_mismatch", true);
        add("oom", true);
        add("committed", true);
        add("descriptor_hex", false);
        add("descriptor_sha256", false);
        add("selector_structural_sha256", false);
        static const char* const roles[] = {
            "raw", "repaired", "dispatch", "wh1"
        };
        static const char* const role_suffixes[] = {
            "encode_result",
            "decode_construct_result",
            "feed_result",
            "recover_result",
            "recovery_ok",
            "encoded_symbols",
            "received_symbols",
            "overhead",
            "payload_sha256",
            "recovered_sha256"
        };
        for (const char* role : roles)
        {
            for (unsigned i = 0u;
                 i < sizeof(role_suffixes) /
                        sizeof(role_suffixes[0]); ++i)
            {
                add(
                    std::string(role) + "_" + role_suffixes[i],
                    i < 8u);
            }
        }
        add("forced_a_result", true);
        add("forced_b_result", true);
        add("forced_a_payload_sha256", false);
        add("forced_b_payload_sha256", false);
        add("forced_equal", true);
        add("direct_fixed_prefix_symbols", true);
        add("repair_intermediate_bytes", true);
        add("repair_intermediate_sha256", false);
        add("repair_direct_executed", true);
        add("repair_direct_result", true);
        add("repair_direct_intermediate_bytes", true);
        add("repair_direct_intermediate_sha256", false);
        add("repair_direct_witness_equal", true);
        add("dispatch_intermediate_bytes", true);
        add("dispatch_intermediate_sha256", false);
        add("dispatch_direct_executed", true);
        add("dispatch_direct_result", true);
        add("dispatch_direct_intermediate_bytes", true);
        add("dispatch_direct_intermediate_sha256", false);
        add("dispatch_direct_witness_equal", true);

        static const char* const stats_names[] = {
            "packet_rows",
            "peeled_columns",
            "inactivated_columns",
            "residual_rows",
            "residual_rank",
            "binary_residual_rank",
            "binary_row_references",
            "binary_row_storage_bytes",
            "binary_adjacency_storage_bytes",
            "binary_row_storage_allocations",
            "binary_adjacency_storage_allocations",
            "block_xors",
            "block_muladds",
            "block_copies",
            "block_zero_fills",
            "block_add_sets",
            "block_add_set_sources",
            "peel_adjacency_visits",
            "peel_row_scan_steps",
            "peel_heap_operations",
            "projection_word_xors",
            "residual_coeff_word_xors",
            "residual_coeff_byte_ops",
            "build_ns",
            "peel_ns",
            "project_ns",
            "residual_ns",
            "backsub_ns",
            "packet_seed_attempt",
            "tiny_mixed_fastpath_acceptances",
            "mixed_joint_source_xors",
            "mixed_joint_marginal_xors",
            "mixed_joint_marginal_copies",
            "mixed_joint_scratch_bytes",
            "mixed_joint_active_deltas",
            "mixed_dual_source_columns"
        };
        add("attempt", true);
        add("probe_executed", true);
        add("probe_result", true);
        add("probe_stats_available", true);
        for (const char* name : stats_names) {
            add(std::string("probe_") + name, true);
        }
        add("real_executed", true);
        add("real_result", true);
        add("real_stats_available", true);
        for (const char* name : stats_names) {
            add(std::string("real_") + name, true);
        }

        add("timing_panel", false);
        add("timing_panel_index", true);
        add("timing_slot", true);
        add("timing_pair", true);
        add("timing_label", false);
        add("timing_role", false);
        add("timing_scope", false);
        add("timing_eligible", true);
        add("timing_executed", true);
        add("timing_censor_reason", false);
        add("timing_construct_result", true);
        add("timing_result", true);
        add("timing_recover_result", true);
        add("timing_recovery_ok", true);
        add("timing_stable", true);
        add("timing_elapsed_ns", true);
        add("timing_inner_reps", true);
        add("timing_saturated", true);
        add("timing_cpu_before", true);
        add("timing_cpu_after", true);
        add("timing_cpu_migrated", true);
        add("timing_minflt", true);
        add("timing_majflt", true);
        add("timing_fault_contaminated", true);
        add("timing_stats_available", true);
        for (const char* name : stats_names) {
            add(std::string("timing_") + name, true);
        }
        add("timing_fixed_prefix_symbols", true);
        add("timing_intermediate_bytes", true);
        add("timing_intermediate_sha256", false);
        return fields;
    }();
    return kFields;
}

struct RepairTimingRow
{
    std::vector<std::string> Values;

    explicit RepairTimingRow(
        const std::vector<RepairTimingField>& schema)
    {
        Values.reserve(schema.size());
        for (const RepairTimingField& field : schema) {
            Values.push_back(
                field.Numeric ? "-1" : "not_applicable");
        }
    }

    void Set(size_t index, const std::string& value)
    {
        if (index >= Values.size()) {
            throw std::logic_error("repairtiming row index overflow");
        }
        Values[index] = value;
    }

    void Set(size_t index, uint64_t value)
    {
        Set(index, std::to_string(value));
    }

    void SetSigned(size_t index, int64_t value)
    {
        Set(index, std::to_string(value));
    }

    std::string Text() const
    {
        std::ostringstream text;
        for (size_t i = 0u; i < Values.size(); ++i)
        {
            if (i != 0u) text << ',';
            text << Values[i];
        }
        text << '\n';
        const std::string result = text.str();
        if (result.size() > 32768u) {
            throw std::logic_error("repairtiming row exceeds limit");
        }
        return result;
    }
};

static const size_t kRepairTimingCommonFields = 7u;
static const size_t kRepairTimingSummaryStart =
    kRepairTimingCommonFields;
static const size_t kRepairTimingAttemptStart =
    kRepairTimingSummaryStart + 73u;
static const size_t kRepairTimingTimingStart =
    kRepairTimingAttemptStart + 79u;
static const size_t kRepairTimingFieldCount =
    kRepairTimingTimingStart + 64u;

const std::vector<RepairTimingField>& RepairTimingTypedSchema(
    const char* kind)
{
    static const std::vector<RepairTimingField> kSummary = []() {
        const std::vector<RepairTimingField>& all =
            RepairTimingSchema();
        return std::vector<RepairTimingField>(
            all.begin(), all.begin() + kRepairTimingAttemptStart);
    }();
    static const std::vector<RepairTimingField> kAttempt = []() {
        const std::vector<RepairTimingField>& all =
            RepairTimingSchema();
        std::vector<RepairTimingField> result(
            all.begin(), all.begin() + kRepairTimingCommonFields);
        result.insert(
            result.end(),
            all.begin() + kRepairTimingAttemptStart,
            all.begin() + kRepairTimingTimingStart);
        return result;
    }();
    static const std::vector<RepairTimingField> kTiming = []() {
        const std::vector<RepairTimingField>& all =
            RepairTimingSchema();
        std::vector<RepairTimingField> result(
            all.begin(), all.begin() + kRepairTimingCommonFields);
        result.insert(
            result.end(),
            all.begin() + kRepairTimingTimingStart,
            all.end());
        return result;
    }();
    if (!std::strcmp(kind, "summary")) return kSummary;
    if (!std::strcmp(kind, "attempt")) return kAttempt;
    if (!std::strcmp(kind, "timing")) return kTiming;
    throw std::logic_error("repairtiming unknown typed schema");
}

void InitializeRepairTimingCommonRow(
    RepairTimingRow& row,
    const RepairTimingCell& cell,
    const RepairTimingOptions& options,
    const char* kind)
{
    row.Set(0u, kind);
    row.Set(1u, cell.Index);
    row.Set(
        2u,
        cell.Index >= options.WarmupReplicates ? 1u : 0u);
    row.Set(3u, cell.ConstructionRoot);
    row.Set(4u, cell.LossSeed);
    row.Set(5u, cell.TraceSha256);
    row.Set(6u, cell.MessageSha256);
}

void SetRepairTimingStatsRow(
    RepairTimingRow& row,
    size_t& index,
    const wirehair_v2::PrecodeSolveStats* stats,
    bool available)
{
    // Frozen repairtiming v3 36-field row subset.  Do not append later
    // PrecodeSolveStats members without versioning the native and parser
    // schemas together.
    if (!available || !stats)
    {
        index += 36u;
        return;
    }
    const uint64_t values[] = {
        stats->PacketRows,
        stats->PeeledColumns,
        stats->InactivatedColumns,
        stats->ResidualRows,
        stats->ResidualRank,
        stats->BinaryResidualRank,
        stats->BinaryRowReferences,
        stats->BinaryRowStorageBytes,
        stats->BinaryAdjacencyStorageBytes,
        stats->BinaryRowStorageAllocations,
        stats->BinaryAdjacencyStorageAllocations,
        stats->BlockXors,
        stats->BlockMulAdds,
        stats->BlockCopies,
        stats->BlockZeroFills,
        stats->BlockAddSets,
        stats->BlockAddSetSources,
        stats->PeelAdjacencyVisits,
        stats->PeelRowScanSteps,
        stats->PeelHeapOperations,
        stats->ProjectionWordXors,
        stats->ResidualCoeffWordXors,
        stats->ResidualCoeffByteOps,
        stats->BuildNanoseconds,
        stats->PeelNanoseconds,
        stats->ProjectNanoseconds,
        stats->ResidualNanoseconds,
        stats->BackSubNanoseconds,
        stats->PacketSeedAttempt,
        stats->TinyMixedFastPathAcceptances,
        stats->MixedJointSourceXors,
        stats->MixedJointMarginalXors,
        stats->MixedJointMarginalCopies,
        stats->MixedJointScratchBytes,
        stats->MixedJointActiveDeltas,
        stats->MixedDualSourceColumns
    };
    static_assert(
        sizeof(values) / sizeof(values[0]) == 36u,
        "repairtiming v3 stats row must retain 36 fields");
    for (uint64_t value : values) {
        row.Set(index++, value);
    }
}

void SetRepairTimingRecoveryGroup(
    RepairTimingRow& row,
    size_t& index,
    const BandTimingPayload& payload,
    const BandTimingRecovery& recovery)
{
    row.SetSigned(index++, (int)payload.Result);
    if (payload.Result == Wirehair_Success)
    {
        row.SetSigned(index, (int)recovery.ConstructResult);
    }
    ++index;
    if (payload.Result == Wirehair_Success &&
        recovery.ConstructResult == Wirehair_Success)
    {
        row.SetSigned(index, (int)recovery.Result);
    }
    ++index;
    if (payload.Result == Wirehair_Success &&
        recovery.ConstructResult == Wirehair_Success &&
        recovery.Result == Wirehair_Success)
    {
        row.SetSigned(index, (int)recovery.RecoveryResult);
    }
    ++index;
    row.Set(index++, recovery.RecoveryOk ? 1u : 0u);
    row.Set(
        index++,
        payload.Result == Wirehair_Success ?
            (uint64_t)payload.DataBytes.size() : 0u);
    row.Set(index++, recovery.ReceivedSymbols);
    row.SetSigned(index++, recovery.Overhead);
    if (IsPeelTimingLowerHexSha256(payload.Sha256)) {
        row.Set(index, payload.Sha256);
    }
    ++index;
    if (IsPeelTimingLowerHexSha256(recovery.RecoveredSha256)) {
        row.Set(index, recovery.RecoveredSha256);
    }
    ++index;
}

std::string MakeRepairTimingSummaryRow(
    const RepairTimingOptions& options,
    const RepairTimingCell& cell)
{
    RepairTimingRow row(RepairTimingTypedSchema("summary"));
    InitializeRepairTimingCommonRow(row, cell, options, "summary");
    size_t index = kRepairTimingSummaryStart;
    row.SetSigned(index++, (int)cell.Selector.Result);
    row.SetSigned(
        index++,
        cell.Selector.SelectedAttempt ==
                wirehair_v2::kRepairV1NoAttempt ?
            -1 : (int64_t)cell.Selector.SelectedAttempt);
    row.Set(index++, cell.Selector.AttemptsExecuted);
    row.Set(index++, cell.Selector.CallsExecuted);
    row.Set(index++, cell.Selector.RealConfigurationCalls);
    row.Set(index++, cell.Selector.StructuralProbeCalls);
    row.Set(index++, cell.Selector.CapExhausted ? 1u : 0u);
    row.Set(
        index++,
        cell.Selector.FatalAttemptZeroMismatch ? 1u : 0u);
    row.Set(index++, cell.Selector.Oom ? 1u : 0u);
    row.Set(index++, cell.Selector.Committed ? 1u : 0u);
    if (cell.HaveDescriptor)
    {
        row.Set(index, cell.Descriptor.Hex);
        row.Set(index + 1u, cell.Descriptor.Sha256);
    }
    index += 2u;
    row.Set(index++, cell.StructuralSha256);

    SetRepairTimingRecoveryGroup(
        row, index, cell.RawPayload, cell.RawRecovery);
    SetRepairTimingRecoveryGroup(
        row, index, cell.RepairedPayload, cell.RepairedRecovery);
    SetRepairTimingRecoveryGroup(
        row, index, cell.DispatchPayload, cell.DispatchRecovery);
    SetRepairTimingRecoveryGroup(
        row, index, cell.Wh1Payload, cell.Wh1Recovery);

    if (cell.HaveDescriptor)
    {
        row.SetSigned(index, (int)cell.ForcedAPayload.Result);
        row.SetSigned(index + 1u, (int)cell.ForcedBPayload.Result);
        if (IsPeelTimingLowerHexSha256(cell.ForcedAPayload.Sha256)) {
            row.Set(index + 2u, cell.ForcedAPayload.Sha256);
        }
        if (IsPeelTimingLowerHexSha256(cell.ForcedBPayload.Sha256)) {
            row.Set(index + 3u, cell.ForcedBPayload.Sha256);
        }
    }
    index += 4u;
    const bool forced_equal =
        cell.HaveDescriptor &&
        cell.ForcedAPayload.Result ==
            cell.RepairedPayload.Result &&
        cell.ForcedBPayload.Result ==
            cell.RepairedPayload.Result &&
        cell.ForcedAPayload.Sha256 ==
            cell.RepairedPayload.Sha256 &&
        cell.ForcedBPayload.Sha256 ==
            cell.RepairedPayload.Sha256 &&
        cell.ForcedAPayload.IntermediateBytes ==
            cell.RepairedPayload.IntermediateBytes &&
        cell.ForcedBPayload.IntermediateBytes ==
            cell.RepairedPayload.IntermediateBytes &&
        cell.ForcedAPayload.IntermediateSha256 ==
            cell.RepairedPayload.IntermediateSha256 &&
        cell.ForcedBPayload.IntermediateSha256 ==
            cell.RepairedPayload.IntermediateSha256;
    row.Set(index++, forced_equal ? 1u : 0u);

    row.Set(index++, cell.DirectPrefixSymbols);
    if (cell.RepairedPayload.Result == Wirehair_Success)
    {
        row.Set(index, cell.RepairedPayload.IntermediateBytes);
        if (IsPeelTimingLowerHexSha256(
                cell.RepairedPayload.IntermediateSha256))
        {
            row.Set(index + 1u,
                cell.RepairedPayload.IntermediateSha256);
        }
    }
    index += 2u;
    row.Set(index++, cell.RepairDirectExecuted ? 1u : 0u);
    if (cell.RepairDirectExecuted) {
        row.SetSigned(index, (int)cell.RepairDirect.Result);
    }
    if (cell.RepairDirectExecuted &&
        cell.RepairDirect.Result == Wirehair_Success)
    {
        row.Set(index + 1u, cell.RepairDirect.IntermediateBytes);
        if (IsPeelTimingLowerHexSha256(
                cell.RepairDirect.IntermediateSha256))
        {
            row.Set(index + 2u,
                cell.RepairDirect.IntermediateSha256);
        }
    }
    index += 3u;
    row.Set(
        index++,
        cell.RepairDirectExecuted &&
        cell.RepairDirect.Result == Wirehair_Success &&
        cell.RepairDirect.IntermediateBytes ==
            cell.RepairedPayload.IntermediateBytes &&
        cell.RepairDirect.IntermediateSha256 ==
            cell.RepairedPayload.IntermediateSha256 ? 1u : 0u);

    if (cell.DispatchPayload.Result == Wirehair_Success)
    {
        row.Set(index, cell.DispatchPayload.IntermediateBytes);
        if (IsPeelTimingLowerHexSha256(
                cell.DispatchPayload.IntermediateSha256))
        {
            row.Set(index + 1u,
                cell.DispatchPayload.IntermediateSha256);
        }
    }
    index += 2u;
    row.Set(index++, cell.DispatchDirectExecuted ? 1u : 0u);
    if (cell.DispatchDirectExecuted) {
        row.SetSigned(index, (int)cell.DispatchDirect.Result);
    }
    if (cell.DispatchDirectExecuted &&
        cell.DispatchDirect.Result == Wirehair_Success)
    {
        row.Set(index + 1u, cell.DispatchDirect.IntermediateBytes);
        if (IsPeelTimingLowerHexSha256(
                cell.DispatchDirect.IntermediateSha256))
        {
            row.Set(index + 2u,
                cell.DispatchDirect.IntermediateSha256);
        }
    }
    index += 3u;
    row.Set(
        index++,
        cell.DispatchDirectExecuted &&
        cell.DispatchDirect.Result == Wirehair_Success &&
        cell.DispatchDirect.IntermediateBytes ==
            cell.DispatchPayload.IntermediateBytes &&
        cell.DispatchDirect.IntermediateSha256 ==
            cell.DispatchPayload.IntermediateSha256 ? 1u : 0u);
    if (index != kRepairTimingAttemptStart ||
        row.Values.size() != 80u)
    {
        throw std::logic_error("repairtiming summary width mismatch");
    }
    return row.Text();
}

std::string MakeRepairTimingAttemptRow(
    const RepairTimingOptions& options,
    const RepairTimingCell& cell,
    uint32_t attempt)
{
    RepairTimingRow row(RepairTimingTypedSchema("attempt"));
    InitializeRepairTimingCommonRow(row, cell, options, "attempt");
    const wirehair_v2::RepairV1AttemptTelemetry& observation =
        cell.Selector.Attempts[attempt];
    size_t index = kRepairTimingCommonFields;
    row.Set(index++, attempt);
    row.Set(index++, observation.ProbeExecuted ? 1u : 0u);
    if (observation.ProbeExecuted) {
        row.SetSigned(index, (int)observation.ProbeResult);
    }
    ++index;
    row.Set(index++, observation.ProbeStatsAvailable ? 1u : 0u);
    SetRepairTimingStatsRow(
        row, index,
        observation.ProbeStatsAvailable ?
            &observation.ProbeStats : nullptr,
        observation.ProbeStatsAvailable);
    row.Set(index++, observation.RealExecuted ? 1u : 0u);
    if (observation.RealExecuted) {
        row.SetSigned(index, (int)observation.RealResult);
    }
    ++index;
    row.Set(index++, observation.RealStatsAvailable ? 1u : 0u);
    SetRepairTimingStatsRow(
        row, index,
        observation.RealStatsAvailable ?
            &observation.RealStats : nullptr,
        observation.RealStatsAvailable);
    if (index != 86u || row.Values.size() != 86u) {
        throw std::logic_error("repairtiming attempt width mismatch");
    }
    return row.Text();
}

bool AccumulateRepairTimingMeasurement(
    const RepairTimingOptions& options,
    const RepairTimingCell& cell,
    RepairTimingRole role,
    const RepairTimingOperation& expected,
    std::vector<uint8_t>& eviction,
    RepairTimingOperation& measured,
    bool& stable,
    bool& saturated)
{
    measured = RepairTimingOperation();
    stable = true;
    saturated = false;
    for (uint32_t inner = 0u; inner < options.InnerReps; ++inner)
    {
        RepairTimingOperation one;
        if (!RunRepairTimingOperation(
                options, cell, role, eviction, true, one))
        {
            return false;
        }
        stable = stable &&
            SameRepairTimingOperation(one, expected);
        if (inner == 0u)
        {
            measured = one;
        }
        else
        {
            BandTimingAddU64(
                measured.ElapsedNanoseconds,
                one.ElapsedNanoseconds, saturated);
            if (measured.CpuBefore != one.CpuBefore ||
                measured.CpuAfter != one.CpuAfter)
            {
                stable = false;
            }
            if (one.MinorFaults >
                    std::numeric_limits<int64_t>::max() -
                        measured.MinorFaults ||
                one.MajorFaults >
                    std::numeric_limits<int64_t>::max() -
                        measured.MajorFaults)
            {
                saturated = true;
            }
            else
            {
                measured.MinorFaults += one.MinorFaults;
                measured.MajorFaults += one.MajorFaults;
            }
            if (measured.StatsAvailable && one.StatsAvailable)
            {
                BandTimingAddU64(
                    measured.Stats.BuildNanoseconds,
                    one.Stats.BuildNanoseconds, saturated);
                BandTimingAddU64(
                    measured.Stats.PeelNanoseconds,
                    one.Stats.PeelNanoseconds, saturated);
                BandTimingAddU64(
                    measured.Stats.ProjectNanoseconds,
                    one.Stats.ProjectNanoseconds, saturated);
                BandTimingAddU64(
                    measured.Stats.ResidualNanoseconds,
                    one.Stats.ResidualNanoseconds, saturated);
                BandTimingAddU64(
                    measured.Stats.BackSubNanoseconds,
                    one.Stats.BackSubNanoseconds, saturated);
            }
        }
    }
    return true;
}

std::string MakeRepairTimingTimingRow(
    const RepairTimingOptions& options,
    const RepairTimingCell& cell,
    const RepairTimingPanel& panel,
    uint32_t panel_index,
    uint32_t slot,
    char label,
    bool eligible,
    const char* censor_reason,
    const RepairTimingOperation* measured,
    bool stable,
    bool saturated)
{
    RepairTimingRow row(RepairTimingTypedSchema("timing"));
    InitializeRepairTimingCommonRow(row, cell, options, "timing");
    const RepairTimingRole role =
        label == 'A' ? panel.A : panel.B;
    size_t index = kRepairTimingCommonFields;
    row.Set(index++, panel.Name);
    row.Set(index++, panel_index);
    row.Set(index++, slot);
    row.Set(index++, slot / 2u);
    row.Set(index++, std::string(1u, label));
    row.Set(index++, RepairTimingRoleName(role));
    row.Set(index++, RepairTimingScopeName(role));
    row.Set(index++, eligible ? 1u : 0u);
    row.Set(index++, measured ? 1u : 0u);
    row.Set(index++, censor_reason);
    if (measured)
    {
        row.SetSigned(index++, (int)measured->ConstructResult);
        row.SetSigned(index++, (int)measured->Result);
        if ((role == RepairTimingRole::RepairDecoderFull ||
             role == RepairTimingRole::Wh1DecoderFull) &&
            measured->Result == Wirehair_Success)
        {
            row.SetSigned(index, (int)measured->RecoveryResult);
        }
        ++index;
        row.Set(index++, measured->RecoveryOk ? 1u : 0u);
        row.Set(index++, stable ? 1u : 0u);
        row.Set(index++, measured->ElapsedNanoseconds);
        row.Set(index++, options.InnerReps);
        row.Set(index++, saturated ? 1u : 0u);
        row.SetSigned(index++, measured->CpuBefore);
        row.SetSigned(index++, measured->CpuAfter);
        row.Set(
            index++,
            PeelTimingCpuMigrated(
                measured->CpuBefore, measured->CpuAfter));
        row.SetSigned(index++, measured->MinorFaults);
        row.SetSigned(index++, measured->MajorFaults);
        row.Set(
            index++,
            PeelTimingFaultContaminated(
                measured->MinorFaults,
                measured->MajorFaults));
        row.Set(index++, measured->StatsAvailable ? 1u : 0u);
        SetRepairTimingStatsRow(
            row, index,
            measured->StatsAvailable ? &measured->Stats : nullptr,
            measured->StatsAvailable);
        if (role == RepairTimingRole::RepairDirect ||
            role == RepairTimingRole::DispatchDirect)
        {
            row.Set(index, cell.DirectPrefixSymbols);
        }
        ++index;
        const bool witness_role =
            measured->Result == Wirehair_Success &&
            (
                role == RepairTimingRole::RepairSelectorEncoder ||
                role == RepairTimingRole::RepairForcedEncoder ||
                role == RepairTimingRole::RepairDirect ||
                role == RepairTimingRole::DispatchDirect
            );
        if (witness_role)
        {
            row.Set(index, measured->IntermediateBytes);
            if (IsPeelTimingLowerHexSha256(
                    measured->IntermediateSha256))
            {
                row.Set(index + 1u,
                    measured->IntermediateSha256);
            }
        }
        index += 2u;
    }
    else {
        index += 15u + 36u + 3u;
    }
    if (index != 71u || row.Values.size() != 71u ||
        RepairTimingSchema().size() != kRepairTimingFieldCount)
    {
        throw std::logic_error(
            std::string("repairtiming timing width mismatch index=") +
            std::to_string(index) + " typed=" +
            std::to_string(row.Values.size()) + " all=" +
            std::to_string(RepairTimingSchema().size()) +
            " expected=" +
            std::to_string(kRepairTimingFieldCount));
    }
    return row.Text();
}

BandTimingOptions RepairTimingBandOptions(
    const RepairTimingOptions& options)
{
    BandTimingOptions band;
    band.BlockCount = options.BlockCount;
    band.BlockBytes = options.BlockBytes;
    band.ConstructionSeedBase = options.ConstructionSeedBase;
    band.Loss = options.Loss;
    band.LossSeedBase = options.LossSeedBase;
    band.Schedule = options.Schedule;
    band.WarmupReplicates = options.WarmupReplicates;
    band.Replicates = options.Replicates;
    band.InnerReps = options.InnerReps;
    band.MaxOverhead = options.MaxOverhead;
    band.CacheState = options.CacheState;
    band.SystematicCache = options.SystematicCache;
    band.EvictBytes = options.EvictBytes;
    band.ContextSha256 = options.ContextSha256;
    band.RequiredMargin = options.RequiredMargin;
    band.DispatchContract = options.DispatchContract;
    return band;
}

std::string RepairTimingTraceDigest(
    const RepairTimingOptions& options,
    uint32_t construction_root,
    uint64_t loss_seed,
    const std::vector<uint32_t>& ids)
{
    std::ostringstream trace;
    trace
        << "wirehair-wh2-repairtiming-loss-trace-v1\n"
        << "K=" << options.BlockCount << "\n"
        << "block_bytes=" << options.BlockBytes << "\n"
        << "construction_root=" << construction_root << "\n"
        << "loss_seed=" << loss_seed << "\n"
        << "loss=" << std::setprecision(17) << options.Loss << "\n"
        << "schedule=" << PacketScheduleName(options.Schedule) << "\n"
        << "message_seed_policy="
            "replicate-loss-seed-partial-final-v1\n"
        << "ids=";
    for (size_t i = 0u; i < ids.size(); ++i)
    {
        if (i != 0u) trace << ',';
        trace << ids[i];
    }
    trace << '\n';
    return Sha256Hex(trace.str());
}

bool SameRepairTimingDescriptor(
    const RepairTimingDescriptor& left,
    const RepairTimingDescriptor& right)
{
    return left.Contract == right.Contract &&
        left.ConstructionRoot == right.ConstructionRoot &&
        left.SelectedAttempt == right.SelectedAttempt &&
        left.Bytes == right.Bytes &&
        left.Hex == right.Hex &&
        left.Sha256 == right.Sha256;
}

bool ExerciseRepairTimingDescriptorFault(
    const RepairTimingDescriptor& valid,
    const std::string& fault,
    const RepairTimingOptions& options,
    const RepairTimingCell& cell)
{
    RepairTimingDescriptor destination = valid;
    const RepairTimingDescriptor original = destination;
    std::vector<uint8_t> malformed = valid.Bytes;
    if (fault == "bad-magic") {
        malformed[0] ^= 1u;
    }
    else if (fault == "bad-id") {
        malformed[4] ^= 1u;
    }
    else if (fault == "bad-attempt") {
        malformed[16] = (uint8_t)wirehair_v2::kRepairV1AttemptCount;
    }
    else if (fault == "trailing") {
        malformed.push_back(0u);
    }
    else {
        return false;
    }

    wirehair_v2::MessagePrecodeDecoder decoder;
    if (decoder.InitializeRepairV1ResultForTesting(
            cell.Message.Size(), options.BlockBytes,
            *valid.Contract, valid.ConstructionRoot,
            valid.SelectedAttempt, options.SystematicCache) !=
            Wirehair_Success)
    {
        return false;
    }
    const uint32_t received_before = decoder.ReceivedCount();
    const bool parsed =
        ParseRepairTimingDescriptor(malformed, destination);
    const wirehair_v2::RepairV1Contract* view_contract =
        valid.Contract;
    const uint32_t view_root =
        valid.ConstructionRoot ^ UINT32_C(0xa5a5a5a5);
    const uint32_t view_attempt = wirehair_v2::kRepairV1NoAttempt;
    const wirehair_v2::RepairV1Contract* parsed_view_contract =
        view_contract;
    uint32_t parsed_view_root = view_root;
    uint32_t parsed_view_attempt = view_attempt;
    const bool view_parsed =
        ParseRepairTimingDescriptorView(
            malformed, parsed_view_contract,
            parsed_view_root, parsed_view_attempt);
    return !parsed &&
        !view_parsed &&
        SameRepairTimingDescriptor(destination, original) &&
        parsed_view_contract == view_contract &&
        parsed_view_root == view_root &&
        parsed_view_attempt == view_attempt &&
        decoder.ReceivedCount() == received_before &&
        !decoder.IsDecoded();
}

bool BuildRepairTimingCell(
    const RepairTimingOptions& options,
    uint32_t cell_index,
    RepairTimingCell& cell,
    bool& descriptor_fault_exercised)
{
    cell = RepairTimingCell();
    cell.Index = cell_index;
    cell.ConstructionRoot =
        RepairTimingConstructionRoot(options, cell_index);
    cell.LossSeed =
        BandTimingDerivedLossSeed(options.LossSeedBase, cell_index);
    cell.Ids = BuildPacketSchedule(
        options.BlockCount,
        options.BlockCount + options.MaxOverhead,
        options.Loss, cell.LossSeed, options.Schedule);
    if (cell.Ids.size() !=
        (size_t)options.BlockCount + options.MaxOverhead)
    {
        return false;
    }
    cell.FirstKIds.resize(options.BlockCount);
    for (uint32_t id = 0u; id < options.BlockCount; ++id) {
        cell.FirstKIds[id] = id;
    }
    const uint64_t message_bytes =
        ((uint64_t)options.BlockCount - 1u) *
            options.BlockBytes + 1u;
    if (message_bytes >
            (uint64_t)std::numeric_limits<size_t>::max() ||
        !FillBandTimingMessage(
            cell.Message, (size_t)message_bytes, cell.LossSeed))
    {
        return false;
    }
    if (std::find_if(
            cell.Message.Data(),
            cell.Message.Data() + cell.Message.Size(),
            [](uint8_t value) { return value != 0u; }) ==
        cell.Message.Data() + cell.Message.Size())
    {
        return false;
    }
    cell.MessageSha256 =
        Sha256HexBytes(cell.Message.Data(), cell.Message.Size());
    cell.TraceSha256 = RepairTimingTraceDigest(
        options, cell.ConstructionRoot, cell.LossSeed, cell.Ids);

    wirehair_v2::PrecodeSolveStats ignored_stats;
    (void)BuildRepairTimingPayload(
        options, cell.ConstructionRoot, cell.Message, cell.Ids,
        RepairTimingEncoderKind::RawAttemptZero, 0u,
        cell.RawPayload, nullptr, &ignored_stats);
    (void)BuildRepairTimingPayload(
        options, cell.ConstructionRoot, cell.Message, cell.Ids,
        RepairTimingEncoderKind::LazySelector, 0u,
        cell.RepairedPayload, &cell.Selector, &ignored_stats);

    if (cell.Selector.Committed)
    {
        if (cell.Selector.Result != Wirehair_Success ||
            cell.Selector.SelectedAttempt >=
                wirehair_v2::kRepairV1AttemptCount ||
            !BuildRepairTimingDescriptor(
                *options.Contract, cell.ConstructionRoot,
                cell.Selector.SelectedAttempt, cell.Descriptor))
        {
            return false;
        }
        cell.HaveDescriptor = true;
        RepairTimingDescriptor parsed;
        if (!ParseRepairTimingDescriptor(
                cell.Descriptor.Bytes, parsed) ||
            !SameRepairTimingDescriptor(parsed, cell.Descriptor))
        {
            return false;
        }
        if (!options.DescriptorFault.empty() &&
            !descriptor_fault_exercised)
        {
            if (!ExerciseRepairTimingDescriptorFault(
                    cell.Descriptor, options.DescriptorFault,
                    options, cell))
            {
                return false;
            }
            descriptor_fault_exercised = true;
            return false;
        }
        (void)BuildRepairTimingPayload(
            options, cell.ConstructionRoot, cell.Message, cell.Ids,
            RepairTimingEncoderKind::ForcedSelected,
            cell.Descriptor.SelectedAttempt,
            cell.ForcedAPayload, nullptr, &ignored_stats);
        (void)BuildRepairTimingPayload(
            options, cell.ConstructionRoot, cell.Message, cell.Ids,
            RepairTimingEncoderKind::ForcedSelected,
            cell.Descriptor.SelectedAttempt,
            cell.ForcedBPayload, nullptr, &ignored_stats);
        if (options.ForcedWitnessFault)
        {
            if (!IsPeelTimingLowerHexSha256(
                    cell.ForcedBPayload.IntermediateSha256))
            {
                return false;
            }
            char& first =
                cell.ForcedBPayload.IntermediateSha256[0];
            first = first == '0' ? '1' : '0';
        }
        if (cell.ForcedAPayload.Result !=
                cell.RepairedPayload.Result ||
            cell.ForcedBPayload.Result !=
                cell.RepairedPayload.Result ||
            cell.ForcedAPayload.Sha256 !=
                cell.RepairedPayload.Sha256 ||
            cell.ForcedBPayload.Sha256 !=
                cell.RepairedPayload.Sha256 ||
            cell.ForcedAPayload.IntermediateBytes !=
                cell.RepairedPayload.IntermediateBytes ||
            cell.ForcedBPayload.IntermediateBytes !=
                cell.RepairedPayload.IntermediateBytes ||
            cell.ForcedAPayload.IntermediateSha256 !=
                cell.RepairedPayload.IntermediateSha256 ||
            cell.ForcedBPayload.IntermediateSha256 !=
                cell.RepairedPayload.IntermediateSha256)
        {
            return false;
        }
        cell.RepairedRecovery = RunRepairTimingRecovery(
            options, cell.Descriptor, cell.Message,
            cell.Ids, cell.RepairedPayload);
        cell.HaveRepairExplicit =
            BuildRepairTimingExplicitDescriptor(
                options, cell.Descriptor, cell.RepairExplicit);
        if (!cell.HaveRepairExplicit) return false;
    }

    if (!cell.Selector.Attempts[0].RealExecuted ||
        cell.RawPayload.Result !=
            cell.Selector.Attempts[0].RealResult ||
        cell.RepairedPayload.Result != cell.Selector.Result ||
        cell.HaveDescriptor != cell.Selector.Committed)
    {
        return false;
    }

    RepairTimingDescriptor raw_descriptor;
    if (cell.RawPayload.Result == Wirehair_Success)
    {
        if (!BuildRepairTimingDescriptor(
                *options.Contract, cell.ConstructionRoot, 0u,
                raw_descriptor))
        {
            return false;
        }
        cell.RawRecovery = RunRepairTimingRecovery(
            options, raw_descriptor, cell.Message,
            cell.Ids, cell.RawPayload);
    }
    else
    {
        cell.RawRecovery.ConstructResult =
            cell.RawPayload.ConstructResult;
        cell.RawRecovery.Result = cell.RawPayload.Result;
    }

    const BandTimingOptions band = RepairTimingBandOptions(options);
    (void)BuildBandTimingRealDispatchPayload(
        band, cell.ConstructionRoot, cell.Message,
        cell.Ids, cell.DispatchPayload);
    cell.DispatchRecovery = RunBandTimingRealDispatchRecovery(
        band, cell.ConstructionRoot, cell.Message,
        cell.Ids, cell.DispatchPayload);
    (void)BuildBandTimingPayload(
        BandTimingRole::Wh1, band, nullptr,
        cell.ConstructionRoot, cell.Message,
        cell.Ids, cell.Wh1Payload);
    cell.Wh1Recovery = RunBandTimingRecovery(
        BandTimingRole::Wh1, band, nullptr,
        cell.ConstructionRoot, cell.Message,
        cell.Ids, cell.Wh1Payload);

    cell.HaveDispatchExplicit = BuildBandTimingDescriptor(
        BandTimingRole::Dispatch, band,
        cell.ConstructionRoot, cell.DispatchExplicit);
    if (!cell.HaveDispatchExplicit ||
        cell.DispatchExplicit.GF256Rows !=
            wirehair_v2::kMixedGF256Rows ||
        cell.DispatchExplicit.GF16Rows !=
            wirehair_v2::kMixedGF16Rows ||
        cell.DispatchExplicit.Explicit.Params.DenseRows != 4u ||
        cell.DispatchExplicit.Explicit.Params.HeavyRows !=
            wirehair_v2::kMixedGF256Rows +
                wirehair_v2::kMixedGF16Rows)
    {
        return false;
    }

    if (cell.RepairedRecovery.RecoveryOk &&
        cell.DispatchRecovery.RecoveryOk)
    {
        cell.DirectPrefixSymbols = std::max(
            cell.RepairedRecovery.ReceivedSymbols,
            cell.DispatchRecovery.ReceivedSymbols);
    }
    else
    {
        cell.DirectPrefixSymbols =
            options.BlockCount + 64u;
    }
    const uint64_t delivered =
        (uint64_t)options.BlockCount + options.MaxOverhead;
    if (cell.DirectPrefixSymbols < options.BlockCount ||
        cell.DirectPrefixSymbols > options.BlockCount + 64u ||
        cell.DirectPrefixSymbols > delivered)
    {
        return false;
    }

    if (cell.HaveRepairExplicit &&
        cell.RepairedPayload.Result == Wirehair_Success)
    {
        cell.RepairDirectExecuted = true;
        cell.RepairDirect = RunBandTimingDirect(
            band, cell.RepairExplicit, cell.Ids,
            cell.RepairedPayload, cell.DirectPrefixSymbols, true);
    }
    if (cell.DispatchPayload.Result == Wirehair_Success)
    {
        cell.DispatchDirectExecuted = true;
        cell.DispatchDirect = RunBandTimingDirect(
            band, cell.DispatchExplicit, cell.Ids,
            cell.DispatchPayload, cell.DirectPrefixSymbols, true);
    }
    if (cell.RepairDirectExecuted &&
        !BandTimingDirectWitnessValid(
            band, cell.RepairExplicit,
            cell.RepairDirect, cell.RepairedPayload))
    {
        return false;
    }
    if (cell.DispatchDirectExecuted &&
        !BandTimingDirectWitnessValid(
            band, cell.DispatchExplicit,
            cell.DispatchDirect, cell.DispatchPayload))
    {
        return false;
    }
    cell.StructuralSha256 = RepairTimingStructuralDigest(
        options, cell.ConstructionRoot, cell.Selector,
        cell.HaveDescriptor ? &cell.Descriptor : nullptr);
    return IsPeelTimingLowerHexSha256(cell.TraceSha256) &&
        IsPeelTimingLowerHexSha256(cell.MessageSha256) &&
        IsPeelTimingLowerHexSha256(cell.StructuralSha256);
}

std::string RepairTimingHeader(
    const std::vector<RepairTimingField>& schema)
{
    std::ostringstream text;
    for (size_t i = 0u; i < schema.size(); ++i)
    {
        if (i != 0u) text << ',';
        text << schema[i].Name;
    }
    const std::string header = text.str();
    if (header.size() > 65536u) {
        throw std::logic_error("repairtiming header exceeds limit");
    }
    return header;
}

bool RepairTimingRoleRequiresSelector(RepairTimingRole role)
{
    return
        role == RepairTimingRole::RepairSelectorEncoder ||
        role == RepairTimingRole::RepairForcedEncoder ||
        role == RepairTimingRole::RepairDecoderFeed ||
        role == RepairTimingRole::RepairDecoderFull ||
        role == RepairTimingRole::RepairDirect;
}

int CmdRepairTiming(int argc, char** argv)
{
    RepairTimingOptions options;
    if (!ParseRepairTimingOptions(argc, argv, options)) return 1;
    if (!ValidateRepairTimingEnvironment()) return 1;
    BandTimingHookCleanup hook_cleanup;
#if !defined(__linux__)
    std::fprintf(stderr,
        "repairtiming requires Linux CPU and resource-usage telemetry\n");
    return 1;
#endif
    if (!PinPeelTimingToFirstAllowedCpu())
    {
        std::fprintf(stderr,
            "repairtiming could not pin its first allowed CPU\n");
        return 2;
    }
    const uint64_t message_bytes =
        ((uint64_t)options.BlockCount - 1u) *
            options.BlockBytes + 1u;
    const uint64_t cells =
        (uint64_t)options.WarmupReplicates + options.Replicates;
    const uint64_t expected_summary_rows = cells;
    const uint64_t expected_attempt_rows = cells * 8u;
    const uint64_t expected_timing_rows = cells * 112u;
    const uint64_t expected_rows = cells * 121u;
    const uint64_t working_blocks =
        32u * ((uint64_t)options.BlockCount + 64u) +
        48u * options.BlockCount + 4096u;
    if (message_bytes >
            (uint64_t)std::numeric_limits<size_t>::max() ||
        working_blocks >
            (UINT64_MAX - options.EvictBytes) /
                options.BlockBytes ||
        working_blocks * options.BlockBytes +
            options.EvictBytes > kPeelTimingWorkingBytesMax)
    {
        std::fprintf(stderr,
            "repairtiming working set exceeds protocol limits\n");
        return 1;
    }

    uint64_t started_monotonic_ns = 0u;
    if (!ReadPeelTimingMonotonicNanoseconds(
            started_monotonic_ns))
    {
        std::fprintf(stderr,
            "repairtiming could not read monotonic start clock\n");
        return 2;
    }

    std::ostringstream panels;
    for (size_t i = 0u;
         i < sizeof(kRepairTimingPanels) /
                sizeof(kRepairTimingPanels[0]); ++i)
    {
        if (i != 0u) panels << '+';
        panels << kRepairTimingPanels[i].Name;
    }
    std::ostringstream contract_id;
    contract_id << std::hex << std::setfill('0')
        << std::setw(16) << options.Contract->ProvisionalId;
    std::ostringstream dispatch_id;
    dispatch_id << std::hex << std::setfill('0')
        << std::setw(16)
        << options.DispatchContract->ProfileId;

    std::ostringstream receipt;
    receipt << std::setprecision(17);
    receipt
        << "# repairtiming"
        << ",protocol=wirehair-v2-bench:repairtiming:repair-v1:v3"
        << ",schema=wirehair.wh2.repairtiming.v3"
        << ",repair_arm=" << options.Contract->Name
        << ",repair_contract_id=" << contract_id.str()
        << ",repair_contract_sha256="
            << options.Contract->CanonicalSha256
        << ",repair_policy_sha256="
            << wirehair_v2::kRepairV1PolicySha256
        << ",repair_attempt_cap="
            << wirehair_v2::kRepairV1AttemptCount
        << ",descriptor_magic=W2R3"
        << ",descriptor_encoding="
            "ascii-W2R3-u64le-contract-id-u32le-root-u8-attempt-v1"
        << ",descriptor_bytes=17"
        << ",descriptor_policy="
            "exact-length-two-known-ids-attempt-lt8-transactional-v1"
        << ",K=" << options.BlockCount
        << ",bb=" << options.BlockBytes
        << ",message_bytes=" << message_bytes
        << ",message_tail_bytes=1"
        << ",dispatch_profile=dispatch-v1"
        << ",dispatch_contract_id=" << dispatch_id.str()
        << ",construction_seed_base="
            << options.ConstructionSeedBase
        << ",construction_seed_derivation="
            << RepairTimingRootDerivationName(
                options.RootDerivation)
        << ",loss=" << options.Loss
        << ",loss_seed_base=" << options.LossSeedBase
        << ",loss_seed_derivation="
            "base_xor_9e3779b97f4a7c15_times_rep_plus_1_v1"
        << ",schedule=" << PacketScheduleName(options.Schedule)
        << ",loss_model=packet-schedule-v1"
        << ",trace_encoding="
            "wirehair-wh2-repairtiming-loss-trace-v1"
        << ",message_seed_policy="
            "replicate-loss-seed-partial-final-v1"
        << ",warmup_replicates="
            << options.WarmupReplicates
        << ",replicates=" << options.Replicates
        << ",cells=" << cells
        << ",summary_rows_per_cell=1"
        << ",attempt_rows_per_cell=8"
        << ",timing_panels=" << panels.str()
        << ",timing_panels_per_cell=14"
        << ",timing_slots_per_panel=8"
        << ",timing_rows_per_cell=112"
        << ",timing_order=ABBABAAB"
        << ",timing_label_policy=fixed-logical-role-v1"
        << ",inner_reps=" << options.InnerReps
        << ",max_overhead=" << options.MaxOverhead
        << ",cache_state="
            << PeelTimingCacheStateName(options.CacheState)
        << ",systematic_cache="
            << (options.SystematicCache ? "on" : "off")
        << ",evict_bytes=" << options.EvictBytes
        << ",payload_alignment=64"
        << ",prefault=1"
        << ",cpu_affinity_policy="
            "first-allowed-affinity-v1"
        << ",encoder_selector_scope="
            "lazy-selector-init-first-k-encode-w2r3-serialize-v1"
        << ",encoder_forced_scope="
            "forced-selected-init-first-k-encode-w2r3-serialize-v1"
        << ",decoder_feed_scope="
            "exact-selected-descriptor-bound-feed-through-first-success-v1"
        << ",decoder_full_scope="
            "w2r3-parse-create-feed-through-first-success-recover-v1"
        << ",direct_scope="
            "selected-attempt-pair-local-fixed-prefix-"
            "encoder-intermediate-witnessed-solve-v1"
        << ",codec_reuse=none-fresh-object-every-inner-v1"
        << ",scratch_reporting="
            "complete-precode-solve-stats-per-call-v1"
        << ",timing_stats_policy="
            "first-inner-structural-counters-summed-phase-ns-v1"
        << ",selector_timing_stats_policy="
            "all-executed-selector-subcalls-summed-per-inner-v1"
        << ",structural_digest_encoding="
            "wirehair-wh2-repair-v1-selector-structural-v1"
        << ",context_sha256=" << options.ContextSha256
        << ",uncertainty=paired-log-ratio-t95/v1"
        << ",required_margin=" << options.RequiredMargin
        << ",margin_rule="
            "upper-log-cost-lt-negative-required-margin-and-aa-floors-v1"
        << ",clock_domain=" << PeelTimingClockDomain()
        << ",stream_hash_scope=body-plus-done-prefix-v1"
        << ",started_monotonic_ns=" << started_monotonic_ns
        << ",expected_summary_rows="
            << expected_summary_rows
        << ",expected_attempt_rows="
            << expected_attempt_rows
        << ",expected_timing_rows="
            << expected_timing_rows
        << ",expected_rows=" << expected_rows << '\n';
    receipt
        << "# repairtiming_summary_columns,"
        << RepairTimingHeader(
            RepairTimingTypedSchema("summary")) << '\n'
        << "# repairtiming_attempt_columns,"
        << RepairTimingHeader(
            RepairTimingTypedSchema("attempt")) << '\n'
        << "# repairtiming_timing_columns,"
        << RepairTimingHeader(
            RepairTimingTypedSchema("timing")) << '\n';

    std::vector<uint8_t> eviction(
        (size_t)options.EvictBytes, uint8_t{0});
    EvictPeelTimingCache(eviction);
    static const char kOrder[] = {
        'A', 'B', 'B', 'A', 'B', 'A', 'A', 'B'
    };
    uint64_t summary_rows = 0u;
    uint64_t attempt_rows = 0u;
    uint64_t timing_rows = 0u;
    bool descriptor_fault_exercised = false;

    for (uint32_t cell_index = 0u;
         cell_index < cells; ++cell_index)
    {
        RepairTimingCell cell;
        if (!BuildRepairTimingCell(
                options, cell_index, cell,
                descriptor_fault_exercised))
        {
            if (!options.DescriptorFault.empty() &&
                descriptor_fault_exercised)
            {
                std::fprintf(stderr,
                    "repairtiming injected malformed W2R3 "
                    "descriptor was rejected transactionally\n");
            }
            else
            {
                std::fprintf(stderr,
                    "repairtiming semantic preflight failed "
                    "cell=%u\n", cell_index);
            }
            return 2;
        }

        receipt << MakeRepairTimingSummaryRow(options, cell);
        ++summary_rows;
        for (uint32_t attempt = 0u;
             attempt < wirehair_v2::kRepairV1AttemptCount;
             ++attempt)
        {
            receipt << MakeRepairTimingAttemptRow(
                options, cell, attempt);
            ++attempt_rows;
        }

        static const size_t kRoleCount = 9u;
        RepairTimingOperation expected[kRoleCount];
        for (size_t role_index = 0u;
             role_index < kRoleCount; ++role_index)
        {
            const RepairTimingRole role =
                (RepairTimingRole)role_index;
            if (!RunRepairTimingOperation(
                    options, cell, role, eviction,
                    false, expected[role_index]))
            {
                std::fprintf(stderr,
                    "repairtiming role preflight failed cell=%u "
                    "role=%s\n", cell_index,
                    RepairTimingRoleName(role));
                return 2;
            }
        }
        if (cell.HaveDescriptor &&
            expected[(size_t)
                RepairTimingRole::RepairSelectorEncoder].
                    SelectorStructuralSha256 !=
                cell.StructuralSha256)
        {
            std::fprintf(stderr,
                "repairtiming selector path witness changed "
                "cell=%u\n", cell_index);
            return 2;
        }

        for (uint32_t panel_index = 0u;
             panel_index <
                sizeof(kRepairTimingPanels) /
                    sizeof(kRepairTimingPanels[0]);
             ++panel_index)
        {
            const RepairTimingPanel& panel =
                kRepairTimingPanels[panel_index];
            const RepairTimingOperation& expected_a =
                expected[(size_t)panel.A];
            const RepairTimingOperation& expected_b =
                expected[(size_t)panel.B];
            const bool selector_missing =
                !cell.HaveDescriptor &&
                (
                    RepairTimingRoleRequiresSelector(panel.A) ||
                    RepairTimingRoleRequiresSelector(panel.B)
                );
            const bool preflight_eligible =
                !selector_missing &&
                expected_a.Available &&
                expected_b.Available;
            const char* initial_reason =
                selector_missing ?
                    "selector_uncommitted" :
                    (preflight_eligible ?
                        "none" : "role_preflight_failed");

            RepairTimingOperation measured[8];
            bool stable[8] = {};
            bool saturated[8] = {};
            if (preflight_eligible)
            {
                for (uint32_t slot = 0u; slot < 8u; ++slot)
                {
                    const RepairTimingRole role =
                        kOrder[slot] == 'A' ?
                            panel.A : panel.B;
                    const RepairTimingOperation& expected_role =
                        expected[(size_t)role];
                    if (!AccumulateRepairTimingMeasurement(
                            options, cell, role, expected_role,
                            eviction, measured[slot],
                            stable[slot], saturated[slot]))
                    {
                        std::fprintf(stderr,
                            "repairtiming timed operation failed "
                            "cell=%u panel=%u slot=%u\n",
                            cell_index, panel_index, slot);
                        return 2;
                    }
                }
            }

            const char* final_reason = initial_reason;
            bool timing_eligible = preflight_eligible;
            if (preflight_eligible)
            {
                bool migrated = false;
                bool faulted = false;
                bool any_saturated = false;
                bool unstable = false;
                for (uint32_t slot = 0u; slot < 8u; ++slot)
                {
                    migrated = migrated ||
                        PeelTimingCpuMigrated(
                            measured[slot].CpuBefore,
                            measured[slot].CpuAfter) != 0;
                    faulted = faulted ||
                        PeelTimingFaultContaminated(
                            measured[slot].MinorFaults,
                            measured[slot].MajorFaults) != 0;
                    any_saturated =
                        any_saturated || saturated[slot];
                    unstable = unstable || !stable[slot];
                }
                if (migrated) {
                    final_reason = "cpu_migrated";
                }
                else if (faulted) {
                    final_reason = "fault_contaminated";
                }
                else if (any_saturated) {
                    final_reason = "saturated";
                }
                else if (unstable) {
                    final_reason = "outcome_unstable";
                }
                timing_eligible =
                    !migrated && !faulted &&
                    !any_saturated && !unstable;
            }

            for (uint32_t slot = 0u; slot < 8u; ++slot)
            {
                receipt << MakeRepairTimingTimingRow(
                    options, cell, panel, panel_index, slot,
                    kOrder[slot], timing_eligible, final_reason,
                    preflight_eligible ? &measured[slot] : nullptr,
                    preflight_eligible ? stable[slot] : false,
                    preflight_eligible ?
                        saturated[slot] : false);
                ++timing_rows;
            }
        }
    }

    if (!options.DescriptorFault.empty())
    {
        std::fprintf(stderr,
            "repairtiming requested descriptor fault was not "
            "exercised\n");
        return 2;
    }
    const uint64_t emitted_rows =
        summary_rows + attempt_rows + timing_rows;
    if (summary_rows != expected_summary_rows ||
        attempt_rows != expected_attempt_rows ||
        timing_rows != expected_timing_rows ||
        emitted_rows != expected_rows)
    {
        std::fprintf(stderr,
            "repairtiming internal row count mismatch\n");
        return 2;
    }
    uint64_t finished_monotonic_ns = 0u;
    if (!ReadPeelTimingMonotonicNanoseconds(
            finished_monotonic_ns) ||
        finished_monotonic_ns <= started_monotonic_ns)
    {
        std::fprintf(stderr,
            "repairtiming could not read monotonic finish clock\n");
        return 2;
    }
    const std::string body = receipt.str();
    std::ostringstream done_prefix_stream;
    done_prefix_stream
        << "# repairtiming_done"
        << ",complete=1"
        << ",cells=" << cells
        << ",summary_rows=" << summary_rows
        << ",attempt_rows=" << attempt_rows
        << ",timing_rows=" << timing_rows
        << ",rows=" << emitted_rows
        << ",finished_monotonic_ns="
            << finished_monotonic_ns
        << ",stream_sha256=";
    const std::string done_prefix = done_prefix_stream.str();
    const std::string stream_sha256 =
        Sha256Hex(body + done_prefix);
    const std::string done =
        done_prefix + stream_sha256 + '\n';
    if (std::fwrite(
            body.data(), 1u, body.size(), stdout) != body.size() ||
        std::fwrite(
            done.data(), 1u, done.size(), stdout) != done.size() ||
        std::fflush(stdout) != 0)
    {
        std::fprintf(stderr,
            "repairtiming could not write complete receipt\n");
        return 2;
    }
    return 0;
}
