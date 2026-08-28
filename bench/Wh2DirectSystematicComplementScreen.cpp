#include "Wh2FrozenTrace.h"
#include "Wh2NativeCodec.h"
#include "Wh2NativePanel.h"

#include <algorithm>
#include <array>
#include <cerrno>
#include <chrono>
#include <csignal>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#if defined(__linux__)
#include <sched.h>
#endif

#ifndef WIREHAIR_WH2_SOURCE_GIT_COMMIT
#error "direct systematic complement screen requires an exact source commit receipt"
#endif

namespace {

using wirehair_wh2_bench::NativeArm;
using wirehair_wh2_bench::NativeArmSpec;
using wirehair_wh2_bench::NativeEncoderFixture;
using wirehair_wh2_bench::NativeEncoderStorageReceipt;
using wirehair_wh2_bench::NativePanelArm;
using wirehair_wh2_bench::NativePanelDisposition;
using wirehair_wh2_bench::NativePanelInvocation;
using wirehair_wh2_bench::NativePanelInvocationResult;
using wirehair_wh2_bench::NativePanelOrder;
using wirehair_wh2_bench::NativePanelResult;
using wirehair_wh2_bench::NativePanelSide;
using wirehair_wh2_bench::NativePanelStatus;
using wirehair_wh2_bench::NativeSystematicEmission;
using wirehair_wh2_bench::ResolvedNativeWh2Configuration;

static const char kCampaign[] =
    "wh2-retained-direct-systematic-complement-v2-r0";
static const char kConfigSchema[] =
    "wirehair.wh2.direct-systematic-complement-screen.config.v2";
static const char kPanelSchema[] =
    "wirehair.wh2.direct-systematic-complement-screen.panel.v2";
static const char kPanelKeySchema[] =
    "wirehair.wh2.direct-systematic-complement-screen.panel-key.v2";
static const char kTerminalSchema[] =
    "wirehair.wh2.direct-systematic-complement-screen.terminal.v2";
static const char kDescriptorSchema[] =
    "wirehair.wh2.direct-systematic-complement-arm.v2";
static const char kScope[] = "fresh-systematic-total-v2";
static const char kGoldenPanelKeySha256[] =
    "a0a0933ae824bacc6d11eb00d14d63efa6a55ec4e2bfe5f2f1e0d590830e1b4e";
static const char kSourceStorage[] = "native-arm-owned-kxb-v1";
static const char kSourceAcquisition[] =
    "fixture-copy-before-clock-move-v1";
static const char kTimedWork[] =
    "fresh-eager-init-plus-systematic-ids-0-through-k-minus-1-v2";
static const char kSolvedStateScope[] =
    "identity-system-rows-packet-runtime-configuration-intermediate-bytes-"
    "non-timing-solve-statistics-v1";

static const uint32_t kPanelReplicates = 12u;
static const int kTargetCpu = 120;
static const uint32_t kInvocationBudget = 65536u;
static const uint32_t kMinimumInvocations = 24u;
static const uint32_t kInternalDeadlineSeconds = 105u;
static const uint32_t kExpectedPanelCount = 396u;
static const uint32_t kExpectedRecordCount = 398u;
static const uint64_t kExpectedMeasuredInvocationCount = UINT64_C(43224);
static const uint64_t kExpectedWarmupInvocationCount = UINT64_C(792);
static const uint64_t kExpectedInvocationCount = UINT64_C(44016);
static const uint64_t kSourceSeedBase = UINT64_C(0x547000f8f0b19596);

struct CellShape
{
    uint32_t K;
    uint32_t BlockBytes;
};

// This is the exact complement of the invalid broad-r0 and valid initial
// systematic screen inside the frozen 24-cell development grid.  Its order is
// part of the JSONL roster and is duplicated in the adjudicator.
static const CellShape kCellShapes[] = {
    { 32u, 1280u },
    { 100u, 64u },
    { 100u, 1280u },
    { 1000u, 1280u },
    { 2048u, 64u },
    { 8192u, 64u },
    { 8192u, 1280u },
    { 20000u, 64u },
    { 20000u, 1280u },
    { 32768u, 64u },
    { 32768u, 1280u }
};

static const uint32_t kExpectedInvocationsPerSlot[11][12] = {
    { 171u, 171u, 171u, 171u, 171u, 171u,
      171u, 171u, 170u, 170u, 170u, 170u },
    { 55u, 55u, 55u, 55u, 55u, 55u,
      55u, 55u, 54u, 54u, 54u, 54u },
    { 55u, 55u, 55u, 55u, 55u, 55u,
      55u, 55u, 54u, 54u, 54u, 54u },
    { 6u, 6u, 6u, 6u, 6u, 6u,
      5u, 5u, 5u, 5u, 5u, 5u },
    { 3u, 3u, 3u, 3u, 3u, 3u,
      3u, 3u, 2u, 2u, 2u, 2u },
    { 2u, 2u, 2u, 2u, 2u, 2u,
      2u, 2u, 2u, 2u, 2u, 2u },
    { 2u, 2u, 2u, 2u, 2u, 2u,
      2u, 2u, 2u, 2u, 2u, 2u },
    { 2u, 2u, 2u, 2u, 2u, 2u,
      2u, 2u, 2u, 2u, 2u, 2u },
    { 2u, 2u, 2u, 2u, 2u, 2u,
      2u, 2u, 2u, 2u, 2u, 2u },
    { 2u, 2u, 2u, 2u, 2u, 2u,
      2u, 2u, 2u, 2u, 2u, 2u },
    { 2u, 2u, 2u, 2u, 2u, 2u,
      2u, 2u, 2u, 2u, 2u, 2u }
};

struct FrozenCellSeals
{
    const char* Configuration;
    const char* Source;
    const char* FullRepair;
    const char* FirstRepair;
    const char* HighId;
    const char* DirectHit;
    const char* SolvedState;
};

// Frozen from the sole deterministic receipt-generation selftest.  Each row
// is configuration, source/systematic, full repair, first repair, high-ID,
// direct-hit, and solved-state-equivalence SHA-256 in kCellShapes order.
static const FrozenCellSeals kExpectedCellSeals[] = {
    { "b283539ee89370374952a920b7f4605de7930fa22a428d8edac7d40f9f305971",
      "1e296fa19681d28489ea4c245f975e8f61c63085c9ef19897d4ff7f1b2315d71",
      "da9ba77df030059b0e68f78d544590b5b0cbf5faaa673356d63c6d391cefa062",
      "e4ec3da26bc910aab4a7f1d5143496526e7b040a019aef5bb30a98e6cbde7e11",
      "006ba6b059f6b86dc5a92d7b5c53b8dcdf29e4d6e9a83276ea4e6397606b4fcb",
      "57fe8cc94797c30786446d7d21feefd78523fe0a42f55106532530a54956ab68",
      "4cf3c08f22a4fd9055c894fd8c7d253c746c1089d1483d5a84f112c69b7b47da" },
    { "58b35a698480884fb8632a5be4ecd2b12252fed394eaa410d439258c2db61ba4",
      "a2e57558b1c8027ae2a5c4a0634903a105a0f89935bb72748c503e55a65bd655",
      "d5ce6769b52e553aabf28f17e65628250c0c6d7ead7a83e46ec2879131060325",
      "91ab0bb101199ff7143cc8ac3e6ee460fce3cd3264d1678b8824c2cf45fc3a3a",
      "072925662f69d22a448aeffb7d9ec48d05ec3a8eaf3e96b85e3f80944ccb5ad5",
      "c49f184f939634d268b9ac30f452e172a6acf6cfcd5952ee51fc9c58b0b32dfc",
      "dbb815ef7b18adb456c552044eb1e3843830dddb310667e1694de2792faf1ab4" },
    { "07cbf4890e7b9e9650c6d16db922550f562fef1ad7d85cc91ce105edd546e215",
      "17f6f7bc90fb70cd21e0968809b57aa627991388537105ce5749df3a7d11736c",
      "04c2233417e2797c59bdd22d05344d0128b3c7120d9c43f20e694e00dbf5888f",
      "2032d39d529b666847e477d34e3727606f3fa1261b2e12d5ffe0cd9caac72b97",
      "b8af8ef80d2f64a62503afcfd4aa1a390f2f8270f663a5a1a7955a001f1a4083",
      "c49f184f939634d268b9ac30f452e172a6acf6cfcd5952ee51fc9c58b0b32dfc",
      "ee95c5f56c1ac3bc5e749884e8e0639dd274a32691d48328d07cf6e81e19f2ee" },
    { "92cf189ef939a4aba39c68a701aa0378ead1a04c6714c613c366e67847bec4ed",
      "ff645fe4936b314fd4c471c626e7b2e020e533ae1206e1875d2d46924f16d7ca",
      "857cbb1384675e86abd40a7fca1560892e3f191b40c1d03f377d9b68e85f929d",
      "695f23f5df60b3bd7f04225f53c4963fe2af815aede65128701791a8661680ea",
      "48d051490a1702234fba2216ed74c249e834c60e37f5db60e800e2b334f1f2a8",
      "a6f9cb15de11c8cb55c30f134843006627d11bfa4a196886392c8a9616801139",
      "4f4a4af87c304a13810f251582c3cbceeab0cc704ac07c29da9f1cd31af2c7ac" },
    { "53b2410721b558e8705c5512c74268233cff1353220a8981ba0f89f594dec147",
      "9c250263d5b2e75f8f82d28a3abe72cc644e0347cc014caed4164f1e6ee0724a",
      "c039209bc65ab8557abace8e1a4931dcde135fb19e54c7f5f532a6bb95752be1",
      "96db9bc74ffe6ef938c84fa678361f4a5c7afa10b799e74a3ec312a386f056c5",
      "817a3c100489f2e8a838faedb71facbb4d23558c70c090b1d5676e1db3ec3e83",
      "88364e26249a7f897da55030719c17543675f4ddf040dc5239cd75fc82389795",
      "69c492fa7e3b5c3332d7b92273f06ed7928fbd0ba4cc354fc799e734f9dbb12f" },
    { "14ce58ca623b635a98c9aa802edfb0a4f24095d1cedfc9069e5fdca88b16daf2",
      "4c1a75ed90a9df9aee57aec7649656147ce8d2ea4b8cd75f0bd6bfc46c062952",
      "9aa31b1e1feb78010cf81fbb3c2127b7fdaa469612a78bb9ed4f8bc347517800",
      "7d6337d1b55f1e59715df8c825a82288178f1432361cfaec49241f15da257179",
      "b566c35130def561985c04e42607753600035d69c462e88b1bcfd5058df63192",
      "61ed57e6b1b924da3ed34a053a840beca82faf648b958f85eb849001d2ad7181",
      "34d482bff01ee09cb458bb912b490161cf1ffa723bbf3fc9cd4784de3de32862" },
    { "c3c69374cec038a6f38127707e78cce35e9293f9f8e607625a00bda608b07511",
      "4ca30eb14e3e1b516390933fc18a8e6773ee52b5265082d23ef42e9566ed673f",
      "d4917f59db9e3a57baa26edd9f0f14f7aeead9372e6d5678749efb565c5bf307",
      "71a876791a4e443544e2730e2a9a1eb361f709b7c591eff980b437e9f0ec25b1",
      "4e3e79768942bfa306d5198840c44bd90e7a923dbdf2575bd4bc48d4ff0d6ae8",
      "61ed57e6b1b924da3ed34a053a840beca82faf648b958f85eb849001d2ad7181",
      "29448aec3574f98177cddbc61d4a8a3bb754359ffb0ea7faf351268fd34387c8" },
    { "f877ec4bb14d2c19023620c58d7438c0b8c6ffeebe1a71671eef4e980fc93a79",
      "3db8c777ae42ad27e188d074cd30ce5f993f1a9e682f33d08d549f39d6b76d5b",
      "d701db59cab3dae48e8bc5ea33ec11b9d3afcf62b4266ef2e5b29e370a79c4d3",
      "e64980ab7053b4b2acd93b2f33b5cba875429144c58b30c3e7642ce1e6a4cb64",
      "656b06de7ef753f0999d89f176b298510592a9574f6c0969029fb02a6675d0a6",
      "4486b30e27d851608c266d583e55ebf9ccbabc706638a2a214c12499fb3074e7",
      "92622eb136713cc4330ac7304c71e02783bcc0cf7079dcba8e66050c481d605e" },
    { "ccbfc46c15380255317dda067e6a8cb814c43fe2f28a6b347a804482b37cfb68",
      "0b7004dd9c291f91a1aef3219c77862ca35d792ec2a6fc276695f3e47044d3a4",
      "502588e72e1ac409690d4897a063ec81e5e0a45bcd3a32da47432e0cd27a4d7e",
      "dc8c2ed8ed9a8947ef73376222c23eba0af5da221748218bddbe22292898d169",
      "64d6afc42f2c61acee4384143191e5d432f5f7f3ee5967bebcee4aa06bc49923",
      "4486b30e27d851608c266d583e55ebf9ccbabc706638a2a214c12499fb3074e7",
      "1c32565b969f3e443d5ac1943fcf3f3dffbd1b1a24398f1462be244227bafd81" },
    { "e361af19658f217f663436822d1df28b42fc7644622b2b384e5081934b710626",
      "3628fea7aada74f0dedf60f34ffcc7e76cb634b23d9f28738a1ab9a72461a15e",
      "80720243cc518ea1b559c9f23ea46b6fc87dd9a4f276350bcfca72f1d929d982",
      "bb73eca6a8929f39be371badaa3485fd45efbdb0303c2a94c9bdf16a1a3493ac",
      "c7989e346ba0a3bb693ec5f909bcdc93545a3c43893c17af740d4923ee8fca45",
      "4338dea9a6928832fbfb1e45e23f51313a13b4ad4bde7ccf87b9a9e0c3737565",
      "5ed759119604192b2f8f038825f09783113b362218f6e59ac6c9389923af53af" },
    { "4d9978d65eacfd17cb7b17d5dbc8cae28c0bfd989e4aed3f5f4ec9c24ee85fed",
      "13aaa0e514bcdf73cc24645cfa67155df23fe645d8ba4161601779668f93f955",
      "9cfec479624a860ecc6d3feffdb4f73f47dfac77a84ffce7287032b84cfa3769",
      "40d3fe7e775ad2993747cafe0fcd8663f5ba91cf92db129355305533790f1903",
      "9e6b0375d9b4c8d9ffa9c48e4d957962f5239a4916f718754add9a39d2c7594e",
      "4338dea9a6928832fbfb1e45e23f51313a13b4ad4bde7ccf87b9a9e0c3737565",
      "5e00fc2f973da725376fb6d606e4f04c0e467e692483f245cf1f7f6618788f49" }
};

struct ScreenCell
{
    CellShape Shape = {};
    uint64_t SourceSeed = 0u;
    std::vector<uint8_t> Source;
    std::vector<uint8_t> ExpectedRepairs;
    std::string EquationConfigurationJson;
    std::string EquationConfigurationSha256;
    std::string SourceSha256;
    std::string SystematicOracleSha256;
    std::string FullRepairSha256;
    std::string FirstRepairSha256;
    std::string HighIdOracleSha256;
    std::string DirectHitOracleSha256;
    std::string SolvedStateEquivalenceReceiptSha256;
    bool ConstructionEquivalent = false;
    bool SystematicOracleVerified = false;
    bool FullRepairOracleVerified = false;
    bool HighIdOracleVerified = false;
    bool DirectHitOracleVerified = false;
};

struct Comparison
{
    const char* Name;
    NativeSystematicEmission Left;
    NativeSystematicEmission Right;
};

static const Comparison kComparisons[] = {
    { "baseline-aa",
      NativeSystematicEmission::EquationEvaluation,
      NativeSystematicEmission::EquationEvaluation },
    { "candidate-aa",
      NativeSystematicEmission::RetainedSourceDirect,
      NativeSystematicEmission::RetainedSourceDirect },
    { "candidate-over-baseline",
      NativeSystematicEmission::RetainedSourceDirect,
      NativeSystematicEmission::EquationEvaluation }
};

const char* SideName(NativePanelSide side)
{
    return side == NativePanelSide::Left ? "left" : "right";
}

const char* OrderName(NativePanelOrder order)
{
    return order == NativePanelOrder::ABBA ? "ABBA" : "BAAB";
}

const char* EmissionMode(NativeSystematicEmission emission)
{
    if (emission == NativeSystematicEmission::EquationEvaluation) {
        return "equation";
    }
    if (emission == NativeSystematicEmission::RetainedSourceDirect) {
        return "retained";
    }
    return nullptr;
}

bool ParseTargetCpu(const char* text, int& cpu_out)
{
    if (!text || !*text) return false;
    errno = 0;
    char* end = nullptr;
    const long value = std::strtol(text, &end, 10);
    if (errno != 0 || !end || *end != '\0' || value != kTargetCpu)
    {
        return false;
    }
    cpu_out = static_cast<int>(value);
    return true;
}

#if defined(__linux__)
bool IsExactSingletonCpuSet(const cpu_set_t& set, int target_cpu)
{
    if (target_cpu < 0 || target_cpu >= CPU_SETSIZE ||
        !CPU_ISSET(target_cpu, &set))
    {
        return false;
    }
    unsigned selected_count = 0u;
    for (int cpu = 0; cpu < CPU_SETSIZE; ++cpu) {
        selected_count += CPU_ISSET(cpu, &set) ? 1u : 0u;
    }
    return selected_count == 1u;
}
#endif

bool TestSingletonCpuSetMatcher()
{
#if defined(__linux__)
    static const int first_cpu = 0;
    static const int second_cpu = 1;
    if (CPU_SETSIZE <= second_cpu) return false;
    cpu_set_t set;
    CPU_ZERO(&set);
    if (IsExactSingletonCpuSet(set, first_cpu)) return false;
    CPU_SET(first_cpu, &set);
    if (!IsExactSingletonCpuSet(set, first_cpu) ||
        IsExactSingletonCpuSet(set, second_cpu))
    {
        return false;
    }
    CPU_SET(second_cpu, &set);
    if (IsExactSingletonCpuSet(set, first_cpu) ||
        IsExactSingletonCpuSet(set, second_cpu))
    {
        return false;
    }
    CPU_CLR(first_cpu, &set);
    return IsExactSingletonCpuSet(set, second_cpu) &&
        !IsExactSingletonCpuSet(set, -1) &&
        !IsExactSingletonCpuSet(set, CPU_SETSIZE);
#else
    return true;
#endif
}

bool PinAndVerifyTargetCpu(int target_cpu, std::string& diagnostic)
{
    diagnostic.clear();
#if defined(__linux__)
    if (target_cpu != kTargetCpu || target_cpu < 0 ||
        target_cpu >= CPU_SETSIZE)
    {
        diagnostic = "target CPU is outside the frozen singleton domain";
        return false;
    }
    cpu_set_t requested;
    CPU_ZERO(&requested);
    CPU_SET(target_cpu, &requested);
    errno = 0;
    if (sched_setaffinity(0, sizeof(requested), &requested) != 0)
    {
        const int error = errno;
        diagnostic = std::string("sched_setaffinity failed: ") +
            std::strerror(error);
        return false;
    }
    cpu_set_t observed;
    CPU_ZERO(&observed);
    errno = 0;
    if (sched_getaffinity(0, sizeof(observed), &observed) != 0)
    {
        const int error = errno;
        diagnostic = std::string("sched_getaffinity failed: ") +
            std::strerror(error);
        return false;
    }
    if (!IsExactSingletonCpuSet(observed, target_cpu))
    {
        diagnostic = "worker affinity is not the frozen singleton CPU";
        return false;
    }
    errno = 0;
    const int observed_cpu = sched_getcpu();
    if (observed_cpu < 0)
    {
        const int error = errno;
        diagnostic = std::string("sched_getcpu failed: ") +
            std::strerror(error);
        return false;
    }
    if (observed_cpu != target_cpu)
    {
        std::ostringstream message;
        message << "worker is running on CPU " << observed_cpu
                << " instead of frozen CPU " << target_cpu;
        diagnostic = message.str();
        return false;
    }
    return true;
#else
    (void)target_cpu;
    diagnostic = "direct systematic complement requires Linux affinity";
    return false;
#endif
}

uint64_t SourceSeed(const CellShape& shape)
{
    return kSourceSeedBase ^
        (static_cast<uint64_t>(shape.K) << 17) ^ shape.BlockBytes;
}

uint32_t TotalInvocations(uint32_t K)
{
    const uint32_t budget = (kInvocationBudget + K - 1u) / K;
    return std::max(kMinimumInvocations, budget);
}

uint32_t InvocationsPerSlot(uint32_t K, uint32_t replicate)
{
    const uint32_t total = TotalInvocations(K);
    const uint32_t quotient = total / kPanelReplicates;
    const uint32_t remainder = total % kPanelReplicates;
    return quotient + (replicate < remainder ? 1u : 0u);
}

std::string PanelKeySha256(
    const CellShape& shape,
    const char* comparison)
{
    if (!comparison || !*comparison) return std::string();
    std::ostringstream json;
    json << "{\"K\":" << shape.K
         << ",\"block_bytes\":" << shape.BlockBytes
         << ",\"campaign\":\"" << kCampaign
         << "\",\"comparison\":\"" << comparison
         << "\",\"schema\":\"" << kPanelKeySchema
         << "\",\"scope\":\"" << kScope << "\"}";
    return wirehair::wh2_benchmark::Sha256Hex(json.str());
}

NativePanelOrder PanelOrder(
    const CellShape& shape,
    const char* comparison,
    uint32_t replicate)
{
    const std::string digest = PanelKeySha256(shape, comparison);
    if (digest.size() != 64u) {
        return static_cast<NativePanelOrder>(UINT32_MAX);
    }
    const char last = digest.back();
    const uint32_t nibble = last >= '0' && last <= '9' ?
        static_cast<uint32_t>(last - '0') :
        static_cast<uint32_t>(last - 'a' + 10);
    return ((replicate & 1u) == (nibble & 1u)) ?
        NativePanelOrder::ABBA : NativePanelOrder::BAAB;
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
        left.DenseAnchors == right.DenseAnchors &&
        left.HeavyFamily == right.HeavyFamily &&
        left.Seed == right.Seed;
}

bool SameResolved(
    const ResolvedNativeWh2Configuration& left,
    const ResolvedNativeWh2Configuration& right)
{
    return SamePrecodeParams(left.Params, right.Params) &&
        left.PacketConfig.PeelSeed == right.PacketConfig.PeelSeed &&
        left.PacketConfig.MixCount == right.PacketConfig.MixCount &&
        left.PrecodeAttempt == right.PrecodeAttempt &&
        left.PacketAttempt == right.PacketAttempt;
}

NativeArmSpec ArmSpec(NativeSystematicEmission emission)
{
    NativeArmSpec spec = wirehair_wh2_bench::MakeCertifiedWh2Arm(0u);
    spec.SystematicEmission = emission;
    return spec;
}

std::string DescriptorJson(NativeSystematicEmission emission)
{
    const char* const emission_name =
        wirehair_wh2_bench::NativeSystematicEmissionName(emission);
    if (!emission_name || !EmissionMode(emission)) return std::string();
    std::string json;
    json.reserve(420u);
    json += "{\"arm\":\"wirehair2_head\",\"codec\":\"wirehair2_certified";
    json += "\",\"equation_transform\":\"none\",\"metric\":\"";
    json += kScope;
    json += "\",\"schema\":\"";
    json += kDescriptorSchema;
    json += "\",\"source_acquisition\":\"";
    json += kSourceAcquisition;
    json += "\",\"source_storage\":\"";
    json += kSourceStorage;
    json += "\",\"systematic_emission\":\"";
    json += emission_name;
    json += "\",\"timed_work\":\"";
    json += kTimedWork;
    json += "\"}";
    return json;
}

std::string ConfigurationJson(
    const CellShape& shape,
    const ResolvedNativeWh2Configuration& resolved)
{
    const wirehair_v2::PrecodeParams& params = resolved.Params;
    std::ostringstream json;
    json << "{\"K\":" << shape.K
         << ",\"block_bytes\":" << shape.BlockBytes
         << ",\"dense_anchor_layout\":"
         << static_cast<uint32_t>(params.DenseAnchors)
         << ",\"dense_identity_corner\":"
         << (params.DenseIdentityCorner ? "true" : "false")
         << ",\"dense_rows\":" << params.DenseRows
         << ",\"heavy_family\":"
         << static_cast<uint32_t>(params.HeavyFamily)
         << ",\"heavy_rows\":" << params.HeavyRows
         << ",\"mix_count\":" << resolved.PacketConfig.MixCount
         << ",\"packet_attempt\":" << resolved.PacketAttempt
         << ",\"packet_peel_seed\":" << resolved.PacketConfig.PeelSeed
         << ",\"precode_attempt\":" << resolved.PrecodeAttempt
         << ",\"precode_seed\":" << params.Seed
         << ",\"source_hits\":" << params.SourceHits
         << ",\"staircase\":" << params.Staircase << "}";
    return json.str();
}

std::string DirectHitReceiptJson(const CellShape& shape)
{
    std::ostringstream json;
    json << "{\"equation_direct_systematic_packets\":0"
         << ",\"retained_direct_systematic_packets\":" << shape.K
         << ",\"verified_ids\":[0," << shape.K - 1u << ',' << shape.K
         << ',' << UINT32_MAX << "]}";
    return json.str();
}

std::string SolvedStateEquivalenceReceiptJson(const ScreenCell& cell)
{
    std::string json;
    json.reserve(760u);
    json += "{\"configuration_sha256\":\"";
    json += cell.EquationConfigurationSha256;
    json += "\",\"direct_hit_oracle_sha256\":\"";
    json += cell.DirectHitOracleSha256;
    json += "\",\"first_repair_sha256\":\"";
    json += cell.FirstRepairSha256;
    json += "\",\"full_repair_sha256\":\"";
    json += cell.FullRepairSha256;
    json += "\",\"high_id_oracle_sha256\":\"";
    json += cell.HighIdOracleSha256;
    json += "\",\"source_sha256\":\"";
    json += cell.SourceSha256;
    json += "\",\"systematic_oracle_sha256\":\"";
    json += cell.SystematicOracleSha256;
    json += "\",\"verified_state_scope\":\"";
    json += kSolvedStateScope;
    json += "\"}";
    return json;
}

bool ValidStorageReceipt(
    const NativeEncoderStorageReceipt& receipt,
    NativeSystematicEmission expected)
{
    return receipt.SourceStorage && receipt.SourceAcquisition &&
        std::strcmp(receipt.SourceStorage, kSourceStorage) == 0 &&
        std::strcmp(receipt.SourceAcquisition, kSourceAcquisition) == 0 &&
        receipt.SystematicEmission == expected;
}

bool ExactPolicyRoster(
    const NativeArm& equation_arm,
    const NativeArm& retained_arm,
    uint32_t K)
{
    return !equation_arm.UsesRetainedSourceFor(0u) &&
        !equation_arm.UsesRetainedSourceFor(K - 1u) &&
        !equation_arm.UsesRetainedSourceFor(K) &&
        !equation_arm.UsesRetainedSourceFor(UINT32_MAX) &&
        retained_arm.UsesRetainedSourceFor(0u) &&
        retained_arm.UsesRetainedSourceFor(K - 1u) &&
        !retained_arm.UsesRetainedSourceFor(K) &&
        !retained_arm.UsesRetainedSourceFor(UINT32_MAX);
}

bool BuildCell(const CellShape& shape, ScreenCell& cell_out)
{
    ScreenCell cell;
    cell.Shape = shape;
    cell.SourceSeed = SourceSeed(shape);
    if (!wirehair_wh2_bench::MakeDeterministicSource(
            shape.K, shape.BlockBytes, cell.SourceSeed, cell.Source))
    {
        return false;
    }

    const NativeArmSpec equation =
        ArmSpec(NativeSystematicEmission::EquationEvaluation);
    const NativeArmSpec retained =
        ArmSpec(NativeSystematicEmission::RetainedSourceDirect);
    ResolvedNativeWh2Configuration equation_config;
    ResolvedNativeWh2Configuration retained_config;
    if (!wirehair_wh2_bench::ResolveNativeWh2Configuration(
            equation, shape.K, shape.BlockBytes, equation_config) ||
        !wirehair_wh2_bench::ResolveNativeWh2Configuration(
            retained, shape.K, shape.BlockBytes, retained_config) ||
        !SameResolved(equation_config, retained_config))
    {
        return false;
    }
    cell.EquationConfigurationJson =
        ConfigurationJson(shape, equation_config);
    cell.EquationConfigurationSha256 =
        wirehair::wh2_benchmark::Sha256Hex(
            cell.EquationConfigurationJson);
    cell.SourceSha256 = wirehair::wh2_benchmark::Sha256Hex(
        cell.Source.data(), cell.Source.size());
    cell.SystematicOracleSha256 = cell.SourceSha256;

    NativeArm equation_arm;
    NativeArm retained_arm;
    if (equation_arm.Initialize(
            equation, shape.K, shape.BlockBytes, cell.Source) !=
            Wirehair_Success ||
        retained_arm.Initialize(
            retained, shape.K, shape.BlockBytes, cell.Source) !=
            Wirehair_Success ||
        !equation_arm.HasSameWh2SolvedState(retained_arm))
    {
        return false;
    }
    cell.ConstructionEquivalent = true;
    if (!ExactPolicyRoster(equation_arm, retained_arm, shape.K)) {
        return false;
    }
    cell.DirectHitOracleSha256 = wirehair::wh2_benchmark::Sha256Hex(
        DirectHitReceiptJson(shape));
    cell.DirectHitOracleVerified =
        cell.DirectHitOracleSha256.size() == 64u;

    std::vector<uint8_t> equation_packet(shape.BlockBytes, 0u);
    std::vector<uint8_t> retained_packet(shape.BlockBytes, 0u);
    for (uint32_t id = 0u; id < shape.K; ++id)
    {
        if (equation_arm.EncodeInto(
                id, equation_packet.data(), shape.BlockBytes) !=
                Wirehair_Success ||
            retained_arm.EncodeInto(
                id, retained_packet.data(), shape.BlockBytes) !=
                Wirehair_Success ||
            equation_packet != retained_packet ||
            std::memcmp(
                equation_packet.data(),
                cell.Source.data() +
                    static_cast<size_t>(id) * shape.BlockBytes,
                shape.BlockBytes) != 0)
        {
            return false;
        }
    }
    cell.SystematicOracleVerified = true;

    try {
        cell.ExpectedRepairs.assign(
            static_cast<size_t>(shape.K) * shape.BlockBytes, 0u);
    }
    catch (...) {
        return false;
    }
    for (uint32_t offset = 0u; offset < shape.K; ++offset)
    {
        uint8_t* const expected = cell.ExpectedRepairs.data() +
            static_cast<size_t>(offset) * shape.BlockBytes;
        if (equation_arm.EncodeInto(
                shape.K + offset, expected, shape.BlockBytes) !=
                Wirehair_Success ||
            retained_arm.EncodeInto(
                shape.K + offset,
                retained_packet.data(),
                shape.BlockBytes) != Wirehair_Success ||
            std::memcmp(
                expected, retained_packet.data(), shape.BlockBytes) != 0)
        {
            return false;
        }
    }
    cell.FullRepairSha256 = wirehair::wh2_benchmark::Sha256Hex(
        cell.ExpectedRepairs.data(), cell.ExpectedRepairs.size());
    cell.FirstRepairSha256 = wirehair::wh2_benchmark::Sha256Hex(
        cell.ExpectedRepairs.data(), shape.BlockBytes);
    cell.FullRepairOracleVerified = true;

    std::vector<uint8_t> high_outputs;
    try {
        high_outputs.assign(static_cast<size_t>(3u) * shape.BlockBytes, 0u);
    }
    catch (...) {
        return false;
    }
    const uint32_t high_ids[] = {
        shape.K + 7u, UINT32_C(0x80000000), UINT32_MAX
    };
    for (size_t index = 0u; index < 3u; ++index)
    {
        uint8_t* const output = high_outputs.data() +
            index * shape.BlockBytes;
        if (equation_arm.EncodeInto(
                high_ids[index], output, shape.BlockBytes) !=
                Wirehair_Success ||
            retained_arm.EncodeInto(
                high_ids[index], retained_packet.data(), shape.BlockBytes) !=
                Wirehair_Success ||
            std::memcmp(
                output, retained_packet.data(), shape.BlockBytes) != 0)
        {
            return false;
        }
    }
    cell.HighIdOracleSha256 = wirehair::wh2_benchmark::Sha256Hex(
        high_outputs.data(), high_outputs.size());
    cell.HighIdOracleVerified = true;

    cell.SolvedStateEquivalenceReceiptSha256 =
        wirehair::wh2_benchmark::Sha256Hex(
            SolvedStateEquivalenceReceiptJson(cell));
    const bool hashes_valid =
        cell.EquationConfigurationSha256.size() == 64u &&
        cell.SourceSha256.size() == 64u &&
        cell.SystematicOracleSha256.size() == 64u &&
        cell.FullRepairSha256.size() == 64u &&
        cell.FirstRepairSha256.size() == 64u &&
        cell.HighIdOracleSha256.size() == 64u &&
        cell.DirectHitOracleSha256.size() == 64u &&
        cell.SolvedStateEquivalenceReceiptSha256.size() == 64u;
    if (!hashes_valid || !cell.ConstructionEquivalent ||
        !cell.SystematicOracleVerified ||
        !cell.FullRepairOracleVerified ||
        !cell.HighIdOracleVerified ||
        !cell.DirectHitOracleVerified)
    {
        return false;
    }
    // Only the digest is evidence for this fresh-systematic worker.  Release
    // the large untimed repair oracle before cells are retained for timing so
    // it cannot perturb the campaign's resident set or cache state.
    std::vector<uint8_t>().swap(cell.ExpectedRepairs);
    cell_out = std::move(cell);
    return true;
}

class FreshSystematicInvocation : public NativePanelInvocation
{
public:
    FreshSystematicInvocation(
        const std::shared_ptr<const ScreenCell>& cell,
        NativeSystematicEmission emission,
        const std::string& identity)
        : Cell(cell)
        , Emission(emission)
        , IdentityValue(identity)
    {
        const NativeArmSpec spec = ArmSpec(emission);
        Ready = Cell && Fixture.Initialize(
            spec,
            Cell->Shape.K,
            Cell->Shape.BlockBytes,
            Cell->Source) == Wirehair_Success;
        NativeEncoderStorageReceipt receipt;
        Ready = Ready && Fixture.GetStorageReceipt(receipt) &&
            ValidStorageReceipt(receipt, emission);
    }

    std::string Identity() const override { return IdentityValue; }

    NativePanelInvocationResult Invoke() override
    {
        if (!Ready) {
            return NativePanelInvocationResult(
                NativePanelDisposition::Fatal,
                Wirehair_Error,
                false,
                0u,
                0u);
        }
        const wirehair_wh2_bench::TimedArmResult timed = Fixture.Run();
        const uint32_t expected_direct =
            Emission == NativeSystematicEmission::RetainedSourceDirect ?
                Cell->Shape.K : 0u;
        if (timed.Result != Wirehair_Success || !timed.BytesVerified ||
            timed.ElapsedNanoseconds == 0u ||
            timed.DirectSystematicPackets != expected_direct)
        {
            return NativePanelInvocationResult(
                NativePanelDisposition::PreflightFailure,
                timed.Result,
                true,
                0u,
                0u);
        }
        return NativePanelInvocationResult(
            NativePanelDisposition::Success,
            timed.Result,
            true,
            expected_direct,
            timed.ElapsedNanoseconds);
    }

private:
    std::shared_ptr<const ScreenCell> Cell;
    NativeSystematicEmission Emission;
    std::string IdentityValue;
    NativeEncoderFixture Fixture;
    bool Ready = false;
};

struct DescriptorSeal
{
    NativeSystematicEmission Emission;
    std::string Json;
    std::string Sha256;
};

const DescriptorSeal* FindDescriptor(
    const std::vector<DescriptorSeal>& seals,
    NativeSystematicEmission emission)
{
    for (const DescriptorSeal& seal : seals) {
        if (seal.Emission == emission) return &seal;
    }
    return nullptr;
}

bool BuildDescriptorSeals(std::vector<DescriptorSeal>& seals_out)
{
    std::vector<DescriptorSeal> seals;
    for (NativeSystematicEmission emission : {
             NativeSystematicEmission::EquationEvaluation,
             NativeSystematicEmission::RetainedSourceDirect })
    {
        DescriptorSeal seal;
        seal.Emission = emission;
        seal.Json = DescriptorJson(emission);
        seal.Sha256 = wirehair::wh2_benchmark::Sha256Hex(seal.Json);
        if (seal.Json.empty() || seal.Sha256.size() != 64u) return false;
        seals.push_back(std::move(seal));
    }
    seals_out.swap(seals);
    return true;
}

NativePanelArm MakePanelArm(
    const std::shared_ptr<const ScreenCell>& cell,
    NativeSystematicEmission emission,
    const std::string& descriptor_sha256)
{
    if (!cell || descriptor_sha256.size() != 64u ||
        !EmissionMode(emission))
    {
        return NativePanelArm();
    }
    std::ostringstream identity;
    identity << kScope << ':' << cell->Shape.K << ':'
             << cell->Shape.BlockBytes << ':' << descriptor_sha256;
    const std::string identity_value = identity.str();
    return NativePanelArm(
        identity_value,
        [cell, emission, identity_value]() {
            return std::unique_ptr<NativePanelInvocation>(
                new FreshSystematicInvocation(
                    cell, emission, identity_value));
        });
}

bool EmitConfig(
    std::ostream& output,
    int target_cpu,
    const std::vector<std::shared_ptr<const ScreenCell> >& cells,
    const std::vector<DescriptorSeal>& descriptors)
{
    output << "{\"campaign\":\"" << kCampaign << "\",\"cells\":[";
    for (size_t i = 0u; i < cells.size(); ++i)
    {
        if (i != 0u) output << ',';
        const ScreenCell& cell = *cells[i];
        output << "{\"K\":" << cell.Shape.K
                  << ",\"block_bytes\":" << cell.Shape.BlockBytes
                  << ",\"construction_equivalent\":"
                  << (cell.ConstructionEquivalent ? "true" : "false")
                  << ",\"direct_hit_oracle_sha256\":\""
                  << cell.DirectHitOracleSha256
                  << "\",\"direct_hit_oracle_verified\":"
                  << (cell.DirectHitOracleVerified ? "true" : "false")
                  << ",\"equation_configuration\":"
                  << cell.EquationConfigurationJson
                  << ",\"equation_configuration_sha256\":\""
                  << cell.EquationConfigurationSha256
                  << "\",\"first_repair_sha256\":\""
                  << cell.FirstRepairSha256
                  << "\",\"full_repair_oracle_verified\":"
                  << (cell.FullRepairOracleVerified ? "true" : "false")
                  << ",\"full_repair_sha256\":\""
                  << cell.FullRepairSha256
                  << "\",\"high_id_oracle_sha256\":\""
                  << cell.HighIdOracleSha256
                  << "\",\"high_id_oracle_verified\":"
                  << (cell.HighIdOracleVerified ? "true" : "false")
                  << ",\"invocations_by_replicate\":[";
        for (uint32_t replicate = 0u;
             replicate < kPanelReplicates;
             ++replicate)
        {
            if (replicate != 0u) output << ',';
            output << InvocationsPerSlot(cell.Shape.K, replicate);
        }
        output << ']'
                  << ",\"solved_state_equivalence_receipt_sha256\":\""
                  << cell.SolvedStateEquivalenceReceiptSha256
                  << "\",\"solved_state_equivalent\":true"
                  << ",\"solved_state_scope\":\"" << kSolvedStateScope
                  << "\""
                  << ",\"source_seed\":\"0x" << std::hex
                  << cell.SourceSeed << std::dec
                  << "\",\"source_sha256\":\"" << cell.SourceSha256
                  << "\",\"systematic_oracle_sha256\":\""
                  << cell.SystematicOracleSha256
                  << "\",\"systematic_oracle_verified\":"
                  << (cell.SystematicOracleVerified ? "true" : "false")
                  << "}";
    }
    output << "],\"comparisons\":[";
    for (size_t i = 0u;
         i < sizeof(kComparisons) / sizeof(kComparisons[0]);
         ++i)
    {
        if (i != 0u) output << ',';
        output << "{\"left_mode\":\""
                  << EmissionMode(kComparisons[i].Left)
                  << "\",\"name\":\"" << kComparisons[i].Name
                  << "\",\"right_mode\":\""
                  << EmissionMode(kComparisons[i].Right) << "\"}";
    }
    output << "],\"descriptors\":[";
    for (size_t i = 0u; i < descriptors.size(); ++i)
    {
        if (i != 0u) output << ',';
        output << "{\"descriptor\":" << descriptors[i].Json
                  << ",\"descriptor_sha256\":\""
                  << descriptors[i].Sha256
                  << "\",\"mode\":\""
                  << EmissionMode(descriptors[i].Emission) << "\"}";
    }
    output << "],\"expected_invocations\":"
              << kExpectedInvocationCount
              << ",\"expected_measured_invocations\":"
              << kExpectedMeasuredInvocationCount
              << ",\"expected_panels\":" << kExpectedPanelCount
              << ",\"expected_records\":" << kExpectedRecordCount
              << ",\"expected_warmup_invocations\":"
              << kExpectedWarmupInvocationCount
              << ",\"internal_deadline_seconds\":"
              << kInternalDeadlineSeconds
              << ",\"invocation_budget\":" << kInvocationBudget
              << ",\"minimum_invocations\":" << kMinimumInvocations
              << ",\"panel_key_schema\":\"" << kPanelKeySchema
              << "\",\"panel_replicates\":" << kPanelReplicates
              << ",\"schema\":\"" << kConfigSchema
              << "\",\"scope\":\"" << kScope
              << "\",\"source_git_commit\":\""
              << WIREHAIR_WH2_SOURCE_GIT_COMMIT
              << "\",\"source_seed_base\":\"0x" << std::hex
              << kSourceSeedBase << std::dec
              << "\",\"target_cpu\":" << target_cpu << "}\n";
    output.flush();
    return output.good();
}

bool EmitPanel(
    const ScreenCell& cell,
    const char* comparison,
    uint32_t invocations_per_slot,
    const std::string& left_sha256,
    NativePanelOrder order,
    const std::string& panel_key_sha256,
    uint32_t replicate,
    const std::string& right_sha256,
    const NativePanelResult& panel)
{
    std::cout << "{\"K\":" << cell.Shape.K
              << ",\"block_bytes\":" << cell.Shape.BlockBytes
              << ",\"comparison\":\"" << comparison
              << "\",\"invocations_per_slot\":" << invocations_per_slot
              << ",\"left_descriptor_sha256\":\"" << left_sha256
              << "\",\"left_direct_systematic_packets\":"
              << panel.LeftPreflight.DecodedExtra
              << ",\"left_outcome_code\":"
              << panel.LeftPreflight.OutcomeCode
              << ",\"order\":\"" << OrderName(order)
              << "\",\"panel_key_sha256\":\"" << panel_key_sha256
              << "\",\"primary_metric\":\"" << kScope
              << "\",\"replicate\":" << replicate
              << ",\"right_descriptor_sha256\":\"" << right_sha256
              << "\",\"right_direct_systematic_packets\":"
              << panel.RightPreflight.DecodedExtra
              << ",\"right_outcome_code\":"
              << panel.RightPreflight.OutcomeCode
              << ",\"schema\":\"" << kPanelSchema
              << "\",\"scope\":\"" << kScope << "\",\"slots\":[";
    for (size_t i = 0u; i < panel.Slots.size(); ++i)
    {
        if (i != 0u) std::cout << ',';
        const uint32_t slot_invocations = i < 4u ?
            invocations_per_slot / 2u + invocations_per_slot % 2u :
            invocations_per_slot / 2u;
        std::cout << "{\"elapsed_ns\":"
                  << panel.Slots[i].ElapsedNanoseconds
                  << ",\"invocation_count\":" << slot_invocations
                  << ",\"side\":\"" << SideName(panel.Slots[i].Side)
                  << "\"}";
    }
    std::cout << "],\"target_cpu\":" << panel.TargetCpu << "}\n";
    std::cout.flush();
    return std::cout.good();
}

bool CellSealsMatch(size_t index, const ScreenCell& cell)
{
    if (index >= sizeof(kExpectedCellSeals) /
                     sizeof(kExpectedCellSeals[0]))
    {
        return false;
    }
    const FrozenCellSeals& expected = kExpectedCellSeals[index];
    return expected.Configuration[0] != '\0' &&
        cell.EquationConfigurationSha256 == expected.Configuration &&
        cell.SourceSha256 == expected.Source &&
        cell.FullRepairSha256 == expected.FullRepair &&
        cell.FirstRepairSha256 == expected.FirstRepair &&
        cell.HighIdOracleSha256 == expected.HighId &&
        cell.DirectHitOracleSha256 == expected.DirectHit &&
        cell.SolvedStateEquivalenceReceiptSha256 == expected.SolvedState;
}

void PrintCellSeals(const ScreenCell& cell)
{
    std::cerr << "{ \"" << cell.EquationConfigurationSha256
              << "\", \"" << cell.SourceSha256
              << "\", \"" << cell.FullRepairSha256
              << "\", \"" << cell.FirstRepairSha256
              << "\", \"" << cell.HighIdOracleSha256
              << "\", \"" << cell.DirectHitOracleSha256
              << "\", \""
              << cell.SolvedStateEquivalenceReceiptSha256 << "\" }, // K="
              << cell.Shape.K << " B=" << cell.Shape.BlockBytes << '\n';
}

bool TestRosterArithmetic()
{
    uint64_t invocation_count = 0u;
    uint64_t measured_invocation_count = 0u;
    uint64_t warmup_invocation_count = 0u;
    uint32_t panel_count = 0u;
    std::set<std::string> panel_keys;
    for (size_t cell_index = 0u;
         cell_index < sizeof(kCellShapes) / sizeof(kCellShapes[0]);
         ++cell_index)
    {
        const CellShape& shape = kCellShapes[cell_index];
        uint32_t total = 0u;
        for (uint32_t replicate = 0u;
             replicate < kPanelReplicates;
             ++replicate)
        {
            const uint32_t n = InvocationsPerSlot(shape.K, replicate);
            if (n < 2u ||
                n != kExpectedInvocationsPerSlot[cell_index][replicate])
            {
                return false;
            }
            total += n;
            for (const Comparison& comparison : kComparisons)
            {
                const NativePanelOrder order = PanelOrder(
                    shape, comparison.Name, replicate);
                if (order != NativePanelOrder::ABBA &&
                    order != NativePanelOrder::BAAB)
                {
                    return false;
                }
                measured_invocation_count += UINT64_C(4) * n;
                warmup_invocation_count += 2u;
                invocation_count += UINT64_C(4) * n + 2u;
                ++panel_count;
            }
        }
        if (total != TotalInvocations(shape.K)) return false;
        for (const Comparison& comparison : kComparisons)
        {
            const std::string key = PanelKeySha256(shape, comparison.Name);
            if (key.size() != 64u || !panel_keys.insert(key).second) {
                return false;
            }
            uint32_t abba = 0u;
            uint32_t baab = 0u;
            for (uint32_t replicate = 0u;
                 replicate < kPanelReplicates;
                 ++replicate)
            {
                if (PanelOrder(shape, comparison.Name, replicate) ==
                    NativePanelOrder::ABBA)
                {
                    ++abba;
                }
                else {
                    ++baab;
                }
            }
            if (abba != 6u || baab != 6u) return false;
        }
    }
    return panel_keys.size() == 33u &&
        panel_count == kExpectedPanelCount &&
        measured_invocation_count == kExpectedMeasuredInvocationCount &&
        warmup_invocation_count == kExpectedWarmupInvocationCount &&
        invocation_count == kExpectedInvocationCount &&
        kExpectedRecordCount == kExpectedPanelCount + 2u;
}

bool SelfTest()
{
    int parsed_cpu = -1;
    if (!TestSingletonCpuSetMatcher() ||
        !TestRosterArithmetic() ||
        sizeof(kCellShapes) / sizeof(kCellShapes[0]) != 11u ||
        sizeof(kExpectedCellSeals) / sizeof(kExpectedCellSeals[0]) != 11u ||
        DescriptorJson(NativeSystematicEmission::Invalid).size() != 0u ||
        EmissionMode(NativeSystematicEmission::Invalid) != nullptr ||
        PanelKeySha256(kCellShapes[0], "").size() != 0u ||
        PanelKeySha256(kCellShapes[0], "baseline-aa") !=
            kGoldenPanelKeySha256 ||
        !ParseTargetCpu("120", parsed_cpu) || parsed_cpu != kTargetCpu)
    {
        return false;
    }
    parsed_cpu = -1;
    if (ParseTargetCpu("119", parsed_cpu) || parsed_cpu != -1 ||
        ParseTargetCpu("121", parsed_cpu) || parsed_cpu != -1 ||
        ParseTargetCpu("120x", parsed_cpu) || parsed_cpu != -1 ||
        ParseTargetCpu(nullptr, parsed_cpu) || parsed_cpu != -1)
    {
        return false;
    }
    std::vector<DescriptorSeal> descriptors;
    if (!BuildDescriptorSeals(descriptors) || descriptors.size() != 2u ||
        descriptors[0].Json == descriptors[1].Json ||
        descriptors[0].Sha256 == descriptors[1].Sha256)
    {
        return false;
    }
    bool seals_match = true;
    std::vector<std::shared_ptr<const ScreenCell> > cells;
    for (size_t i = 0u;
         i < sizeof(kCellShapes) / sizeof(kCellShapes[0]);
         ++i)
    {
        std::shared_ptr<ScreenCell> cell(new ScreenCell);
        if (!BuildCell(kCellShapes[i], *cell)) return false;
        if (!CellSealsMatch(i, *cell)) {
            PrintCellSeals(*cell);
            seals_match = false;
        }
        cells.push_back(cell);
    }
    std::ostringstream config;
    if (!EmitConfig(config, kTargetCpu, cells, descriptors)) return false;
    const std::string config_text = config.str();
    const std::string expected_config_boundary =
        std::string("\"panel_key_schema\":\"") + kPanelKeySchema +
        "\",\"panel_replicates\":12,\"schema\":\"" + kConfigSchema + "\"";
    const size_t boundary_offset = config_text.find(expected_config_boundary);
    return seals_match && !config_text.empty() &&
        config_text.front() == '{' && config_text.back() == '\n' &&
        boundary_offset != std::string::npos &&
        config_text.find(expected_config_boundary, boundary_offset + 1u) ==
            std::string::npos;
}

bool RunScreen(int target_cpu)
{
    std::string affinity_diagnostic;
    if (!PinAndVerifyTargetCpu(target_cpu, affinity_diagnostic)) {
        std::cerr << "direct systematic complement affinity failed: "
                  << affinity_diagnostic << '\n';
        return false;
    }
    const std::chrono::steady_clock::time_point deadline =
        std::chrono::steady_clock::now() +
        std::chrono::seconds(kInternalDeadlineSeconds);

    std::vector<std::shared_ptr<const ScreenCell> > cells;
    for (size_t index = 0u;
         index < sizeof(kCellShapes) / sizeof(kCellShapes[0]);
         ++index)
    {
        std::shared_ptr<ScreenCell> cell(new ScreenCell);
        if (!BuildCell(kCellShapes[index], *cell) ||
            !CellSealsMatch(index, *cell))
        {
            std::cerr <<
                "cannot build exact direct systematic complement cell\n";
            return false;
        }
        cells.push_back(cell);
    }
    std::vector<DescriptorSeal> descriptors;
    if (!BuildDescriptorSeals(descriptors)) {
        std::cerr << "invalid direct systematic complement descriptor seal\n";
        return false;
    }
    if (!EmitConfig(std::cout, target_cpu, cells, descriptors)) {
        std::cerr << "direct systematic complement config write failed\n";
        return false;
    }

    uint32_t panel_count = 0u;
    uint64_t measured_invocation_count = 0u;
    uint64_t warmup_invocation_count = 0u;
    uint64_t invocation_count = 0u;
    for (uint32_t replicate = 0u;
         replicate < kPanelReplicates;
         ++replicate)
    {
        for (const std::shared_ptr<const ScreenCell>& cell : cells)
        {
            const uint32_t invocations = InvocationsPerSlot(
                cell->Shape.K, replicate);
            for (const Comparison& comparison : kComparisons)
            {
                const DescriptorSeal* const left_descriptor =
                    FindDescriptor(descriptors, comparison.Left);
                const DescriptorSeal* const right_descriptor =
                    FindDescriptor(descriptors, comparison.Right);
                const std::string panel_key_sha256 = PanelKeySha256(
                    cell->Shape, comparison.Name);
                const NativePanelOrder order = PanelOrder(
                    cell->Shape, comparison.Name, replicate);
                if (!left_descriptor || !right_descriptor ||
                    panel_key_sha256.size() != 64u ||
                    (order != NativePanelOrder::ABBA &&
                     order != NativePanelOrder::BAAB))
                {
                    std::cerr <<
                        "invalid direct systematic complement panel key\n";
                    return false;
                }
                if (std::chrono::steady_clock::now() >= deadline) {
                    std::cerr <<
                        "direct systematic complement internal deadline\n";
                    return false;
                }
                const NativePanelArm left = MakePanelArm(
                    cell,
                    comparison.Left,
                    left_descriptor->Sha256);
                const NativePanelArm right = MakePanelArm(
                    cell,
                    comparison.Right,
                    right_descriptor->Sha256);
                const NativePanelResult panel =
                    wirehair_wh2_bench::ExecuteNativeTimingPanel(
                        target_cpu,
                        order,
                        invocations,
                        left,
                        right);
                if (panel.Status != NativePanelStatus::Complete ||
                    !panel.PanelComparable || panel.IsFatal() ||
                    !panel.HasLeftPreflight ||
                    !panel.HasRightPreflight ||
                    panel.LeftPreflight.Disposition !=
                        NativePanelDisposition::Success ||
                    panel.RightPreflight.Disposition !=
                        NativePanelDisposition::Success ||
                    panel.LeftPreflight.OutcomeCode != Wirehair_Success ||
                    panel.RightPreflight.OutcomeCode != Wirehair_Success)
                {
                    std::cerr <<
                        "direct systematic complement panel failed: "
                              << wirehair_wh2_bench::NativePanelStatusName(
                                     panel.Status)
                              << ' ' << panel.Diagnostic << '\n';
                    return false;
                }
                const uint32_t expected_left_direct =
                    comparison.Left ==
                        NativeSystematicEmission::RetainedSourceDirect ?
                        cell->Shape.K : 0u;
                const uint32_t expected_right_direct =
                    comparison.Right ==
                        NativeSystematicEmission::RetainedSourceDirect ?
                        cell->Shape.K : 0u;
                if (!panel.LeftPreflight.HasDecodedExtra ||
                    !panel.RightPreflight.HasDecodedExtra ||
                    panel.LeftPreflight.DecodedExtra !=
                        expected_left_direct ||
                    panel.RightPreflight.DecodedExtra !=
                        expected_right_direct)
                {
                    std::cerr <<
                        "direct systematic complement preflight receipt drifted\n";
                    return false;
                }
                for (const wirehair_wh2_bench::NativePanelSlot& slot :
                     panel.Slots)
                {
                    if (!slot.HasElapsedNanoseconds ||
                        slot.ElapsedNanoseconds == 0u)
                    {
                        std::cerr <<
                            "direct systematic complement panel has null time\n";
                        return false;
                    }
                }
                if (!EmitPanel(
                        *cell,
                        comparison.Name,
                        invocations,
                        left_descriptor->Sha256,
                        order,
                        panel_key_sha256,
                        replicate,
                        right_descriptor->Sha256,
                        panel))
                {
                    std::cerr <<
                        "direct systematic complement panel write failed\n";
                    return false;
                }
                measured_invocation_count += UINT64_C(4) * invocations;
                warmup_invocation_count += 2u;
                invocation_count += UINT64_C(4) * invocations + 2u;
                ++panel_count;
            }
        }
    }
    if (std::chrono::steady_clock::now() >= deadline) {
        std::cerr << "direct systematic complement internal deadline\n";
        return false;
    }
    if (panel_count != kExpectedPanelCount ||
        measured_invocation_count != kExpectedMeasuredInvocationCount ||
        warmup_invocation_count != kExpectedWarmupInvocationCount ||
        invocation_count != kExpectedInvocationCount)
    {
        std::cerr << "direct systematic complement roster count drifted\n";
        return false;
    }
    std::cout << "{\"invocation_count\":" << invocation_count
              << ",\"measured_invocation_count\":"
              << measured_invocation_count
              << ",\"panel_count\":" << panel_count
              << ",\"record_count\":" << kExpectedRecordCount
              << ",\"schema\":\"" << kTerminalSchema
              << "\",\"status\":\"complete\""
              << ",\"warmup_invocation_count\":"
              << warmup_invocation_count << "}\n";
    std::cout.flush();
    if (!std::cout.good()) {
        std::cerr << "direct systematic complement terminal write failed\n";
        return false;
    }
    if (std::chrono::steady_clock::now() >= deadline) {
        std::cerr << "direct systematic complement internal deadline\n";
        return false;
    }
    return true;
}

} // namespace

int main(int argc, char** argv)
{
#if defined(SIGPIPE)
    if (std::signal(SIGPIPE, SIG_IGN) == SIG_ERR) {
        std::cerr << "cannot ignore SIGPIPE for transactional output\n";
        return 1;
    }
#endif
    if (argc == 2 && std::strcmp(argv[1], "--selftest") == 0)
    {
        if (!SelfTest()) {
            std::cerr <<
                "WH2 direct systematic complement selftest failed\n";
            return 1;
        }
        std::cout <<
            "WH2 direct systematic complement selftest passed\n" <<
            std::flush;
        return std::cout.good() ? 0 : 1;
    }
    int target_cpu = -1;
    if (argc != 4 || std::strcmp(argv[1], "--run") != 0 ||
        std::strcmp(argv[2], "--cpu") != 0 ||
        !ParseTargetCpu(argv[3], target_cpu))
    {
        std::cerr << "usage: " << argv[0]
                  << " --selftest | --run --cpu 120\n";
        return 2;
    }
    return RunScreen(target_cpu) ? 0 : 1;
}
