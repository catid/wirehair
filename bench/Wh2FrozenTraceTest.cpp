#include "Wh2FrozenTrace.h"

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <string>
#include <vector>

namespace {

using wirehair::wh2_benchmark::CanonicalRecoveryCellJson;
using wirehair::wh2_benchmark::CopyNestedPrefix;
using wirehair::wh2_benchmark::DevelopmentRecoveryDomainSha256;
using wirehair::wh2_benchmark::EnumerateDevelopmentRecoveryCells;
using wirehair::wh2_benchmark::FrozenCandidateLimit;
using wirehair::wh2_benchmark::FrozenPacketTrace;
using wirehair::wh2_benchmark::FrozenRecoveryCell;
using wirehair::wh2_benchmark::FrozenSchedule;
using wirehair::wh2_benchmark::FrozenTraceSeed;
using wirehair::wh2_benchmark::FrozenTraceStatus;
using wirehair::wh2_benchmark::GenerateFrozenPacketTrace;
using wirehair::wh2_benchmark::PacketIdsSha256;
using wirehair::wh2_benchmark::RecoveryCellSha256;
using wirehair::wh2_benchmark::Sha256Hex;

unsigned g_failures = 0u;

bool Check(bool condition, const char* message)
{
    if (condition) {
        return true;
    }
    std::fprintf(stderr, "WH2 frozen trace test failed: %s\n", message);
    ++g_failures;
    return false;
}

bool CheckTraceGolden(
    uint32_t K,
    uint32_t loss_ppm,
    FrozenSchedule schedule,
    const std::vector<uint32_t>& expected_ids,
    uint64_t expected_candidates,
    const char* expected_sha256)
{
    FrozenPacketTrace trace;
    const FrozenTraceStatus status = GenerateFrozenPacketTrace(
        K,
        2u,
        loss_ppm,
        UINT64_C(0xd1b54a32d192ed03),
        schedule,
        trace);
    return
        Check(status == FrozenTraceStatus::Complete,
            "golden trace did not complete") &&
        Check(trace.delivered_ids == expected_ids,
            "golden delivered packet IDs changed") &&
        Check(trace.attempted_candidates == expected_candidates,
            "golden attempted-candidate count changed") &&
        Check(trace.trace_sha256 == expected_sha256,
            "golden trace SHA-256 changed") &&
        Check(PacketIdsSha256(trace.delivered_ids) == expected_sha256,
            "packet-ID SHA-256 disagrees with trace receipt") &&
        Check(trace.candidate_limit ==
                UINT64_C(256) * (static_cast<uint64_t>(K) + 4u) + 65536u,
            "golden candidate limit changed");
}

void TestSha256()
{
    Check(
        Sha256Hex(static_cast<const void*>(NULL), 0u) ==
            "e3b0c44298fc1c149afbf4c8996fb924"
            "27ae41e4649b934ca495991b7852b855",
        "empty SHA-256 golden changed");
    Check(
        Sha256Hex(std::string("abc")) ==
            "ba7816bf8f01cfea414140de5dae2223"
            "b00361a396177a9cb410ff61f20015ad",
        "abc SHA-256 golden changed");
    Check(Sha256Hex(static_cast<const void*>(NULL), 1u).empty(),
        "non-empty null SHA input was accepted");

    const std::vector<uint32_t> words = {
        UINT32_C(0), UINT32_C(0x12345678), UINT32_C(0xffffffff)
    };
    Check(
        PacketIdsSha256(words) ==
            "cdf516d8944dafbd0e8a427e7517d767"
            "a63667051927d23e2d5380d6a6d86921",
        "packet IDs are not hashed as explicit little-endian words");
}

void TestDefaultRecoveryCellIsSafelyInvalid()
{
    const FrozenRecoveryCell cell;
    Check(CanonicalRecoveryCellJson(cell).empty() &&
            RecoveryCellSha256(cell).empty(),
        "default recovery cell was not safely rejected");
}

void TestDevelopmentEnumeration()
{
    static const uint32_t short_K[] = {
        2u, 3u, 4u, 5u, 6u, 8u, 16u, 32u, 64u, 100u,
        101u, 128u, 256u, 512u, 513u, 1000u,
        1001u, 2048u, 4096u, 5000u,
        5001u, 8192u, 10000u,
        10001u, 16384u, 20000u,
        20001u, 32768u, 49152u, 64000u
    };
    static const uint64_t roots[] = {
        UINT64_C(0xd1b54a32d192ed03),
        UINT64_C(0x94d049bb133111eb),
        UINT64_C(0x8538ecb5bd456ea3)
    };
    static const FrozenSchedule schedules[] = {
        FrozenSchedule::Iid,
        FrozenSchedule::Burst,
        FrozenSchedule::Adversarial,
        FrozenSchedule::RepairOnly
    };
    static const uint32_t losses[] = {
        100000u, 500000u, 500000u, 500000u
    };

    const std::vector<FrozenRecoveryCell> cells =
        EnumerateDevelopmentRecoveryCells();
    Check(cells.size() == 360u, "development cell cardinality is not 360");
    for (std::size_t ordinal = 0u; ordinal < cells.size(); ++ordinal)
    {
        const FrozenRecoveryCell& cell = cells[ordinal];
        const uint32_t trial = static_cast<uint32_t>(ordinal / 120u);
        const uint32_t stratum =
            static_cast<uint32_t>((ordinal % 120u) / 30u);
        const std::size_t k_index = ordinal % 30u;
        Check(cell.ordinal == ordinal, "development ordinal changed");
        Check(cell.trial == trial, "development trial ordering changed");
        Check(cell.base_seed_attempt == trial,
            "development public attempt pairing changed");
        Check(cell.loss_seed == roots[trial],
            "development loss-root pairing changed");
        Check(cell.stratum_index == stratum,
            "development stratum ordering changed");
        Check(cell.schedule == schedules[stratum],
            "development schedule ordering changed");
        Check(cell.loss_ppm == losses[stratum],
            "development loss stratum changed");
        Check(cell.K == short_K[k_index],
            "development K ordering changed");
        Check(cell.phase == "development" && cell.block_bytes == 2u &&
                cell.overhead_cap == 4u,
            "development fixed cell field changed");
    }

    const std::string first_json =
        "{\"K\":2,\"band\":\"2-100\",\"base_seed_attempt\":0,"
        "\"block_bytes\":2,\"loss_ppm\":100000,"
        "\"loss_seed\":\"0xd1b54a32d192ed03\",\"overhead_cap\":4,"
        "\"phase\":\"development\",\"schedule\":\"iid\",\"trial\":0}";
    Check(CanonicalRecoveryCellJson(cells[0]) == first_json,
        "first canonical development cell changed");
    Check(
        RecoveryCellSha256(cells[0]) ==
            "e0bbe2f2afee6477e2091d25ee4b52f6"
            "be79031c7861862f308eb62b9864068c",
        "first canonical development cell hash changed");
    Check(
        RecoveryCellSha256(cells[30]) ==
            "b96d81b0c2995b41d16c6fbb6337c8f"
            "8d8dae104269a0683a707969c03b550a6",
        "first burst cell hash changed");
    Check(
        RecoveryCellSha256(cells[359]) ==
            "22dde87355124189d330a4777ad070533"
            "fd7ef0c35d2ef511319e48065387600",
        "last development cell hash changed");
    Check(
        DevelopmentRecoveryDomainSha256() ==
            "f97f28c211428cd77aed97160073b192"
            "d93014cb4a61a844bc7d76375ac61b77",
        "development recovery-domain SHA-256 changed");

    FrozenPacketTrace valid;
    Check(GenerateFrozenPacketTrace(cells[0], valid) ==
            FrozenTraceStatus::Complete,
        "enumerated development cell was rejected by trace generator");
    FrozenRecoveryCell relabeled = cells[0];
    relabeled.phase = "final_raw";
    FrozenPacketTrace invalid;
    Check(GenerateFrozenPacketTrace(relabeled, invalid) ==
            FrozenTraceStatus::InvalidInput,
        "relabeled development cell was accepted");
    relabeled = cells[0];
    relabeled.ordinal = 1u;
    Check(CanonicalRecoveryCellJson(relabeled).empty() &&
            GenerateFrozenPacketTrace(relabeled, invalid) ==
                FrozenTraceStatus::InvalidInput,
        "wrong-ordinal development cell was accepted");
}

void TestFirstCellTraceGoldens()
{
    Check(
        FrozenTraceSeed(2u, 2u, UINT64_C(0xd1b54a32d192ed03)) ==
            UINT64_C(0x936b379a16cfde5b),
        "first-cell trace seed changed");

    CheckTraceGolden(
        2u,
        100000u,
        FrozenSchedule::Iid,
        std::vector<uint32_t>{0u, 1u, 3u, 4u, 5u, 6u},
        7u,
        "37391e3501b8cb36bf90f089771535ab"
        "6824a2ac518e0c00ac7c600746b5c3ab");
    CheckTraceGolden(
        2u,
        500000u,
        FrozenSchedule::Burst,
        std::vector<uint32_t>{0u, 1u, 2u, 3u, 4u, 5u},
        6u,
        "cd9a54ed1f18bf97db08914e280ea734"
        "9e11ca2c4885a4d8052552ceba84208d");
    CheckTraceGolden(
        2u,
        500000u,
        FrozenSchedule::Adversarial,
        std::vector<uint32_t>{
            UINT32_C(4294967295), UINT32_C(4294967291),
            UINT32_C(4294967287), UINT32_C(4294967283),
            UINT32_C(4294967281), UINT32_C(4294967271)
        },
        13u,
        "e5ce1cd184f9ff951dcfcf196f1656ba"
        "a15fadb745b4dd94f05fb5d1e10375ae");
    CheckTraceGolden(
        2u,
        500000u,
        FrozenSchedule::RepairOnly,
        std::vector<uint32_t>{2u, 4u, 6u, 8u, 9u, 14u},
        13u,
        "ff58d29f143b0d1d52f6a6f25b0c2e"
        "f53ed4a123898639918cf72091ea6d3b5a");
}

void TestK8EvidenceGoldens()
{
    Check(
        FrozenTraceSeed(8u, 2u, UINT64_C(0xd1b54a32d192ed03)) ==
            UINT64_C(0x5ebe09231208c6d9),
        "K=8 trace seed changed");

    static const FrozenSchedule schedules[] = {
        FrozenSchedule::Iid,
        FrozenSchedule::Burst,
        FrozenSchedule::Adversarial,
        FrozenSchedule::RepairOnly
    };
    static const char* hashes[] = {
        "d0077759638cb99b6084690f4d063d18"
        "a22cf3ab4286495f544698de72878e55",
        "38ded04446783ead8f40d1307a5abfb1"
        "30d8f00f3320d455e2be447886874ebd",
        "9252d166ff6476991af38963e2052ed0"
        "aad5d708a77a0cc53d937739c5dc0419",
        "52cbd184e73ad47ae1f99caaca72a555"
        "7c9025d238570e3022007c300e88fb72"
    };
    for (unsigned i = 0u; i < 4u; ++i)
    {
        FrozenPacketTrace trace;
        Check(GenerateFrozenPacketTrace(
                8u,
                2u,
                500000u,
                UINT64_C(0xd1b54a32d192ed03),
                schedules[i],
                trace) == FrozenTraceStatus::Complete,
            "K=8 evidence trace did not complete");
        Check(trace.trace_sha256 == hashes[i],
            "K=8 evidence trace hash changed");
    }
}

void TestDeterminismAndNestedPrefixes()
{
    static const FrozenSchedule schedules[] = {
        FrozenSchedule::Iid,
        FrozenSchedule::Burst,
        FrozenSchedule::Adversarial,
        FrozenSchedule::RepairOnly
    };
    for (unsigned i = 0u; i < 4u; ++i)
    {
        FrozenPacketTrace first;
        FrozenPacketTrace second;
        const uint32_t loss_ppm =
            schedules[i] == FrozenSchedule::Iid ? 100000u : 500000u;
        Check(GenerateFrozenPacketTrace(
                128u,
                1280u,
                loss_ppm,
                UINT64_C(0x94d049bb133111eb),
                schedules[i],
                first) == FrozenTraceStatus::Complete,
            "first deterministic trace did not complete");
        Check(GenerateFrozenPacketTrace(
                128u,
                1280u,
                loss_ppm,
                UINT64_C(0x94d049bb133111eb),
                schedules[i],
                second) == FrozenTraceStatus::Complete,
            "second deterministic trace did not complete");
        Check(first.delivered_ids == second.delivered_ids &&
                first.attempted_candidates == second.attempted_candidates &&
                first.trace_sha256 == second.trace_sha256,
            "identical trace inputs were not deterministic");

        std::vector<uint32_t> previous;
        for (uint32_t overhead = 0u; overhead <= 4u; ++overhead)
        {
            std::vector<uint32_t> prefix;
            Check(CopyNestedPrefix(first, overhead, prefix),
                "valid nested prefix was rejected");
            Check(prefix.size() == 128u + overhead,
                "nested prefix has wrong K+overhead size");
            Check(std::equal(
                    prefix.begin(), prefix.end(), first.delivered_ids.begin()),
                "nested prefix changed the frozen delivered-ID stream");
            Check(previous.empty() || std::equal(
                    previous.begin(), previous.end(), prefix.begin()),
                "nested thresholds are not exact prefixes");
            previous.swap(prefix);
        }

        std::vector<uint32_t> rejected(1u, 123u);
        Check(!CopyNestedPrefix(first, 5u, rejected) && rejected.empty(),
            "overhead above the frozen cap was accepted");

        FrozenPacketTrace corrupted = first;
        corrupted.delivered_ids[0] ^= 1u;
        Check(!CopyNestedPrefix(corrupted, 0u, rejected) && rejected.empty(),
            "trace mutated after hashing produced a nested prefix");
    }
}

void TestFullDevelopmentTraceAggregate()
{
    const std::vector<FrozenRecoveryCell> cells =
        EnumerateDevelopmentRecoveryCells();
    std::string aggregate;
    aggregate.reserve(25200u);
    uint64_t total_attempts = 0u;
    uint64_t maximum_attempts = 0u;
    for (std::size_t i = 0u; i < cells.size(); ++i)
    {
        FrozenPacketTrace trace;
        Check(GenerateFrozenPacketTrace(cells[i], trace) ==
                FrozenTraceStatus::Complete,
            "development trace aggregate contains an incomplete cell");
        Check(trace.delivered_ids.size() ==
                static_cast<std::size_t>(cells[i].K) + 4u,
            "development trace aggregate cell has wrong cardinality");
        aggregate += trace.trace_sha256;
        aggregate += ':';
        aggregate += std::to_string(trace.attempted_candidates);
        aggregate += '\n';
        total_attempts += trace.attempted_candidates;
        maximum_attempts = std::max(
            maximum_attempts, trace.attempted_candidates);
    }
    Check(total_attempts == UINT64_C(5345254),
        "development trace aggregate candidate count changed");
    Check(maximum_attempts == UINT64_C(128334),
        "development trace aggregate maximum candidate count changed");
    Check(
        Sha256Hex(aggregate) ==
            "ec6f5ed3976e5e9664bb5db927c525a7"
            "3e8ab952d080492aaecd99c366db41a9",
        "development trace aggregate SHA-256 changed");
}

void TestCandidateCapAndInvalidInputs()
{
    uint64_t candidate_limit = 123u;
    Check(FrozenCandidateLimit(6u, candidate_limit) &&
            candidate_limit == 67072u,
        "first-cell candidate cap changed");

    const uint64_t maximum = std::numeric_limits<uint64_t>::max();
    const uint64_t largest_safe = (maximum - 65536u) / 256u;
    Check(FrozenCandidateLimit(largest_safe, candidate_limit) &&
            candidate_limit == largest_safe * 256u + 65536u,
        "largest safe candidate-cap input was rejected");
    Check(!FrozenCandidateLimit(largest_safe + 1u, candidate_limit) &&
            candidate_limit == 0u,
        "overflowing candidate-cap input was accepted");

    FrozenPacketTrace capped;
    Check(GenerateFrozenPacketTrace(
            2u,
            2u,
            1000000u,
            UINT64_C(0xd1b54a32d192ed03),
            FrozenSchedule::Iid,
            capped) == FrozenTraceStatus::CandidateLimitReached,
        "full-loss trace did not stop at the candidate cap");
    Check(capped.attempted_candidates == capped.candidate_limit &&
            capped.candidate_limit == 67072u &&
            capped.delivered_ids.empty() && capped.trace_sha256.empty(),
        "candidate-cap failure receipt is inconsistent");

    FrozenPacketTrace boundary;
    Check(GenerateFrozenPacketTrace(
            64000u,
            4096u,
            0u,
            UINT64_C(0x8538ecb5bd456ea3),
            FrozenSchedule::RepairOnly,
            boundary) == FrozenTraceStatus::Complete,
        "maximum frozen K trace did not complete");
    Check(boundary.delivered_ids.size() == 64004u &&
            boundary.attempted_candidates == 64004u &&
            boundary.delivered_ids.front() == 64000u &&
            boundary.delivered_ids.back() == 128003u &&
            boundary.candidate_limit == 16450560u,
        "maximum frozen K trace arithmetic changed");

    FrozenPacketTrace invalid;
    const struct InvalidCase {
        uint32_t K;
        uint32_t block_bytes;
        uint32_t loss_ppm;
        FrozenSchedule schedule;
    } cases[] = {
        {1u, 2u, 0u, FrozenSchedule::Iid},
        {64001u, 2u, 0u, FrozenSchedule::Iid},
        {2u, 0u, 0u, FrozenSchedule::Iid},
        {2u, 2u, 1000001u, FrozenSchedule::Iid},
        {2u, 2u, 0u, static_cast<FrozenSchedule>(99)}
    };
    for (std::size_t i = 0u; i < sizeof(cases) / sizeof(cases[0]); ++i)
    {
        invalid = boundary;
        Check(GenerateFrozenPacketTrace(
                cases[i].K,
                cases[i].block_bytes,
                cases[i].loss_ppm,
                UINT64_C(1),
                cases[i].schedule,
                invalid) == FrozenTraceStatus::InvalidInput,
            "invalid frozen-trace input was accepted");
        Check(invalid.delivered_ids.empty() &&
                invalid.trace_sha256.empty() &&
                invalid.attempted_candidates == 0u,
            "invalid input retained prior trace state");
    }

    std::vector<uint32_t> prefix;
    Check(!CopyNestedPrefix(capped, 0u, prefix),
        "incomplete trace produced a nested prefix");
}

} // namespace

int main()
{
    TestSha256();
    TestDefaultRecoveryCellIsSafelyInvalid();
    TestDevelopmentEnumeration();
    TestFirstCellTraceGoldens();
    TestK8EvidenceGoldens();
    TestDeterminismAndNestedPrefixes();
    TestFullDevelopmentTraceAggregate();
    TestCandidateCapAndInvalidInputs();

    if (g_failures != 0u) {
        std::fprintf(stderr, "WH2 frozen trace test: %u failure(s)\n", g_failures);
        return 1;
    }
    std::printf("WH2 frozen trace test: PASS\n");
    return 0;
}
