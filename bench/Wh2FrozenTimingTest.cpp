#include "Wh2FrozenTiming.h"

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

namespace {

using namespace wirehair::wh2_benchmark;

unsigned g_failures = 0u;

bool Check(bool condition, const char* message)
{
    if (condition) {
        return true;
    }
    std::fprintf(stderr, "WH2 frozen timing test failed: %s\n", message);
    ++g_failures;
    return false;
}

void TestTimingSeeds()
{
    static const uint32_t attempts[] = {
        0u, 0u, 0u, 0u, 0u, 0u, 0u, 0u,
        0u, 0u, 0u, 0u, 0u, 0u, 0u, 0u
    };
    static const uint64_t seeds[] = {
        UINT64_C(0x2d0f28c7e7e786b2),
        UINT64_C(0x49c18e491a66eb93),
        UINT64_C(0xe93dcff7eed78181),
        UINT64_C(0xe9e7bf14d7c14038),
        UINT64_C(0x8886c61661bd0c88),
        UINT64_C(0x4694f3b073a0c6c7),
        UINT64_C(0x7d771d8bc3d7af5f),
        UINT64_C(0xe4db579c79a43708),
        UINT64_C(0xc34be295f3f62201),
        UINT64_C(0x8df810ac40206963),
        UINT64_C(0xc1ba1dde6ccd5c85),
        UINT64_C(0xc5b7018f93a18591),
        UINT64_C(0x554ea6750448cc4c),
        UINT64_C(0x5a1602ce508795aa),
        UINT64_C(0x6e670f6df8b37a57),
        UINT64_C(0x2302e830b01627be)
    };
    for (uint32_t replicate = 0u; replicate < 16u; ++replicate)
    {
        uint32_t attempt = UINT32_MAX;
        uint64_t seed = UINT64_MAX;
        Check(DevelopmentTimingSeed(replicate, attempt, seed),
            "valid development timing replicate was rejected");
        Check(attempt == attempts[replicate],
            "fixed production construction attempt changed");
        Check(seed == seeds[replicate],
            "paired SplitMix64 timing loss seed changed");
    }

    uint32_t attempt = 9u;
    uint64_t seed = 9u;
    for (uint32_t replicate = 16u; replicate < 96u; ++replicate)
    {
        uint32_t repeated_attempt = UINT32_MAX;
        uint64_t repeated_seed = UINT64_MAX;
        Check(DevelopmentTimingSeed(replicate, attempt, seed) &&
                DevelopmentTimingSeed(
                    replicate, repeated_attempt, repeated_seed) &&
                attempt == 0u && repeated_attempt == 0u &&
                seed == repeated_seed,
            "extended timing seed domain is not deterministic");
    }
    Check(seed == UINT64_C(0xecc8a5db602b9a94),
        "last extended timing loss-seed golden changed");

    Check(!DevelopmentTimingSeed(96u, attempt, seed) &&
            attempt == 0u && seed == 0u,
        "out-of-domain timing replicate was accepted");
}

void TestInvocationsPerSlotFormula()
{
    Check(DevelopmentTimingInvocationsPerSlot(0u) == 0u,
        "zero K did not produce an invalid zero batch count");
    Check(DevelopmentTimingInvocationsPerSlot(1u) == 65536u,
        "K=1 batch count changed");
    Check(DevelopmentTimingInvocationsPerSlot(2u) == 32768u,
        "exactly divisible batch count changed");
    Check(DevelopmentTimingInvocationsPerSlot(3u) == 21846u,
        "batch count did not round division upward");
    Check(DevelopmentTimingInvocationsPerSlot(65535u) == 2u &&
            DevelopmentTimingInvocationsPerSlot(65536u) == 2u &&
            DevelopmentTimingInvocationsPerSlot(65537u) == 2u &&
            DevelopmentTimingInvocationsPerSlot(UINT32_MAX) == 2u,
        "batch-count boundary or overflow handling changed");
}

void TestExecutionGeometry()
{
    Check(DevelopmentTimingWorkerCount() == 8u &&
            DevelopmentTimingCohortCount() == 264u,
        "frozen timing worker/cohort cardinality changed");
    const std::vector<FrozenTimingBaseCell> cells =
        EnumerateDevelopmentTimingBaseCells();
    const std::vector<FrozenTimingPanel> panels =
        EnumerateOneCandidateTimingPanels("candidate");
    Check(cells.size() == 2304u && panels.size() == 11u &&
            cells.size() / 96u * panels.size() ==
                DevelopmentTimingCohortCount(),
        "timing cohort count does not match cells and panels");

    for (std::size_t cohort = 0u;
        cohort < DevelopmentTimingCohortCount();
        ++cohort)
    {
        for (uint32_t round = 0u; round < 12u; ++round)
        {
            bool seen[8] = {};
            for (uint32_t local = 0u; local < 8u; ++local)
            {
                const uint32_t replicate = round * 8u + local;
                uint32_t slot = UINT32_MAX;
                const uint32_t expected = static_cast<uint32_t>(
                    (local + cohort + round) % 8u);
                Check(DevelopmentTimingWorkerSlot(
                        cohort, replicate, slot) &&
                        slot == expected && !seen[slot],
                    "timing wave is not an exact Latin worker rotation");
                if (slot < 8u) {
                    seen[slot] = true;
                }
            }
            Check(std::count(seen, seen + 8, true) == 8,
                "timing wave does not use every worker exactly once");
        }
    }

    for (uint32_t replicate = 0u; replicate < 96u; ++replicate)
    {
        uint32_t counts[8] = {};
        for (std::size_t cohort = 0u;
            cohort < DevelopmentTimingCohortCount();
            ++cohort)
        {
            uint32_t slot = UINT32_MAX;
            if (DevelopmentTimingWorkerSlot(cohort, replicate, slot) &&
                slot < 8u)
            {
                ++counts[slot];
            }
        }
        Check(std::count(counts, counts + 8, 33u) == 8,
            "one timing replicate is not balanced over worker slots");
    }

    uint32_t slot = UINT32_MAX;
    Check(!DevelopmentTimingWorkerSlot(264u, 0u, slot) && slot == 0u,
        "out-of-domain timing cohort was accepted");
    slot = UINT32_MAX;
    Check(!DevelopmentTimingWorkerSlot(0u, 96u, slot) && slot == 0u,
        "out-of-domain timing replicate received a worker slot");
}

void TestDefaultTimingValuesAreSafelyInvalid()
{
    const FrozenTimingBaseCell base_cell;
    Check(CanonicalTimingBaseCellJson(base_cell).empty() &&
            TimingBaseCellSha256(base_cell).empty(),
        "default timing base cell was not safely rejected");
    const FrozenTimingCell cell;
    Check(CanonicalTimingCellJson(cell).empty() &&
            TimingCellSha256(cell).empty() &&
            CanonicalTimingSourceJson(cell).empty() &&
            TimingSourceIdentitySha256(cell).empty(),
        "default timing cell was not safely rejected");
    const FrozenTimingPanel panel;
    Check(CanonicalTimingPanelJson(panel).empty() &&
            TimingPanelSha256(panel).empty() &&
            TimingPanelOrder(panel, 0u) == FrozenTimingOrder::Invalid,
        "default timing panel was not safely rejected");
}

void TestTimingBaseCellEnumeration()
{
    static const uint32_t K_values[] = {
        8u, 32u, 100u, 128u, 512u, 1000u,
        2048u, 5000u, 8192u, 20000u, 32768u, 64000u
    };
    static const uint32_t widths[] = { 64u, 1280u };

    const std::vector<FrozenTimingBaseCell> cells =
        EnumerateDevelopmentTimingBaseCells();
    if (!Check(cells.size() == 2304u,
            "development timing base-cell cardinality is not 2304"))
    {
        return;
    }
    for (std::size_t ordinal = 0u; ordinal < cells.size(); ++ordinal)
    {
        const FrozenTimingBaseCell& cell = cells[ordinal];
        const uint32_t round = static_cast<uint32_t>(ordinal / 192u);
        const std::size_t round_ordinal = ordinal % 192u;
        const std::size_t width_index = round_ordinal / 96u;
        const std::size_t k_index = (round_ordinal % 96u) / 8u;
        const uint32_t lane = static_cast<uint32_t>(round_ordinal % 8u);
        const uint32_t replicate = round * 8u + lane;
        uint32_t expected_attempt = 0u;
        uint64_t expected_seed = 0u;
        Check(DevelopmentTimingSeed(
                replicate, expected_attempt, expected_seed),
            "enumerated timing replicate has no paired seed");
        Check(cell.ordinal == ordinal,
            "development timing base ordinal changed");
        Check(cell.block_bytes == widths[width_index],
            "development timing width ordering changed");
        Check(cell.replicate == replicate,
            "development timing replicate ordering changed");
        Check(cell.K == K_values[k_index],
            "development timing K ordering changed");
        Check(cell.base_seed_attempt == expected_attempt &&
                cell.base_loss_seed == expected_seed,
            "development timing base-seed pairing changed");
        Check(cell.phase == "development" &&
                cell.loss_ppm == 100000u &&
                cell.schedule == FrozenSchedule::Iid &&
                cell.fixed_received_overhead == 4u &&
                cell.receive_overhead_cap == 256u &&
                cell.interleave_policy ==
                    "self-counterbalanced-repeat-major-v1" &&
                cell.invocations_per_slot ==
                    DevelopmentTimingInvocationsPerSlot(cell.K),
            "development timing fixed field changed");
    }
    Check(cells[7].K == 8u && cells[7].replicate == 7u &&
            cells[8].K == 32u && cells[8].replicate == 0u &&
            cells[191].block_bytes == 1280u &&
            cells[191].K == 64000u && cells[191].replicate == 7u &&
            cells[192].block_bytes == 64u &&
            cells[192].K == 8u && cells[192].replicate == 8u,
        "development timing bases are not frozen in round-major order");

    const std::string first_json =
        "{\"K\":8,\"band\":\"2-100\","
        "\"base_loss_seed\":\"0x2d0f28c7e7e786b2\","
        "\"base_seed_attempt\":0,\"block_bytes\":64,"
        "\"fixed_received_overhead\":4,"
        "\"interleave_policy\":\"self-counterbalanced-repeat-major-v1\","
        "\"invocations_per_slot\":8192,\"loss_ppm\":100000,"
        "\"phase\":\"development\",\"receive_overhead_cap\":256,"
        "\"replicate\":0,\"schedule\":\"iid\"}";
    Check(CanonicalTimingBaseCellJson(cells[0]) == first_json,
        "first canonical timing base cell changed");
    Check(
        TimingBaseCellSha256(cells[0]) ==
            "158931044b779851eca21bd3112bce798"
            "77e6e59f5896aebf698466789c47814",
        "first timing base-cell SHA-256 changed");
    Check(
        TimingBaseCellSha256(cells[96]) ==
            "014c217757d7fbc7d0f41826edde37627"
            "ec2f915050e21b483d7d21695e99a29",
        "first wide timing base-cell SHA-256 changed");
    Check(
        TimingBaseCellSha256(cells[2303]) ==
            "430f5f12c551c5eddbbc6364e7586d56"
            "8334a858cfab5f6169fd1fd037bd7b4b",
        "last timing base-cell SHA-256 changed");
    Check(
        DevelopmentTimingBaseDomainSha256() ==
            "eab25fe96642d8dd12d4b64e91fe00b5"
            "33d464f1bd7487c84e670ce888e2d164",
        "development timing base-domain SHA-256 changed");

    const std::vector<FrozenTimingBaseCell> repeat =
        EnumerateDevelopmentTimingBaseCells();
    Check(repeat.size() == cells.size(),
        "repeated timing base enumeration changed cardinality");
    for (std::size_t i = 0u; i < cells.size() && i < repeat.size(); ++i) {
        Check(CanonicalTimingBaseCellJson(cells[i]) ==
                CanonicalTimingBaseCellJson(repeat[i]),
            "repeated timing base enumeration changed canonical bytes");
    }

    FrozenTimingBaseCell mutant = cells[0];
    mutant.ordinal = 1u;
    Check(CanonicalTimingBaseCellJson(mutant).empty(),
        "wrong timing-base ordinal was accepted");
    mutant = cells[0];
    mutant.schedule = FrozenSchedule::Burst;
    Check(CanonicalTimingBaseCellJson(mutant).empty(),
        "non-IID development timing base was accepted");
    mutant = cells[0];
    mutant.base_loss_seed ^= 1u;
    Check(CanonicalTimingBaseCellJson(mutant).empty(),
        "wrong paired timing base loss seed was accepted");
    mutant = cells[0];
    --mutant.invocations_per_slot;
    Check(CanonicalTimingBaseCellJson(mutant).empty(),
        "wrong timing invocation batch was accepted");
    mutant = cells[0];
    mutant.interleave_policy = "self-counterbalanced-slot-major-v1";
    Check(CanonicalTimingBaseCellJson(mutant).empty(),
        "wrong timing interleave policy was accepted");
    mutant = cells[0];
    mutant.receive_overhead_cap = 255u;
    Check(CanonicalTimingBaseCellJson(mutant).empty(),
        "wrong timing receive-overhead cap was accepted");
    mutant = cells[0];
    mutant.base_seed_attempt = 1u;
    Check(CanonicalTimingBaseCellJson(mutant).empty(),
        "non-production development timing attempt was accepted");
}

void TestTimingQualificationAndSourceIdentity()
{
    const std::vector<FrozenTimingBaseCell> bases =
        EnumerateDevelopmentTimingBaseCells();
    if (!Check(bases.size() == 2304u,
            "cannot test qualification without the frozen base domain"))
    {
        return;
    }

    std::vector<uint32_t> retry_offsets(2304u, 0u);
    const std::vector<FrozenTimingCell> cells =
        EnumerateDevelopmentTimingCells(retry_offsets);
    if (!Check(cells.size() == bases.size(),
            "explicit retry-zero map did not qualify every timing cell"))
    {
        return;
    }
    Check(EnumerateDevelopmentTimingCells(
            std::vector<uint32_t>()).empty() &&
            EnumerateDevelopmentTimingCells(
                std::vector<uint32_t>(2303u, 0u)).empty() &&
            EnumerateDevelopmentTimingCells(
                std::vector<uint32_t>(2305u, 0u)).empty() &&
            DevelopmentTimingDomainSha256(
                std::vector<uint32_t>()).empty() &&
            DevelopmentTimingDomainSha256(
                std::vector<uint32_t>(2303u, 0u)).empty() &&
            DevelopmentTimingDomainSha256(
                std::vector<uint32_t>(2305u, 0u)).empty(),
        "missing or oversized timing qualification map was accepted");
    std::vector<uint32_t> invalid_offsets(2304u, 0u);
    invalid_offsets[137u] = 256u;
    Check(EnumerateDevelopmentTimingCells(invalid_offsets).empty() &&
            DevelopmentTimingDomainSha256(invalid_offsets).empty(),
        "non-uint8 timing retry offset was accepted");

    const std::string first_json =
        "{\"K\":8,\"band\":\"2-100\","
        "\"base_cell_sha256\":\"158931044b779851eca21bd3112bce798"
        "77e6e59f5896aebf698466789c47814\","
        "\"base_loss_seed\":\"0x2d0f28c7e7e786b2\","
        "\"base_seed_attempt\":0,\"block_bytes\":64,"
        "\"fixed_received_overhead\":4,"
        "\"interleave_policy\":\"self-counterbalanced-repeat-major-v1\","
        "\"invocations_per_slot\":8192,\"loss_ppm\":100000,"
        "\"loss_retry_offset\":0,"
        "\"loss_seed\":\"0x2d0f28c7e7e786b2\","
        "\"phase\":\"development\",\"receive_overhead_cap\":256,"
        "\"replicate\":0,\"schedule\":\"iid\"}";
    Check(CanonicalTimingCellJson(cells[0]) == first_json,
        "first canonical qualified timing cell changed");
    Check(TimingCellSha256(cells[0]) ==
            "e85871b9c09ac946f6c8f6e0e077bd75"
            "8d4862bfb5dddd73f6ee05e3ba35f2d7",
        "retry-zero timing cell SHA-256 changed");
    Check(cells[0].base_cell_sha256 == TimingBaseCellSha256(bases[0]) &&
            cells[0].loss_retry_offset == 0u &&
            cells[0].loss_seed == cells[0].base_loss_seed,
        "retry-zero timing qualification fields are inconsistent");

    FrozenTimingCell retry_one;
    FrozenTimingCell retry_two;
    FrozenTimingCell retry_255;
    Check(QualifyDevelopmentTimingCell(bases[0], 1u, retry_one) &&
            QualifyDevelopmentTimingCell(bases[0], 2u, retry_two) &&
            QualifyDevelopmentTimingCell(bases[0], 255u, retry_255),
        "valid timing retry offsets were rejected");
    Check(retry_one.loss_seed == UINT64_C(0xcb46a281673202c7) &&
            retry_two.loss_seed == UINT64_C(0x697e1c3ae67c7edc) &&
            retry_255.loss_seed == UINT64_C(0xc651688db3191f9d),
        "timing retry seed stride or uint64 wrap changed");
    Check(TimingCellSha256(retry_one) ==
            "4c216d6ca9416585f777980fc1ffd0a9f"
            "162644606d5f0468accec353470cf37" &&
            TimingCellSha256(retry_two) ==
            "3e3aaf22b4854b9395c43aec5d50ca16"
            "524fe6259245c4ced377b0b12e2cd9be" &&
            TimingCellSha256(retry_255) ==
            "22968aeb7c50a96c6eab4ce75b7653cd"
            "82b2a2df0ac8d93182d6d06f47b46ae1",
        "qualified timing retry identity golden changed");

    uint64_t wrapped_seed = UINT64_MAX;
    Check(DevelopmentTimingLossSeed(
            UINT64_MAX, 1u, wrapped_seed) &&
            wrapped_seed == UINT64_C(0x9e3779b97f4a7c14),
        "explicit timing retry did not wrap modulo 2^64");
    wrapped_seed = UINT64_MAX;
    Check(!DevelopmentTimingLossSeed(
            bases[0].base_loss_seed, 256u, wrapped_seed) &&
            wrapped_seed == 0u,
        "out-of-domain timing retry retained an output seed");
    FrozenTimingCell rejected = retry_255;
    Check(!QualifyDevelopmentTimingCell(bases[0], 256u, rejected) &&
            CanonicalTimingCellJson(rejected).empty() &&
            rejected.base_cell_sha256.empty() && rejected.loss_seed == 0u,
        "out-of-domain timing retry retained a qualified cell");
    FrozenTimingBaseCell invalid_base = bases[0];
    invalid_base.ordinal = 1u;
    rejected = retry_255;
    Check(!QualifyDevelopmentTimingCell(invalid_base, 0u, rejected) &&
            rejected.base_cell_sha256.empty(),
        "invalid timing base retained a qualified cell");

    const std::string source_json = CanonicalTimingBaseCellJson(bases[0]);
    const std::string source_sha256 = TimingBaseCellSha256(bases[0]);
    Check(CanonicalTimingSourceJson(cells[0]) == source_json &&
            CanonicalTimingSourceJson(retry_255) == source_json &&
            TimingSourceIdentitySha256(cells[0]) == source_sha256 &&
            TimingSourceIdentitySha256(retry_255) == source_sha256,
        "loss retry changed the deterministic timing source identity");
    Check(CanonicalTimingCellJson(cells[0]) !=
            CanonicalTimingCellJson(retry_255) &&
            TimingCellSha256(cells[0]) != TimingCellSha256(retry_255),
        "loss retry did not change the qualified timing identity");

    Check(DevelopmentTimingDomainSha256(retry_offsets) ==
            "4d0a2423b4fe54783e034860c0044c27"
            "1ae26cc0157300e208984032a3773726",
        "retry-zero qualified timing-domain SHA-256 changed");
    retry_offsets[0] = 2u;
    Check(DevelopmentTimingDomainSha256(retry_offsets) ==
            "7e622a31d82ccd80b93b2ed4118ef019"
            "3af10742d6c8bcaa98508b6ca4e509e7",
        "single-retry qualified timing-domain SHA-256 changed");
    retry_offsets.assign(2304u, 0u);
    retry_offsets[0] = 255u;
    retry_offsets[2303] = 1u;
    Check(DevelopmentTimingDomainSha256(retry_offsets) ==
            "cdfe34f7626c997567ea14c28b76c324"
            "a99068337096c0d223d1e5c28294e1bb",
        "wrapped-retry qualified timing-domain SHA-256 changed");

    FrozenTimingCell mutant = cells[0];
    mutant.base_cell_sha256[0] = '0';
    Check(CanonicalTimingCellJson(mutant).empty() &&
            CanonicalTimingSourceJson(mutant).empty(),
        "forged timing base-cell identity was accepted");
    mutant = cells[0];
    mutant.loss_retry_offset = 1u;
    Check(CanonicalTimingCellJson(mutant).empty(),
        "retry offset inconsistent with its loss seed was accepted");
    mutant = cells[0];
    mutant.loss_seed ^= 1u;
    Check(CanonicalTimingCellJson(mutant).empty(),
        "wrong realized timing loss seed was accepted");
    mutant = cells[0];
    mutant.base_loss_seed ^= 1u;
    Check(CanonicalTimingCellJson(mutant).empty() &&
            TimingSourceIdentitySha256(mutant).empty(),
        "qualified cell with a forged base loss seed was accepted");
}

void TestAllTimingTracesAndReceipts()
{
    const std::vector<uint32_t> retry_offsets(2304u, 0u);
    const std::vector<FrozenTimingCell> cells =
        EnumerateDevelopmentTimingCells(retry_offsets);
    if (!Check(cells.size() == 2304u,
            "cannot test traces without a qualified timing domain"))
    {
        return;
    }
    std::string trace_aggregate;
    std::string pair_aggregate;
    trace_aggregate.reserve(170000u);
    pair_aggregate.reserve(310000u);
    uint64_t total_attempts = 0u;
    uint64_t maximum_attempts = 0u;

    for (std::size_t i = 0u; i < cells.size(); ++i)
    {
        FrozenPacketTrace trace;
        FrozenTimingTraceReceipt receipt;
        Check(GenerateDevelopmentTimingTrace(cells[i], trace, receipt) ==
                FrozenTraceStatus::Complete,
            "frozen development timing trace did not complete");
        Check(trace.schedule == FrozenSchedule::Iid &&
                trace.delivered_ids.size() ==
                    static_cast<std::size_t>(cells[i].K) + 256u,
            "timing trace did not use one IID K+256 prefix");
        Check(receipt.ordinal == i && receipt.K == cells[i].K &&
                receipt.cell_sha256 == TimingCellSha256(cells[i]) &&
                receipt.trace_sha256 == trace.trace_sha256 &&
                receipt.attempted_candidates ==
                    trace.attempted_candidates &&
                receipt.candidate_limit == trace.candidate_limit,
            "timing trace receipt does not bind native trace fields");
        Check(!CanonicalTimingTraceManifestRow(cells[i], receipt).empty(),
            "valid timing trace manifest row was rejected");

        trace_aggregate += trace.trace_sha256;
        trace_aggregate += ':';
        trace_aggregate += std::to_string(trace.attempted_candidates);
        trace_aggregate += '\n';
        pair_aggregate += receipt.cell_sha256;
        pair_aggregate += ':';
        pair_aggregate += receipt.trace_sha256;
        pair_aggregate += '\n';
        total_attempts += trace.attempted_candidates;
        maximum_attempts = std::max(
            maximum_attempts, trace.attempted_candidates);

        if (i == 0u)
        {
            const std::vector<uint32_t> expected_prefix = {
                0u, 1u, 3u, 4u, 5u, 7u,
                9u, 10u, 11u, 12u, 13u, 14u
            };
            Check(trace.trace_seed == UINT64_C(0x0aa53e4b248d085a) &&
                    trace.delivered_ids.size() == 264u &&
                    std::equal(
                        expected_prefix.begin(), expected_prefix.end(),
                        trace.delivered_ids.begin()) &&
                    trace.attempted_candidates == 305u &&
                    trace.trace_sha256 ==
                        "768c7e52d3bfc315b9d9eb44287766f9"
                        "401904062a5595be439f8ba6b20c90aa",
                "first native timing trace golden changed");
            Check(
                CanonicalTimingTraceManifestRow(cells[i], receipt) ==
                    "{\"cell_sha256\":\"e85871b9c09ac946f6c8f6e0e077bd75"
                    "8d4862bfb5dddd73f6ee05e3ba35f2d7\",\"ordinal\":0,"
                    "\"trace_sha256\":\"768c7e52d3bfc315b9d9eb44287766f9"
                    "401904062a5595be439f8ba6b20c90aa\"}",
                "first canonical timing trace manifest row changed");
        }
        if (i == 2303u) {
            Check(trace.trace_sha256 ==
                    "78aa1d6d24567b3d8787a1de516842c8"
                    "baa8a482a8606934892bd9d711763516" &&
                    trace.attempted_candidates == 71235u,
                "last native timing trace golden changed");
        }
    }

    Check(total_attempts == UINT64_C(29197317),
        "all-timing-trace candidate total changed");
    Check(maximum_attempts == UINT64_C(71646),
        "all-timing-trace maximum candidate count changed");
    Check(
        Sha256Hex(trace_aggregate) ==
            "ee77fce6265c7610865ebe2b3a6e2812"
            "461216cf98ab7ddea077022bfc284be4",
        "all 2304 timing trace hashes disagree with Python reference");
    Check(
        Sha256Hex(pair_aggregate) ==
            "f863671a31aa63ca3b6123d8d1f1f08b"
            "9f864e0db328baa288c61d52d23a7a46",
        "timing cell-to-trace mapping disagrees with Python reference");

    FrozenPacketTrace first_trace;
    FrozenTimingTraceReceipt first_receipt;
    FrozenPacketTrace repeated_trace;
    FrozenTimingTraceReceipt repeated_receipt;
    Check(GenerateDevelopmentTimingTrace(
            cells[0], first_trace, first_receipt) ==
                FrozenTraceStatus::Complete &&
            GenerateDevelopmentTimingTrace(
                cells[0], repeated_trace, repeated_receipt) ==
                FrozenTraceStatus::Complete &&
            first_trace.delivered_ids == repeated_trace.delivered_ids &&
            CanonicalTimingTraceManifestRow(cells[0], first_receipt) ==
                CanonicalTimingTraceManifestRow(
                    cells[0], repeated_receipt),
        "repeated timing trace/receipt generation was not deterministic");

    const std::vector<FrozenTimingBaseCell> bases =
        EnumerateDevelopmentTimingBaseCells();
    FrozenTimingCell retry_two;
    FrozenPacketTrace retry_trace;
    FrozenTimingTraceReceipt retry_receipt;
    Check(bases.size() == 2304u &&
            QualifyDevelopmentTimingCell(bases[0], 2u, retry_two) &&
            GenerateDevelopmentTimingTrace(
                retry_two, retry_trace, retry_receipt) ==
                    FrozenTraceStatus::Complete &&
            retry_trace.delivered_ids != first_trace.delivered_ids &&
            retry_trace.trace_sha256 != first_trace.trace_sha256 &&
            retry_trace.trace_seed == UINT64_C(0x4ed40ab62516f034) &&
            retry_trace.attempted_candidates == 292u &&
            retry_trace.trace_sha256 ==
                "55c8e148ab0140d568efd6a64df3fb59"
                "157163420c089f5853d3acccaf02ee49" &&
            !CanonicalTimingTraceManifestRow(
                retry_two, retry_receipt).empty() &&
            TimingSourceIdentitySha256(retry_two) ==
                TimingSourceIdentitySha256(cells[0]),
        "retry-qualified trace did not vary independently of source bytes");

    FrozenPacketTrace recovery_compatible;
    FrozenPacketTrace explicit_four;
    Check(GenerateFrozenPacketTrace(
            cells[0].K, cells[0].block_bytes, cells[0].loss_ppm,
            cells[0].loss_seed, FrozenSchedule::Iid,
            recovery_compatible) == FrozenTraceStatus::Complete &&
            GenerateFrozenPacketTrace(
                cells[0].K, cells[0].block_bytes, cells[0].loss_ppm,
                cells[0].loss_seed, FrozenSchedule::Iid, 4u,
                explicit_four) == FrozenTraceStatus::Complete &&
            recovery_compatible.delivered_ids ==
                explicit_four.delivered_ids &&
            recovery_compatible.trace_sha256 == explicit_four.trace_sha256,
        "recovery-compatible K+4 trace overload changed bytes");
    FrozenPacketTrace excessive;
    Check(GenerateFrozenPacketTrace(
            cells[0].K, cells[0].block_bytes, cells[0].loss_ppm,
            cells[0].loss_seed, FrozenSchedule::Iid, 65537u,
            excessive) == FrozenTraceStatus::InvalidInput &&
            excessive.delivered_ids.empty(),
        "excessive explicit trace overhead was accepted");
    std::vector<uint32_t> nested;
    Check(!CopyNestedPrefix(first_trace, 4u, nested) && nested.empty(),
        "long timing trace was accepted as a K+4 recovery trace");

    FrozenTimingCell invalid = cells[0];
    invalid.replicate = 96u;
    Check(GenerateDevelopmentTimingTrace(
            invalid, repeated_trace, repeated_receipt) ==
                FrozenTraceStatus::InvalidInput &&
            repeated_trace.delivered_ids.empty() &&
            repeated_receipt.cell_sha256.empty(),
        "invalid timing cell retained prior trace or receipt state");

    FrozenTimingTraceReceipt forged = first_receipt;
    forged.candidate_limit += 1u;
    Check(CanonicalTimingTraceManifestRow(cells[0], forged).empty(),
        "wrong timing trace candidate cap was accepted");
    forged = first_receipt;
    forged.trace_sha256[0] = 'A';
    Check(CanonicalTimingTraceManifestRow(cells[0], forged).empty(),
        "noncanonical timing trace SHA-256 was accepted");
    forged = first_receipt;
    forged.K = 32u;
    FrozenCandidateLimit(288u, forged.candidate_limit);
    Check(CanonicalTimingTraceManifestRow(cells[0], forged).empty(),
        "timing trace K inconsistent with its ordinal was accepted");
    forged = first_receipt;
    forged.attempted_candidates = 1u;
    Check(CanonicalTimingTraceManifestRow(cells[0], forged).empty(),
        "forged timing trace candidate count was accepted");

    FrozenTimingCell wrong_cell = cells[1];
    Check(CanonicalTimingTraceManifestRow(wrong_cell, first_receipt).empty(),
        "timing trace receipt was accepted for a different qualified cell");
}

void TestOneCandidatePanelsAndOrders()
{
    const std::vector<FrozenTimingPanel> panels =
        EnumerateOneCandidateTimingPanels("candidate");
    Check(panels.size() == 11u,
        "one-candidate timing panel count is not 11");

    static const FrozenPanelKind kinds[] = {
        FrozenPanelKind::AA, FrozenPanelKind::AA,
        FrozenPanelKind::AA, FrozenPanelKind::AA,
        FrozenPanelKind::AA, FrozenPanelKind::AA, FrozenPanelKind::AA,
        FrozenPanelKind::AB, FrozenPanelKind::AB,
        FrozenPanelKind::AB, FrozenPanelKind::AB
    };
    static const FrozenTimingScope scopes[] = {
        FrozenTimingScope::DecoderSolve,
        FrozenTimingScope::ReceiveToSuccess,
        FrozenTimingScope::EncoderInitPlusFirstKSymbols,
        FrozenTimingScope::EncoderInitPlusFirstKSymbols,
        FrozenTimingScope::DecoderSolve,
        FrozenTimingScope::ReceiveToSuccess,
        FrozenTimingScope::EncoderInitPlusFirstKSymbols,
        FrozenTimingScope::DecoderSolve,
        FrozenTimingScope::ReceiveToSuccess,
        FrozenTimingScope::EncoderInitPlusFirstKSymbols,
        FrozenTimingScope::EncoderInitPlusFirstKSymbols
    };
    static const char* hashes[] = {
        "c9352024e96b76166e65de682d3b8fa3"
        "589cd3a8d8c78b5203a4b0f7b568af7f",
        "f29e96d72671cf1a81f40ee6bca585c8"
        "9c57a4dbf3f1372e147933104ad9a77b",
        "69863a0140f2a8f152597cecfd8e56f9"
        "8c664ccb5a8aa5f7b09649aee6474ef8",
        "7284f7f8daceea2ed821bf572a95ba54"
        "bea0315f13c414e2ed45554f5bb3b662",
        "9d122b393d66bdc4ebaacfad2385157a"
        "928fb20478a7e1e41320ec7246325dc0",
        "7f019413efc18c82623b2daa9f0c80bb"
        "23af074cdec48e221f56770217d2b8d6",
        "b5126326e270c4893ca2854cb7b8f30f"
        "cda2aa9dd5ce2d1bb158b51eef1070ed",
        "996031fe4f0403f2eafbdf603a2d77f5"
        "9f47f602d45628d62ea17a7777ee01e4",
        "2bee7f4d4926bd723933978315ff4b6e7"
        "e9d76fc502e70807b2174dd3e0fa06f",
        "dac88e2090f4c2b8ad2ae4d25addbc6"
        "3b74d416dfc23bb703c4512939989ca60",
        "00d0fc48ad9e784a089ab97a954d55c0"
        "5a2195e6a0f292915c9c99e15a4a24f5"
    };
    static const FrozenTimingOrder even_orders[] = {
        FrozenTimingOrder::BAAB, FrozenTimingOrder::BAAB,
        FrozenTimingOrder::ABBA, FrozenTimingOrder::ABBA,
        FrozenTimingOrder::ABBA, FrozenTimingOrder::ABBA,
        FrozenTimingOrder::BAAB, FrozenTimingOrder::ABBA,
        FrozenTimingOrder::BAAB, FrozenTimingOrder::ABBA,
        FrozenTimingOrder::BAAB
    };

    std::vector<std::string> canonical_panels;
    for (std::size_t i = 0u; i < panels.size(); ++i)
    {
        Check(panels[i].ordinal == i &&
                panels[i].panel_kind == kinds[i] &&
                panels[i].scope == scopes[i],
            "one-candidate panel ordering changed");
        const std::string canonical = CanonicalTimingPanelJson(panels[i]);
        Check(!canonical.empty(), "valid timing panel was rejected");
        Check(std::find(
                canonical_panels.begin(), canonical_panels.end(), canonical) ==
                    canonical_panels.end(),
            "duplicate canonical timing panel was generated");
        canonical_panels.push_back(canonical);
        Check(TimingPanelSha256(panels[i]) == hashes[i],
            "canonical timing panel SHA-256 changed");
        Check(TimingPanelOrder(panels[i], 0u) == even_orders[i],
            "replicate-zero ABBA/BAAB assignment changed");
        const FrozenTimingOrder odd_order =
            even_orders[i] == FrozenTimingOrder::ABBA ?
                FrozenTimingOrder::BAAB : FrozenTimingOrder::ABBA;
        Check(TimingPanelOrder(panels[i], 1u) == odd_order &&
                TimingPanelOrder(panels[i], 2u) == even_orders[i],
            "timing order does not alternate solely by replicate parity");
    }

    Check(
        CanonicalTimingPanelJson(panels[0]) ==
            "{\"left_arm\":\"wirehair2_head\",\"panel_kind\":\"AA\","
            "\"right_arm\":\"wirehair2_head\",\"scope\":\"decoder_solve\"}",
        "first canonical control panel changed");
    Check(
        CanonicalTimingPanelJson(panels[4]) ==
            "{\"left_arm\":\"candidate\",\"panel_kind\":\"AA\","
            "\"right_arm\":\"candidate\",\"scope\":\"decoder_solve\"}",
        "first canonical candidate panel changed");

    const std::vector<FrozenTimingBaseCell> cells =
        EnumerateDevelopmentTimingBaseCells();
    uint32_t abba = 0u;
    uint32_t baab = 0u;
    for (std::size_t cell = 0u; cell < cells.size(); ++cell) {
        for (std::size_t panel = 0u; panel < panels.size(); ++panel) {
            if (TimingPanelOrder(panels[panel], cells[cell].replicate) ==
                FrozenTimingOrder::ABBA)
            {
                ++abba;
            }
            else {
                ++baab;
            }
        }
    }
    Check(abba == 12672u && baab == 12672u,
        "25,344 timing rows are not parity-counterbalanced");

    FrozenTimingPanel mutant = panels[0];
    mutant.ordinal = 1u;
    Check(CanonicalTimingPanelJson(mutant).empty() &&
            TimingPanelOrder(mutant, 0u) == FrozenTimingOrder::Invalid,
        "relabeled timing panel was accepted");
    Check(EnumerateOneCandidateTimingPanels("").empty() &&
            EnumerateOneCandidateTimingPanels("wirehair1").empty() &&
            EnumerateOneCandidateTimingPanels("wirehair2_head").empty(),
        "empty or reserved candidate arm was accepted");

    const std::string unicode_candidate("caf\xc3\xa9");
    const std::vector<FrozenTimingPanel> unicode =
        EnumerateOneCandidateTimingPanels(unicode_candidate);
    Check(unicode.size() == 11u,
        "valid UTF-8 candidate arm was rejected");
    if (unicode.size() == 11u) {
        Check(
            CanonicalTimingPanelJson(unicode[4]) ==
                "{\"left_arm\":\"caf\\u00e9\",\"panel_kind\":\"AA\","
                "\"right_arm\":\"caf\\u00e9\",\"scope\":\"decoder_solve\"}" &&
            TimingPanelSha256(unicode[4]) ==
                "1bf5e33e6da6f261ff00749e5e721cb"
                "3714f76f1fd014c6e1315ed0b7ed879fb",
            "UTF-8 candidate differs from Python ensure_ascii JSON");
    }
    const std::string emoji_candidate("emoji\xf0\x9f\x98\x80");
    const std::vector<FrozenTimingPanel> emoji =
        EnumerateOneCandidateTimingPanels(emoji_candidate);
    Check(emoji.size() == 11u &&
            CanonicalTimingPanelJson(emoji[4]) ==
                "{\"left_arm\":\"emoji\\ud83d\\ude00\","
                "\"panel_kind\":\"AA\",\"right_arm\":"
                "\"emoji\\ud83d\\ude00\",\"scope\":\"decoder_solve\"}",
        "non-BMP UTF-8 candidate did not use Python surrogate escaping");
    std::string malformed_utf8("candidate");
    malformed_utf8 += static_cast<char>(0xc3);
    Check(EnumerateOneCandidateTimingPanels(malformed_utf8).empty(),
        "malformed UTF-8 candidate arm was accepted");

    std::string escaped_candidate("cand\"\\");
    escaped_candidate += '\n';
    escaped_candidate += static_cast<char>(1);
    const std::vector<FrozenTimingPanel> escaped =
        EnumerateOneCandidateTimingPanels(escaped_candidate);
    Check(escaped.size() == 11u,
        "ASCII candidate requiring JSON escaping was rejected");
    if (escaped.size() == 11u) {
        Check(
            CanonicalTimingPanelJson(escaped[4]) ==
                "{\"left_arm\":\"cand\\\"\\\\\\n\\u0001\","
                "\"panel_kind\":\"AA\","
                "\"right_arm\":\"cand\\\"\\\\\\n\\u0001\","
                "\"scope\":\"decoder_solve\"}" &&
            TimingPanelSha256(escaped[4]) ==
                "e38f07f1d5673ed3edfbacafc600890"
                "6d2e556e4b5e7b0bc11a342e422a82ebb",
            "escaped candidate panel differs from Python canonical JSON");
    }
}

} // namespace

int main()
{
    TestTimingSeeds();
    TestInvocationsPerSlotFormula();
    TestExecutionGeometry();
    TestDefaultTimingValuesAreSafelyInvalid();
    TestTimingBaseCellEnumeration();
    TestTimingQualificationAndSourceIdentity();
    TestAllTimingTracesAndReceipts();
    TestOneCandidatePanelsAndOrders();

    if (g_failures != 0u) {
        std::fprintf(stderr,
            "WH2 frozen timing test: %u failure(s)\n", g_failures);
        return 1;
    }
    std::printf("WH2 frozen timing test: PASS\n");
    return 0;
}
