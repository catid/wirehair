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
    Check(!DevelopmentTimingSeed(16u, attempt, seed) &&
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
            DevelopmentTimingInvocationsPerSlot(65536u) == 1u &&
            DevelopmentTimingInvocationsPerSlot(65537u) == 1u &&
            DevelopmentTimingInvocationsPerSlot(UINT32_MAX) == 1u,
        "batch-count boundary or overflow handling changed");
}

void TestExecutionGeometry()
{
    Check(DevelopmentTimingWorkerCount() == 8u &&
            DevelopmentTimingCohortCount() == 264u,
        "frozen timing worker/cohort cardinality changed");
    const std::vector<FrozenTimingCell> cells =
        EnumerateDevelopmentTimingCells();
    const std::vector<FrozenTimingPanel> panels =
        EnumerateOneCandidateTimingPanels("candidate");
    Check(cells.size() == 384u && panels.size() == 11u &&
            cells.size() / 16u * panels.size() ==
                DevelopmentTimingCohortCount(),
        "timing cohort count does not match cells and panels");

    for (std::size_t cohort = 0u;
        cohort < DevelopmentTimingCohortCount();
        ++cohort)
    {
        for (uint32_t wave = 0u; wave < 2u; ++wave)
        {
            bool seen[8] = {};
            for (uint32_t local = 0u; local < 8u; ++local)
            {
                const uint32_t replicate = wave * 8u + local;
                uint32_t slot = UINT32_MAX;
                const uint32_t expected = static_cast<uint32_t>(
                    (local + cohort + wave) % 8u);
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

    for (uint32_t replicate = 0u; replicate < 16u; ++replicate)
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
    Check(!DevelopmentTimingWorkerSlot(0u, 16u, slot) && slot == 0u,
        "out-of-domain timing replicate received a worker slot");
}

void TestDefaultTimingValuesAreSafelyInvalid()
{
    const FrozenTimingCell cell;
    Check(CanonicalTimingCellJson(cell).empty() && TimingCellSha256(cell).empty(),
        "default timing cell was not safely rejected");
    const FrozenTimingPanel panel;
    Check(CanonicalTimingPanelJson(panel).empty() &&
            TimingPanelSha256(panel).empty() &&
            TimingPanelOrder(panel, 0u) == FrozenTimingOrder::Invalid,
        "default timing panel was not safely rejected");
}

void TestTimingCellEnumeration()
{
    static const uint32_t K_values[] = {
        8u, 32u, 100u, 128u, 512u, 1000u,
        2048u, 5000u, 8192u, 20000u, 32768u, 64000u
    };
    static const uint32_t widths[] = { 64u, 1280u };

    const std::vector<FrozenTimingCell> cells =
        EnumerateDevelopmentTimingCells();
    Check(cells.size() == 384u,
        "development timing cell cardinality is not 384");
    for (std::size_t ordinal = 0u; ordinal < cells.size(); ++ordinal)
    {
        const FrozenTimingCell& cell = cells[ordinal];
        const std::size_t width_index = ordinal / 192u;
        const uint32_t replicate =
            static_cast<uint32_t>((ordinal % 192u) / 12u);
        const std::size_t k_index = ordinal % 12u;
        uint32_t expected_attempt = 0u;
        uint64_t expected_seed = 0u;
        Check(DevelopmentTimingSeed(
                replicate, expected_attempt, expected_seed),
            "enumerated timing replicate has no paired seed");
        Check(cell.ordinal == ordinal,
            "development timing ordinal changed");
        Check(cell.block_bytes == widths[width_index],
            "development timing width ordering changed");
        Check(cell.replicate == replicate,
            "development timing replicate ordering changed");
        Check(cell.K == K_values[k_index],
            "development timing K ordering changed");
        Check(cell.base_seed_attempt == expected_attempt &&
                cell.loss_seed == expected_seed,
            "development timing seed pairing changed");
        Check(cell.phase == "development" &&
                cell.loss_ppm == 100000u &&
                cell.schedule == FrozenSchedule::Iid &&
                cell.fixed_received_overhead == 4u &&
                cell.invocations_per_slot ==
                    DevelopmentTimingInvocationsPerSlot(cell.K),
            "development timing fixed field changed");
    }

    const std::string first_json =
        "{\"K\":8,\"band\":\"2-100\",\"base_seed_attempt\":0,"
        "\"block_bytes\":64,\"fixed_received_overhead\":4,"
        "\"invocations_per_slot\":8192,"
        "\"loss_ppm\":100000,\"loss_seed\":\"0x2d0f28c7e7e786b2\","
        "\"phase\":\"development\",\"replicate\":0,"
        "\"schedule\":\"iid\"}";
    Check(CanonicalTimingCellJson(cells[0]) == first_json,
        "first canonical timing cell changed");
    Check(
        TimingCellSha256(cells[0]) ==
            "1a059880adb0a19f5c0d1d457df25447"
            "9688a6fa89677503e7d12cc4736f0a9f",
        "first timing cell SHA-256 changed");
    Check(
        TimingCellSha256(cells[192]) ==
            "66f9fe3a15beb8a6284e3609a79c4743"
            "c73fc34f698cf3720389633956104ee9",
        "first wide timing cell SHA-256 changed");
    Check(
        TimingCellSha256(cells[383]) ==
            "08d5ab5c41e778574fdc1c57d8284780"
            "b65242c2446a57808c72245987be94f5",
        "last timing cell SHA-256 changed");
    Check(
        DevelopmentTimingDomainSha256() ==
            "4d094a6d351d3b57cbc2654bf269b1da"
            "41886a68c7bebd14803b2892d1cd89d7",
        "development timing domain SHA-256 changed");

    const std::vector<FrozenTimingCell> repeat =
        EnumerateDevelopmentTimingCells();
    Check(repeat.size() == cells.size(),
        "repeated timing enumeration changed cardinality");
    for (std::size_t i = 0u; i < cells.size() && i < repeat.size(); ++i) {
        Check(CanonicalTimingCellJson(cells[i]) ==
                CanonicalTimingCellJson(repeat[i]),
            "repeated timing enumeration changed canonical bytes");
    }

    FrozenTimingCell mutant = cells[0];
    mutant.ordinal = 1u;
    Check(CanonicalTimingCellJson(mutant).empty(),
        "wrong timing-cell ordinal was accepted");
    mutant = cells[0];
    mutant.schedule = FrozenSchedule::Burst;
    Check(CanonicalTimingCellJson(mutant).empty(),
        "non-IID development timing cell was accepted");
    mutant = cells[0];
    mutant.loss_seed ^= 1u;
    Check(CanonicalTimingCellJson(mutant).empty(),
        "wrong paired timing loss seed was accepted");
    mutant = cells[0];
    --mutant.invocations_per_slot;
    Check(CanonicalTimingCellJson(mutant).empty(),
        "wrong timing invocation batch was accepted");
    mutant = cells[0];
    mutant.base_seed_attempt = 1u;
    Check(CanonicalTimingCellJson(mutant).empty(),
        "non-production development timing attempt was accepted");
}

void TestAllTimingTracesAndReceipts()
{
    const std::vector<FrozenTimingCell> cells =
        EnumerateDevelopmentTimingCells();
    std::string trace_aggregate;
    std::string pair_aggregate;
    trace_aggregate.reserve(26800u);
    pair_aggregate.reserve(50000u);
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
                    static_cast<std::size_t>(cells[i].K) + 4u,
            "timing trace did not use one IID K+4 prefix");
        Check(receipt.ordinal == i && receipt.K == cells[i].K &&
                receipt.cell_sha256 == TimingCellSha256(cells[i]) &&
                receipt.trace_sha256 == trace.trace_sha256 &&
                receipt.attempted_candidates ==
                    trace.attempted_candidates &&
                receipt.candidate_limit == trace.candidate_limit,
            "timing trace receipt does not bind native trace fields");
        Check(!CanonicalTimingTraceManifestRow(receipt).empty(),
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
            const std::vector<uint32_t> expected_ids = {
                0u, 1u, 3u, 4u, 5u, 7u,
                9u, 10u, 11u, 12u, 13u, 14u
            };
            Check(trace.trace_seed == UINT64_C(0x0aa53e4b248d085a) &&
                    trace.delivered_ids == expected_ids &&
                    trace.attempted_candidates == 15u &&
                    trace.trace_sha256 ==
                        "0e773013a48e8a4aaa2858ca7eedec25"
                        "4bbef05c82c93d6e0e75bcf482e36e94",
                "first native timing trace golden changed");
            Check(
                CanonicalTimingTraceManifestRow(receipt) ==
                    "{\"cell_sha256\":\"1a059880adb0a19f5c0d1d457df25447"
                    "9688a6fa89677503e7d12cc4736f0a9f\",\"ordinal\":0,"
                    "\"trace_sha256\":\"0e773013a48e8a4aaa2858ca7eedec25"
                    "4bbef05c82c93d6e0e75bcf482e36e94\"}",
                "first canonical timing trace manifest row changed");
        }
        if (i == 383u) {
            Check(trace.trace_sha256 ==
                    "24c43a92dcc1705e46a772abbfc2bd09e"
                    "d491f1f9829972025f0dec9397fb123" &&
                    trace.attempted_candidates == 71131u,
                "last native timing trace golden changed");
        }
    }

    Check(total_attempts == UINT64_C(4759082),
        "all-timing-trace candidate total changed");
    Check(maximum_attempts == UINT64_C(71360),
        "all-timing-trace maximum candidate count changed");
    Check(
        Sha256Hex(trace_aggregate) ==
            "f326fec5818e4f6e15ec31fec161623a0"
            "96bce8c352657d119989ef6a9902161",
        "all 384 timing trace hashes disagree with Python reference");
    Check(
        Sha256Hex(pair_aggregate) ==
            "0fdb16de07899ee0ab4bf8e116efa120"
            "666856d77283c557af0f8e15460582d8",
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
            CanonicalTimingTraceManifestRow(first_receipt) ==
                CanonicalTimingTraceManifestRow(repeated_receipt),
        "repeated timing trace/receipt generation was not deterministic");

    FrozenTimingCell invalid = cells[0];
    invalid.replicate = 16u;
    Check(GenerateDevelopmentTimingTrace(
            invalid, repeated_trace, repeated_receipt) ==
                FrozenTraceStatus::InvalidInput &&
            repeated_trace.delivered_ids.empty() &&
            repeated_receipt.cell_sha256.empty(),
        "invalid timing cell retained prior trace or receipt state");

    FrozenTimingTraceReceipt forged = first_receipt;
    forged.candidate_limit += 1u;
    Check(CanonicalTimingTraceManifestRow(forged).empty(),
        "wrong timing trace candidate cap was accepted");
    forged = first_receipt;
    forged.trace_sha256[0] = 'A';
    Check(CanonicalTimingTraceManifestRow(forged).empty(),
        "noncanonical timing trace SHA-256 was accepted");
    forged = first_receipt;
    forged.K = 32u;
    FrozenCandidateLimit(36u, forged.candidate_limit);
    Check(CanonicalTimingTraceManifestRow(forged).empty(),
        "timing trace K inconsistent with its ordinal was accepted");
    forged = first_receipt;
    forged.attempted_candidates = 1u;
    Check(CanonicalTimingTraceManifestRow(forged).empty(),
        "forged timing trace candidate count was accepted");
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

    const std::vector<FrozenTimingCell> cells =
        EnumerateDevelopmentTimingCells();
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
    Check(abba == 2112u && baab == 2112u,
        "4,224 timing rows are not parity-counterbalanced");

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
    TestTimingCellEnumeration();
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
