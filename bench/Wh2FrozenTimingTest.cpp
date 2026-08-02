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
        0u, 1u, 2u, 0u, 1u, 2u, 0u, 1u
    };
    static const uint64_t seeds[] = {
        UINT64_C(0x2d0f28c7e7e786b2),
        UINT64_C(0x49c18e491a66eb93),
        UINT64_C(0xe93dcff7eed78181),
        UINT64_C(0xe9e7bf14d7c14038),
        UINT64_C(0x8886c61661bd0c88),
        UINT64_C(0x4694f3b073a0c6c7),
        UINT64_C(0x7d771d8bc3d7af5f),
        UINT64_C(0xe4db579c79a43708)
    };
    for (uint32_t replicate = 0u; replicate < 8u; ++replicate)
    {
        uint32_t attempt = UINT32_MAX;
        uint64_t seed = UINT64_MAX;
        Check(DevelopmentTimingSeed(replicate, attempt, seed),
            "valid development timing replicate was rejected");
        Check(attempt == attempts[replicate],
            "paired public construction attempt changed");
        Check(seed == seeds[replicate],
            "paired SplitMix64 timing loss seed changed");
    }

    uint32_t attempt = 9u;
    uint64_t seed = 9u;
    Check(!DevelopmentTimingSeed(8u, attempt, seed) &&
            attempt == 0u && seed == 0u,
        "out-of-domain timing replicate was accepted");
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
    Check(cells.size() == 192u,
        "development timing cell cardinality is not 192");
    for (std::size_t ordinal = 0u; ordinal < cells.size(); ++ordinal)
    {
        const FrozenTimingCell& cell = cells[ordinal];
        const std::size_t width_index = ordinal / 96u;
        const uint32_t replicate =
            static_cast<uint32_t>((ordinal % 96u) / 12u);
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
                cell.fixed_received_overhead == 4u,
            "development timing fixed field changed");
    }

    const std::string first_json =
        "{\"K\":8,\"band\":\"2-100\",\"base_seed_attempt\":0,"
        "\"block_bytes\":64,\"fixed_received_overhead\":4,"
        "\"loss_ppm\":100000,\"loss_seed\":\"0x2d0f28c7e7e786b2\","
        "\"phase\":\"development\",\"replicate\":0,"
        "\"schedule\":\"iid\"}";
    Check(CanonicalTimingCellJson(cells[0]) == first_json,
        "first canonical timing cell changed");
    Check(
        TimingCellSha256(cells[0]) ==
            "aa5c4c40257eacb6d01d44ef5ef9ac45"
            "e0e7c43a454b8b16ffa858364140a0a2",
        "first timing cell SHA-256 changed");
    Check(
        TimingCellSha256(cells[96]) ==
            "d663c2956a16c515da0b33ad684fe84d"
            "b5f67988b87c89cbfe5fde3d5f935bc7",
        "first wide timing cell SHA-256 changed");
    Check(
        TimingCellSha256(cells[191]) ==
            "56ea96402d6b7d39d5211fe0774491bb"
            "57f3a32a204271f159a37c0d6d726764",
        "last timing cell SHA-256 changed");
    Check(
        DevelopmentTimingDomainSha256() ==
            "1a15ee48e893280013b74f448d81607e"
            "45630fbf631d4e73960b3a2194e204c3",
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
}

void TestAllTimingTracesAndReceipts()
{
    const std::vector<FrozenTimingCell> cells =
        EnumerateDevelopmentTimingCells();
    std::string trace_aggregate;
    std::string pair_aggregate;
    trace_aggregate.reserve(13400u);
    pair_aggregate.reserve(25000u);
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
                    "{\"cell_sha256\":\"aa5c4c40257eacb6d01d44ef5ef9ac45"
                    "e0e7c43a454b8b16ffa858364140a0a2\",\"ordinal\":0,"
                    "\"trace_sha256\":\"0e773013a48e8a4aaa2858ca7eedec25"
                    "4bbef05c82c93d6e0e75bcf482e36e94\"}",
                "first canonical timing trace manifest row changed");
        }
        if (i == 191u) {
            Check(trace.trace_sha256 ==
                    "8829d12cd3ef85df09051bef9cad285e"
                    "965e1c2d45f6b51c43aacc551620dd20" &&
                    trace.attempted_candidates == 71181u,
                "last native timing trace golden changed");
        }
    }

    Check(total_attempts == UINT64_C(2378674),
        "all-timing-trace candidate total changed");
    Check(maximum_attempts == UINT64_C(71238),
        "all-timing-trace maximum candidate count changed");
    Check(
        Sha256Hex(trace_aggregate) ==
            "62dd9de7f685c7b4b373d8dea7603a45"
            "8525ba6aa99f7c8226b526ec245ec26d",
        "all 192 timing trace hashes disagree with Python reference");
    Check(
        Sha256Hex(pair_aggregate) ==
            "2f4f9fb10c2747715b053c3d1e7334e"
            "2aee42a24ddf7f0a40930b90334ac9fe8",
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
    invalid.replicate = 8u;
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
    Check(abba == 1056u && baab == 1056u,
        "2,112 timing rows are not parity-counterbalanced");

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
