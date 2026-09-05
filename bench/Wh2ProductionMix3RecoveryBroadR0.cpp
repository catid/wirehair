// Additive fresh-root pilot: reuse the reviewed public API/oracle implementation.
// The included frozen r0 entry point is never reachable from this executable.
#define main wh2_recovery_r0_unused_main
#include "Wh2ProductionMix3RecoveryR0.cpp"
#undef main

namespace {
const char kBroadProtocol[] = "wirehair.wh2.production-mix3-k3k6-broad-r0";

std::vector<uint64_t> BroadRoots(const std::string& phase)
{
    Require(phase == "train" || phase == "holdout", "root phase");
    std::vector<uint64_t> roots;
    const unsigned count = phase == "train" ? 16u : 64u;
    for (unsigned i = 0; i < count; ++i) {
        const std::string hash = Sha256Hex(std::string(kBroadProtocol) + ":" +
            phase + "/" + std::to_string(i));
        Require(hash.size() == 64u, "root hash");
        roots.push_back(std::stoull(hash.substr(0u, 16u), nullptr, 16));
    }
    return roots;
}
void BroadBegin(const std::string& phase)
{
    std::cout << "{\"type\":\"begin\",\"protocol\":" << Quoted(kBroadProtocol)
        << ",\"phase\":" << Quoted(phase) << ",\"source_commit\":"
        << Quoted(WIREHAIR_WH2_SOURCE_GIT_COMMIT) << ",\"library_path\":"
        << Quoted(LibraryPath()) << "}\n";
}
void BroadSelftest()
{
    // Input derivation only: none of the fresh roster is decoded by this test.
    const std::vector<uint64_t> train = BroadRoots("train"), holdout = BroadRoots("holdout");
    Require(train.size() == 16u && train.front() == UINT64_C(0x77e5349d4b9f2c07) &&
        train.back() == UINT64_C(0x7d0de1a454c28717) && holdout.size() == 64u &&
        holdout.front() == UINT64_C(0x25b1e0cb2091ddb3) &&
        holdout.back() == UINT64_C(0x752bfbe6effcc2f2), "broad root selftest");
    Selftest();
}
void BroadRun(const std::string& phase, uint32_t seconds,
    const std::array<uint32_t, 2>& map)
{
    const Clock::time_point deadline = Clock::now() + std::chrono::seconds(seconds);
    Require(wirehair_init() == Wirehair_Success, "initialization");
    const std::vector<uint64_t> roots = BroadRoots(phase);
    BroadBegin(phase);
    uint32_t rows = 0u;
    std::array<int, 2> selected = {{-1, -1}};
    for (size_t k = 0; k < kKs.size(); ++k) {
        uint32_t actual = 0u;
        EmitArm(phase, "baseline", kKs[k], -1, roots, deadline, rows, actual);
        if (phase == "train") {
            for (uint32_t attempt = 0u; attempt < 256u; ++attempt) {
                if (EmitArm(phase, "candidate", kKs[k], static_cast<int>(attempt),
                        roots, deadline, rows, actual)) {
                    selected[k] = static_cast<int>(attempt);
                    break;
                }
            }
        } else {
            selected[k] = static_cast<int>(map[k]);
            EmitArm(phase, "candidate", kKs[k], selected[k], roots, deadline, rows, actual);
        }
    }
    Guard(deadline);
    std::cout << "{\"type\":\"terminal\",\"phase\":" << Quoted(phase)
        << ",\"rows\":" << rows << ",\"selected_attempts\":[" << selected[0] << ','
        << selected[1] << "]}\n";
    std::cout.flush();
    Require(std::cout.good(), "terminal stream");
}
} // namespace

int main(int argc, char** argv)
{
    try {
        if (argc == 2 && std::string(argv[1]) == "--describe") { BroadBegin("describe"); return 0; }
        if (argc == 2 && std::string(argv[1]) == "--selftest") { BroadSelftest(); return 0; }
        if (argc == 3 && std::string(argv[1]) == "--train") {
            const uint32_t seconds = Number(argv[2], 60u);
            Require(seconds != 0u, "zero deadline");
            BroadRun("train", seconds, {{0u, 0u}}); return 0;
        }
        if (argc == 5 && std::string(argv[1]) == "--holdout") {
            const uint32_t seconds = Number(argv[2], 60u);
            Require(seconds != 0u, "zero deadline");
            BroadRun("holdout", seconds, {{Number(argv[3], 255u), Number(argv[4], 255u)}}); return 0;
        }
        throw std::runtime_error("usage: --describe | --selftest | --train seconds | --holdout seconds a3 a6");
    } catch (const std::exception& error) {
        std::cerr << "production MIX3 recovery broad r0: " << error.what() << '\n';
        return 2;
    }
}
