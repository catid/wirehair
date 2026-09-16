struct Row {
    unsigned rep, cell, pair, order, position, arm;
    uint64_t ns;
    unsigned steps;
    std::uintptr_t source, output, descriptor, descriptor_out, low, distant;
};
std::vector<Row> rows;
bool published = false;
void Publish()
{
    if (published) return;
    published = true;
    std::puts("rep,cell,pair,order,position,arm,ns,steps,source,output,descriptor,descriptor_out,low,distant");
    for (const auto& r : rows)
        std::printf("%u,%u,%u,%u,%u,%u,%llu,%u,%llu,%llu,%llu,%llu,%llu,%llu\n",
            r.rep,r.cell,r.pair,r.order,r.position,r.arm,static_cast<unsigned long long>(r.ns),r.steps,
            static_cast<unsigned long long>(r.source),static_cast<unsigned long long>(r.output),
            static_cast<unsigned long long>(r.descriptor),static_cast<unsigned long long>(r.descriptor_out),
            static_cast<unsigned long long>(r.low),static_cast<unsigned long long>(r.distant));
    Check(!std::ferror(stdout) && std::fflush(stdout) == 0,"publish");
}
} // namespace

int main(int argc, char** argv)
{
    try {
        Check(argc == 4,"arguments");
        const std::string mode = argv[3];
        Check(mode == "neutral" || mode == "neutral-reverse" || mode == "run" || mode == "run-reverse","mode");
        const bool timing = mode == "run" || mode == "run-reverse";
        const bool reverse = mode == "run-reverse" || mode == "neutral-reverse";
        const auto start = Now();
        if (timing) {
            cpu_set_t mask; CPU_ZERO(&mask); CPU_SET(50,&mask);
            Check(sched_setaffinity(0,sizeof(mask),&mask) == 0 && sched_getaffinity(0,sizeof(mask),&mask) == 0 &&
                  CPU_COUNT(&mask) == 1 && CPU_ISSET(50,&mask) && sched_getcpu() == 50,"CPU50 singleton");
        }
        Library first(argv[reverse ? 2 : 1]), second(argv[reverse ? 1 : 2]);
        Check(first.lib != second.lib,"distinct libraries");
        Library* baseline = reverse ? &second : &first;
        Library* candidate = reverse ? &first : &second;
        Api apis[3] = {{baseline,false},{candidate,false},{baseline,true}};
        std::vector<std::array<Fixture,3>> fixtures;
        for (unsigned route = 0; route < 6; ++route) for (unsigned variant = 0; variant < 5; ++variant)
            for (unsigned policy = 0; policy < 2; ++policy) {
                const unsigned ks[6] = {3,3,5,8,3,16}, widths[5] = {2,2,64,1280,1280};
                const unsigned k = ks[route], width = widths[variant];
                Shape s = {k,width,variant == 1 || variant == 4 ? 1u : width,policy,route};
                std::array<Fixture,3> f;
                for (unsigned arm = 0; arm < 3; ++arm) {
                    f[arm] = Prepare(apis[arm],s);
                    if (arm < 2) {
                        WirehairV2Profile profile = {};
                        Check(apis[arm].lib->parse(f[arm].profile.data(),32,&profile) == 0 &&
                              profile.profile_id == (route >= 4 ? WIREHAIR_V2_PROFILE_CERTIFIED_2026_07 : ProfileId(k)) &&
                              profile.message_bytes == s.Message() && profile.block_bytes == width,"exact selected route");
                    }
                }
                Check(f[0].profile == f[1].profile && f[0].source == f[1].source && f[0].source == f[2].source &&
                      f[0].packets == f[1].packets && f[0].repairs == f[1].repairs &&
                      f[0].steps[0] == f[1].steps[0] && f[0].steps[1] == f[1].steps[1],"A/B equation parity");
                fixtures.push_back(f);
            }
        Check(fixtures.size() == 60,"fixture roster");
        const unsigned pairs[5][2] = {{0,0},{1,1},{2,2},{0,1},{2,1}};
        const unsigned sides[18] = {0,1,0,1,1,0,1,0,0,1,1,0,0,1,0,1,1,0};
        if (timing) rows.reserve(Reps*180*5*2*18);
        for (unsigned rep = 0; rep < (timing ? Reps : 1u); ++rep)
            for (unsigned slot = 0; slot < 180; ++slot) {
                const unsigned cell = (slot+rep*17)%180, metric = cell%3+1;
                const auto& fs = fixtures[cell/3];
                Fixture matched = fs[0];
                Profile produced = {};
                const std::size_t bytes = Batch*std::size_t(3*matched.shape.k)*matched.shape.b;
                std::vector<Byte> output(bytes+32,173);
                for (unsigned ps = 0; ps < (timing ? 5u : 1u); ++ps)
                    for (unsigned order = 0; order < (timing ? 2u : 1u); ++order) {
                        const unsigned pair = (ps+rep+cell)%5;
                        for (unsigned pos = 0; pos < (timing ? 18u : 3u); ++pos) {
                            const unsigned side = sides[pos]^order, arm = timing ? pairs[pair][side] : pos;
                            matched.profile = fs[arm].profile;
                            for (unsigned family = 0; family < 2; ++family) {
                                Check(matched.repairs[family].size() == fs[arm].repairs[family].size(),"staging extent");
                                std::copy(fs[arm].repairs[family].begin(),fs[arm].repairs[family].end(),matched.repairs[family].begin());
                                matched.steps[family] = fs[arm].steps[family];
                            }
                            std::fill(output.begin()+16,output.end()-16,0);
                            const auto begin = Now();
                            Work(apis[arm],matched,metric,nullptr,output.data()+16,produced);
                            const auto end = Now();
                            if (timing) rows.push_back(Row{rep,cell,pair,order,pos,arm,end-begin,
                                metric == 1 ? 0u : matched.steps[metric-2],
                                reinterpret_cast<std::uintptr_t>(matched.source.data()),
                                reinterpret_cast<std::uintptr_t>(output.data()+16),
                                reinterpret_cast<std::uintptr_t>(matched.profile.data()),
                                reinterpret_cast<std::uintptr_t>(produced.data()),
                                reinterpret_cast<std::uintptr_t>(matched.repairs[0].data()),
                                reinterpret_cast<std::uintptr_t>(matched.repairs[1].data())});
                            Verify(fs[arm],metric,output);
                            Check(Now()-start < UINT64_C(240000000000),"240-second worker cap");
                        }
                    }
            }
        if (!timing) { std::puts("PASS 60 fixtures x 3 lifecycle metrics x 3 APIs; descriptor/packet/first-success parity"); return 0; }
        Check(rows.size() == Reps*180*5*2*18,"complete cohort");
        Publish();
    } catch (const std::exception& e) {
        std::fprintf(stderr,"FAIL: %s\n",e.what());
        try { if (!rows.empty()) Publish(); } catch (...) {}
        return 1;
    }
}
