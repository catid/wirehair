int main(int argc, char** argv)
{
    try {
        Check(argc == 4,"arguments");
        const std::string mode = argv[3];
        Check(mode == "neutral" || mode == "neutral-reverse" || mode == "run" || mode == "run-reverse","mode");
        const bool timing = mode == "run" || mode == "run-reverse";
        const bool reverse = mode == "run-reverse" || mode == "neutral-reverse";
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
        for (unsigned k : {3u,5u,8u,16u}) for (unsigned variant = 0; variant < (k == 16 ? 2u : 7u); ++variant)
            for (unsigned policy = 0; policy < 2; ++policy) {
                const unsigned widths[7] = {2,2,64,256,257,1280,1280};
                const unsigned width = k == 16 ? (variant ? 1280 : 64) : widths[variant];
                const unsigned tail = k != 16 && (variant == 1 || variant == 6) ? 1u : width;
                Shape s = {k,width,tail,policy};
                std::array<Fixture,3> f;
                for (unsigned arm = 0; arm < 3; ++arm) f[arm] = Prepare(apis[arm],s);
                Check(f[0].profile == f[1].profile && f[0].source == f[1].source && f[0].source == f[2].source &&
                      f[0].packets == f[1].packets && f[0].repairs == f[1].repairs &&
                      f[0].steps[0] == f[1].steps[0] && f[0].steps[1] == f[1].steps[1],"A/B equation parity");
                fixtures.push_back(f);
            }
        Check(fixtures.size() == 46,"fixture roster");
        const unsigned pairs[5][2] = {{0,0},{1,1},{2,2},{0,1},{2,1}};
        const unsigned sides[18] = {0,1,0,1,1,0,1,0,0,1,1,0,0,1,0,1,1,0};
        rows.reserve(Reps*138*5*2*18);
        const auto start = Now();
        for (unsigned rep = 0; rep < (timing ? Reps : 1u); ++rep)
            for (unsigned slot = 0; slot < 138; ++slot) {
                const unsigned cell = (slot+rep*17)%138, metric = cell%3+1;
                const auto& fs = fixtures[cell/3];
                Fixture matched = fs[0];
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
                            const auto begin = Now(); Work(apis[arm],matched,metric,nullptr,output.data()+16); const auto end = Now();
                            if (timing) rows.push_back(Row{rep,cell,pair,order,pos,arm,end-begin,
                                reinterpret_cast<std::uintptr_t>(matched.source.data()),
                                reinterpret_cast<std::uintptr_t>(output.data()+16),
                                reinterpret_cast<std::uintptr_t>(matched.profile.data()),
                                reinterpret_cast<std::uintptr_t>(matched.repairs[0].data()),
                                reinterpret_cast<std::uintptr_t>(matched.repairs[1].data())});
                            Verify(fs[arm],metric,output);
                            Check(Now()-start < UINT64_C(120000000000),"120-second lifecycle cap");
                        }
                    }
            }
        if (!timing) { std::puts("PASS 46 fixtures x 3 lifecycle metrics x 3 APIs; shared workspace parity"); return 0; }
        Check(rows.size() == Reps*138*5*2*18,"complete cohort");
        Publish();
    } catch (const std::exception& e) {
        std::fprintf(stderr,"FAIL: %s\n",e.what());
        try { if (!rows.empty()) Publish(); } catch (...) {}
        return 1;
    }
}
