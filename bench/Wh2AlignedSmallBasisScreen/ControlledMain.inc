// The preceding source prefix supplies only allocator/loader/oracle helpers.
unsigned active_arm = 0;
unsigned Skip() { return active_arm && arena::offset ? 64-arena::offset : 0; }
void CheckGuards(Shape s)
{
    UNPOISON(arena::basis,sizeof(arena::basis));
    UNPOISON(arena::object,sizeof(arena::object));
    UNPOISON(arena::evaluator,sizeof(arena::evaluator));
    const unsigned begin = 4096+arena::offset+Skip();
    for (unsigned slot = 0; slot < 2; ++slot)
        for (unsigned i = 0; i < arena::CarrierBytes; ++i) {
            const auto wanted = slot == arena::carrier && i >= begin && i < begin+s.k*s.bytes ? source[i-begin] : 173;
            Check(arena::basis[slot][i] == wanted,"aligned basis copy and all guards");
        }
    for (unsigned i = 0; i < 512; ++i) {
        if (i >= arena::extents[0]) Check(arena::object[i] == 173,"handle guard");
        if (i >= arena::extents[2]) Check(arena::evaluator[i] == 173,"evaluator guard");
    }
    POISON(arena::basis,sizeof(arena::basis));
    POISON(arena::object,sizeof(arena::object));
    POISON(arena::evaluator,sizeof(arena::evaluator));
    for (unsigned i = 0; i < 3; ++i) if (arena::live[i]) UNPOISON(arena::live[i],arena::extents[i]);
}
__attribute__((noinline)) void Work(Library& lib, WirehairV2Codec h, Shape s)
{
    for (unsigned batch = 0; batch < Batch; ++batch) for (unsigned j = 0; j < 2*s.k; ++j) {
        std::uint32_t written = 0;
        Check(lib.encode(h,Packet(s.k,j),output+64+(batch*2*s.k+j)*s.bytes,s.bytes,&written) == WirehairV2_Success &&
              written == s.bytes,"repair work");
    }
}
struct Row { unsigned rep,cell,pair,order,pos,side,arm,carrier,offset; std::uint64_t ns; };
std::vector<Row> rows;
bool published = false;
void Publish()
{
    if (published) return;
    published = true;
    std::puts("rep,cell,pair,order,position,side,arm,carrier,offset,ns,handle,evaluator,raw_basis,basis,source,output");
    for (const auto& r : rows) {
        const auto raw = reinterpret_cast<std::uintptr_t>(arena::basis[r.carrier]+4096+r.offset);
        std::printf("%u,%u,%u,%u,%u,%u,%u,%u,%u,%llu,%llu,%llu,%llu,%llu,%llu,%llu\n",
            r.rep,r.cell,r.pair,r.order,r.pos,r.side,r.arm,r.carrier,r.offset,
            static_cast<unsigned long long>(r.ns),
            static_cast<unsigned long long>(reinterpret_cast<std::uintptr_t>(arena::object)),
            static_cast<unsigned long long>(reinterpret_cast<std::uintptr_t>(arena::evaluator)),
            static_cast<unsigned long long>(raw),static_cast<unsigned long long>(raw+(r.arm && r.offset ? 64-r.offset : 0)),
            static_cast<unsigned long long>(reinterpret_cast<std::uintptr_t>(source)),
            static_cast<unsigned long long>(reinterpret_cast<std::uintptr_t>(output+64)));
    }
    Check(!std::ferror(stdout) && std::fflush(stdout) == 0,"output");
}
void Placement(unsigned pair, unsigned side, unsigned& arm, unsigned& carrier, unsigned& offset)
{
    if (pair < 8) { arm = pair/4; carrier = side; offset = pair%4*16; }
    else { arm = side; carrier = (pair-8)/4; offset = (pair-8)%4*16; }
}
void Run(Library& baseline, Library& candidate, bool timing)
{
    Library* libraries[2] = {&baseline,&candidate};
    for (unsigned i = 0; i < MaxMessage; ++i) source[i] = static_cast<std::uint8_t>(i*37+i/11);
    rows.reserve(Reps*6*16*2*18);
    const auto start = Now();
    for (unsigned rep = 0; rep < (timing ? Reps : 1u); ++rep)
        for (unsigned slot = 0; slot < 6; ++slot) {
            const unsigned cell = (slot+rep*5)%6; const auto s = ShapeOf(cell);
            PrepareExpected(baseline,s);
            for (unsigned ps = 0; ps < 16; ++ps) for (unsigned order = 0; order < 2; ++order) {
                const unsigned pair = (ps+rep+cell)%16;
                for (unsigned pos = 0; pos < (timing ? 18u : 2u); ++pos) {
                    const unsigned side = Sides[pos]^order;
                    unsigned arm,carrier,offset; Placement(pair,side,arm,carrier,offset);
                    active_arm = arm;
                    PrepareArena(s,carrier,offset);
                    arena::message_bytes += arm ? 63u : 0u;
                    Check(4096+offset+arena::message_bytes < arena::CarrierBytes,"owned extent");
                    arena::active = true;
                    auto h = Create(*libraries[arm],s);
                    arena::active = false;
                    Check(!arena::bad && arena::calls == 3 && static_cast<void*>(h) == arena::object &&
                          arena::live[1] == arena::basis[carrier]+4096+offset && arena::live[2] == arena::evaluator,"exact allocations");
                    CheckGuards(s);
                    arena::working = true;
                    const auto begin = Now(); Work(*libraries[arm],h,s); const auto end = Now();
                    arena::working = false;
                    const unsigned bytes = 2*s.k*s.bytes;
                    for (unsigned batch = 0; batch < Batch; ++batch)
                        Check(std::memcmp(output+64+batch*bytes,expected,bytes) == 0,"baseline packet parity");
                    Check(std::all_of(output,output+64,[](std::uint8_t x){return x==173;}) &&
                          std::all_of(output+64+Batch*bytes,output+sizeof(output),[](std::uint8_t x){return x==173;}),"output guards");
                    CheckGuards(s);
                    Check(arena::work_allocations == 0,"allocation-free WORK");
                    libraries[arm]->free(h);
                    Check(!arena::bad && !arena::live[0] && !arena::live[1] && !arena::live[2],"exact original-pointer frees");
                    CheckGuards(s);
                    if (timing) rows.push_back(Row{rep,cell,pair,order,pos,side,arm,carrier,offset,end-begin});
                    Check(Now()-start < UINT64_C(120000000000),"120-second controlled cap");
                }
            }
        }
}
}

int main(int argc, char** argv)
{
    try {
        Check(argc == 4,"arguments");
        const bool timing = !std::strcmp(argv[3],"run") || !std::strcmp(argv[3],"run-reverse");
        const bool reverse = !std::strcmp(argv[3],"neutral-reverse") || !std::strcmp(argv[3],"run-reverse");
        Check(timing || !std::strcmp(argv[3],"neutral") || !std::strcmp(argv[3],"neutral-reverse"),"mode");
        cpu_set_t mask; CPU_ZERO(&mask); CPU_SET(50,&mask);
        Check(sched_setaffinity(0,sizeof(mask),&mask) == 0 && sched_getaffinity(0,sizeof(mask),&mask) == 0 &&
              CPU_COUNT(&mask) == 1 && CPU_ISSET(50,&mask) && sched_getcpu() == 50,"CPU50 singleton");
        Library first(argv[reverse ? 2 : 1]), second(argv[reverse ? 1 : 2]);
        Check(first.lib != second.lib,"distinct libraries");
        Run(reverse ? second : first,reverse ? first : second,timing);
        if (!timing) { std::puts("PASS 384 controlled baseline/candidate constructions; copies, packets, guards, frees"); return 0; }
        Check(rows.size() == Reps*6*16*2*18,"complete cohort");
        Publish();
    } catch (const std::exception& e) {
        arena::active = false; arena::working = false;
        std::fprintf(stderr,"FAIL: %s\n",e.what());
        try { if (!rows.empty()) Publish(); } catch (...) {}
        return 1;
    }
}
