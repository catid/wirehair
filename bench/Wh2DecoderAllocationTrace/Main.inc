namespace {
bool Disjoint(std::uintptr_t a, std::size_t na, std::uintptr_t b, std::size_t nb)
{
    return a <= b ? na <= b-a : nb <= a-b;
}
void Validate(const Fixture& f, unsigned family)
{
    using namespace allocation_trace;
    const Shape& s = f.shape;
    Check(!overflow && count == Batch*6,"exact small-decoder event count");
    const unsigned feeds = f.steps[family], stride = feeds+3;
    Check(call_count == Batch*stride,"exact public call count");
    const std::size_t object = s.k == 3 ? 104 : s.k == 5 ? 136 : 200;
    for (unsigned cycle = 0; cycle < Batch; ++cycle) {
        const Event* e = events+cycle*6;
        const std::size_t sizes[3] = {296,object,std::size_t(s.k+1)*s.b};
        for (unsigned j = 0; j < 3; ++j) {
            Check(e[j].allocate && e[j].phase == 0 && e[j].array == (j == 2) && e[j].bytes == sizes[j],
                "handle/decoder/slab allocation order and size");
            Check(!e[5-j].allocate && e[5-j].phase == 3 && e[5-j].array == e[j].array &&
                  e[5-j].pointer == e[j].pointer,"reverse original-pointer frees");
            for (unsigned i = 0; i < j; ++i)
                Check(Disjoint(e[i].pointer,e[i].bytes,e[j].pointer,e[j].bytes),"disjoint live allocations");
        }
        const Invocation* c = calls+cycle*stride;
        Check(c[0].phase == 0 && !c[0].handle && c[0].bytes == 32 && c[0].first && c[0].second,"constructor arguments");
        Check(c[0].first == reinterpret_cast<std::uintptr_t>(f.profile.data()) &&
              Disjoint(c[0].first,32,c[0].second,sizeof(void*)),"staged descriptor and public pointer-out");
        for (unsigned j = 1; j <= feeds; ++j)
            Check(c[j].phase == 1 && c[j].handle == e[0].pointer && c[j].first && !c[j].second &&
                  c[j].bytes == s.b && c[j].packet == (family ? UINT32_MAX-2*(j-1) : s.k+j-1) &&
                  c[j].first == reinterpret_cast<std::uintptr_t>(f.repairs[family].data()+std::size_t(j-1)*s.b),"feed arguments");
        Check(c[feeds+1].phase == 2 && c[feeds+1].handle == e[0].pointer && c[feeds+1].first &&
              c[feeds+1].second && c[feeds+1].bytes == s.Message(),"recover arguments");
        Check(Disjoint(c[feeds+1].first,s.Message(),c[feeds+1].second,sizeof(uint64_t)) &&
              c[feeds+1].first == calls[feeds+1].first + std::size_t(cycle)*3*s.k*s.b,"recover public outputs/stride");
        Check(c[feeds+2].phase == 3 && c[feeds+2].handle == e[0].pointer &&
              !c[feeds+2].first && !c[feeds+2].second && !c[feeds+2].bytes,"free argument");
        for (unsigned j = 0; j < stride-1; ++j) for (unsigned a = 0; a < 3; ++a) {
            Check(Disjoint(c[j].first,c[j].bytes,e[a].pointer,e[a].bytes),"public/private buffer disjointness");
            if (c[j].second) Check(Disjoint(c[j].second,j == 0 ? sizeof(void*) : sizeof(uint64_t),
                e[a].pointer,e[a].bytes),"public pointer/count output disjointness");
        }
    }
}
void Print(unsigned fixture, const Shape& s, unsigned family, unsigned history,
           unsigned order, unsigned position, unsigned arm, const Library& library)
{
    Dl_info owner = {}, observer = {};
    Check(dladdr(dlsym(library.lib,"wirehair_init_"), &owner),"library owner");
    Check(dladdr(reinterpret_cast<void*>(&allocation_trace::Start), &observer),"observer owner");
    for (unsigned i = 0; i < allocation_trace::count; ++i) {
        const auto& e = allocation_trace::events[i];
        Dl_info caller = {};
        Check(dladdr(reinterpret_cast<void*>(e.caller),&caller),"resolved allocation/free caller");
        // The final handle delete is a tail call, so its return address is
        // in Work, not in the DSO whose destructor selected that delete.
        const bool from_observer = caller.dli_fbase == observer.dli_fbase;
        Check(caller.dli_fbase == owner.dli_fbase || (from_observer && i%6 == 5),"caller owner");
        std::printf("%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,0,%u,%u,%u,%zu,%llu,%u,%llu,%llu,%llu,%u,0,0,0\n",
            fixture,s.k,s.b,s.tail,s.policy,family,history,order,position,arm,i,
            unsigned(e.allocate),unsigned(e.array),e.bytes,
            static_cast<unsigned long long>(e.pointer),
            unsigned(from_observer),
            static_cast<unsigned long long>(e.caller-reinterpret_cast<std::uintptr_t>(caller.dli_fbase)),
            static_cast<unsigned long long>(e.stack),
            static_cast<unsigned long long>(reinterpret_cast<std::uintptr_t>(owner.dli_fbase)),e.phase);
    }
    for (unsigned i = 0; i < allocation_trace::call_count; ++i) {
        const auto& c = allocation_trace::calls[i];
        std::printf("%u,%u,%u,%u,%u,%u,%u,%u,%u,%u,1,%u,0,0,%zu,%llu,0,0,0,%llu,%u,%llu,%llu,%u\n",
            fixture,s.k,s.b,s.tail,s.policy,family,history,order,position,arm,i,c.bytes,
            static_cast<unsigned long long>(c.handle),
            static_cast<unsigned long long>(reinterpret_cast<std::uintptr_t>(owner.dli_fbase)),c.phase,
            static_cast<unsigned long long>(c.first),static_cast<unsigned long long>(c.second),c.packet);
    }
}
}

int main(int argc, char** argv)
{
    try {
        Check(argc == 4 && (std::strcmp(argv[3],"normal") == 0 ||
              std::strcmp(argv[3],"reverse") == 0),"arguments");
        const bool reverse = std::strcmp(argv[3],"reverse") == 0;
        Library first(argv[reverse ? 2 : 1]), second(argv[reverse ? 1 : 2]);
        Check(first.lib != second.lib,"distinct libraries");
        Api apis[2] = {{reverse ? &second : &first,false},{reverse ? &first : &second,false}};
        std::vector<std::array<Fixture,2>> fixtures;
        for (unsigned k : {3u,5u,8u}) for (unsigned variant = 0; variant < 7; ++variant)
            for (unsigned policy = 0; policy < 2; ++policy) {
                const unsigned widths[7] = {2,2,64,256,257,1280,1280};
                const unsigned width = widths[variant], tail = variant == 1 || variant == 6 ? 1u : width;
                Shape s = {k,width,tail,policy};
                std::array<Fixture,2> f = {{Prepare(apis[0],s),Prepare(apis[1],s)}};
                Check(f[0].profile == f[1].profile && f[0].source == f[1].source &&
                      f[0].packets == f[1].packets && f[0].repairs == f[1].repairs &&
                      f[0].steps[0] == f[1].steps[0] && f[0].steps[1] == f[1].steps[1],"fixture parity");
                fixtures.push_back(std::move(f));
            }
        Check(fixtures.size() == 42,"roster");
        std::puts("fixture,k,width,tail,policy,family,history,order,position,arm,kind,event,allocate,array,bytes,pointer,caller_owner,caller_offset,stack,dso_base,phase,arg1,arg2,packet");
        unsigned observations = 0;
        for (unsigned fixture = 0; fixture < fixtures.size(); ++fixture) {
            Fixture matched = fixtures[fixture][0];
            const Shape& s = matched.shape;
            std::vector<Byte> output(Batch*std::size_t(3*s.k)*s.b+32,173);
            for (unsigned family = 0; family < 2; ++family)
                for (unsigned history = 0; history < 3; ++history)
                    for (unsigned order = 0; order < 2; ++order)
                        for (unsigned position = 0; position < 4; ++position) {
                            const unsigned sides[4] = {0,1,1,0};
                            const unsigned arm = sides[position]^order;
                            // History is a non-timing allocator precondition, not a
                            // reproduction of the original timing process's heap.
                            if (history) {
                                std::fill(output.begin()+16,output.end()-16,0);
                                Work(apis[history-1],matched,1,nullptr,output.data()+16);
                                Verify(matched,1,output);
                            }
                            std::fill(output.begin()+16,output.end()-16,0);
                            allocation_trace::Start();
                            Work(apis[arm],matched,family+2,nullptr,output.data()+16);
                            allocation_trace::Stop();
                            Validate(matched,family);
                            Verify(matched,family+2,output);
                            Print(fixture,s,family,history,order,position,arm,*apis[arm].lib);
                            ++observations;
                        }
        }
        Check(observations == 2016 && !std::ferror(stdout) && std::fflush(stdout) == 0,"complete trace");
        std::fprintf(stderr,"PASS 2016 observations, 64512 decoders, 387072 C++ allocation/free events; no timing\n");
    } catch (const std::exception& e) {
        allocation_trace::Stop();
        std::fprintf(stderr,"FAIL: %s\n",e.what());
        return 1;
    }
}
