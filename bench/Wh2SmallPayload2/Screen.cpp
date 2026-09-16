// Bounded screening only, not a production/all-K qualification. DSO calls
// prevent loop-invariant hoisting; every output is checked outside WORK.
#include <wirehair/wirehair.h>
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <dlfcn.h>
#include <sched.h>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using Byte = uint8_t;
using Profile = std::array<Byte, 32>;
const unsigned Batch = 32, Reps = 12;
void Check(bool ok, const char* why) { if (!ok) throw std::runtime_error(why); }
template<class F> F Symbol(void* library, const char* name, void* owner)
{
    void* p = dlsym(library, name);
    Dl_info info = {};
    Check(p && dladdr(p, &info) && info.dli_fbase == owner, "DSO symbol owner");
    F f; static_assert(sizeof(f) == sizeof(p), "Linux function pointer size");
    std::memcpy(&f, &p, sizeof(f));
    return f;
}
struct Library {
    void* lib;
    decltype(&wirehair_v2_encoder_create_profile_id_with_options) create;
    decltype(&wirehair_v2_encode) encode;
    decltype(&wirehair_v2_decoder_create) decoder;
    decltype(&wirehair_v2_decode) feed;
    decltype(&wirehair_v2_recover) recover;
    decltype(&wirehair_v2_free) free;
    decltype(&wirehair_encoder_create_ex) wcreate, wowned;
    decltype(&wirehair_encode) wencode;
    decltype(&wirehair_decoder_create_ex) wdecoder;
    decltype(&wirehair_decode) wfeed;
    decltype(&wirehair_recover) wrecover;
    decltype(&wirehair_free) wfree;
    explicit Library(const char* path): lib(dlopen(path, RTLD_NOW | RTLD_LOCAL))
    {
        Check(lib, "dlopen");
        Dl_info owner = {};
        Check(dladdr(dlsym(lib, "wirehair_init_"), &owner), "DSO base");
        Check(Symbol<decltype(&wirehair_init_)>(lib, "wirehair_init_", owner.dli_fbase)(WIREHAIR_VERSION) == Wirehair_Success, "init");
#define LOAD(FIELD, NAME) FIELD = Symbol<decltype(FIELD)>(lib, NAME, owner.dli_fbase)
        LOAD(create, "wirehair_v2_encoder_create_profile_id_with_options");
        LOAD(encode, "wirehair_v2_encode"); LOAD(decoder, "wirehair_v2_decoder_create");
        LOAD(feed, "wirehair_v2_decode"); LOAD(recover, "wirehair_v2_recover"); LOAD(free, "wirehair_v2_free");
        LOAD(wcreate, "wirehair_encoder_create_ex"); LOAD(wowned, "wirehair_encoder_create_owned_ex");
        LOAD(wencode, "wirehair_encode"); LOAD(wdecoder, "wirehair_decoder_create_ex");
        LOAD(wfeed, "wirehair_decode"); LOAD(wrecover, "wirehair_recover"); LOAD(wfree, "wirehair_free");
#undef LOAD
    }
    ~Library() { dlclose(lib); }
    Library(const Library&) = delete;
    Library& operator=(const Library&) = delete;
};
uint64_t ProfileId(unsigned k)
{
    return k == 3 ? WIREHAIR_V2_PROFILE_SMALL_K3_2026_09 :
        k == 5 ? WIREHAIR_V2_PROFILE_SMALL_K5_2026_09 : WIREHAIR_V2_PROFILE_SMALL_K8_2026_09;
}
struct Shape {
    unsigned k, b, tail, policy;
    size_t Message() const { return size_t(k-1)*b + tail; }
};
struct Api {
    Library* lib; bool wh1;
    int Create(const Shape& s, const void* source, Profile& p, void*& h) const
    {
        if (wh1) {
            WirehairCodec c = nullptr;
            const int r = (s.policy ? lib->wcreate : lib->wowned)(nullptr, source, s.Message(), s.b, &c);
            h = c; return r;
        }
        WirehairV2EncoderOptions opts = WIREHAIR_V2_ENCODER_OPTIONS_INIT;
        opts.source_policy = s.policy ? WirehairV2EncoderSource_BorrowedImmutable : WirehairV2EncoderSource_Independent;
        WirehairV2Codec c = nullptr; uint32_t n = 0;
        const int r = lib->create(ProfileId(s.k), source, s.Message(), s.b, &opts, p.data(), 32, &n, &c);
        h = c; return r == 0 && n != 32 ? -1 : r;
    }
    int Encode(void* h, uint32_t id, void* out, unsigned b, uint32_t* n) const
    {
        return wh1 ? int(lib->wencode(static_cast<WirehairCodec>(h), id, out, b, n)) :
            int(lib->encode(static_cast<WirehairV2Codec>(h), id, out, b, n));
    }
    int Decoder(const Shape& s, const Profile& p, void*& h) const
    {
        if (wh1) {
            WirehairCodec c = nullptr; int r = lib->wdecoder(nullptr, s.Message(), s.b, &c); h = c; return r;
        }
        WirehairV2Codec c = nullptr; int r = lib->decoder(p.data(), 32, &c); h = c; return r;
    }
    int Feed(void* h, uint32_t id, const void* in, unsigned n) const
    { return wh1 ? int(lib->wfeed(static_cast<WirehairCodec>(h), id, in, n)) : int(lib->feed(static_cast<WirehairV2Codec>(h), id, in, n)); }
    int Recover(void* h, void* out, size_t bytes) const
    {
        if (wh1) return lib->wrecover(static_cast<WirehairCodec>(h), out, bytes);
        uint64_t n = 0;
        int r = lib->recover(static_cast<WirehairV2Codec>(h), out, bytes, &n);
        return r == 0 && n != bytes ? -1 : r;
    }
    void Free(void* h) const
    { if (wh1) lib->wfree(static_cast<WirehairCodec>(h)); else lib->free(static_cast<WirehairV2Codec>(h)); }
};
struct Handle {
    Api api; void* value = nullptr;
    explicit Handle(Api a): api(a) {}
    ~Handle() { if (value) api.Free(value); }
};
uint32_t Packet(const Shape& s, unsigned slot)
{ return slot < 2*s.k ? slot : UINT32_MAX - 2*(slot-2*s.k); }
struct Fixture {
    Shape shape; std::vector<Byte> source, packets;
    std::array<std::vector<Byte>, 2> repairs;
    Profile profile = {};
    unsigned steps[2] = {};
};
Fixture Prepare(Api api, Shape s)
{
    Fixture f; f.shape = s; f.source.resize(s.Message()); f.packets.resize(size_t(3*s.k)*s.b, 0);
    for (size_t i = 0; i < f.source.size(); ++i) f.source[i] = static_cast<Byte>(i*37+i/11+s.k);
    Handle encoder(api);
    Check(api.Create(s, f.source.data(), f.profile, encoder.value) == 0 && encoder.value, "fixture encoder");
    for (unsigned i = 0; i < 3*s.k; ++i) {
        uint32_t n = 0;
        Check(api.Encode(encoder.value, Packet(s,i), f.packets.data()+size_t(i)*s.b, s.b, &n) == 0 &&
            n == (i == s.k-1 ? s.tail : s.b), "fixture packet");
    }
    for (unsigned family = 0; family < 2; ++family) {
        f.repairs[family].resize(size_t(32)*s.b);
        Handle decoder(api);
        Check(api.Decoder(s, f.profile, decoder.value) == 0 && decoder.value, "fixture decoder");
        for (unsigned i = 0; i < 32; ++i) {
            uint32_t id = family ? UINT32_MAX-2*i : s.k+i, n = 0;
            Byte* out = f.repairs[family].data()+size_t(i)*s.b;
            Check(api.Encode(encoder.value, id, out, s.b, &n) == 0 && n == s.b, "fixture repair");
            int r = api.Feed(decoder.value, id, out, n);
            Check(r == 0 || r == 1, "fixture feed");
            if (r == 0) { f.steps[family] = i+1; break; }
        }
        Check(f.steps[family] >= s.k, "fixture success bound");
        std::vector<Byte> recovered(s.Message());
        Check(api.Recover(decoder.value, recovered.data(), recovered.size()) == 0 && recovered == f.source, "fixture recovery");
    }
    return f;
}

// metric 0: prebuilt repair; 1: full encoder; 2/3: repair-only decoder, own first success.
__attribute__((noinline)) void Work(Api api, const Fixture& f, unsigned metric, void* prebuilt, Byte* output)
{
    const Shape& s = f.shape;
    const size_t stride = size_t(3*s.k)*s.b;
    for (unsigned cycle = 0; cycle < Batch; ++cycle) {
        Byte* out = output + cycle*stride;
        Handle h(api); Profile p = {};
        if (metric == 0) {
            for (unsigned j = s.k; j < 3*s.k; ++j) {
                uint32_t n = 0;
                Check(api.Encode(prebuilt, Packet(s,j), out+size_t(j)*s.b, s.b, &n) == 0 && n == s.b, "repair work");
            }
        } else if (metric == 1) {
            Check(api.Create(s, f.source.data(), p, h.value) == 0 && h.value && p == f.profile, "encoder work");
            for (unsigned j = 0; j < 3*s.k; ++j) {
                uint32_t n = 0;
                Check(api.Encode(h.value, Packet(s,j), out+size_t(j)*s.b, s.b, &n) == 0 &&
                    n == (j == s.k-1 ? s.tail : s.b), "encode work");
            }
        } else {
            const unsigned family = metric-2;
            Check(api.Decoder(s, f.profile, h.value) == 0 && h.value, "decoder work");
            for (unsigned j = 0; j < f.steps[family]; ++j) {
                const uint32_t id = family ? UINT32_MAX-2*j : s.k+j;
                Check(api.Feed(h.value, id, f.repairs[family].data()+size_t(j)*s.b, s.b) ==
                    (j+1 == f.steps[family] ? 0 : 1), "feed work");
            }
            Check(api.Recover(h.value, out, s.Message()) == 0, "recover work");
        }
    }
}
void Verify(const Fixture& f, unsigned metric, const std::vector<Byte>& output)
{
    const Shape& s = f.shape; const size_t stride = size_t(3*s.k)*s.b;
    for (unsigned cycle = 0; cycle < Batch; ++cycle) {
        const Byte* out = output.data()+16+cycle*stride;
        for (size_t i = 0; i < stride; ++i) {
            Byte expected = 0;
            if (metric >= 2) expected = i < s.Message() ? f.source[i] : 0;
            else if (metric != 0 || i >= size_t(s.k)*s.b) expected = f.packets[i];
            Check(out[i] == expected, "work bytes including padding");
        }
    }
    Check(std::all_of(output.begin(), output.begin()+16, [](Byte b){return b==173;}) &&
        std::all_of(output.end()-16, output.end(), [](Byte b){return b==173;}), "work guard");
}
uint64_t Now() { return std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now().time_since_epoch()).count(); }
struct Row { unsigned rep, cell, pair, order, position, arm; uint64_t ns; };
std::vector<Row> rows;
bool published = false;
void Publish()
{
    if (published) return;
    published = true;
    std::puts("rep,cell,pair,order,position,arm,ns");
    for (const auto& r : rows) std::printf("%u,%u,%u,%u,%u,%u,%llu\n",r.rep,r.cell,r.pair,r.order,r.position,r.arm,static_cast<unsigned long long>(r.ns));
    Check(!std::ferror(stdout) && std::fflush(stdout) == 0, "publish");
}
} // namespace

int main(int argc, char** argv)
{
    try {
        Check(argc == 4, "arguments");
        const std::string mode = argv[3];
        Check(mode == "neutral" || mode == "neutral-reverse" || mode == "run" || mode == "run-reverse", "mode");
        const bool timing = mode == "run" || mode == "run-reverse";
        const bool reverse = mode == "run-reverse" || mode == "neutral-reverse";
        if (timing) {
            cpu_set_t cpus; CPU_ZERO(&cpus); CPU_SET(50, &cpus);
            Check(sched_setaffinity(0,sizeof(cpus),&cpus) == 0 && sched_getcpu() == 50, "pin CPU50");
        }
        Library first(argv[reverse ? 2 : 1]), second(argv[reverse ? 1 : 2]);
        Check(first.lib != second.lib, "distinct baseline/candidate loader handles");
        Library* baseline = reverse ? &second : &first;
        Library* candidate = reverse ? &first : &second;
        Api apis[3] = {{baseline,false},{candidate,false},{baseline,true}};
        std::vector<std::array<Fixture,3>> cells;
        for (unsigned k : {3u,5u,8u}) for (unsigned shape = 0; shape < 4; ++shape)
            for (unsigned policy = 0; policy < 2; ++policy) {
                const unsigned b = shape < 2 ? 2 : shape == 2 ? 64 : 1280;
                Shape s = {k,b,shape == 1 ? 1 : b,policy};
                std::array<Fixture,3> f;
                for (unsigned a = 0; a < 3; ++a) f[a] = Prepare(apis[a],s);
                Check(f[0].profile == f[1].profile && f[0].packets == f[1].packets &&
                    f[0].repairs == f[1].repairs && f[0].steps[0] == f[1].steps[0] && f[0].steps[1] == f[1].steps[1], "A/B equation parity");
                cells.push_back(f);
            }
        const unsigned pairs[5][2] = {{0,0},{1,1},{2,2},{0,1},{2,1}};
        const unsigned sides[18] = {0,1,0,1,1,0,1,0,0,1,1,0,0,1,0,1,1,0};
        rows.reserve(Reps*cells.size()*4*5*2*18);
        const uint64_t start = Now();
        for (unsigned rep = 0; rep < (timing ? Reps : 1); ++rep)
            for (unsigned slot = 0; slot < cells.size()*4; ++slot) {
                unsigned cell = (slot + rep*17) % (cells.size()*4), metric = cell%4;
                const auto& fs = cells[cell/4];
                Handle h0(apis[0]), h1(apis[1]), h2(apis[2]);
                Handle* hs[] = {&h0,&h1,&h2};
                if (metric == 0) for (unsigned a = 0; a < 3; ++a) {
                    Profile p = {}; Check(apis[a].Create(fs[a].shape,fs[a].source.data(),p,hs[a]->value) == 0 &&
                        hs[a]->value && p == fs[a].profile, "prebuilt create");
                }
                const size_t size = Batch*size_t(3*fs[0].shape.k)*fs[0].shape.b;
                std::vector<Byte> output[2] = {std::vector<Byte>(size+32,173),std::vector<Byte>(size+32,173)};
                for (unsigned ps = 0; ps < (timing ? 5u : 1u); ++ps) for (unsigned order = 0; order < (timing ? 2u : 1u); ++order) {
                    unsigned pair = (ps+rep+cell)%5;
                    for (unsigned pos = 0; pos < (timing ? 18u : 3u); ++pos) {
                        unsigned side = sides[pos]^order, arm = timing ? pairs[pair][side] : pos;
                        std::fill(output[side].begin()+16,output[side].end()-16,0);
                        const uint64_t begin = Now();
                        Work(apis[arm], fs[arm], metric, hs[arm]->value, output[side].data()+16);
                        const uint64_t end = Now();
                        if (timing) rows.push_back(Row{rep,cell,pair,order,pos,arm,end-begin});
                        Verify(fs[arm],metric,output[side]);
                    }
                }
                Check(!timing || Now()-start < UINT64_C(180000000000), "180-second screen cap");
            }
        if (!timing) { std::puts("PASS 24 fixtures x 4 metrics x 3 APIs, complete output parity"); return 0; }
        Publish();
        return 0;
    } catch (const std::exception& e) {
        std::fprintf(stderr,"FAIL: %s\n",e.what());
        try { if (!rows.empty()) Publish(); } catch (...) {}
        return 1;
    }
}
