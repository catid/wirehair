// Count only scalar coefficient multiplies in an isolated copy of the core.
// Public-library templates and payload kernels are never instrumented.
#include "gf256.h"
#include <cstdint>
namespace {
std::uint64_t counted_multiplies = 0;
__attribute__((noinline)) std::uint8_t CountedMul(std::uint8_t a, std::uint8_t b)
{
    ++counted_multiplies;
    return gf256_mul(a,b);
}
}
#define gf256_mul CountedMul
#define wirehair_small_core diagonal_probe_core
#include "WirehairSmallCore.h"
#undef wirehair_small_core
#undef gf256_mul
#include <algorithm>
#include <array>
#include <cstdio>
#include <cstdlib>
#include <vector>

namespace {
namespace S = diagonal_probe_core;
using Byte = std::uint8_t;
std::uint64_t cases = 0, feeds = 0, recoveries = 0, eliminations = 0, insertions = 0;
std::uint64_t full_echelon_checks = 0, nonidentity_echelon_checks = 0;
std::uint64_t baseline_multiplies = 0, candidate_multiplies = 0, actual_multiplies = 0;
void Check(bool ok, const char* why)
{
    if (!ok) { std::fprintf(stderr,"FAIL: %s\n",why); std::abort(); }
}
Byte Mul(Byte a, Byte b)
{
    unsigned product = 0;
    for (unsigned bit = 0; bit < 8; ++bit) if (b & (1u << bit)) product ^= unsigned(a) << bit;
    for (int bit = 14; bit >= 8; --bit)
        if (product & (1u << bit)) product ^= 0x14du << (bit-8);
    return static_cast<Byte>(product);
}
Byte Inv(Byte value)
{
    Check(value != 0,"nonzero oracle pivot");
    Byte result = 1;
    for (unsigned exponent = 254; exponent; exponent >>= 1) {
        if (exponent & 1u) result = Mul(result,value);
        value = Mul(value,value);
    }
    return result;
}
template<unsigned K> using Row = std::array<Byte,K>;

template<unsigned K> struct Table {
    using G = S::detail::Geometry<K>;
    std::vector<Byte> bytes;
    Table(): bytes(G::LookupBytes,0)
    {
        for (unsigned id = 0; id < 1024; ++id) {
            Row<K> row = {};
            for (unsigned c = 0; c < K; ++c)
                row[c] = static_cast<Byte>((id*37+c*19+id*c*7+11)&255);
            if (id < K) { row.fill(0); row[id] = 1; }
            Set(id,row);
        }
        // Both low phases agree, and every middle/high operator is identity.
        // Thus the independent map is simply low[id&1023], even at high IDs.
        for (unsigned phase = 0; phase < 2; ++phase) for (unsigned entry = 0; entry < 128; ++entry)
            for (unsigned d = 0; d < K; ++d) {
                bytes[G::Middle10+phase*G::MiddlePhase+entry*G::MatrixBytes+d*K+d] = 1;
                bytes[G::Middle17+phase*G::MiddlePhase+entry*G::MatrixBytes+d*K+d] = 1;
            }
        for (unsigned entry = 0; entry < 256; ++entry) for (unsigned d = 0; d < K; ++d)
            bytes[G::High24+entry*G::MatrixBytes+d*K+d] = 1;
        Row<K> zero = {}; Set(96,zero);
    }
    void Set(unsigned id, const Row<K>& row)
    {
        Check(id < 1024,"test table row bound");
        for (unsigned phase = 0; phase < 2; ++phase)
            std::copy(row.begin(),row.end(),bytes.begin()+phase*G::LowPhase+id*K);
    }
    Row<K> Coefficients(std::uint32_t id) const
    {
        Row<K> row = {};
        std::copy_n(bytes.data()+(id&1023u)*K,K,row.data());
        return row;
    }
    S::Lookup View() const { return {bytes.data(),bytes.size()}; }
};

struct Cost { unsigned baseline = 0, saved = 0, eliminated = 0, inserted = 0; };
template<unsigned K> struct Oracle {
    std::array<Row<K>,K> rows = {};
    std::array<std::vector<Byte>,K> rhs;
    unsigned rank = 0;
    explicit Oracle(unsigned width)
    { for (auto& value : rhs) value.resize(width,0); }
    S::Status Feed(Row<K> row, const std::vector<Byte>& packet, Cost& cost)
    {
        std::vector<Byte> value(rhs[0].size(),0);
        std::copy(packet.begin(),packet.end(),value.begin());
        for (unsigned p = 0; p < K; ++p) if (rows[p][p] && row[p]) {
            Check(rows[p][p] == 1,"oracle retained unit diagonal");
            const Byte factor = row[p];
            // Oracle deliberately multiplies all columns, including zeros.
            for (unsigned c = 0; c < K; ++c) row[c] ^= Mul(rows[p][c],factor);
            for (std::size_t j = 0; j < value.size(); ++j) value[j] ^= Mul(rhs[p][j],factor);
            cost.baseline += K-p; ++cost.saved; ++cost.eliminated;
        }
        unsigned pivot = 0;
        while (pivot < K && !row[pivot]) ++pivot;
        if (pivot == K) {
            if (std::any_of(value.begin(),value.end(),[](Byte b){ return b != 0; })) return S::Status::Conflict;
            return rank == K ? S::Status::Success : S::Status::NeedMore;
        }
        Check(!rows[pivot][pivot],"oracle new pivot unoccupied");
        const Byte inverse = Inv(row[pivot]);
        for (Byte& b : row) b = Mul(b,inverse);
        for (Byte& b : value) b = Mul(b,inverse);
        cost.baseline += K-pivot; ++cost.saved; ++cost.inserted;
        rows[pivot] = row; rhs[pivot] = value; ++rank;
        return rank == K ? S::Status::Success : S::Status::NeedMore;
    }
    void Solve()
    {
        Check(rank == K,"oracle solve rank");
        // Full Gauss-Jordan, not the production triangular back-substitution.
        for (unsigned p = 0; p < K; ++p) for (unsigned r = 0; r < K; ++r) if (r != p) {
            const Byte factor = rows[r][p];
            for (unsigned c = 0; c < K; ++c) rows[r][c] ^= Mul(rows[p][c],factor);
            for (std::size_t j = 0; j < rhs[r].size(); ++j) rhs[r][j] ^= Mul(rhs[p][j],factor);
        }
    }
};

template<unsigned K> void Exercise(unsigned width, unsigned tail, unsigned pattern, Byte diagonal)
{
    using D = S::Decoder<K>;
    Table<K> table;
    for (unsigned p = 0; p < K; ++p) {
        Row<K> row = {}; row[p] = diagonal;
        for (unsigned c = p+1; c < K; ++c) row[c] = static_cast<Byte>((17*p+29*c+3)&255);
        table.Set(32+p,row);
    }
    const std::size_t message = std::size_t(K-1)*width+tail;
    std::vector<Byte> source(message), output(message+32,0xa5);
    for (std::size_t j = 0; j < message; ++j) source[j] = static_cast<Byte>((j*47+j/13+K)&255);
    std::unique_ptr<D> decoder;
    Check(D::Create(table.View(),message,width,decoder) == S::Status::Success && decoder,"decoder create");
    Oracle<K> oracle(width);
    auto submit = [&](std::uint32_t id, bool corrupt) {
        const Row<K> row = table.Coefficients(id);
        const unsigned required = id == K-1 ? tail : width;
        std::vector<Byte> packet(required,0);
        for (unsigned j = 0; j < required; ++j) for (unsigned c = 0; c < K; ++c) {
            const std::size_t index = std::size_t(c)*width+j;
            if (index < message) packet[j] ^= Mul(row[c],source[index]);
        }
        if (corrupt) packet[0] ^= 1;
        Cost cost;
        const S::Status expected = oracle.Feed(row,packet,cost);
        if (corrupt) Check(expected == S::Status::Conflict,"only corrupt dependent test rows");
        if (id >= 1024) cost.baseline += K*K*unsigned(
            ((id >> 24) != 0) + (((id >> 17)&127u) != 0) + (((id >> 10)&127u) != 0));
        std::vector<Byte> guarded(required+2,0x73);
        std::copy(packet.begin(),packet.end(),guarded.begin()+1);
        counted_multiplies = 0;
        const auto answer = decoder->Feed(id,guarded.data()+1,required);
        const std::uint64_t measured = counted_multiplies;
        Check(answer.status == expected && answer.bytes_required == required && !answer.bytes_written &&
              decoder->Rank() == oracle.rank,"feed status/extent/rank parity");
#ifdef UNIT_DIAGONAL_CANDIDATE
        Check(measured == cost.baseline-cost.saved,"candidate exact scalar work reduction");
#else
        Check(measured == cost.baseline,"baseline exact scalar work");
#endif
        Check(guarded.front() == 0x73 && guarded.back() == 0x73 &&
              std::equal(packet.begin(),packet.end(),guarded.begin()+1),"feed input unchanged");
        baseline_multiplies += cost.baseline; candidate_multiplies += cost.baseline-cost.saved;
        actual_multiplies += measured; eliminations += cost.eliminated; insertions += cost.inserted; ++feeds;
    };
    auto recover = [&]() {
        std::fill(output.begin(),output.end(),0xa5);
        counted_multiplies = 0;
        Check(decoder->Recover(output.data()+16,message-1).status == S::Status::BufferTooSmall &&
              counted_multiplies == 0,"undersized recovery rejects before coefficient work");
        Check(std::all_of(output.begin(),output.end(),[](Byte b){ return b == 0xa5; }),"rejected recovery preserves output");
        const auto answer = decoder->Recover(output.data()+16,message);
        Check(answer.status == (oracle.rank == K ? S::Status::Success : S::Status::NeedMore) &&
              answer.bytes_required == message && answer.bytes_written == (oracle.rank == K ? message : 0) &&
              counted_multiplies == 0,"recovery status/extent and unchanged scalar work");
        if (oracle.rank == K) {
            oracle.Solve();
            for (unsigned c = 0; c < K; ++c) for (unsigned j = 0; j < (c == K-1 ? tail : width); ++j)
                Check(oracle.rhs[c][j] == source[std::size_t(c)*width+j],"independent solved bytes");
            Check(std::equal(source.begin(),source.end(),output.begin()+16),"exact recovered message");
        } else Check(std::all_of(output.begin(),output.end(),[](Byte b){ return b == 0xa5; }),"incomplete output unchanged");
        Check(std::all_of(output.begin(),output.begin()+16,[](Byte b){ return b == 0xa5; }) &&
              std::all_of(output.end()-16,output.end(),[](Byte b){ return b == 0xa5; }),"recovery guards");
        ++recoveries;
    };
    bool checked_first_full_rank = false;
    auto triplet = [&](std::uint32_t id) {
        submit(id,false); submit(id,false); submit(id,true);
        if (oracle.rank == K && !checked_first_full_rank) {
            bool nonidentity = false;
            for (unsigned r = 0; r < K; ++r) for (unsigned c = r+1; c < K; ++c)
                nonidentity = nonidentity || oracle.rows[r][c] != 0;
            if (pattern == 0) Check(nonidentity,"forward full echelon is not yet identity");
            // A novel dense equation must be dependent before the first Solve,
            // even when the retained full-rank basis still has off-diagonals.
            submit(UINT32_C(0xfedc0040),false);
            submit(UINT32_C(0xfedc0040),true);
            ++full_echelon_checks;
            nonidentity_echelon_checks += nonidentity;
            checked_first_full_rank = true;
        }
        recover(); recover();
    };
    counted_multiplies = 0;
    Check(decoder->Feed(0,nullptr,width).status == S::Status::InvalidInput &&
          decoder->Feed(0,source.data(),0).status == S::Status::InvalidInput &&
          decoder->Rank() == 0 && counted_multiplies == 0,"invalid feed preserves empty basis");
    submit(96,false); submit(96,true); recover();
    if (pattern < 2) for (unsigned slot = 0; slot < K; ++slot)
        triplet(32+(pattern == 0 ? slot : K-1-slot));
    else {
        const std::uint32_t high[] = {0u,1024u,1u<<17,1u<<24,0xfffffc00u};
        for (unsigned slot = 0; slot < 2*K; ++slot) triplet(high[slot%5] | (64+slot));
    }
    for (unsigned slot = K; slot-- > 0;) triplet(slot);
    Check(decoder->Rank() == K && checked_first_full_rank,"systematic completion and full-echelon coverage");
    // Feed after Solve, including high-ID identity operators and zero rows.
    triplet(UINT32_C(0xfffffc20)); triplet(96); recover();
    ++cases;
}

template<unsigned K> void Dimension()
{
    for (unsigned width : {1u,2u,3u,17u,33u,64u,1280u}) for (unsigned tail : {1u,width})
        for (unsigned pattern = 0; pattern < 3; ++pattern) Exercise<K>(width,tail,pattern,Byte(173));
    for (unsigned value = 1; value < 256; ++value) for (unsigned pattern = 0; pattern < 2; ++pattern)
        Exercise<K>(2,2,pattern,static_cast<Byte>(value));
}
}

int main()
{
    Check(gf256_init() == 0,"GF initialization");
    for (unsigned a = 0; a < 256; ++a) for (unsigned b = 0; b < 256; ++b)
        Check(Mul(static_cast<Byte>(a),static_cast<Byte>(b)) == gf256_mul(static_cast<Byte>(a),static_cast<Byte>(b)),
              "exhaustive polynomial product oracle");
    for (unsigned a = 1; a < 256; ++a) Check(Mul(static_cast<Byte>(a),Inv(static_cast<Byte>(a))) == 1,"all nonzero inverses");
    Dimension<2>(); Dimension<3>(); Dimension<4>(); Dimension<5>(); Dimension<6>(); Dimension<8>();
    Check(cases == 3312 && full_echelon_checks == cases && nonidentity_echelon_checks >= 6*255 &&
          baseline_multiplies-candidate_multiplies == eliminations+insertions &&
          candidate_multiplies < baseline_multiplies,"complete six-dimension bounded mechanism roster");
    std::printf("cases=%llu feeds=%llu recoveries=%llu eliminations=%llu insertions=%llu baseline_mul=%llu candidate_mul=%llu actual_mul=%llu full_echelons=%llu nonidentity_echelons=%llu\n",
        static_cast<unsigned long long>(cases),static_cast<unsigned long long>(feeds),
        static_cast<unsigned long long>(recoveries),static_cast<unsigned long long>(eliminations),
        static_cast<unsigned long long>(insertions),static_cast<unsigned long long>(baseline_multiplies),
        static_cast<unsigned long long>(candidate_multiplies),static_cast<unsigned long long>(actual_multiplies),
        static_cast<unsigned long long>(full_echelon_checks),static_cast<unsigned long long>(nonidentity_echelon_checks));
    return 0;
}
