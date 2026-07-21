#include "WirehairV2Fingerprint.h"

#include "WirehairV2Plan.h"
#include "WirehairV2Precode.h"
#include "WirehairV2PrecodeEncode.h"
#include "WirehairV2Seeds.h"
#include "WirehairV2Solve.h"

#include "../WirehairTools.h"
#include "../gf256.h"

#include <wirehair/wirehair.h>

#include <stddef.h>

#include <vector>

namespace wirehair_v2 {
namespace {

static_assert(kEquationFingerprintMinBlockCount == CAT_WIREHAIR_MIN_N &&
    kEquationFingerprintMaxBlockCount == CAT_WIREHAIR_MAX_N,
    "the fingerprint domain must cover every supported block count");
static_assert(kEquationFingerprintHexChars ==
    2u * kEquationFingerprintBytes,
    "hex encoding must cover the complete digest");

//------------------------------------------------------------------------------
// SHA-256 (FIPS 180-4)
//
// The fingerprint is a compatibility constant shared with external release
// records, so it uses standard SHA-256 rather than a project-local hash.
// This is a self-contained portable implementation; performance is bounded
// by row generation, not hashing.

const uint32_t kSha256K[64] = {
    0x428a2f98u, 0x71374491u, 0xb5c0fbcfu, 0xe9b5dba5u,
    0x3956c25bu, 0x59f111f1u, 0x923f82a4u, 0xab1c5ed5u,
    0xd807aa98u, 0x12835b01u, 0x243185beu, 0x550c7dc3u,
    0x72be5d74u, 0x80deb1feu, 0x9bdc06a7u, 0xc19bf174u,
    0xe49b69c1u, 0xefbe4786u, 0x0fc19dc6u, 0x240ca1ccu,
    0x2de92c6fu, 0x4a7484aau, 0x5cb0a9dcu, 0x76f988dau,
    0x983e5152u, 0xa831c66du, 0xb00327c8u, 0xbf597fc7u,
    0xc6e00bf3u, 0xd5a79147u, 0x06ca6351u, 0x14292967u,
    0x27b70a85u, 0x2e1b2138u, 0x4d2c6dfcu, 0x53380d13u,
    0x650a7354u, 0x766a0abbu, 0x81c2c92eu, 0x92722c85u,
    0xa2bfe8a1u, 0xa81a664bu, 0xc24b8b70u, 0xc76c51a3u,
    0xd192e819u, 0xd6990624u, 0xf40e3585u, 0x106aa070u,
    0x19a4c116u, 0x1e376c08u, 0x2748774cu, 0x34b0bcb5u,
    0x391c0cb3u, 0x4ed8aa4au, 0x5b9cca4fu, 0x682e6ff3u,
    0x748f82eeu, 0x78a5636fu, 0x84c87814u, 0x8cc70208u,
    0x90befffau, 0xa4506cebu, 0xbef9a3f7u, 0xc67178f2u
};

inline uint32_t Rotr32(uint32_t x, unsigned bits)
{
    return (x >> bits) | (x << (32u - bits));
}

class Sha256
{
public:
    Sha256()
    {
        State[0] = 0x6a09e667u;
        State[1] = 0xbb67ae85u;
        State[2] = 0x3c6ef372u;
        State[3] = 0xa54ff53au;
        State[4] = 0x510e527fu;
        State[5] = 0x9b05688cu;
        State[6] = 0x1f83d9abu;
        State[7] = 0x5be0cd19u;
        TotalBytes = 0u;
        BufferBytes = 0u;
    }

    void Update(const void* data, size_t bytes)
    {
        const uint8_t* input = reinterpret_cast<const uint8_t*>(data);
        TotalBytes += bytes;
        if (BufferBytes != 0u)
        {
            while (bytes != 0u && BufferBytes < 64u) {
                Buffer[BufferBytes++] = *input++;
                --bytes;
            }
            if (BufferBytes < 64u) {
                return;
            }
            ProcessBlock(Buffer);
            BufferBytes = 0u;
        }
        while (bytes >= 64u) {
            ProcessBlock(input);
            input += 64u;
            bytes -= 64u;
        }
        while (bytes != 0u) {
            Buffer[BufferBytes++] = *input++;
            --bytes;
        }
    }

    void Finalize(uint8_t* digest_out)
    {
        const uint64_t total_bits = TotalBytes * 8u;
        const uint8_t pad_byte = 0x80u;
        Update(&pad_byte, 1u);
        const uint8_t zero = 0u;
        while (BufferBytes != 56u) {
            Update(&zero, 1u);
        }
        uint8_t length_bytes[8];
        for (unsigned i = 0; i < 8u; ++i) {
            length_bytes[i] = (uint8_t)(total_bits >> (56u - 8u * i));
        }
        Update(length_bytes, sizeof(length_bytes));
        for (unsigned i = 0; i < 8u; ++i)
        {
            digest_out[4u * i] = (uint8_t)(State[i] >> 24);
            digest_out[4u * i + 1u] = (uint8_t)(State[i] >> 16);
            digest_out[4u * i + 2u] = (uint8_t)(State[i] >> 8);
            digest_out[4u * i + 3u] = (uint8_t)State[i];
        }
    }

private:
    void ProcessBlock(const uint8_t* block)
    {
        uint32_t w[64];
        for (unsigned i = 0; i < 16u; ++i)
        {
            w[i] = ((uint32_t)block[4u * i] << 24) |
                ((uint32_t)block[4u * i + 1u] << 16) |
                ((uint32_t)block[4u * i + 2u] << 8) |
                (uint32_t)block[4u * i + 3u];
        }
        for (unsigned i = 16u; i < 64u; ++i)
        {
            const uint32_t s0 = Rotr32(w[i - 15u], 7u) ^
                Rotr32(w[i - 15u], 18u) ^ (w[i - 15u] >> 3);
            const uint32_t s1 = Rotr32(w[i - 2u], 17u) ^
                Rotr32(w[i - 2u], 19u) ^ (w[i - 2u] >> 10);
            w[i] = w[i - 16u] + s0 + w[i - 7u] + s1;
        }

        uint32_t a = State[0], b = State[1], c = State[2], d = State[3];
        uint32_t e = State[4], f = State[5], g = State[6], h = State[7];
        for (unsigned i = 0; i < 64u; ++i)
        {
            const uint32_t s1 = Rotr32(e, 6u) ^ Rotr32(e, 11u) ^
                Rotr32(e, 25u);
            const uint32_t ch = (e & f) ^ (~e & g);
            const uint32_t temp1 = h + s1 + ch + kSha256K[i] + w[i];
            const uint32_t s0 = Rotr32(a, 2u) ^ Rotr32(a, 13u) ^
                Rotr32(a, 22u);
            const uint32_t maj = (a & b) ^ (a & c) ^ (b & c);
            const uint32_t temp2 = s0 + maj;
            h = g;
            g = f;
            f = e;
            e = d + temp1;
            d = c;
            c = b;
            b = a;
            a = temp1 + temp2;
        }
        State[0] += a; State[1] += b; State[2] += c; State[3] += d;
        State[4] += e; State[5] += f; State[6] += g; State[7] += h;
    }

    uint32_t State[8];
    uint64_t TotalBytes;
    uint8_t Buffer[64];
    size_t BufferBytes;
};

//------------------------------------------------------------------------------
// Canonical little-endian field emission

struct StreamHasher
{
    Sha256 Digest;

    void Bytes(const void* data, size_t bytes)
    {
        Digest.Update(data, bytes);
    }
    void U8(uint8_t value)
    {
        Digest.Update(&value, 1u);
    }
    void U16(uint16_t value)
    {
        const uint8_t bytes[2] = {
            (uint8_t)value, (uint8_t)(value >> 8)
        };
        Digest.Update(bytes, sizeof(bytes));
    }
    void U32(uint32_t value)
    {
        const uint8_t bytes[4] = {
            (uint8_t)value, (uint8_t)(value >> 8),
            (uint8_t)(value >> 16), (uint8_t)(value >> 24)
        };
        Digest.Update(bytes, sizeof(bytes));
    }
    void U64(uint64_t value)
    {
        U32((uint32_t)value);
        U32((uint32_t)(value >> 32));
    }
    void Row(const std::vector<uint32_t>& columns)
    {
        U32((uint32_t)columns.size());
        for (uint32_t column : columns) {
            U32(column);
        }
    }
};

/// Exactly the policy fields MatrixSeedFromProfile() folds into the matrix
/// seed.  The remaining PeelPolicy fields never reach a public V2 equation.
void HashSeedFoldedPolicy(StreamHasher& hasher, const PeelPolicy& policy)
{
    hasher.U32(policy.Codec.MinDegree);
    hasher.U32(policy.Codec.MaxDegree);
    hasher.U32(policy.Codec.SolverCandidateLimit);
    hasher.U32((uint32_t)policy.Structure);
    hasher.U32((uint32_t)policy.Solver);
}

/// Block-byte probe values spanning every peel-policy byte class and both
/// class boundaries (12 KiB and 320 KiB), plus the extreme legal sizes.
/// Mixed profiles accept only even block sizes, so their probes shift the
/// odd values to the adjacent even ones.
const uint32_t kProbeBlockBytesGF256[] = {
    1u, 2u, 16u, 1280u, 12287u, 12288u,
    327679u, 327680u, 1048576u, 2147483647u
};
const uint32_t kProbeBlockBytesMixed[] = {
    2u, 4u, 16u, 1280u, 12286u, 12288u,
    327678u, 327680u, 1048576u, 2147483646u
};

/// Attempts pinned in addition to attempt zero.  The stepping rules are
/// linear in the attempt, so the first step plus the last supported attempt
/// freeze the complete deterministic attempt domain.
const uint32_t kSteppedAttempts[] = { 1u, 255u };

/// Two full mod-244 coefficient periods, so the window wrap rule itself is
/// part of the digest rather than an assumption.
const uint32_t kCoefficientWindowPeriods = 2u;

const EquationFingerprintContract kContracts[] = {
    {
        WIREHAIR_V2_PROFILE_CERTIFIED_2026_07,
        "certified_2026_07",
        CompletionField::GF256,
        kCertifiedPacketMixCount
    },
    {
        WIREHAIR_V2_PROFILE_MIXED_2026_07,
        "mixed_2026_07",
        CompletionField::MixedGF256GF16,
        kCertifiedPacketMixCount
    },
    {
        WIREHAIR_V2_PROFILE_MIXED_MIX2_2026_07,
        "mixed_mix2_2026_07",
        CompletionField::MixedGF256GF16,
        2u
    }
};

bool HashCompletionCoefficients(
    StreamHasher& hasher,
    const EquationFingerprintContract& contract)
{
    if (contract.Completion == CompletionField::GF256)
    {
        const uint32_t heavy_rows = 12u;
        const uint32_t window =
            kCoefficientWindowPeriods * (256u - heavy_rows);
        hasher.U32(heavy_rows);
        hasher.U32(window);
        for (uint32_t row = 0; row < heavy_rows; ++row) {
            for (uint32_t column = 0; column < window; ++column) {
                hasher.U8(HeavyCoefficient(row, column, heavy_rows));
            }
        }
        return true;
    }

    const MixedCoefficientRows* rows = GetMixedCoefficientRows();
    if (!rows) {
        return false;
    }
    const uint32_t gf256_rows = ActiveMixedGF256Rows();
    const uint32_t gf16_rows = ActiveMixedGF16Rows();
    const uint32_t period = ActiveMixedCoefficientPeriod();
    if (gf256_rows > kMixedGF256RowsMax || gf16_rows > kMixedGF16RowsMax ||
        period == 0u || period > kMixedCoefficientPeriod)
    {
        return false;
    }
    hasher.U32(gf256_rows);
    hasher.U32(gf16_rows);
    hasher.U32(period);
    hasher.U32((uint32_t)ActiveMixedCoefficientGeometry());
    hasher.U32((uint32_t)ActiveMixedResidueSchedule());
    hasher.U32(ActiveMixedResidueSkew());
    hasher.U32(ActiveMixedResidueHashSeed());
    hasher.U8(ActiveMixedResiduesRotated() ? 1u : 0u);
    hasher.U8(ActiveMixedIndependentExtensionResidues() ? 1u : 0u);
    for (uint32_t row = 0; row < gf256_rows; ++row) {
        for (uint32_t residue = 0; residue < period; ++residue) {
            hasher.U8(rows->Subfield[row][residue]);
        }
    }
    for (uint32_t row = 0; row < gf16_rows; ++row) {
        for (uint32_t residue = 0; residue < period; ++residue) {
            hasher.U16(rows->Extension[row][residue]);
        }
    }
    const uint32_t residue_window = kCoefficientWindowPeriods * period;
    for (uint32_t column = 0; column < residue_window; ++column)
    {
        hasher.U32(ActiveMixedCoefficientResidue(column));
        hasher.U32(ActiveMixedExtensionCoefficientResidue(column));
    }
    for (uint32_t block = 0; block < 4u; ++block)
    {
        hasher.U32(ActiveMixedResidueBlockShift(block));
        hasher.U32(ActiveMixedExtensionResidueBlockShift(block));
    }
    return true;
}

bool HashBlockCount(
    StreamHasher& hasher,
    const EquationFingerprintContract& contract,
    uint32_t block_count,
    const uint32_t* probes,
    uint32_t probe_count)
{
    const bool mixed =
        contract.Completion == CompletionField::MixedGF256GF16;

    // Mirror ExpandProfile() in WirehairV2Profile.cpp statement for
    // statement: legacy table selection, seed derivation, contract binding,
    // and attempt stepping must digest exactly what the public path derives.
    const SeedProfile base = SelectSeedProfile(
        block_count, kEquationFingerprintCanonicalBlockBytes);
    hasher.U32(block_count);
    hasher.U32(base.DenseCount);
    hasher.U32(base.PeelSeed);
    hasher.U32(base.DenseSeed);
    HashSeedFoldedPolicy(hasher, base.Policy);

    const uint64_t matrix_seed = MatrixSeedFromProfile(
        base, 0u, kMessagePrecodeSeedSalt);
    PrecodeParams params = mixed ?
        MakeMixedParams(block_count, matrix_seed) :
        MakeCertifiedParams(block_count, matrix_seed);
    params.Staircase = base.DenseCount;
    params.DenseIdentityCorner = false;

    PacketRowConfig packet;
    packet.PeelSeed = PacketPeelSeedFromProfile(
        base, kMessageRecoveryRowSeedSalt);
    packet.MixCount = contract.RecoveryMixCount;

    const PrecodeParams params0 = PrecodeParamsForAttempt(params, 0u);
    const PacketRowConfig packet0 = PacketConfigForAttempt(packet, 0u);
    hasher.U32(params0.Staircase);
    hasher.U32(params0.DenseRows);
    hasher.U32(params0.HeavyRows);
    hasher.U32(params0.SourceHits);
    hasher.U32((uint32_t)params0.Field);
    hasher.U32((uint32_t)params0.HeavyFamily);
    hasher.U8(params0.DenseIdentityCorner ? 1u : 0u);
    hasher.U64(params0.Seed);
    hasher.U32(packet0.PeelSeed);
    hasher.U32(packet0.MixCount);

    for (uint32_t attempt : kSteppedAttempts)
    {
        hasher.U32(attempt);
        hasher.U64(PrecodeParamsForAttempt(params, attempt).Seed);
        hasher.U32(PacketConfigForAttempt(packet, attempt).PeelSeed);
    }

    PrecodeSystem system;
    if (!BuildPrecodeSystem(params0, system)) {
        return false;
    }
    hasher.U32((uint32_t)system.StaircaseRows.size());
    for (const std::vector<uint32_t>& row : system.StaircaseRows) {
        hasher.Row(row);
    }
    hasher.U32((uint32_t)system.DenseRowColumns.size());
    for (const std::vector<uint32_t>& row : system.DenseRowColumns) {
        hasher.Row(row);
    }

    const uint64_t precode_count_wide = (uint64_t)params0.Staircase +
        params0.DenseRows + params0.HeavyRows;
    if (precode_count_wide > UINT32_MAX) {
        return false;
    }
    const uint32_t precode_count = (uint32_t)precode_count_wide;
    PacketRowRuntime runtime;
    if (!runtime.Initialize(block_count, precode_count, packet0.MixCount)) {
        return false;
    }
    const uint32_t packet_ids[] = {
        0u, 1u, block_count - 1u, block_count, block_count + 1u,
        block_count + 7u, block_count + 255u,
        2u * block_count + 12345u, 0x1ffffu, 0xfffffffeu
    };
    hasher.U32((uint32_t)(sizeof(packet_ids) / sizeof(packet_ids[0])));
    for (uint32_t packet_id : packet_ids)
    {
        const std::vector<uint32_t> row = GeneratePacketMatrixRowWithRuntime(
            block_count, precode_count, packet_id, packet0, runtime);
        if (row.empty()) {
            return false;
        }
        hasher.U32(packet_id);
        hasher.Row(row);
    }

    for (uint32_t probe_index = 0; probe_index < probe_count; ++probe_index)
    {
        const uint32_t probe_bytes = probes[probe_index];
        const SeedProfile probe_profile =
            SelectSeedProfile(block_count, probe_bytes);
        hasher.U32(probe_bytes);
        HashSeedFoldedPolicy(hasher, probe_profile.Policy);
        hasher.U64(MatrixSeedFromProfile(
            probe_profile, 0u, kMessagePrecodeSeedSalt));
        hasher.U32(PacketPeelSeedFromProfile(
            probe_profile, kMessageRecoveryRowSeedSalt));
    }
    return true;
}

} // namespace

const EquationFingerprintContract* EquationFingerprintContracts(
    uint32_t& count_out)
{
    count_out = (uint32_t)(sizeof(kContracts) / sizeof(kContracts[0]));
    return kContracts;
}

bool ComputeEquationFingerprint(
    const EquationFingerprintContract& contract,
    uint32_t min_block_count,
    uint32_t max_block_count,
    uint8_t* digest_out,
    EquationFingerprintProgress progress,
    void* progress_context)
{
    if (!digest_out ||
        min_block_count < kEquationFingerprintMinBlockCount ||
        max_block_count > kEquationFingerprintMaxBlockCount ||
        min_block_count > max_block_count ||
        (contract.Completion != CompletionField::GF256 &&
         contract.Completion != CompletionField::MixedGF256GF16))
    {
        return false;
    }
    if (gf256_init() != 0) {
        return false;
    }

    const bool mixed =
        contract.Completion == CompletionField::MixedGF256GF16;
    const uint32_t* probes =
        mixed ? kProbeBlockBytesMixed : kProbeBlockBytesGF256;
    const uint32_t probe_count = mixed ?
        (uint32_t)(sizeof(kProbeBlockBytesMixed) /
            sizeof(kProbeBlockBytesMixed[0])) :
        (uint32_t)(sizeof(kProbeBlockBytesGF256) /
            sizeof(kProbeBlockBytesGF256[0]));

    StreamHasher hasher;
    static const char kStreamTag[8] = {
        'W', 'H', '2', 'E', 'Q', 'F', 'P', '1'
    };
    hasher.Bytes(kStreamTag, sizeof(kStreamTag));
    hasher.U64(contract.ProfileId);
    hasher.U32((uint32_t)contract.Completion);
    hasher.U32(contract.RecoveryMixCount);
    hasher.U32(PrecodeContractVersion(contract.Completion));
    hasher.U32(kPacketRowContractVersion);
    hasher.U64(kMessagePrecodeSeedSalt);
    hasher.U64(kMessageRecoveryRowSeedSalt);
    hasher.U32(kMaxPacketSeedAttempts);
    hasher.U32(min_block_count);
    hasher.U32(max_block_count);
    hasher.U32(kEquationFingerprintCanonicalBlockBytes);
    hasher.U32(probe_count);
    for (uint32_t probe_index = 0; probe_index < probe_count; ++probe_index) {
        hasher.U32(probes[probe_index]);
    }
    if (!HashCompletionCoefficients(hasher, contract)) {
        return false;
    }

    for (uint32_t block_count = min_block_count; ; ++block_count)
    {
        if (!HashBlockCount(
                hasher, contract, block_count, probes, probe_count))
        {
            return false;
        }
        if (progress) {
            progress(progress_context, block_count);
        }
        if (block_count == max_block_count) {
            break;
        }
    }

    hasher.Digest.Finalize(digest_out);
    return true;
}

void FormatEquationFingerprintHex(const uint8_t* digest, char* hex_out)
{
    static const char kHexDigits[] = "0123456789abcdef";
    for (uint32_t i = 0; i < kEquationFingerprintBytes; ++i)
    {
        hex_out[2u * i] = kHexDigits[digest[i] >> 4];
        hex_out[2u * i + 1u] = kHexDigits[digest[i] & 0xfu];
    }
    hex_out[kEquationFingerprintHexChars] = '\0';
}

} // namespace wirehair_v2
