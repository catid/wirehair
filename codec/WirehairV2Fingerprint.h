#pragma once

#include "WirehairV2GF16.h"

#include <stdint.h>

/*
    All-K equation fingerprints for the public serialized V2 profiles.

    A public WIREHAIR_V2_PROFILE_* ID promises byte-identical equations
    forever, but the profile expansion transitively consumes the legacy
    dense-count/dense-seed/peel-seed tables, the block-byte peel-policy
    classification, the matrix/packet seed derivation, the certified precode
    row construction, the completion coefficient tables, and the version-4
    packet-row iterators.  The existing self-rebuild guard in
    WirehairV2Profile.cpp only compares two expansions of the SAME
    implementation, so a coherent future change to any of those inputs would
    silently rewrite every published profile ID.

    This module folds the complete equation-affecting expansion for EVERY
    supported block count K into one SHA-256 per public profile ID.  The
    digests are compared against checked-in golden constants by
    V2FingerprintTest.cpp; any drift under an existing profile ID is a
    compatibility bug, not a golden to update.

    Stream layout (fingerprint version 1).  Every integer is unsigned
    little-endian; SHA-256 output bytes follow the FIPS 180-4 big-endian word
    convention.  All counts precede their arrays, so the stream is
    unambiguous without field tags:

    - preamble: ASCII tag "WH2EQFP1"; profile ID (u64); completion field
      (u32); recovery mix count (u32); precode contract version (u32);
      packet-row contract version (u32); precode seed salt (u64);
      recovery-row seed salt (u64); max packet seed attempts (u32); min K
      (u32); max K (u32); canonical block bytes (u32); probe count (u32) and
      each probe block-byte value (u32).
    - completion coefficients (K-independent):
      - GF256 profiles: heavy row count (u32); coefficient window = 488
        columns spanning two mod-244 periods (u32); HeavyCoefficient(r, c,
        12) bytes row-major over that window.
      - mixed profiles: GF(256) row count, GF(2^16) row count, and
        coefficient period (u32 each); coefficient geometry, residue
        schedule, residue skew, and residue hash seed (u32 each); rotated and
        independent-extension flags (u8 each); the subfield coefficient
        period bytes row-major; the extension coefficient period (u16 each)
        row-major; the subfield and extension residues (u32 pairs) for
        columns [0, 2 * period); the residue block shifts and extension block
        shifts (u32 pairs) for period blocks [0, 4).
    - for each K from min K to max K ascending:
      - K (u32).
      - legacy table outputs at the canonical block bytes: dense count, peel
        seed, dense seed (u32 each).
      - the peel-policy fields folded into MatrixSeedFromProfile(): min
        degree, max degree, solver candidate limit, structure, solver (u32
        each).
      - attempt-0 precode parameters: staircase, dense rows, heavy rows,
        source hits, completion field, heavy family (u32 each); dense
        identity corner (u8); precode seed (u64).
      - attempt-0 packet config: packet peel seed and mix count (u32 each).
      - attempt stepping probes: for attempts 1 and 255, the attempt (u32),
        stepped precode seed (u64), and stepped packet peel seed (u32).
      - the complete attempt-0 BuildPrecodeSystem() structure: staircase row
        count (u32), then each row as length (u32) and its sorted columns
        (u32 each); dense row count (u32), then each row likewise.
      - attempt-0 packet rows: row count (u32), then for each probe block id
        the id (u32), row length (u32), and every generated column (u32).
        The probe ids are 0, 1, K-1, K, K+1, K+7, K+255, 2K+12345, 0x1ffff,
        and 0xfffffffe, covering systematic ids, the repair boundary, and
        deep repair ids.
      - block-byte probes: for each probe value, the probe (u32), the five
        seed-folded policy fields (u32 each), the precode matrix seed (u64),
        and the packet peel seed (u32).  Block bytes reach the equations only
        through these values, so probing them across every byte-class
        boundary pins the full BlockBytes dependency without expanding rows
        at every block size.

    The digest deliberately contains no solver, runtime, tuning, fixup
    diagnostic, or cached-prime state: those may legitimately change without
    changing wire equations.  Every hashed value above feeds the published
    equations directly, so none of them can change under a frozen profile ID.
*/

namespace wirehair_v2 {

/// SHA-256 digest size and its lowercase hex encoding (plus terminator).
static const uint32_t kEquationFingerprintBytes = 32u;
static const uint32_t kEquationFingerprintHexChars = 64u;

/// Complete supported block-count domain (CAT_WIREHAIR_MIN_N/MAX_N).
static const uint32_t kEquationFingerprintMinBlockCount = 2u;
static const uint32_t kEquationFingerprintMaxBlockCount = 64000u;

/// Canonical block bytes for the full per-K row expansion.  Even, so the
/// same expansion domain is valid for the GF(256) and mixed profiles.
static const uint32_t kEquationFingerprintCanonicalBlockBytes = 1280u;

/// One public serialized profile ID and the equation configuration it binds.
/// This mirrors the contract table in WirehairV2Profile.cpp; the fingerprint
/// test cross-checks the IDs against the public header constants.
struct EquationFingerprintContract
{
    uint64_t ProfileId;
    const char* Name;
    CompletionField Completion;
    uint32_t RecoveryMixCount;
};

/// Returns the frozen public profile contracts and their count.
const EquationFingerprintContract* EquationFingerprintContracts(
    uint32_t& count_out);

/// Optional per-K progress callback for long full-range computations.
typedef void (*EquationFingerprintProgress)(
    void* context,
    uint32_t block_count);

/**
    Compute the equation fingerprint for one public profile contract over
    block counts [min_block_count, max_block_count].

    Only the full range produces the frozen golden value; smaller ranges are
    a debugging aid.  Returns false for invalid arguments, arithmetic table
    initialization failure, or any expansion failure (which would itself be
    an equation-contract bug for a supported K).  digest_out receives
    kEquationFingerprintBytes bytes on success and is unchanged on failure.
*/
bool ComputeEquationFingerprint(
    const EquationFingerprintContract& contract,
    uint32_t min_block_count,
    uint32_t max_block_count,
    uint8_t* digest_out,
    EquationFingerprintProgress progress = nullptr,
    void* progress_context = nullptr);

/**
    Format a digest as lowercase hex.  hex_out must hold
    kEquationFingerprintHexChars + 1 bytes and receives a terminated string.
*/
void FormatEquationFingerprintHex(const uint8_t* digest, char* hex_out);

} // namespace wirehair_v2
