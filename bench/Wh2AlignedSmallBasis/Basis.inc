    // Keep tiny/irregular allocations unchanged. Whole cache-line strides
    // preserve this alignment for every complete source block.
    const size_t extra = profile.block_bytes >= 256u &&
        profile.block_bytes % 64u == 0u ? 63u : 0u;
    if (profile.message_bytes > std::numeric_limits<size_t>::max() - extra) {
        return WirehairV2_UnsupportedPlatform;
    }
    codec->SmallBasis.reset(new (std::nothrow) uint8_t[(size_t)profile.message_bytes + extra]);
    if (!codec->SmallBasis) {
        return WirehairV2_OOM;
    }
    // Ownership stays with the original new[] pointer. The existing evaluator
    // retains only the interior view; its partial-tail padding is unchanged.
    uint8_t* basis = codec->SmallBasis.get();
    const size_t skew = reinterpret_cast<uintptr_t>(basis) % 64u;
    if (extra && skew) {
        basis += 64u - skew;
    }
    std::memcpy(basis, message, (size_t)profile.message_bytes);
    const WirehairV2Result result = MapSmallResult(
        wirehair_small_core::Encoder<K>::Create(lookup, basis,
            profile.message_bytes, profile.block_bytes, codec->SmallEncoder));
