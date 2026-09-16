    explicit PublicCodec(uint8_t small_k = 0) : SmallK(small_k)
    {
        // Small profiles never use the certified core. Preserve its storage
        // and every header offset, but begin its lifetime only when needed.
        if (!SmallK) {
            new (&Impl) wirehair_v2::Codec;
        }
    }
    ~PublicCodec()
    {
        if (!SmallK) {
            Impl.~Codec();
        }
    }
    PublicCodec(const PublicCodec&) = delete;
    PublicCodec& operator=(const PublicCodec&) = delete;
    PublicCodec(PublicCodec&&) = delete;
    PublicCodec& operator=(PublicCodec&&) = delete;

    // SmallK is immutable; an inactive core is never read, destroyed or reused.
    union { wirehair_v2::Codec Impl; };
