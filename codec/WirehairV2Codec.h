#pragma once

#include "WirehairV2PrecodeEncode.h"
#include "WirehairV2Seeds.h"

#include "../WirehairCodec.h"

#include <memory>

namespace wirehair_v2 {

class Codec
{
public:
    Codec();
    ~Codec();

    WirehairResult InitializeEncoder(
        const void* message,
        uint64_t message_bytes,
        uint32_t block_bytes,
        const SeedProfile* seed_override = 0);

    WirehairResult InitializePrecodeEncoder(
        const void* message,
        uint64_t message_bytes,
        uint32_t block_bytes,
        const SeedProfile* seed_override = 0,
        const MessagePrecodeEncoderOptions* options = 0);

    WirehairResult InitializePrecodeDecoder(
        uint64_t message_bytes,
        uint32_t block_bytes,
        const SeedProfile* seed_override = 0,
        const MessagePrecodeEncoderOptions* options = 0);

    WirehairResult InitializeDecoder(
        uint64_t message_bytes,
        uint32_t block_bytes,
        const SeedProfile* seed_override = 0);

    WirehairResult Encode(
        uint32_t block_id,
        void* block_out,
        uint32_t out_bytes,
        uint32_t* data_bytes_out);

    WirehairResult Decode(
        uint32_t block_id,
        const void* block_in,
        uint32_t block_bytes);

    WirehairResult Recover(void* message_out, uint64_t message_bytes);

    const SeedProfile& Profile() const;

private:
    wirehair::Codec Impl;
    std::unique_ptr<MessagePrecodeEncoder> PrecodeImpl;
    std::unique_ptr<MessagePrecodeDecoder> PrecodeDecoderImpl;
    SeedProfile CurrentProfile = {};
};

} // namespace wirehair_v2
