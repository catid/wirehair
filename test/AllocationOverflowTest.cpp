#include "../WirehairCodec.h"
#include "../WirehairTools.h"

#include <cinttypes>
#include <cstdio>
#include <cstdlib>

namespace {

unsigned CallocCalls = 0;
size_t LastCallocCount = 0;
size_t LastCallocSize = 0;
int CallocFailureCountdown = -1;

void* RejectingCalloc(size_t count, size_t size)
{
    ++CallocCalls;
    LastCallocCount = count;
    LastCallocSize = size;
    return nullptr;
}

void* ControlledCalloc(size_t count, size_t size)
{
    ++CallocCalls;
    LastCallocCount = count;
    LastCallocSize = size;
    if (CallocFailureCountdown == 0) {
        return nullptr;
    }
    if (CallocFailureCountdown > 0) {
        --CallocFailureCountdown;
    }
    return std::calloc(count, size);
}

class CallocHookGuard
{
public:
    explicit CallocHookGuard(wirehair::detail::SIMDSafeCallocForTesting callback)
    {
        wirehair::detail::SetSIMDSafeCallocForTesting(callback);
    }

    ~CallocHookGuard()
    {
        wirehair::detail::SetSIMDSafeCallocForTesting(nullptr);
    }

private:
    CallocHookGuard(const CallocHookGuard&);
    CallocHookGuard& operator=(const CallocHookGuard&);
};

bool Check(bool condition, const char* description)
{
    if (condition) {
        return true;
    }
    std::fprintf(stderr, "allocation overflow test failed: %s\n", description);
    return false;
}

bool TestCheckedBoundary()
{
    size_t padding = 0;
    if (!Check(
            wirehair::detail::GetSIMDSafeAllocationSize(0, padding),
            "zero-size padding calculation"))
    {
        return false;
    }
    if (!Check(padding > 0, "allocator padding must be non-zero")) {
        return false;
    }

    const size_t largest_safe = SIZE_MAX - padding;
    size_t allocation_size = 0;
    if (!Check(
            wirehair::detail::GetSIMDSafeAllocationSize(
                largest_safe, allocation_size),
            "largest safe request rejected"))
    {
        return false;
    }
    if (!Check(
            allocation_size == SIZE_MAX,
            "largest safe request did not produce SIZE_MAX"))
    {
        return false;
    }

    CallocCalls = 0;
    LastCallocCount = 0;
    LastCallocSize = 0;
    {
        CallocHookGuard hook(RejectingCalloc);
        if (!Check(
                wirehair::SIMDSafeAllocate(largest_safe) == nullptr,
                "largest safe rejected-calloc request returned non-null"))
        {
            return false;
        }
        if (!Check(
                CallocCalls == 1 && LastCallocCount == 1 &&
                    LastCallocSize == SIZE_MAX,
                "largest safe request did not reach calloc exactly once"))
        {
            return false;
        }
    }

    allocation_size = 123;
    CallocCalls = 0;
    {
        CallocHookGuard hook(RejectingCalloc);
        if (!Check(
                !wirehair::detail::GetSIMDSafeAllocationSize(
                    largest_safe + 1, allocation_size),
                "first overflowing request was accepted"))
        {
            return false;
        }
        if (!Check(
                allocation_size == 0,
                "overflowing size calculation did not clear its output"))
        {
            return false;
        }
        if (!Check(
                wirehair::SIMDSafeAllocate(largest_safe + 1) == nullptr,
                "first overflowing allocation returned non-null"))
        {
            return false;
        }
        if (!Check(
                CallocCalls == 0,
                "overflowing allocation called calloc"))
        {
            return false;
        }
    }
    return true;
}

bool TestRealAllocations()
{
    size_t alignment = 0;
    if (!wirehair::detail::GetSIMDSafeAllocationSize(0, alignment)) {
        return Check(false, "failed to obtain allocator alignment");
    }

    const size_t sizes[] = {
        0,
        1,
        alignment - 1,
        alignment,
        alignment + 1,
        4096,
        3u * 1024u * 1024u
    };
    for (size_t case_i = 0;
        case_i < sizeof(sizes) / sizeof(sizes[0]); ++case_i)
    {
        const size_t size = sizes[case_i];
        uint8_t* allocation = wirehair::SIMDSafeAllocate(size);
        if (!Check(allocation != nullptr, "ordinary allocation failed")) {
            return false;
        }
        if (!Check(
                reinterpret_cast<uintptr_t>(allocation) % alignment == 0,
                "ordinary allocation is not aligned"))
        {
            wirehair::SIMDSafeFree(allocation);
            return false;
        }
        for (size_t i = 0; i < size; ++i)
        {
            if (!Check(allocation[i] == 0, "calloc bytes are not zero")) {
                wirehair::SIMDSafeFree(allocation);
                return false;
            }
            allocation[i] = static_cast<uint8_t>(i);
        }
        wirehair::SIMDSafeFree(allocation);
    }

    wirehair::SIMDSafeFree(nullptr);
    return true;
}

void ResetControlledCalloc(int failure_countdown)
{
    CallocCalls = 0;
    LastCallocCount = 0;
    LastCallocSize = 0;
    CallocFailureCountdown = failure_countdown;
}

bool TestAllocationFailureCallers()
{
    if (!Check(gf256_init() == 0, "gf256 initialization failed")) {
        return false;
    }

    uint8_t message[257];
    for (size_t i = 0; i < sizeof(message); ++i) {
        message[i] = static_cast<uint8_t>(i * 29u + 7u);
    }

    // AllocateWorkspace: Failure must leave a reusable, non-operational codec.
    wirehair::Codec workspace_encoder;
    ResetControlledCalloc(0);
    {
        CallocHookGuard hook(ControlledCalloc);
        if (!Check(
                workspace_encoder.InitializeEncoder(sizeof(message), 32) ==
                    Wirehair_OOM && CallocCalls == 1,
                "encoder workspace allocation failure was not transactional"))
        {
            return false;
        }
    }
    if (!Check(
            !workspace_encoder.CanEncode() &&
                workspace_encoder.InitializeEncoder(sizeof(message), 32) ==
                    Wirehair_Success,
            "encoder workspace failure corrupted codec reuse"))
    {
        return false;
    }

    // AllocateInput: Owned encoder input failure must be recoverable.
    wirehair::Codec input_encoder;
    if (!Check(
            input_encoder.InitializeEncoder(sizeof(message), 32) ==
                Wirehair_Success,
            "encoder input test setup failed"))
    {
        return false;
    }
    ResetControlledCalloc(0);
    {
        CallocHookGuard hook(ControlledCalloc);
        if (!Check(
                input_encoder.EncodeFeed(message, true) == Wirehair_OOM &&
                    CallocCalls == 1,
                "encoder input allocation failure was not transactional"))
        {
            return false;
        }
    }
    if (!Check(
            !input_encoder.CanEncode() &&
                input_encoder.InitializeEncoder(sizeof(message), 32) ==
                    Wirehair_Success &&
                input_encoder.EncodeFeed(message, false) == Wirehair_Success,
            "encoder input failure corrupted codec reuse"))
    {
        return false;
    }

    // AllocateMatrix: Solver failure must clear aliases and permit a retry.
    wirehair::Codec matrix_encoder;
    if (!Check(
            matrix_encoder.InitializeEncoder(sizeof(message), 32) ==
                Wirehair_Success,
            "encoder matrix test setup failed"))
    {
        return false;
    }
    ResetControlledCalloc(0);
    {
        CallocHookGuard hook(ControlledCalloc);
        if (!Check(
                matrix_encoder.EncodeFeed(message, false) == Wirehair_OOM &&
                    CallocCalls == 1,
                "matrix allocation failure was not transactional"))
        {
            return false;
        }
    }
    if (!Check(
            !matrix_encoder.CanEncode() &&
                matrix_encoder.InitializeEncoder(sizeof(message), 32) ==
                    Wirehair_Success &&
                matrix_encoder.EncodeFeed(message, false) == Wirehair_Success,
            "matrix allocation failure corrupted codec reuse"))
    {
        return false;
    }

    // Decoder input failure and the later workspace failure exercise both
    // allocation sites and the partial-initialization cleanup/reuse path.
    wirehair::Codec input_decoder;
    ResetControlledCalloc(0);
    {
        CallocHookGuard hook(ControlledCalloc);
        if (!Check(
                input_decoder.InitializeDecoder(sizeof(message), 32) ==
                    Wirehair_OOM && CallocCalls == 1,
                "decoder input allocation failure was not transactional"))
        {
            return false;
        }
    }
    if (!Check(
            !input_decoder.CanDecode() &&
                input_decoder.InitializeDecoder(sizeof(message), 32) ==
                    Wirehair_Success,
            "decoder input failure corrupted codec reuse"))
    {
        return false;
    }

    wirehair::Codec workspace_decoder;
    ResetControlledCalloc(1);
    {
        CallocHookGuard hook(ControlledCalloc);
        if (!Check(
                workspace_decoder.InitializeDecoder(sizeof(message), 32) ==
                    Wirehair_OOM && CallocCalls == 2,
                "decoder workspace allocation failure was not transactional"))
        {
            return false;
        }
    }
    if (!Check(
            !workspace_decoder.CanDecode() &&
                workspace_decoder.InitializeDecoder(sizeof(message), 32) ==
                    Wirehair_Success,
            "decoder workspace failure corrupted codec reuse"))
    {
        return false;
    }

    return true;
}

bool TestI686CodecReachability()
{
#if SIZE_MAX == UINT32_MAX
    static const uint32_t kEncoderN = 77;
    static const uint32_t kEncoderBlockBytes = 7354383;
    static const uint32_t kDecoderN = 1170;
    static const uint32_t kDecoderBlockBytes = 2561032;

    const auto is_expected_size_failure = [](WirehairResult result) {
        return result == Wirehair_OOM || result == Wirehair_InvalidInput;
    };

    wirehair::Codec encoder;
    encoder.OverrideSeeds(500, 0, 0);
    WirehairResult result = encoder.InitializeEncoder(
        static_cast<uint64_t>(kEncoderN) * kEncoderBlockBytes,
        kEncoderBlockBytes);
    if (!Check(
            is_expected_size_failure(result),
            "i686 encoder overflow case did not fail safely"))
    {
        return false;
    }
    if (!Check(
            !encoder.CanEncode() && !encoder.CanDecode(),
            "failed i686 encoder remained usable"))
    {
        return false;
    }
    if (!Check(
            encoder.InitializeEncoder(2, 1) == Wirehair_Success &&
                encoder.CanEncode(),
            "encoder was corrupted by rejected overflow case"))
    {
        return false;
    }

    wirehair::Codec decoder;
    decoder.OverrideSeeds(500, 0, 0);
    result = decoder.InitializeDecoder(
        static_cast<uint64_t>(kDecoderN) * kDecoderBlockBytes,
        kDecoderBlockBytes);
    if (!Check(
            is_expected_size_failure(result),
            "i686 decoder overflow case did not fail safely"))
    {
        return false;
    }
    if (!Check(
            !decoder.CanEncode() && !decoder.CanDecode(),
            "failed i686 decoder remained usable"))
    {
        return false;
    }
    if (!Check(
            decoder.InitializeDecoder(2, 1) == Wirehair_Success &&
                decoder.CanDecode(),
            "decoder was corrupted by rejected overflow case"))
    {
        return false;
    }
#endif
    return true;
}

} // namespace

int main()
{
    if (!TestCheckedBoundary() ||
        !TestRealAllocations() ||
        !TestAllocationFailureCallers() ||
        !TestI686CodecReachability())
    {
        return 1;
    }

    std::printf(
        "allocation overflow boundary and %zu-bit reachability tests passed\n",
        sizeof(size_t) * 8);
    return 0;
}
