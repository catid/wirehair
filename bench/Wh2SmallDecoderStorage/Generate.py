"""Generate one benchmark-only decoder storage change; no codec execution."""
import hashlib
from pathlib import Path

CORE_SHA = '5b0acdd096d24b76351bacd1718c44ca5b37d4df587fe7334822a9f61f1e0b8c'
TEST_SHA = '61041a29d5464e4f824b45cfcd6e26ddb779567a7899a07c300e10fb6159aa8e'
SMALL_TEST_SHA = '11c57d20ffbfbcc68285e1bad650eb780838e648c79ec45a8610955b002f80a1'


def replace_once(source, old, new):
    if source.count(old) != 1:
        raise ValueError('expected exactly one overlay site: ' + old[:80])
    return source.replace(old, new)


def candidate(raw):
    if hashlib.sha256(raw).hexdigest() != CORE_SHA:
        raise ValueError('re-audit storage overlay after core changes')
    text = raw.decode()
    text = replace_once(text, '''        std::unique_ptr<Decoder> decoder(new (std::nothrow) Decoder(lookup, static_cast<std::size_t>(message), block, tail));
        if (!decoder) return Status::OutOfMemory;
        decoder->slab_.reset(new (std::nothrow) std::uint8_t[decoder->slab_bytes_]);
        if (!decoder->slab_) return Status::OutOfMemory;''', '''        const std::size_t slab_bytes = static_cast<std::size_t>(block) * (K + 1);
        if (slab_bytes > std::numeric_limits<std::size_t>::max() - sizeof(Decoder))
            return Status::OutOfMemory;
        void* storage = ::operator new(sizeof(Decoder) + slab_bytes, std::nothrow);
        if (!storage) return Status::OutOfMemory;
        // Constructor is nonthrowing. The trailing byte storage belongs to the
        // same allocation; only the Decoder object is constructed in its prefix.
        std::unique_ptr<Decoder> decoder(::new (storage) Decoder(lookup, static_cast<std::size_t>(message), block, tail));
        decoder->slab_ = ::new (static_cast<std::uint8_t*>(storage) + sizeof(Decoder)) std::uint8_t[slab_bytes];
        if (!decoder->slab_) return Status::OutOfMemory;''')
    text = replace_once(text, '    ~Decoder() = default;', '''    ~Decoder() = default;
    // The allocation includes trailing payload bytes, so a sized delete using
    // sizeof(Decoder) would be incorrect. Match the raw scalar allocation.
    static void operator delete(void* storage) noexcept { ::operator delete(storage); }''')
    text = replace_once(text,
                        'Decoder(Lookup lookup, std::size_t message, std::uint32_t block, std::uint32_t tail)\n',
                        'Decoder(Lookup lookup, std::size_t message, std::uint32_t block, std::uint32_t tail) noexcept\n')
    text = replace_once(text, '    std::unique_ptr<std::uint8_t[]> slab_;',
                        '    std::uint8_t* slab_ = nullptr;')
    if text.count('slab_.get()') != 4:
        raise ValueError('exact decoder slab accesses')
    text = text.replace('slab_.get()', 'slab_')
    return text


def public_test(raw):
    if hashlib.sha256(raw).hexdigest() != TEST_SHA:
        raise ValueError('re-audit allocation assertion after public test changes')
    return replace_once(raw.decode(), 'count == 3, "decoder allocations"',
                        'count == 2, "decoder allocations"')


def small_test(raw):
    if hashlib.sha256(raw).hexdigest() != SMALL_TEST_SHA:
        raise ValueError('re-audit standalone facade test after source changes')
    text = raw.decode()
    for old, new in (
            ('for (size_t failure = 0; failure < 3; ++failure)',
             'for (size_t failure = 0; failure < 2; ++failure)'),
            ('Stop() == 3 && d.status == WirehairSmall_Success',
             'Stop() == 2 && d.status == WirehairSmall_Success'),
            ('for (unsigned i = 0; i < 3; ++i)',
             'for (unsigned i = 0; i < 2; ++i)')):
        text = replace_once(text, old, new)
    return text


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    root, output = args.root.resolve(), args.output.resolve()
    if output == root or root in output.parents:
        raise ValueError('external generated outputs only')
    core = candidate((root/'codec/WirehairSmallCore.h').read_bytes())
    test = public_test((root/'test/V2SmallCodecTest.cpp').read_bytes())
    small = small_test((root/'test/SmallCodecTest.cpp').read_bytes())
    output.mkdir(exist_ok=True)
    (output/'WirehairSmallCore.h').write_text(core)
    (output/'V2SmallCodecTest.cpp').write_text(test)
    (output/'SmallCodecTest.cpp').write_text(small)
