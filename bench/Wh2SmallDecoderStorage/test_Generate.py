import hashlib
from pathlib import Path
import unittest
import Generate as G

ROOT = Path(__file__).resolve().parents[2]


class GeneratorTest(unittest.TestCase):
    def test_candidate_is_decoder_only(self):
        raw = (ROOT/'codec/WirehairSmallCore.h').read_bytes()
        generated = G.candidate(raw)
        prefix = '// Same private core contract'
        self.assertEqual(raw.decode().split(prefix)[0], generated.split(prefix)[0])
        for token in ('::operator new(sizeof(Decoder) + slab_bytes, std::nothrow)',
                      '::new (storage) Decoder(', 'static void operator delete(void* storage) noexcept',
                      'max() - sizeof(Decoder)', 'std::uint8_t* slab_ = nullptr;'):
            self.assertIn(token, generated)
        self.assertNotIn('slab_.get()', generated)
        self.assertNotIn('slab_.reset(', generated)
        self.assertEqual(generated[generated.index('    Result Feed('):generated.index('private:\n    Decoder(')],
                         raw.decode()[raw.decode().index('    Result Feed('):raw.decode().index('private:\n    Decoder(')]
                         .replace('slab_.get()', 'slab_'))

    def test_hash_and_site_fail_closed(self):
        for fn, path in ((G.candidate, 'codec/WirehairSmallCore.h'),
                         (G.public_test, 'test/V2SmallCodecTest.cpp'),
                         (G.small_test, 'test/SmallCodecTest.cpp')):
            with self.assertRaises(ValueError):
                fn((ROOT/path).read_bytes() + b'\n')
        for source in ('none', 'aa'):
            with self.assertRaises(ValueError): G.replace_once(source, 'a', 'b')

    def test_public_test_changes_only_exact_allocation_count(self):
        raw = (ROOT/'test/V2SmallCodecTest.cpp').read_bytes()
        self.assertEqual(hashlib.sha256(raw).hexdigest(), G.TEST_SHA)
        self.assertEqual(G.public_test(raw).replace('count == 2, "decoder allocations"',
                                                  'count == 3, "decoder allocations"'), raw.decode())

    def test_standalone_test_changes_only_decoder_allocations(self):
        raw = (ROOT/'test/SmallCodecTest.cpp').read_bytes()
        expected = raw.decode().replace('for (size_t failure = 0; failure < 3; ++failure)',
                                        'for (size_t failure = 0; failure < 2; ++failure)')
        expected = expected.replace('Stop() == 3 && d.status == WirehairSmall_Success',
                                    'Stop() == 2 && d.status == WirehairSmall_Success')
        expected = expected.replace('for (unsigned i = 0; i < 3; ++i)', 'for (unsigned i = 0; i < 2; ++i)')
        self.assertEqual(G.small_test(raw), expected)


if __name__ == '__main__': unittest.main()
