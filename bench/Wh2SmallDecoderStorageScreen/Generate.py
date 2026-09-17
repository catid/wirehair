"""Build a new matched-workspace observer from immutable reviewed components."""
import argparse
import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def replace_once(text, old, new):
    if text.count(old) != 1:
        raise ValueError('one exact observer site required: ' + old[:80])
    return text.replace(old, new)


def pinned(name, sha):
    raw = (ROOT/name).read_bytes()
    if hashlib.sha256(raw).hexdigest() != sha:
        raise ValueError('review shared source changes: ' + name)
    return raw.decode()


def generate():
    source = pinned('bench/Wh2SmallPayload2/Screen.cpp',
                    '59a3e1d69fc2230bf1119ec0b8c209ef7f669ad1a613c96ee40a12f23b433494')
    source = source[:source.index('struct Row {')]
    edits = (
        ('#include <wirehair/wirehair.h>', '#include <wirehair/wirehair.h>\n#include <wirehair/wirehair_small.h>'),
        ('Batch = 32, Reps = 12', 'Batch = 64, Reps = 12'),
        ('decltype(&wirehair_v2_encode) encode;', '''decltype(&wirehair_v2_encoder_create_with_options) ordinary;
    decltype(&wirehair_v2_profile_deserialize) parse;
    decltype(&wirehair_small_encoder_create) screate;
    decltype(&wirehair_small_decoder_create) sdecoder;
    decltype(&wirehair_small_encode) sencode;
    decltype(&wirehair_small_decode) sfeed;
    decltype(&wirehair_small_recover) srecover;
    decltype(&wirehair_small_free) sfree;
    decltype(&wirehair_v2_encode) encode;'''),
        ('LOAD(create, "wirehair_v2_encoder_create_profile_id_with_options");', '''LOAD(create, "wirehair_v2_encoder_create_profile_id_with_options");
        LOAD(ordinary, "wirehair_v2_encoder_create_with_options");
        LOAD(parse, "wirehair_v2_profile_deserialize");
        LOAD(screate, "wirehair_small_encoder_create"); LOAD(sdecoder, "wirehair_small_decoder_create");
        LOAD(sencode, "wirehair_small_encode"); LOAD(sfeed, "wirehair_small_decode");
        LOAD(srecover, "wirehair_small_recover"); LOAD(sfree, "wirehair_small_free");'''),
        ('unsigned k, b, tail, policy;', 'unsigned k, b, tail, policy, route;'),
        ('Library* lib; bool wh1;', 'Library* lib; bool wh1, standalone;'),
        ('        WirehairV2EncoderOptions opts = WIREHAIR_V2_ENCODER_OPTIONS_INIT;', '''        if (standalone) {
            const auto result = lib->screate(source, s.Message(), s.b,
                s.policy ? WirehairSmall_BorrowedImmutable : WirehairSmall_Independent, p.data(), p.size());
            h = result.codec; return int(result.status);
        }
        WirehairV2EncoderOptions opts = WIREHAIR_V2_ENCODER_OPTIONS_INIT;'''),
        ('const int r = lib->create(ProfileId(s.k), source, s.Message(), s.b, &opts, p.data(), 32, &n, &c);', '''const uint64_t id = s.route == 4 || s.route == 5 ? WIREHAIR_V2_PROFILE_CERTIFIED_2026_07 : ProfileId(s.k);
        const int r = s.route == 0 ? lib->ordinary(source, s.Message(), s.b, &opts, p.data(), 32, &n, &c) :
            lib->create(id, source, s.Message(), s.b, &opts, p.data(), 32, &n, &c);'''),
        ('        return wh1 ? int(lib->wencode', '''        if (standalone) {
            const auto result = lib->sencode(static_cast<WirehairSmallCodec>(h), id, out, b);
            Check(result.bytes_written <= UINT32_MAX, "standalone packet length bound");
            *n = static_cast<uint32_t>(result.bytes_written); return int(result.status);
        }
        return wh1 ? int(lib->wencode'''),
        ('        WirehairV2Codec c = nullptr; int r = lib->decoder', '''        if (standalone) {
            const auto result = lib->sdecoder(p.data(), p.size()); h = result.codec; return int(result.status);
        }
        WirehairV2Codec c = nullptr; int r = lib->decoder'''),
        ('    { return wh1 ? int(lib->wfeed', '''    { if (standalone) return int(lib->sfeed(static_cast<WirehairSmallCodec>(h), id, in, n));
      return wh1 ? int(lib->wfeed'''),
        ('        if (wh1) return lib->wrecover', '''        if (standalone) {
            const auto result = lib->srecover(static_cast<WirehairSmallCodec>(h), out, bytes);
            return result.status == WirehairSmall_Success && result.bytes_written != bytes ? -1 : int(result.status);
        }
        if (wh1) return lib->wrecover'''),
        ('    { if (wh1) lib->wfree', '''    { if (standalone) lib->sfree(static_cast<WirehairSmallCodec>(h));
      else if (wh1) lib->wfree'''),
        ('void* prebuilt, Byte* output)', 'void* prebuilt, Byte* output, Profile& p)'),
        ('Handle h(api); Profile p = {};', 'Handle h(api); p.fill(0);'),
    )
    for old, new in edits:
        source = replace_once(source, old, new)
    # A literal WHK3 descriptor, independently constructed from the declared
    # shape; never reinterpret WHK3 as a WHV2 profile.
    source += '''Profile SmallDescriptor(const Shape& s)
{
    Profile p = {{'W','H','K','3',1,0,32,0}};
    for (unsigned i = 0; i < 8; ++i) {
        p[8+i] = static_cast<Byte>(WIREHAIR_SMALL_K3_PROFILE_ID >> (8*i));
        p[16+i] = static_cast<Byte>(uint64_t(s.Message()) >> (8*i));
    }
    for (unsigned i = 0; i < 4; ++i) p[24+i] = static_cast<Byte>(s.b >> (8*i));
    return p;
}
'''
    main = pinned('bench/Wh2SmallDormantCoreScreen/Main.inc',
                  'a95e2cfbd93fcf15af1871c33c07a8fa2ae3ecd0201d031847bcb4769865f4ac')
    edits = (
        ('Api apis[3] = {{baseline,false},{candidate,false},{baseline,true}};',
         'Api apis[3] = {{baseline,false,false},{candidate,false,false},{baseline,true,false}};'),
        ('route < 6;', 'route < 7;'),
        ('ks[6] = {3,3,5,8,3,16}', 'ks[7] = {3,3,5,8,3,16,3}'),
        ('f[arm] = Prepare(apis[arm],s);', '''Api api = apis[arm]; api.standalone = route == 6 && !api.wh1;
                    f[arm] = Prepare(api,s);'''),
        ('if (arm < 2) {', '''if (arm < 2 && route == 6) {
                        Check(f[arm].profile == SmallDescriptor(s), "exact standalone WHK3 descriptor");
                    } else if (arm < 2) {'''),
        ('Check(fixtures.size() == 60,', 'Check(fixtures.size() == 70,'),
        ('Work(apis[arm],matched,metric,nullptr,output.data()+16,produced);',
         'Work(api,matched,metric,nullptr,output.data()+16,produced);'),
        ('const auto begin = Now();', '''Api api = apis[arm]; api.standalone = matched.shape.route == 6 && !api.wh1;
                            const auto begin = Now();'''),
        ('PASS 60 fixtures x', 'PASS 70 fixtures x'),
        ('        Publish();\n    } catch', '''        Publish();
        Check(Now()-start < UINT64_C(240000000000), "240-second worker cap includes publication");
    } catch'''),
        ('const bool timing = mode == "run" || mode == "run-reverse";', '''const bool timing = mode == "run" || mode == "run-reverse";
#ifdef STORAGE_SCREEN_SANITIZED
        Check(!timing, "sanitized observer is neutral only");
#endif'''),
    )
    for old, new in edits:
        main = replace_once(main, old, new)
    # All four numeric sites are roster extents, not time caps or sample sizes.
    if main.count('180') != 4:
        raise ValueError('four exact cell-roster extents')
    main = main.replace('180', '210')
    return source + main


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    output = args.output.resolve()
    if ROOT in output.parents:
        raise ValueError('external generated observer only')
    output.write_text(generate())
