"""Generate only the proven small-decoder unit-diagonal substitutions."""
import argparse
import hashlib
from pathlib import Path

CORE_SHA = '5b0acdd096d24b76351bacd1718c44ca5b37d4df587fe7334822a9f61f1e0b8c'
EDITS = (
    ('            for (unsigned c = p; c < K; ++c) row[c] ^= gf256_mul(coefficients_[p][c], factor);',
     '''            // Retained pivots have unit diagonal; factor ^ factor is zero.
            row[p] = 0;
            for (unsigned c = p + 1; c < K; ++c) row[c] ^= gf256_mul(coefficients_[p][c], factor);'''),
    ('        for (unsigned c = pivot; c < K; ++c) row[c] = gf256_mul(row[c], inverse);',
     '''        // The nonzero pivot times its already-computed inverse is one.
        row[pivot] = 1;
        for (unsigned c = pivot + 1; c < K; ++c) row[c] = gf256_mul(row[c], inverse);'''),
)


def replace_once(source, old, new):
    if source.count(old) != 1:
        raise ValueError('one exact unit-diagonal site required')
    return source.replace(old,new)


def candidate(raw):
    if hashlib.sha256(raw).hexdigest() != CORE_SHA:
        raise ValueError('re-audit unit-diagonal proof after core changes')
    source = raw.decode()
    for old,new in EDITS:
        source = replace_once(source,old,new)
    return source


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root',type=Path)
    parser.add_argument('output',type=Path)
    args = parser.parse_args()
    root,output = args.root.resolve(),args.output.resolve()
    if output == root or root in output.parents:
        raise ValueError('external generated output only')
    source = candidate((root/'codec/WirehairSmallCore.h').read_bytes())
    output.mkdir(exist_ok=True)
    (output/'WirehairSmallCore.h').write_text(source)
