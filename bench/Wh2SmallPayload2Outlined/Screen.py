#!/usr/bin/env python3
"""Single fresh screen: outlined helper plus matched per-cell workspaces."""
import argparse
import importlib.util
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('small_payload_screen', HERE.parent/'Wh2SmallPayload2'/'Screen.py')
S = importlib.util.module_from_spec(spec)
spec.loader.exec_module(S)
OUTPUT = Path('/var/tmp/wh2-small-payload2-outlined-screen-r0')


def run(build):
    cache = (build/'CMakeCache.txt').read_text()
    S.require('CMAKE_HOME_DIRECTORY:INTERNAL='+str(HERE) in cache.splitlines(), 'outlined experiment build')
    # These are additional diagnostic sanity checks, not a claim of transitive
    # build provenance. The shared runner pins the generated worker and DSOs.
    worker = (build/'candidate'/'Screen.cpp').read_text()
    S.require('output[side]' not in worker and 'Fixture matched = fs[0];' in worker and
              'Work(apis[arm], matched,' in worker and 'fs[a].shape,matched.source.data(),p' in worker,
              'matched workspace worker')
    S.run(build, output=OUTPUT, here=HERE, protocol='wh2-small-payload2-outlined-screen-r0')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--run', type=Path)
    group.add_argument('--analyze', type=Path)
    args = parser.parse_args()
    if args.run:
        run(args.run)
    else:
        print(json.dumps(S.replay(args.analyze), sort_keys=True))
