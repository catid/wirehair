#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
PYTHON=${PYTHON:-python3}
cache=$(mktemp -d /tmp/wh2-thermal-python-cache.XXXXXX)

cleanup()
{
    rm -rf -- "$cache"
}
trap cleanup EXIT HUP INT TERM

export PYTHONDONTWRITEBYTECODE=1
export PYTHONPYCACHEPREFIX=$cache

"$PYTHON" -m py_compile \
    "$ROOT/bench/wirehair_expo_thermal_sampler.py" \
    "$ROOT/bench/test_wirehair_expo_thermal_sampler.py"
"$PYTHON" -m tabnanny \
    "$ROOT/bench/wirehair_expo_thermal_sampler.py" \
    "$ROOT/bench/test_wirehair_expo_thermal_sampler.py"
"$PYTHON" "$ROOT/bench/test_wirehair_expo_thermal_sampler.py"
