#!/usr/bin/env python3
"""Isolated PEP 517, installation, discovery, and uninstall E2E."""

import argparse
from email.parser import BytesParser
import hashlib
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import venv
import zipfile


VERSION = "2.0.0"
DIST_NAME = "wirehair-fec"
SOURCE_DATE_EPOCH = "1700000000"
BUILD_FRONTEND = "build==1.2.2.post1"


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--library", type=Path, required=True)
    parser.add_argument("--native-prefix", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    return parser.parse_args()


def display(command):
    return " ".join(repr(str(item)) for item in command)


def run(command, *, cwd, env=None, capture=False):
    command = [str(item) for item in command]
    print("+ " + display(command), flush=True)
    completed = subprocess.run(
        command,
        cwd=str(cwd),
        env=env,
        check=False,
        stdout=subprocess.PIPE if capture else None,
        stderr=subprocess.STDOUT if capture else None,
        text=True,
        errors="replace",
    )
    if completed.returncode != 0:
        if completed.stdout:
            sys.stdout.write(completed.stdout)
        raise subprocess.CalledProcessError(completed.returncode, command)
    if capture:
        return completed.stdout.strip()
    return completed


def python_in(environment):
    if os.name == "nt":
        return environment / "Scripts" / "python.exe"
    return environment / "bin" / "python"


def make_venv(path):
    if path.exists():
        shutil.rmtree(path)
    venv.EnvBuilder(with_pip=True, symlinks=os.name != "nt").create(path)
    return python_in(path)


def clean_environment(**updates):
    environment = os.environ.copy()
    for name in (
            "PYTHONHOME", "PYTHONPATH", "WIREHAIR_LIBRARY",
            "WIREHAIR_PREFIX"):
        environment.pop(name, None)
    environment.update({name: str(value) for name, value in updates.items()})
    return environment


def copy_distribution_source(source, destination):
    destination.mkdir(parents=True)
    for relative in (".gitignore", "pyproject.toml", "LICENSE", "README.md"):
        target = destination / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source / relative, target)
    for relative in (
            "python/README.md",
            "python/benchmark_buffers.py",
            "python/whirehair.py"):
        target = destination / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source / relative, target)
    shutil.copytree(
        source / "python" / "wirehair",
        destination / "python" / "wirehair",
        ignore=shutil.ignore_patterns(
            "__pycache__", "*.pyc", "*.pyo", "*.so", "*.dylib", "*.dll"),
    )
    # A detached-tree native artifact must not change either pure-Python
    # distribution.  This catches wheel-target exclusion regressions instead
    # of relying on the checkout to happen to be binary-free.
    (destination / "python" / "wirehair" /
     "must-not-ship.so").write_bytes(b"native sentinel")

    timestamp = int(SOURCE_DATE_EPOCH)
    for path in sorted(destination.rglob("*"), reverse=True):
        os.utime(path, (timestamp, timestamp), follow_symlinks=False)
    os.utime(destination, (timestamp, timestamp), follow_symlinks=False)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def single(directory, pattern):
    matches = sorted(directory.glob(pattern))
    if len(matches) != 1:
        raise RuntimeError(
            "expected one %s under %s, found %d" %
            (pattern, directory, len(matches)))
    return matches[0]


def contains_path(directory, path):
    try:
        path.relative_to(directory)
        return True
    except ValueError:
        return False


def validate_artifacts(wheel, sdist, source):
    expected_wheel = "wirehair_fec-%s-py3-none-any.whl" % VERSION
    expected_sdist = "wirehair_fec-%s.tar.gz" % VERSION
    if wheel.name != expected_wheel or sdist.name != expected_sdist:
        raise RuntimeError(
            "unexpected distribution filenames: %s %s" %
            (wheel.name, sdist.name))

    with zipfile.ZipFile(wheel) as archive:
        names = set(archive.namelist())
        expected = {
            "whirehair.py",
            "wirehair/__init__.py",
            "wirehair_fec-%s.dist-info/METADATA" % VERSION,
            "wirehair_fec-%s.dist-info/WHEEL" % VERSION,
            "wirehair_fec-%s.dist-info/licenses/LICENSE" % VERSION,
            "wirehair_fec-%s.dist-info/RECORD" % VERSION,
        }
        if names != expected:
            raise RuntimeError(
                "wheel payload mismatch: missing=%s unexpected=%s" %
                (sorted(expected - names), sorted(names - expected)))
        forbidden_suffixes = (".so", ".dylib", ".dll", ".a", ".lib")
        if any(name.lower().endswith(forbidden_suffixes) for name in names):
            raise RuntimeError("pure wheel unexpectedly bundles a native library")

        metadata_name = "wirehair_fec-%s.dist-info/METADATA" % VERSION
        metadata = BytesParser().parsebytes(archive.read(metadata_name))
        if metadata["Name"] != DIST_NAME or metadata["Version"] != VERSION:
            raise RuntimeError("wheel name/version metadata mismatch")
        if metadata["Requires-Python"] != ">=3.8":
            raise RuntimeError("wheel Requires-Python metadata mismatch")
        if metadata["License-Expression"] != "BSD-3-Clause":
            raise RuntimeError("wheel SPDX license metadata mismatch")
        if "LICENSE" not in metadata.get_all("License-File", []):
            raise RuntimeError("wheel license-file metadata is missing")
        wheel_text = archive.read(
            "wirehair_fec-%s.dist-info/WHEEL" % VERSION).decode("ascii")
        if "Root-Is-Purelib: true" not in wheel_text or "Tag: py3-none-any" not in wheel_text:
            raise RuntimeError("wheel is not tagged as portable pure Python")
        installed_license = archive.read(
            "wirehair_fec-%s.dist-info/licenses/LICENSE" % VERSION)
        if installed_license != (source / "LICENSE").read_bytes():
            raise RuntimeError("wheel license text differs from the project license")

    with tarfile.open(sdist, "r:gz") as archive:
        names = set(archive.getnames())
    root = "wirehair_fec-%s/" % VERSION
    expected = {
        root + "pyproject.toml",
        root + ".gitignore",
        root + "LICENSE",
        root + "README.md",
        root + "python/README.md",
        root + "python/benchmark_buffers.py",
        root + "python/whirehair.py",
        root + "python/wirehair/__init__.py",
        root + "PKG-INFO",
    }
    if names != expected:
        raise RuntimeError(
            "sdist payload mismatch: missing=%s unexpected=%s" %
            (sorted(expected - names), sorted(names - expected)))
    if any(name.lower().endswith((".so", ".dylib", ".dll")) for name in names):
        raise RuntimeError("sdist unexpectedly contains a native binary")


def assert_imports_absent(python, cwd):
    code = """
import importlib.metadata
import importlib.util
assert importlib.util.find_spec('wirehair') is None
assert importlib.util.find_spec('whirehair') is None
try:
    importlib.metadata.version('wirehair-fec')
except importlib.metadata.PackageNotFoundError:
    pass
else:
    raise AssertionError('wirehair-fec metadata survived uninstall')
"""
    run([python, "-I", "-c", code], cwd=cwd, env=clean_environment())


def check_installed_imports(python, cwd, library, prefix=None):
    environment = clean_environment()
    if prefix is None:
        environment["WIREHAIR_LIBRARY"] = str(library)
    else:
        environment["WIREHAIR_PREFIX"] = str(prefix)
    code = """
import importlib.metadata
from pathlib import Path
import sys
import wirehair
import whirehair
assert wirehair.Encoder is whirehair.Encoder
assert wirehair.Decoder is whirehair.Decoder
assert wirehair.__version__ == %r
assert importlib.metadata.version('wirehair-fec') == %r
Path(wirehair.__file__).resolve().relative_to(Path(sys.prefix).resolve())
native = wirehair.initialize()
assert Path(native._name).resolve() == Path(%r).resolve()
""" % (VERSION, VERSION, str(library))
    run([python, "-I", "-c", code], cwd=cwd, env=environment)


def install_and_test_wheel(
        python, wheel, source, library, native_prefix, work):
    run([python, "-m", "pip", "install", "--no-index", "--no-deps", wheel],
        cwd=work)
    check_installed_imports(python, work, library)
    check_installed_imports(python, work, library, prefix=native_prefix)

    module_dir = Path(run([
        python, "-I", "-c",
        "import pathlib, whirehair; print(pathlib.Path(whirehair.__file__).parent)",
    ], cwd=work, env=clean_environment(WIREHAIR_LIBRARY=library),
        capture=True)).resolve()
    run([
        python, source / "ci" / "python_native_test.py",
        "--module-dir", module_dir,
        "--library", library,
    ], cwd=work, env=clean_environment(WIREHAIR_LIBRARY=library))

    missing = work / "missing-native-library"
    code = """
import wirehair
try:
    wirehair.initialize()
except OSError as error:
    text = str(error)
    assert 'WIREHAIR_LIBRARY' in text and %r in text
else:
    raise AssertionError('missing authoritative library unexpectedly loaded')
""" % str(missing)
    run([python, "-I", "-c", code], cwd=work,
        env=clean_environment(WIREHAIR_LIBRARY=missing))

    run([python, "-m", "pip", "uninstall", "-y", DIST_NAME], cwd=work)
    assert_imports_absent(python, work)


def test_source_and_editable(python, source, library, work):
    environment = clean_environment(
        PYTHONPATH=source / "python",
        WIREHAIR_LIBRARY=library,
    )
    code = """
from pathlib import Path
import wirehair
import whirehair
assert wirehair.Encoder is whirehair.Encoder
Path(wirehair.__file__).resolve().relative_to(Path(%r).resolve())
wirehair.initialize()
""" % str(source)
    run([python, "-c", code], cwd=work, env=environment)

    run([python, "-m", "pip", "install", "--no-deps", "-e", source],
        cwd=work, env=clean_environment())
    editable_code = """
import importlib.metadata
from pathlib import Path
import wirehair
import whirehair
assert wirehair.Encoder is whirehair.Encoder
assert wirehair.__version__ == %r
assert importlib.metadata.version('wirehair-fec') == %r
Path(wirehair.__file__).resolve().relative_to(Path(%r).resolve())
wirehair.initialize()
""" % (VERSION, VERSION, str(source))
    run([python, "-I", "-c", editable_code], cwd=work,
        env=clean_environment(WIREHAIR_LIBRARY=library))
    run([python, "-m", "pip", "uninstall", "-y", DIST_NAME], cwd=work)
    assert_imports_absent(python, work)


def main():
    args = parse_args()
    source = args.source_dir.resolve()
    library = args.library.resolve()
    native_prefix = args.native_prefix.resolve()
    work = args.work_dir.resolve()
    if not source.is_dir():
        raise RuntimeError("source directory does not exist: " + str(source))
    if not library.is_file():
        raise RuntimeError("native library does not exist: " + str(library))
    if not native_prefix.is_dir():
        raise RuntimeError("native prefix does not exist: " + str(native_prefix))
    for protected in (source, native_prefix, library):
        if contains_path(work, protected):
            raise RuntimeError(
                "refusing to remove work directory %s because it contains %s" %
                (work, protected))
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)

    frontend_python = make_venv(work / "frontend-venv")
    run([
        frontend_python, "-m", "pip", "install",
        "--disable-pip-version-check", BUILD_FRONTEND,
    ], cwd=work)

    build_environment = clean_environment(
        SOURCE_DATE_EPOCH=SOURCE_DATE_EPOCH,
        PYTHONHASHSEED="0",
        TZ="UTC",
    )
    detached_source = work / "detached-source"
    copy_distribution_source(source, detached_source)
    # Compare a build of the real checkout with a normalized, minimal source
    # copy.  This simultaneously proves byte reproducibility and catches a
    # backend that accidentally reads VCS state or unlisted checkout files.
    artifact_sets = []
    for iteration, build_source in enumerate((source, detached_source), 1):
        output = work / ("dist-%d" % iteration)
        output.mkdir()
        run([
            frontend_python, "-m", "build", "--sdist", "--wheel",
            "--outdir", output, build_source,
        ], cwd=work, env=build_environment)
        artifact_sets.append((
            single(output, "*.whl"),
            single(output, "*.tar.gz"),
        ))

    first_wheel, first_sdist = artifact_sets[0]
    second_wheel, second_sdist = artifact_sets[1]
    if digest(first_wheel) != digest(second_wheel):
        raise RuntimeError("wheel build is not byte-for-byte reproducible")
    if digest(first_sdist) != digest(second_sdist):
        raise RuntimeError("sdist build is not byte-for-byte reproducible")
    validate_artifacts(first_wheel, first_sdist, source)

    wheel_python = make_venv(work / "wheel-venv")
    install_and_test_wheel(
        wheel_python, first_wheel, source, library, native_prefix, work)

    sdist_python = make_venv(work / "sdist-venv")
    run([
        sdist_python, "-m", "pip", "install", "--no-deps", first_sdist,
    ], cwd=work, env=clean_environment())
    check_installed_imports(sdist_python, work, library)
    run([sdist_python, "-m", "pip", "uninstall", "-y", DIST_NAME], cwd=work)
    assert_imports_absent(sdist_python, work)

    editable_source = work / "editable-source"
    copy_distribution_source(source, editable_source)
    editable_python = make_venv(work / "editable-venv")
    test_source_and_editable(
        editable_python, editable_source, library, work)

    print(
        "Python package E2E passed: wheel_sha256=%s sdist_sha256=%s" %
        (digest(first_wheel), digest(first_sdist)),
        flush=True,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
