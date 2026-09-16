#!/usr/bin/env python3
"""Build the K4 cost observer from an authenticated R1 fresh proof.

This is deliberately a thin, explicit adapter around the byte-pinned R0
observer builder.  It supplies only a proof-validated fresh archive/boundary
root and a new producing provenance record; it never launches the worker.
The R0 ``verify_qualified_library`` guard is not changed.
"""
import argparse
import importlib.util
import inspect
from pathlib import Path
import shlex


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def _load(name, filename):
    spec = importlib.util.spec_from_file_location(name, HERE / filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


C = _load("_k4_fresh_cost_builder_r0", "Wh2K4CostBuildR0.py")
Q = _load("_k4_fresh_qualifier_r1", "Wh2K4FreshQualifierR1.py")


def _pin_records(value):
    """Collect every proof pin, retaining the complete handoff verbatim."""
    records = {}

    def visit(node):
        if isinstance(node, dict):
            if set(node) == {"path", "bytes", "sha256"}:
                record = Q.pin_map([node])[node["path"]]
                old = records.get(node["path"])
                if old is not None:
                    Q.exact(old, record, "conflicting R1 proof pin")
                records[node["path"]] = record
            for child in node.values():
                visit(child)
        elif isinstance(node, list):
            for child in node:
                visit(child)

    visit(value)
    return records


def _proof(path, mode, neutral_root):
    path = Path(path)
    proof = Q.read_json(path, 128 * 1024 ** 2)
    Q.exact(set(proof), {"protocol", "mode", "head", "neutral_result", "neutral_source",
                        "source_pins", "commands", "production", "boundary", "toolchain",
                        "artifacts", "producing_source_closure", "scientific_launch", "scope"},
            "R1 proof schema")
    Q.exact((proof["protocol"], proof["mode"], proof["producing_source_closure"],
             proof["scientific_launch"]), (Q.PROTOCOL, mode, True, False), "R1 proof identity")
    Q.exact(proof["head"], Q.current_head(), "R1 proof current HEAD")
    Q.pin(path)
    neutral_root = Path(neutral_root)
    Q.exact(Path(proof["neutral_result"]["path"]), neutral_root / "RESULT.json",
            "R1 neutral result root")
    Q.exact(Path(proof["neutral_source"]["path"]), neutral_root / "SOURCE.json",
            "R1 neutral source root")
    Q.exact(Q.pin(Path(proof["neutral_result"]["path"])), proof["neutral_result"],
            "R1 neutral result pin")
    Q.exact(Q.pin(Path(proof["neutral_source"]["path"])), proof["neutral_source"],
            "R1 neutral source pin")
    source_map = Q.pin_map(proof["source_pins"])
    Q.exact(source_map, Q.pin_map(proof["source_pins"]), "R1 source pins")
    for record in source_map.values():
        Q.exact(Q.pin(Path(record["path"])), record, "R1 source changed")
    pins = _pin_records(proof)
    pins[str(path)] = Q.pin(path)
    adapter_record = Q.pin(Path(__file__).resolve())
    pins[adapter_record["path"]] = adapter_record
    return proof, pins


def build(mode, proof_path, neutral_root, output):
    Q.require(mode in Q.MODES, "backend")
    proof, pins = _proof(proof_path, mode, neutral_root)
    neutral_root = Path(neutral_root)
    Q.exact(neutral_root, Path(proof["neutral_result"]["path"]).parent,
            "R1 neutral mode root")
    # Bind the R0 builder's fresh paths only in this process.  No R0 source or
    # historical constant is modified on disk.
    C.SMALL = neutral_root.parent.parent / neutral_root.parent.name
    C.U.SMALL = C.SMALL
    C.boundary_recipe = _boundary_recipe
    C.qualified_inputs = lambda requested_mode, unused_output: _inputs(
        requested_mode, proof, pins)
    C.verify_qualified_library = lambda provenance: _verify(provenance, proof)
    # Reuse the reviewed R0 build algorithm but swap only the worker source
    # and its authenticated claim namespace. This leaves the byte-pinned R0
    # module unchanged while binding the fresh worker to the R1 namespace.
    source = inspect.getsource(C.build)
    source = source.replace("HERE/'Wh2K4SerializedCostR0.cpp'",
                            "HERE/'Wh2K4SerializedCostR1.cpp'")
    source = source.replace(
        "A.exact((reader.PROTOCOL, reader.MODES), (PROTOCOL, MODES), 'new cost reader contract')",
        "reader.OUTPUT = Path('/var/tmp/wh2-k4-serialized-cost-r1.R18/science')\n"
        "    A.exact((reader.PROTOCOL, reader.MODES), (PROTOCOL, MODES), 'new cost reader contract')")
    namespace = dict(C.__dict__)
    exec(compile(source, str(HERE / 'Wh2K4FreshCostBuildR1.py'), 'exec'), namespace)
    fresh_build = namespace['build']
    return fresh_build(mode, Path(output))


def _inputs(mode, proof, pins):
    Q.exact(mode, proof["mode"], "R1 mode binding")
    small_archive = Path(proof["boundary"]["archive"]["path"])
    production_archive = Path(proof["production"]["archive"]["path"])
    depfile = Path(proof["boundary"]["depfile"]["path"])
    Q.exact(Q.pin(small_archive), proof["boundary"]["archive"], "R1 boundary archive")
    Q.exact(Q.pin(production_archive), proof["production"]["archive"], "R1 production archive")
    Q.exact(Q.pin(depfile), proof["boundary"]["depfile"], "R1 boundary depfile")
    provenance = dict(protocol=Q.PROTOCOL, mode=mode, producing_source_closure=True,
                      source_pins=proof["source_pins"], input_pins=list(pins.values()),
                      historical_snapshots=[], serialized_dependencies=
                      Q._regular(depfile, 2 * 1024 ** 2).decode(),
                      fresh_proof=proof)
    reused = {Path(record["path"]) for record in pins.values()}
    return [small_archive, production_archive], reused, provenance, {}


def _boundary_recipe(database, build, mode):
    """Parse the current neutral boundary recipe, including scalar's CMake flag."""
    Q.exact(len(database), 8, "eight fresh boundary translation units")
    source = C.HERE / "Wh2SmallSerialized.cpp"
    selected = [entry for entry in database
                if entry.get("file") == str(source)]
    Q.exact(len(selected), 1, "one fresh boundary producer")
    entry = selected[0]
    Q.exact(set(entry), {"directory", "command", "file", "output"},
            "fresh boundary compile schema")
    target = "CMakeFiles/wh2_small_serialized.dir" + str(source) + ".o"
    Q.exact((entry["directory"], entry["output"]), (str(build), target),
            "fresh boundary output binding")
    flags = (["-DANDROID"] if mode == "scalar" else []) + ["-DWH2_SMALL_CODEC_K=4"]
    if mode == "scalar":
        flags.append("-DWH2_SMALL_EXPECT_PORTABLE=1")
    flags += ["-I" + str(C.ROOT), "-I" + str(C.ROOT / "include"), "-I" + str(build)]
    if mode == "scalar":
        flags.append("-DANDROID")
    flags += (["-fsanitize=address,undefined", "-fno-omit-frame-pointer",
               "-march=native", "-g"] if mode == "asan" else ["-O3", "-DNDEBUG"])
    flags += ["-std=c++11", "-fPIC", "-Wall", "-Wextra", "-Wpedantic", "-Werror",
              "-fno-strict-aliasing", "-fno-lto"]
    argv = ["/usr/bin/c++"] + flags + ["-o", target, "-c", str(source)]
    Q.exact(shlex.split(entry["command"]), argv, "fresh boundary compiler recipe")
    return build / target, argv


def _verify(provenance, proof):
    Q.exact(provenance.get("protocol"), Q.PROTOCOL, "R1 producing protocol")
    Q.exact(provenance.get("mode"), proof["mode"], "R1 producing mode")
    Q.exact(provenance.get("producing_source_closure"), True, "R1 producing closure")
    Q.exact(provenance.get("source_pins"), proof["source_pins"], "R1 source handoff")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=Q.MODES)
    parser.add_argument("proof", type=Path)
    parser.add_argument("neutral_mode_root", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args(argv)
    manifest = build(args.mode, args.proof, args.neutral_mode_root, args.output)
    print(Q.canonical({"mode": args.mode, "manifest": manifest,
                       "scientific_launch": False}).decode(), end="", flush=True)


if __name__ == "__main__":
    main()
