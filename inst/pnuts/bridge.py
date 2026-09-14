#!/usr/bin/env python3
"""BridgeStan compilation/initialization and native PNUTS process orchestration.

Sampling uses no R/Python target-gradient callback or sampler implementation.
The PNUTS executable is installed separately. This helper requires BridgeStan 2.9.
"""
import concurrent.futures
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time

import bridgestan as bs
from bridgestan.compile import get_bridgestan_path
import numpy as np


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def compile_model(request):
    make = ["STAN_THREADS=true", *request.get("make_args", [])]
    stanc = request.get("stanc_args", [])
    source = Path(get_bridgestan_path()).resolve()
    key_data = dict(code=request["code"], bridgestan=bs.__version__, source=str(source),
                    make=make, stanc=stanc, platform=platform.platform(),
                    compiler=os.environ.get("CXX", ""))
    # Local compiler settings affect the model and therefore the cache key.
    for name in ("make/local", "Makefile"):
        file = source / name
        key_data[name] = hashlib.sha256(file.read_bytes()).hexdigest() if file.exists() else None
    key = hashlib.sha256(json.dumps(key_data, sort_keys=True).encode()).hexdigest()
    folder = Path(request["cache_dir"]) / key
    folder.mkdir(parents=True, exist_ok=True)
    lock = folder / ".compile-lock"
    start = time.monotonic()
    while True:
        try:
            lock.mkdir()
            break
        except FileExistsError:
            if time.monotonic() - start > 600:
                raise RuntimeError(f"Compilation lock did not clear: {lock}")
            time.sleep(0.1)
    try:
        stan = folder / "model.stan"
        stan.write_text(request["code"])
        library = folder / "model_model.so"
        if not library.exists():
            library = Path(bs.compile_model(stan, make_args=make, stanc_args=stanc))
        return dict(library=str(library.resolve()), stan_file=str(stan.resolve()),
                    code=request["code"], bridgestan=bs.__version__, cache_key=key)
    finally:
        lock.rmdir()


def parse_model(request):
    stan = Path(request["stan_file"])
    command = [str(Path(get_bridgestan_path()) / "bin/stanc"), str(stan), "--o=" + os.devnull]
    result = subprocess.run(command, text=True, capture_output=True)
    if result.returncode:
        raise RuntimeError(result.stderr + result.stdout)
    return dict(code=request["code"])


def constrain_dict(model, q):
    groups = {}
    for name, value in zip(model.param_names(), model.param_constrain(q)):
        name, *index = name.split(".")
        groups.setdefault(name, []).append((tuple(int(x) - 1 for x in index), float(value)))
    output = {}
    for name, values in groups.items():
        if not values[0][0]:
            output[name] = values[0][1]
        else:
            dims = tuple(max(index[j] for index, _ in values) + 1 for j in range(len(values[0][0])))
            array = np.zeros(dims)
            for index, value in values:
                array[index] = value
            output[name] = array.tolist()
    return output


def initialize(model, spec, seed, output):
    rng = np.random.default_rng(seed)
    evaluations = 0
    radius = spec.get("radius", 2.0)
    for attempt in range(100 if spec["kind"] == "random" else 1):
        q = np.zeros(model.param_unc_num()) if spec["kind"] == "zero" else rng.uniform(-radius, radius, model.param_unc_num())
        try:
            if spec["kind"] == "constrained":
                initial = constrain_dict(model, q)
                supplied = spec["values"]
                unknown = set(supplied) - set(initial)
                if unknown:
                    raise ValueError("Unknown initialization parameters: " + ", ".join(sorted(unknown)))
                for name, value in supplied.items():
                    reference = np.asarray(initial[name])
                    value = np.asarray(value)
                    # JSON auto-unboxing loses a length-one R vector's shape.
                    if value.shape != reference.shape:
                        if value.size == reference.size == 1:
                            value = value.reshape(reference.shape)
                        else:
                            raise ValueError(f"Initialization shape for {name}: expected {reference.shape}, got {value.shape}")
                    initial[name] = value.tolist()
                q = model.param_unconstrain_json(initial)
            evaluations += 1
            lp, gradient = model.log_density_gradient(q, propto=True, jacobian=True)
            if not np.isfinite(lp) or not np.isfinite(gradient).all():
                raise ValueError("Nonfinite initial target density or gradient")
            np.savetxt(output, q, fmt="%.17g")
            return dict(path=str(output), gradient_evaluations=evaluations, attempts=attempt + 1)
        except (ValueError, RuntimeError) as error:
            if spec["kind"] != "random":
                raise
            last_error = str(error)
    raise RuntimeError("No target-valid random initialization in 100 attempts: " + last_error)


def sample(request):
    os.environ.update(STAN_NUM_THREADS=str(request["threads"]), OMP_NUM_THREADS="1")
    folder = Path(request["output_dir"])
    folder.mkdir(parents=True, exist_ok=True)
    model = bs.StanModel(request["library"], request["data_file"], seed=request["model_seed"])
    if model.param_unc_num() == 0:
        raise ValueError("PNUTS requires at least one unconstrained parameter")
    commands = []
    initializations = []
    for i, (seed, initial) in enumerate(zip(request["seeds"], request["initializations"]), 1):
        path = folder / f"chain-{i}.csv"
        if path.exists() or Path(str(path) + ".partial").exists():
            raise FileExistsError(f"Refusing to overwrite PNUTS chain output: {path}")
        init_path = folder / f"init-{i}.txt"
        initializations.append(initialize(model, initial, seed, init_path))
        command = [request["executable"], "sample", "--model", request["library"],
                   "--data", request["data_file"], "--output", str(path), "--init", str(init_path),
                   "--draws", str(request["draws"]), "--warmup", str(request["warmup"]),
                   "--seed", str(seed), "--model-seed", str(request["model_seed"]),
                   "--gq-seed", str((seed + 1000003) % 2147483647), "--include-tp", "--include-gq",
                   *request["sampler_args"]]
        if request["save_warmup"]:
            command.append("--save-warmup")
        commands.append(command)
    manifest = dict(commands=commands, initializations=initializations, completed=False,
                    bridgestan=bs.__version__, dimension=model.param_unc_num(),
                    executable_sha256=hashlib.sha256(Path(request["executable"]).read_bytes()).hexdigest(),
                    library_sha256=hashlib.sha256(Path(request["library"]).read_bytes()).hexdigest(),
                    statuses=[])
    write_json(folder / "run.json", manifest)

    def run_chain(i):
        log = folder / f"chain-{i+1}.log"
        start = time.monotonic()
        with log.open("w") as stream:
            stream.write(json.dumps(commands[i]) + "\n")
            stream.flush()
            try:
                result = subprocess.run(commands[i], stdout=stream, stderr=subprocess.STDOUT,
                                        timeout=request["chain_timeout"], check=False)
                return dict(chain=i+1, returncode=result.returncode, timed_out=False,
                            seconds=time.monotonic()-start, log=str(log))
            except subprocess.TimeoutExpired:
                return dict(chain=i+1, returncode=None, timed_out=True,
                            seconds=time.monotonic()-start, log=str(log))

    with concurrent.futures.ThreadPoolExecutor(max_workers=request["cores"]) as pool:
        futures = {pool.submit(run_chain, i): i for i in range(len(commands))}
        for future in concurrent.futures.as_completed(futures):
            status = future.result()
            manifest["statuses"].append(status)
            write_json(folder / "run.json", manifest)
            print(f"PNUTS chain {status['chain']}: " + ("complete" if status["returncode"] == 0 else "failed"), flush=True)
    manifest["statuses"].sort(key=lambda x: x["chain"])
    failed = [x for x in manifest["statuses"] if x["returncode"] != 0]
    if failed:
        messages = [f"Chain {x['chain']} ({'timeout' if x['timed_out'] else 'process error'}): {x['log']}\n" +
                    Path(x["log"]).read_text()[-4000:] for x in failed]
        raise RuntimeError("PNUTS fit failed; no partial ensemble was returned.\n" + "\n".join(messages))
    manifest["completed"] = True
    manifest["files"] = [str(folder / f"chain-{i+1}.csv") for i in range(len(commands))]
    write_json(folder / "run.json", manifest)
    return manifest


def main():
    request = json.loads(Path(sys.argv[1]).read_text())
    if request.get("bridgestan_path"):
        bs.set_bridgestan_path(request["bridgestan_path"])
    if not bs.__version__.startswith("2.9."):
        raise RuntimeError("This experimental adapter requires BridgeStan Python 2.9.x")
    result = {"compile": compile_model, "parse": parse_model, "sample": sample}[request["action"]](request)
    write_json(request["response_file"], result)


if __name__ == "__main__":
    main()
