#!/usr/bin/env python3
"""Re-segment existing reference-masked target CNR measurements with CNVkit 0.9.10.

Uses an explicitly supplied native CNVkit Python runtime; installs nothing.
Replays three original hybrid CNRs before processing all 23 target-only subsets.
Original inputs remain unchanged. Large intermediate CNRs and logs stay in tmp.
"""
import argparse
import concurrent.futures
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
csv.field_size_limit(10_000_000)
PARAMETERS = dict(method="cbs", threshold=0.0001, skip_low=True,
                  skip_outliers=10, min_weight=0, processes=1, smooth_cbs=False)
NATIVE_WORKER = r'''
import logging, sys
from cnvlib import segmentation
from cnvlib.cmdutil import read_cna
from skgenome import tabio
logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
cnarr = read_cna(sys.argv[1])
result = segmentation.do_segmentation(cnarr, "cbs", threshold=0.0001,
    skip_low=True, skip_outliers=10, min_weight=0, processes=1,
    smooth_cbs=False, rscript_path=sys.argv[3])
tabio.write(result, sys.argv[2])
'''
RUNTIME_PROBE = r'''
import hashlib, importlib.metadata, inspect, json
from pathlib import Path
from cnvlib import segmentation, smoothing
from cnvlib.segmentation import cbs
from skgenome import gary
files = {}
for module in (segmentation, cbs, smoothing, gary):
    p = Path(inspect.getfile(module))
    files[module.__name__] = dict(path=str(p), sha256=hashlib.sha256(p.read_bytes()).hexdigest())
versions = {}
for name in ("CNVkit", "numpy", "pandas", "scipy", "biopython", "pysam", "pyfaidx"):
    versions[name] = importlib.metadata.version(name)
print(json.dumps(dict(versions=versions, source_files=files)))
'''


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for b in iter(lambda: handle.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def relative(path):
    return str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)


def read_rows(path, delimiter="\t"):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def write_csv(path, rows):
    with path.open("w", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=rows[0], lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def source_path(path):
    # OVCAN_DATA names the archived data root (contents: "cnvkit wes - new", ...).
    p = Path(path)
    data_root = os.environ.get("OVCAN_DATA")
    prefix = Path("judy_archive/data")
    if data_root and p.is_relative_to(prefix):
        return (Path(data_root) / p.relative_to(prefix)).resolve()
    return (ROOT / p).resolve()


def run_native(python, rscript, source, destination, log):
    started = time.monotonic()
    with log.open("w") as handle:
        subprocess.run([str(python), "-c", NATIVE_WORKER, str(source), str(destination), str(rscript)],
                       check=True, stdout=handle, stderr=subprocess.STDOUT)
    values = read_rows(destination)
    assert values and all(math.isfinite(float(r["log2"])) for r in values)
    assert all(int(r["end"]) > int(r["start"]) and int(r["probes"]) > 0 for r in values)
    return values, time.monotonic() - started


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--python", type=Path,
                        default=os.environ.get("OVCAN_CNV_PYTHON", os.environ.get("OVCAN_CNVKIT_PYTHON", sys.executable)),
                        help="Python interpreter with CNVkit exactly 0.9.10 installed")
    parser.add_argument("--rscript", type=Path, default=shutil.which("Rscript"))
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--scratch-dir", type=Path, default=ROOT / "tmp/wes_cnv_target_only_2026-09-06")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "output/wes_cnv_target_only")
    parser.add_argument("--baseline-only", action="store_true")
    args = parser.parse_args()
    if not args.python or not args.rscript:
        parser.error("Supply --python/OVCAN_CNV_PYTHON and an available Rscript")
    # Keep the virtual-environment interpreter's symlink path: resolving it to
    # the base runtime would discard the selected environment/site-packages.
    python, rscript = args.python.absolute(), args.rscript.absolute()
    scratch, out = args.scratch_dir.resolve(), args.output_dir.resolve()
    scratch.mkdir(parents=True, exist_ok=True)
    out.mkdir(parents=True, exist_ok=True)
    runtime = json.loads(subprocess.check_output([str(python), "-c", RUNTIME_PROBE], text=True))
    assert runtime["versions"]["CNVkit"] == "0.9.10", runtime
    r_versions = subprocess.check_output([str(rscript), "--vanilla", "-e",
        'cat(paste0("R=", getRversion(), "\\nDNAcopy=", packageVersion("DNAcopy"), "\\n"))'], text=True).strip()
    models = read_rows(ROOT / "output/wes_recovered_provenance_cnv_support.csv", ",")
    assert len(models) == 23 and len({r["cell_line"] for r in models}) == 23
    runtime.update(python=str(python), rscript=str(rscript), r_versions=r_versions,
                   parameters=PARAMETERS, segmentation_script_sha256=sha(Path(__file__)),
                   original_input_manifest="output/wes_recovered_provenance_cnv_support.csv",
                   original_input_manifest_sha256=sha(ROOT / "output/wes_recovered_provenance_cnv_support.csv"),
                   source_data_override=os.environ.get("OVCAN_DATA"),
                   subset_rule="gene != Antitarget; target rows retained verbatim; native CNVkit applies low-depth, outlier and weight filters",
                   interpretation="Relative target-supported segments; native CBS is rerun. No alignment, reference normalization, absolute CN calling or baseline source is changed.")
    (out / "runtime.json").write_text(json.dumps(runtime, indent=2) + "\n")
    for model in models:
        for key in ("cnr", "cns"):
            path = source_path(model[key + "_source"])
            assert path.exists() and sha(path) == model[key + "_sha256"], path

    def replay(model):
        name = model["cell_line"]
        target = scratch / f"{name}.baseline.cns"
        log = scratch / f"{name}.baseline.log"
        replayed, elapsed = run_native(python, rscript, source_path(model["cnr_source"]), target, log)
        original = read_rows(source_path(model["cns_source"]))
        n = min(len(original), len(replayed))
        coords = sum(all(a[k] == b[k] for k in ("chromosome", "start", "end")) for a, b in zip(original, replayed))
        probes = sum(int(a["probes"]) == int(b["probes"]) for a, b in zip(original, replayed))
        max_log2 = max(abs(float(a["log2"]) - float(b["log2"])) for a, b in zip(original, replayed))
        okay = len(original) == len(replayed) == coords == probes and max_log2 <= 0.000011
        result = dict(cell_line=name, archived_segments=len(original), replay_segments=len(replayed),
                      ordered_matching_coordinates=coords, ordered_matching_probe_counts=probes,
                      max_ordered_log2_difference=max_log2, baseline_replay_passes=okay,
                      max_ordered_weight_difference=max(abs(float(a["weight"]) - float(b["weight"])) for a, b in zip(original, replayed)),
                      max_ordered_depth_difference=max(abs(float(a["depth"]) - float(b["depth"])) for a, b in zip(original, replayed)),
                      archived_cns_sha256=model["cns_sha256"], replay_cns_sha256=sha(target),
                      replay_log=relative(log), elapsed_seconds=elapsed)
        print(json.dumps(result), flush=True)
        return result

    baseline_models = [r for r in models if r["cell_line"] in ("OV2295", "TOV2835EP", "OV1369-R2")]
    assert len(baseline_models) == 3
    with concurrent.futures.ThreadPoolExecutor(max_workers=min(2, args.workers)) as pool:
        baseline = list(pool.map(replay, baseline_models))
    write_csv(out / "baseline_replay.csv", baseline)
    assert all(r["baseline_replay_passes"] for r in baseline), "Baseline replay mismatch: inspect baseline_replay.csv before any target-only segmentation"
    if args.baseline_only:
        print("All three baseline replays passed; target-only processing not requested in --baseline-only mode.")
        return

    def process(model):
        name = model["cell_line"]
        original = source_path(model["cnr_source"])
        subset = scratch / f"{name}.targetonly.cnr"
        target_count, positive_count = 0, 0
        with original.open() as handle, subset.open("w", newline="") as dest:
            header = handle.readline()
            columns = header.rstrip("\r\n").split("\t")
            gene_i, depth_i = columns.index("gene"), columns.index("depth")
            dest.write(header)
            for line in handle:
                fields = line.rstrip("\r\n").split("\t")
                if fields[gene_i] != "Antitarget":
                    dest.write(line)
                    target_count += 1
                    positive_count += float(fields[depth_i]) > 0
        assert target_count == int(model["target_bins_in_cnr"]) == 204706
        assert positive_count == int(model["positive_target_bins_in_cnr"])
        destination = out / f"{name}.cns"
        pending = scratch / f"{name}.targetonly.cns"
        log = scratch / f"{name}.targetonly.log"
        segments, elapsed = run_native(python, rscript, subset, pending, log)
        assert all("Antitarget" not in r["gene"].split(",") for r in segments)
        # Native outlier filtering may remove additional bins; it must not add bins.
        probes = sum(int(r["probes"]) for r in segments)
        assert probes <= positive_count
        shutil.copyfile(pending, destination)
        result = dict(cell_line=name, sample_id=model["sample_id"], source_cnr=model["cnr_source"],
                      resolved_source_cnr=str(original), source_cnr_sha256=model["cnr_sha256"],
                      target_cnr_sha256=sha(subset), target_cnr_scratch_path=relative(subset),
                      cns_path=relative(destination), cns_sha256=sha(destination),
                      target_cnr_bins=target_count, positive_depth_target_bins=positive_count,
                      contributing_target_bins=probes, segments=len(segments),
                      archived_cns_source=model["cns_source"], archived_cns_sha256=model["cns_sha256"],
                      log_scratch_path=relative(log), elapsed_seconds=elapsed)
        print(json.dumps({k: result[k] for k in ("cell_line", "segments", "contributing_target_bins", "elapsed_seconds")}), flush=True)
        return result

    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        manifest = list(pool.map(process, models))
    assert len(manifest) == 23 and len({r["cell_line"] for r in manifest}) == 23
    write_csv(out / "manifest.csv", manifest)
    print(f"Validated and wrote {len(manifest)} target-only CNS files and manifest in {out}", flush=True)


if __name__ == "__main__":
    main()
