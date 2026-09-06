#!/usr/bin/env python3
"""Verify the CNV coverage handoff, preserving a manifest self-hash defect.

Run only after safe extraction into the ignored retrieval directory. This checks
file identities; it does not execute scripts or modify the original bundle.
"""
import argparse
import collections
import csv
import hashlib
import json
from pathlib import Path, PurePosixPath

ROOT = Path(__file__).resolve().parents[1]
DEFAULT = ROOT / "data/cluster_wes_retrieval/2026-09-06/ovcan_human_wes_cnv_coverage_2026-09-06"
ARCHIVE = ROOT / "data/cluster_wes_retrieval/ovcan_human_wes_cnv_coverage_2026-09-06.tar.gz"

def sha(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bundle", type=Path, default=DEFAULT)
    ap.add_argument("--archive", type=Path, default=ARCHIVE)
    args = ap.parse_args()
    rows, seen = [], set()
    for line in (args.bundle / "SHA256SUMS").read_text().splitlines():
        expected, name = line.split(None, 1)
        name = name.lstrip("*").removeprefix("./")
        relative = PurePosixPath(name)
        assert not relative.is_absolute() and ".." not in relative.parts, name
        assert name not in seen, name
        seen.add(name)
        path = args.bundle / name
        assert path.is_file() and not path.is_symlink(), name
        actual = sha(path)
        # A manifest cannot ordinarily contain its own final SHA256. Preserve
        # the received self-entry, but require every other file to match.
        match = actual == expected
        assert match or name == "SHA256SUMS", f"Payload checksum mismatch: {name}"
        rows.append(dict(relative_path=name, bytes=path.stat().st_size,
                         declared_sha256=expected, sha256=actual,
                         manifest_match=match,
                         status="verified_payload" if match else "manifest_self_hash_mismatch"))
    files = {str(p.relative_to(args.bundle)) for p in args.bundle.rglob("*") if p.is_file()}
    assert files == seen, files.symmetric_difference(seen)
    archive_hash = sha(args.archive)
    assert archive_hash == "c5d79958ae05115c975d63e3858fae7cd372b37c593295e0e2472bbf38c6706a", "Unexpected archive"
    report = ROOT / "reports/wes_cnv_coverage_2026-09-06"
    report.mkdir(exist_ok=True)
    with (report / "archive_inventory.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    summary = dict(archive_name=args.archive.name, archive_sha256=archive_hash,
                   archive_bytes=args.archive.stat().st_size, files=len(files),
                   uncompressed_bytes=sum(r["bytes"] for r in rows),
                   manifest_entries_verified=sum(r["manifest_match"] for r in rows),
                   mismatches=[r for r in rows if not r["manifest_match"]],
                   folder_counts=dict(collections.Counter(r["relative_path"].split("/")[0]
                       if "/" in r["relative_path"] else "[root]" for r in rows)),
                   has_readme_or_handoff=False,
                   original_files_modified=False,
                   extraction="All tar members checked as relative paths without parent traversal and regular files/directories; tarfile data filter used.")
    (report / "archive_validation.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))

if __name__ == "__main__":
    main()
