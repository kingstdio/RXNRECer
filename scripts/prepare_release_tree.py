#!/usr/bin/env python3
"""
Prepare a minimal public release tree from the development repository.

The manifest format is intentionally simple so it can be parsed without
external dependencies:

    key: scalar
    list_key:
      - item1
      - item2
"""

from __future__ import annotations

import argparse
import fnmatch
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, List


PROJECT_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_MANIFEST = PROJECT_ROOT / "release_manifest.yaml"
DEFAULT_OUTPUT = PROJECT_ROOT / "build" / "release_tree"


def normalize_relpath(value: str) -> str:
    return value.strip().strip("/").replace("\\", "/")


def load_manifest(path: Path) -> Dict[str, object]:
    data: Dict[str, object] = {}
    current_list_key: str | None = None

    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.split("#", 1)[0].rstrip()
        if not line.strip():
            continue

        if line.startswith("  - "):
            if current_list_key is None:
                raise ValueError(f"List item without active key: {raw_line}")
            value = normalize_relpath(line[4:])
            items = data.setdefault(current_list_key, [])
            assert isinstance(items, list)
            items.append(value)
            continue

        current_list_key = None
        if ":" not in line:
            raise ValueError(f"Unsupported manifest line: {raw_line}")

        key, value = line.split(":", 1)
        key = key.strip()
        value = value.strip()

        if value:
            data[key] = value
        else:
            data[key] = []
            current_list_key = key

    return data


def git_ls_files(repo_root: Path) -> List[str]:
    result = subprocess.run(
        ["git", "ls-files", "--cached", "--others", "--exclude-standard", "-z"],
        cwd=repo_root,
        check=True,
        capture_output=True,
    )
    files = result.stdout.decode("utf-8").split("\0")
    return sorted({normalize_relpath(item) for item in files if item})


def matches_prefix(path: str, candidate: str) -> bool:
    return path == candidate or path.startswith(candidate + "/")


def is_pruned(path: str, prune_globs: List[str]) -> bool:
    for pattern in prune_globs:
        if fnmatch.fnmatch(path, pattern):
            return True
    return False


def select_release_files(
    tracked_files: List[str],
    include: List[str],
    exclude: List[str],
    prune_globs: List[str],
) -> List[str]:
    selected: List[str] = []
    for relpath in tracked_files:
        if not any(matches_prefix(relpath, item) for item in include):
            continue
        if any(matches_prefix(relpath, item) for item in exclude):
            continue
        if is_pruned(relpath, prune_globs):
            continue
        selected.append(relpath)
    return sorted(selected)


def validate_required_files(required: List[str], selected: List[str]) -> List[str]:
    selected_set = set(selected)
    missing = [item for item in required if item not in selected_set]
    return missing


def copy_release_files(repo_root: Path, output_dir: Path, selected: List[str]) -> None:
    if output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    for relpath in selected:
        src = repo_root / relpath
        dst = output_dir / relpath
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)


def write_metadata(
    repo_root: Path,
    output_dir: Path,
    manifest_path: Path,
    selected: List[str],
) -> None:
    commit = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=repo_root,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    metadata = {
        "source_commit": commit,
        "manifest": str(manifest_path.relative_to(repo_root)),
        "file_count": len(selected),
        "files": selected,
    }
    metadata_path = output_dir / ".release-metadata.json"
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Prepare a minimal release tree")
    parser.add_argument(
        "--manifest",
        default=str(DEFAULT_MANIFEST),
        help="Path to release manifest",
    )
    parser.add_argument(
        "--output",
        default=str(DEFAULT_OUTPUT),
        help="Output directory for the generated release tree",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="Validate manifest selection without copying files",
    )
    parser.add_argument(
        "--list-files",
        action="store_true",
        help="Print selected files to stdout",
    )
    args = parser.parse_args()

    repo_root = PROJECT_ROOT
    manifest_path = Path(args.manifest).resolve()
    output_dir = Path(args.output).resolve()

    manifest = load_manifest(manifest_path)
    include = list(manifest.get("include", []))
    exclude = list(manifest.get("exclude", []))
    prune_globs = list(manifest.get("prune_globs", []))
    required_files = list(manifest.get("required_files", []))

    tracked_files = git_ls_files(repo_root)
    selected = select_release_files(tracked_files, include, exclude, prune_globs)
    missing_required = validate_required_files(required_files, selected)

    if missing_required:
        print("Missing required release files:", file=sys.stderr)
        for item in missing_required:
            print(f"  - {item}", file=sys.stderr)
        return 1

    if args.list_files:
        for item in selected:
            print(item)

    if args.check:
        print(f"Manifest OK: {len(selected)} tracked files selected")
        return 0

    copy_release_files(repo_root, output_dir, selected)
    write_metadata(repo_root, output_dir, manifest_path, selected)
    print(f"Release tree created at: {output_dir}")
    print(f"Tracked files copied: {len(selected)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
