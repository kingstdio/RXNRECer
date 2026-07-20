#!/usr/bin/env python3
"""
Version management script for RXNRECer.

This keeps the public-release version references aligned across packaging,
CLI output, and lightweight release metadata files.
"""

from __future__ import annotations

import argparse
import re
import sys
from datetime import datetime
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parent.parent
SEMVER_RE = re.compile(r"^\d+\.\d+\.\d+$")

TEXT_REPLACEMENTS = {
    "version.py": [
        (r'BUILD_DATE = "[^"]+"', 'BUILD_DATE = "{build_date}"'),
    ],
    "rxnrecer/__init__.py": [
        (r'__version__ = "[^"]+"', '__version__ = "{version}"'),
        (r'__version_info__ = \(\d+, \d+, \d+\)', '__version_info__ = ({major}, {minor}, {patch})'),
        (r'return f"\{__version__\} \([^)]+\)"', 'return f"{{__version__}} ({build_date})"'),
    ],
    ".bumpversion.cfg": [
        (r'^current_version = .+$', 'current_version = {version}'),
    ],
    "rxnrecer/cli/predict.py": [
        (r"version='RXNRECer [^']+'", "version='RXNRECer {version}'"),
        (r'print\("RXNRECer v[^"]+ - Enzyme Reaction Prediction"\)', 'print("RXNRECer v{version} - Enzyme Reaction Prediction")'),
    ],
    "README.md": [
        (r'\*\*RXNRECer v\d+\.\d+\.\d+\*\*', '**RXNRECer v{version}**'),
        (r'pip install rxnrecer==\d+\.\d+\.\d+', 'pip install rxnrecer=={version}'),
        (r'https://pypi\.org/project/rxnrecer/\d+\.\d+\.\d+/', 'https://pypi.org/project/rxnrecer/{version}/'),
    ],
    "docs/INSTALL.md": [
        (r'\*\*Version \d+\.\d+\.\d+\*\*', '**Version {version}**'),
        (r'pip install rxnrecer==\d+\.\d+\.\d+', 'pip install rxnrecer=={version}'),
        (r'https://pypi\.org/project/rxnrecer/\d+\.\d+\.\d+/', 'https://pypi.org/project/rxnrecer/{version}/'),
    ],
    "ckpt/README.md": [
        (r'RXNRECer v\d+\.\d+\.\d+', 'RXNRECer v{version}'),
    ],
    "data/README.md": [
        (r'RXNRECer v\d+\.\d+\.\d+', 'RXNRECer v{version}'),
    ],
    "setup_scm_example.py": [
        (r'RXNRECer v\d+\.\d+\.\d+', 'RXNRECer v{version}'),
    ],
}


def validate_version(version: str) -> str:
    if not SEMVER_RE.match(version):
        raise ValueError(f"Invalid version format: {version}. Expected X.Y.Z")
    return version


def version_parts(version: str) -> tuple[str, str, str]:
    return tuple(version.split("."))  # type: ignore[return-value]


def get_current_version() -> str:
    content = (PROJECT_ROOT / "rxnrecer/__init__.py").read_text(encoding="utf-8")
    match = re.search(r'__version__ = "([^"]+)"', content)
    if not match:
        raise ValueError("Could not determine current version from rxnrecer/__init__.py")
    return match.group(1)


def replace_text(path: Path, replacements: list[tuple[str, str]], values: dict[str, str]) -> bool:
    if not path.exists():
        return False

    content = path.read_text(encoding="utf-8")
    updated = content
    changed = False

    for pattern, replacement_template in replacements:
        replacement = replacement_template.format(**values)
        new_text = re.sub(pattern, replacement, updated, flags=re.MULTILINE)
        if new_text != updated:
            changed = True
            updated = new_text

    if changed:
        path.write_text(updated, encoding="utf-8")
    return changed


def update_all_versions(new_version: str) -> None:
    validate_version(new_version)
    major, minor, patch = version_parts(new_version)
    values = {
        "version": new_version,
        "major": major,
        "minor": minor,
        "patch": patch,
        "build_date": datetime.now().strftime("%Y-%m-%d"),
    }

    print(f"Updating version references to {new_version}...")
    for relpath, replacements in TEXT_REPLACEMENTS.items():
        path = PROJECT_ROOT / relpath
        changed = replace_text(path, replacements, values)
        status = "updated" if changed else "checked"
        print(f"  - {status}: {relpath}")


def bump_version(part: str) -> str:
    major, minor, patch = map(int, get_current_version().split("."))
    if part == "major":
        major, minor, patch = major + 1, 0, 0
    elif part == "minor":
        minor, patch = minor + 1, 0
    elif part == "patch":
        patch += 1
    else:
        raise ValueError(f"Invalid version part: {part}")
    return f"{major}.{minor}.{patch}"


def main() -> int:
    parser = argparse.ArgumentParser(description="RXNRECer Version Manager")
    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    subparsers.add_parser("current", help="Show current version")

    set_parser = subparsers.add_parser("set", help="Set specific version")
    set_parser.add_argument("version", help="Version to set, for example 1.4.0")

    bump_parser = subparsers.add_parser("bump", help="Bump version")
    bump_parser.add_argument("part", choices=["major", "minor", "patch"])

    args = parser.parse_args()

    try:
        if args.command == "current":
            print(get_current_version())
            return 0
        if args.command == "set":
            update_all_versions(args.version)
            return 0
        if args.command == "bump":
            new_version = bump_version(args.part)
            update_all_versions(new_version)
            print(new_version)
            return 0
    except ValueError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    parser.print_help()
    return 1


if __name__ == "__main__":
    sys.exit(main())
