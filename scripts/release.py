#!/usr/bin/env python3
"""
Release workflow for RXNRECer.

Default behavior is local pre-release validation only. Publishing requires an
explicit subcommand.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from version import __version__


PYTHON = sys.executable
RELEASE_TREE = PROJECT_ROOT / "build" / "release_tree"
DIST_DIR = RELEASE_TREE / "dist"
TMP_BUILD_TOOLS = Path("/tmp/rxnrecer_release_build")
TMP_SMOKE_DIR = Path("/tmp/rxnrecer_wheel_smoke")


def run_command(cmd: list[str], cwd: Path | None = None, check: bool = True) -> subprocess.CompletedProcess[str]:
    print("🔧 Running:", " ".join(cmd))
    return subprocess.run(cmd, cwd=cwd, check=check, text=True, capture_output=False)


def check_clean_working_directory() -> bool:
    result = subprocess.run(["git", "status", "--porcelain"], cwd=PROJECT_ROOT, capture_output=True, text=True, check=True)
    if result.stdout.strip():
        print("❌ Working directory is not clean:")
        print(result.stdout)
        return False
    print("✅ Working directory is clean")
    return True


def ensure_build_tools() -> None:
    if TMP_BUILD_TOOLS.exists():
        print(f"✅ Reusing packaging tools in {TMP_BUILD_TOOLS}")
        return
    print(f"📦 Installing packaging tools into {TMP_BUILD_TOOLS}")
    run_command(
        [
            PYTHON,
            "-m",
            "pip",
            "install",
            "--target",
            str(TMP_BUILD_TOOLS),
            "build",
            "twine",
        ],
        cwd=PROJECT_ROOT,
    )


def env_with_build_tools() -> dict[str, str]:
    env = dict(**__import__("os").environ)
    existing = env.get("PYTHONPATH", "")
    env["PYTHONPATH"] = f"{TMP_BUILD_TOOLS}:{existing}" if existing else str(TMP_BUILD_TOOLS)
    return env


def run_python_module(module: str, args: list[str], cwd: Path) -> subprocess.CompletedProcess[str]:
    env = env_with_build_tools()
    cmd = [PYTHON, "-m", module, *args]
    print("🔧 Running:", " ".join(cmd))
    return subprocess.run(cmd, cwd=cwd, check=True, text=True, capture_output=False, env=env)


def sync_versions(version: str) -> None:
    run_command([PYTHON, "scripts/version_manager.py", "set", version], cwd=PROJECT_ROOT)


def prepare_release_tree() -> None:
    run_command([PYTHON, "scripts/prepare_release_tree.py", "--check"], cwd=PROJECT_ROOT)
    run_command([PYTHON, "scripts/prepare_release_tree.py", "--output", str(RELEASE_TREE)], cwd=PROJECT_ROOT)


def run_tests() -> bool:
    result = subprocess.run([PYTHON, "-m", "pytest", "tests/", "-v"], cwd=PROJECT_ROOT, text=True)
    return result.returncode == 0


def clean_release_dist() -> None:
    if DIST_DIR.exists():
        shutil.rmtree(DIST_DIR)


def build_release_artifacts() -> None:
    clean_release_dist()
    run_python_module("build", ["--no-isolation"], cwd=RELEASE_TREE)


def twine_check() -> None:
    run_python_module("twine", ["check", "dist/*"], cwd=RELEASE_TREE)


def smoke_install_wheel() -> None:
    if TMP_SMOKE_DIR.exists():
        shutil.rmtree(TMP_SMOKE_DIR)
    TMP_SMOKE_DIR.mkdir(parents=True, exist_ok=True)
    wheel = next(DIST_DIR.glob("*.whl"))
    run_command([PYTHON, "-m", "pip", "install", "--no-deps", "--target", str(TMP_SMOKE_DIR), str(wheel)], cwd=RELEASE_TREE)
    smoke_env = env_with_build_tools()
    smoke_env["PYTHONPATH"] = f"{TMP_SMOKE_DIR}:{smoke_env['PYTHONPATH']}"
    subprocess.run([PYTHON, "-c", "import rxnrecer; print(rxnrecer.__version__); print(rxnrecer.get_full_version())"], cwd=RELEASE_TREE, check=True, env=smoke_env)
    subprocess.run([PYTHON, "-m", "rxnrecer.cli.predict", "--version"], cwd=RELEASE_TREE, check=True, env=smoke_env)


def list_artifacts() -> None:
    print("📦 Built artifacts:")
    for path in sorted(DIST_DIR.glob("*")):
        print(f"  - {path.name} ({path.stat().st_size} bytes)")


def preflight(version: str, run_repo_tests: bool) -> int:
    print("🚀 RXNRECer Release Preflight")
    print(f"Version: {version}")
    print("=" * 50)
    ensure_build_tools()
    sync_versions(version)
    prepare_release_tree()
    if run_repo_tests:
        ok = run_tests()
        if not ok:
            print("⚠️ Repository tests failed; continuing because this is preflight mode")
    build_release_artifacts()
    twine_check()
    smoke_install_wheel()
    list_artifacts()
    print("✅ Preflight completed successfully")
    print(f"Release tree: {RELEASE_TREE}")
    return 0


def publish(version: str, remote: str, branch: str, push_tag: bool, upload_pypi: bool) -> int:
    print("🚀 RXNRECer Publish")
    print(f"Version: {version}")
    print("=" * 50)
    if not check_clean_working_directory():
        print("❌ Please commit or stash your changes first")
        return 1
    if not DIST_DIR.exists() or not any(DIST_DIR.iterdir()):
        print("❌ No built artifacts found. Run preflight first.")
        return 1
    run_command(["git", "add", "."], cwd=PROJECT_ROOT)
    run_command(["git", "commit", "-m", f"chore: prepare release {version}"], cwd=PROJECT_ROOT)
    run_command(["git", "tag", f"v{version}"], cwd=PROJECT_ROOT)
    run_command(["git", "push", remote, branch], cwd=PROJECT_ROOT)
    if push_tag:
        run_command(["git", "push", remote, "--tags"], cwd=PROJECT_ROOT)
    if upload_pypi:
        run_python_module("twine", ["upload", "dist/*"], cwd=RELEASE_TREE)
    print("✅ Publish steps completed")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description="RXNRECer release workflow")
    subparsers = parser.add_subparsers(dest="command")

    preflight_parser = subparsers.add_parser("preflight", help="Run local release preparation and validation")
    preflight_parser.add_argument("--version", default=__version__)
    preflight_parser.add_argument("--skip-tests", action="store_true", help="Skip repository pytest before packaging")

    publish_parser = subparsers.add_parser("publish", help="Commit/tag/push an already validated release")
    publish_parser.add_argument("--version", default=__version__)
    publish_parser.add_argument("--remote", default="origin")
    publish_parser.add_argument("--branch", default="release")
    publish_parser.add_argument("--push-tag", action="store_true", help="Push the version tag after commit")
    publish_parser.add_argument("--upload-pypi", action="store_true", help="Upload dist/* from release tree to PyPI")

    args = parser.parse_args()
    command = args.command or "preflight"

    if command == "preflight":
        return preflight(args.version, run_repo_tests=not args.skip_tests)
    if command == "publish":
        return publish(args.version, args.remote, args.branch, args.push_tag, args.upload_pypi)
    parser.print_help()
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
