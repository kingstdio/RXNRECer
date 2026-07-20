# Release Workflow

This repository uses a single development source tree and a separate public release track.

## Goals

- Keep `main` as the full internal development workspace.
- Generate a minimal public release tree from a whitelist manifest.
- Avoid publishing internal research assets such as papers, benchmarks, queue scripts, sessions, logs, and temporary outputs.

## Branch and Remote Model

- `main`: internal development branch.
- `release`: public-source branch generated from `main`.
- `origin`: internal GitLab remote.
- `github`: public GitHub remote.

Recommended remote layout:

```bash
git remote add github git@github.com:kingstdio/RXNRECer.git
```

## Release Inputs

The public release tree is defined by [release_manifest.yaml](../release_manifest.yaml).

The manifest is whitelist-first:

- `include`: tracked files or directories allowed into the public tree.
- `exclude`: tracked files or directories always removed even if included above.
- `prune_globs`: cache, editor, and transient patterns removed from the release tree.
- `required_files`: files that must be present in every public release.

## Prepare a Release Tree

Validate the manifest:

```bash
python scripts/prepare_release_tree.py --check
```

Preview the selected files:

```bash
python scripts/prepare_release_tree.py --list-files
```

Generate a clean release tree:

```bash
python scripts/prepare_release_tree.py
```

Default output:

```text
build/release_tree/
```

The generated tree also contains `.release-metadata.json` with the source commit and copied file list.

## Recommended Release Sequence

1. Finish feature and documentation updates on `main`.
2. Run the local preflight:

```bash
python scripts/release.py preflight
```

This will:

- synchronize version references
- regenerate `build/release_tree/`
- build `sdist` and `wheel`
- run `twine check`
- perform a lightweight wheel install/import smoke test

3. If needed, enter the generated tree and inspect the package manually:

```bash
cd build/release_tree
ls dist/
```

If the release environment does not yet have packaging tools, install them first:

```bash
python -m pip install build twine
```

4. After reviewing the generated package, publish explicitly:

```bash
python scripts/release.py publish --remote origin --branch release --push-tag
```

Add `--upload-pypi` only when you really intend to upload the package.

## What Stays Internal

The public release must not contain:

- `paper/`
- `benchmarks/`
- `baselines/`
- internal `results/` artifacts (the public tree only keeps empty placeholders such as `results/README.md` and `.gitkeep`)
- `logs/`
- `sessions/`
- `temp/`
- queue or cluster-specific scripts
- office sync helpers
- internal work notes and reviewer-response materials

If a new directory is added to the internal workflow, update `release_manifest.yaml` before the next release.

## Compatibility Notes

The public tree intentionally stays close to the existing GitHub release layout:

- keep top-level directories such as `examples/`, `results/`, and `extools/msa/`
- keep lightweight compatibility files such as `version.py` and `setup_scm_example.py`
- strip internal benchmark outputs, reviewer materials, and large local artifacts from those paths

This keeps the published structure familiar while still letting `main` evolve as the internal working repository.
