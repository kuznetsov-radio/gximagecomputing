# Releasing Binary Wheels

This project is configured to publish precompiled wheels via GitHub Actions + `cibuildwheel`.

## Goal

Users should be able to run:

```bash
pip install gximagecomputing
```

and receive a wheel containing the compiled `RenderGRFF` extension for their platform (instead of building locally).

## One-Time PyPI Setup

1. Create the package on PyPI (first release only).
2. In PyPI project settings, configure **Trusted Publisher** for this repository.
3. (Recommended) Create the package on TestPyPI as well and configure **Trusted Publisher** there.
4. In GitHub, ensure the workflow can run with `id-token: write` (already set in workflow).

## CI Workflow

Workflow file: `.github/workflows/build_wheels.yml`

- Builds wheels on:
  - Linux (`ubuntu-latest`)
  - Windows (`windows-latest`)
  - macOS Intel (`macos-13`)
  - macOS Apple Silicon (`macos-14`)
- Builds source distribution (`sdist`)
- Uploads to PyPI only when pushing a git tag (`refs/tags/*`)
- Supports manual TestPyPI publish via `workflow_dispatch` input:
  - `publish_target = testpypi`

## Pre-Release Checklist

1. Ensure version is updated in `pyproject.toml` (`[project].version`).
2. Ensure local changes are committed.
3. Verify README/examples commands still work.
4. Verify wheel build config in `pyproject.toml` (`[tool.cibuildwheel]`).
5. Verify required binary/library loading paths still resolve at runtime.

## Local Sanity Commands

Run from repository root.

```bash
python -m pip install -e .
gxrender-mw --help
gxrender-map-view --help
```

Optional local wheel smoke test:

```bash
python -m pip install --upgrade build
python -m build
```

## Release Commands

1. Commit and push your changes:

```bash
git add -A
git commit -m "Release vX.Y.Z"
git push
```

2. Create and push a version tag:

```bash
git tag vX.Y.Z
git push origin vX.Y.Z
```

This triggers wheel+sdist build and PyPI upload in GitHub Actions.

## Optional TestPyPI Dry Run (Recommended)

Before pushing a release tag, run the workflow manually:

1. GitHub -> Actions -> `Build And Publish Wheels` -> `Run workflow`
2. Set `publish_target` to `testpypi`
3. Start run and wait for `upload_testpypi` completion

Then verify install from TestPyPI:

```bash
python -m pip install --upgrade pip
python -m pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple gximagecomputing
python -c "import gximagecomputing; print('ok', gximagecomputing.__file__)"
gxrender-mw --help
```

## Post-Release Validation

In a clean environment:

```bash
python -m pip install --upgrade pip
python -m pip install gximagecomputing
python -c "import gximagecomputing; print('ok', gximagecomputing.__file__)"
gxrender-mw --help
```

If a wheel is not available for the platform/Python version, pip may attempt source build.
