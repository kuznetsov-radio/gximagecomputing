from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[3]
_DEFAULT_EXTERNAL_REPO = _REPO_ROOT.parent / "pyGXrender-test-data"
_DEFAULT_EXTERNAL_RAW = _DEFAULT_EXTERNAL_REPO / "raw"
_ROOT_ENV_VARS = ("GXRENDER_TEST_DATA_ROOT", "GXIMAGECOMPUTING_TEST_DATA_ROOT")


def _normalize_root(path: Path) -> Path:
    path = path.expanduser().resolve()
    if (path / "raw").is_dir() and (path / "scripts").is_dir():
        return path / "raw"
    return path


def candidate_test_data_roots() -> list[Path]:
    roots: list[Path] = []
    seen: set[Path] = set()
    for name in _ROOT_ENV_VARS:
        value = os.environ.get(name, "").strip()
        if not value:
            continue
        path = _normalize_root(Path(value))
        if path not in seen:
            roots.append(path)
            seen.add(path)

    for path in (_DEFAULT_EXTERNAL_RAW,):
        path = path.resolve()
        if path not in seen:
            roots.append(path)
            seen.add(path)
    return roots


def existing_test_data_roots() -> list[Path]:
    return [path for path in candidate_test_data_roots() if path.exists()]


def test_data_setup_hint(subject: str | None = None) -> str:
    prefix = f"Could not locate {subject}. " if subject else "Could not locate the gxrender test dataset. "
    return (
        prefix
        + "Clone ../pyGXrender-test-data next to gximagecomputing and run "
        + "../pyGXrender-test-data/scripts/install_dataset.sh, "
        + "or set GXRENDER_TEST_DATA_ROOT to the extracted raw dataset directory."
    )


def _find_glob(patterns: list[str], *, subject: str) -> Path:
    for root in existing_test_data_roots():
        for pattern in patterns:
            matches = sorted(
                (path for path in root.rglob(pattern) if path.is_file()),
                key=lambda path: (len(path.parts), str(path)),
            )
            if matches:
                return matches[0]
    raise FileNotFoundError(test_data_setup_hint(subject))


def find_model_files(suffix: str | None = None) -> list[Path]:
    suffix_norm = suffix.lower() if suffix else None
    matches: list[Path] = []
    seen: set[Path] = set()
    roots = existing_test_data_roots()
    if not roots:
        return matches
    root = roots[0]
    models_root = root / "models"
    search_root = models_root if models_root.exists() else root
    for path in sorted(search_root.rglob("*"), key=lambda p: (len(p.parts), str(p))):
        if not path.is_file():
            continue
        if suffix_norm and path.suffix.lower() != suffix_norm:
            continue
        if path not in seen:
            matches.append(path)
            seen.add(path)
    return matches


def find_default_model_file(suffix: str | None = None) -> Path:
    matches = find_model_files(suffix)
    if matches:
        return matches[0]
    subject = f"{suffix} model fixture" if suffix else "model fixture"
    raise FileNotFoundError(test_data_setup_hint(subject))


def try_find_default_model_file(suffix: str | None = None) -> Path | None:
    try:
        return find_default_model_file(suffix)
    except FileNotFoundError:
        return None


def try_find_model_file(name: str) -> Path | None:
    try:
        return find_model_file(name)
    except FileNotFoundError:
        return None


def find_model_file(name: str) -> Path:
    return _find_glob([name], subject=f"model file {name!r}")


def try_find_response_file(instrument: str = "aia") -> Path | None:
    try:
        return find_response_file(instrument)
    except FileNotFoundError:
        return None


def find_response_file(instrument: str = "aia") -> Path:
    inst = instrument.strip().lower()
    return _find_glob(
        [f"resp_{inst}_*.sav"],
        subject=f"EUV response file for instrument {inst!r}",
    )


def try_find_ebtel_file(name: str = "ebtel.sav") -> Path | None:
    try:
        return find_ebtel_file(name)
    except FileNotFoundError:
        return None


def find_ebtel_file(name: str = "ebtel.sav") -> Path:
    return _find_glob([name], subject=f"EBTEL table {name!r}")


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Resolve gxrender external test-data fixture paths.")
    sub = parser.add_subparsers(dest="command", required=True)

    sub.add_parser("root", help="Print the first existing test-data root.")

    default_model = sub.add_parser("default-model", help="Print the first installed model file.")
    default_model.add_argument("--suffix", default=None, help="Optional suffix filter, e.g. .h5 or .sav")

    model = sub.add_parser("model", help="Print the path to a model file.")
    model.add_argument("name", help="Model filename, e.g. test.chr.h5")

    response = sub.add_parser("response", help="Print the path to an EUV response file.")
    response.add_argument("instrument", help="Instrument key, e.g. aia or stereo-a")

    ebtel = sub.add_parser("ebtel", help="Print the path to an EBTEL table.")
    ebtel.add_argument("name", nargs="?", default="ebtel.sav", help="EBTEL filename")

    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    try:
        if args.command == "root":
            roots = existing_test_data_roots()
            if not roots:
                raise FileNotFoundError(test_data_setup_hint())
            print(roots[0])
        elif args.command == "default-model":
            print(find_default_model_file(args.suffix))
        elif args.command == "model":
            print(find_model_file(args.name))
        elif args.command == "response":
            print(find_response_file(args.instrument))
        elif args.command == "ebtel":
            print(find_ebtel_file(args.name))
        else:
            raise ValueError(f"Unsupported command: {args.command}")
    except FileNotFoundError as exc:
        print(str(exc), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
