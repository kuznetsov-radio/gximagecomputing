from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from gxrender.io.model import load_model_with_metadata
from gxrender.policy.response_providers import find_supporting_euv_response_provider


@dataclass(slots=True)
class ScriptPolicyResult:
    model_observer: str
    has_saved_fov: bool
    effective_observer: str
    normalized_observer: str
    normalized_instrument: str
    implicit_response_mode: str


def _normalize_text(value: Any) -> str:
    return str(value or "").strip()


def normalize_observer_name(value: Any) -> str:
    raw = _normalize_text(value).lower()
    aliases = {
        "earth": "earth",
        "terra": "earth",
        "sdo": "earth",
        "solar dynamics observatory": "earth",
        "stereo-a": "stereo-a",
        "stereo a": "stereo-a",
        "stereoa": "stereo-a",
        "stereo-b": "stereo-b",
        "stereo b": "stereo-b",
        "stereob": "stereo-b",
        "solar orbiter": "solar orbiter",
        "solar-orbiter": "solar orbiter",
        "solo": "solar orbiter",
    }
    return aliases.get(raw, raw)


def _extract_observer_name(model_metadata: dict[str, Any], observer_metadata: Any) -> str:
    observer_name = _normalize_text(model_metadata.get("observer_name") or model_metadata.get("observer")).lower()
    if observer_name:
        return observer_name
    if isinstance(observer_metadata, dict):
        return _normalize_text(observer_metadata.get("name") or observer_metadata.get("label")).lower()
    return ""


def _extract_has_saved_fov(model_metadata: dict[str, Any], observer_metadata: Any) -> bool:
    if isinstance(observer_metadata, dict):
        fov = observer_metadata.get("fov")
        if isinstance(fov, dict):
            if any(fov.get(key) is not None for key in ("xc_arcsec", "yc_arcsec", "xsize_arcsec", "ysize_arcsec")):
                return True
    return False


def _load_model_metadata(model_path: Path) -> tuple[dict[str, Any], Any]:
    _, _, model_metadata, observer_metadata = load_model_with_metadata(str(model_path))
    return dict(model_metadata), observer_metadata


def _default_instrument_for_observer(normalized_observer: str) -> str:
    mapping = {
        "earth": "aia",
        "stereo-a": "stereo-a",
        "stereo-b": "stereo-b",
        "solar orbiter": "solo-fsi",
    }
    return mapping.get(normalized_observer, "")


def _extract_obs_time_iso(model_metadata: dict[str, Any]) -> str:
    for key in ("obs_time", "date_obs", "date-obs", "obstime"):
        value = _normalize_text(model_metadata.get(key))
        if value:
            return value
    return ""


def evaluate_script_policy(
    *,
    model_path: Path,
    instrument: str | None,
    observer: str | None,
    response_sav: str | None,
) -> ScriptPolicyResult:
    model_metadata, observer_metadata = _load_model_metadata(model_path)
    model_observer = _extract_observer_name(model_metadata, observer_metadata)
    has_saved_fov = _extract_has_saved_fov(model_metadata, observer_metadata)

    effective_observer = _normalize_text(observer) or model_observer
    normalized_observer = normalize_observer_name(effective_observer)
    normalized_instrument = _normalize_text(instrument).lower()
    response_sav_value = _normalize_text(response_sav)
    obs_time_iso = _extract_obs_time_iso(model_metadata)
    implicit_instrument = normalized_instrument or _default_instrument_for_observer(normalized_observer)
    provider_available = bool(implicit_instrument) and find_supporting_euv_response_provider(implicit_instrument, obs_time_iso) is not None

    if not response_sav_value:
        if normalized_instrument and normalized_instrument != "aia" and not provider_available:
            raise ValueError(
                f"INSTRUMENT={instrument} requires an explicit RESPONSE_SAV because time-dependent non-AIA responses are not implemented."
            )
        if not normalized_instrument and normalized_observer and normalized_observer != "earth" and not provider_available:
            raise ValueError(
                f"Effective observer '{effective_observer}' is non-Earth and no RESPONSE_SAV was provided. "
                "Provide RESPONSE_SAV explicitly for non-AIA/non-Earth rendering paths."
            )

    implicit_mode = "none"
    if provider_available:
        implicit_mode = "provider_time_dependent"
    elif not response_sav_value and not normalized_instrument and (not normalized_observer or normalized_observer == "earth"):
        implicit_mode = "earth_aia_time_dependent"

    return ScriptPolicyResult(
        model_observer=model_observer,
        has_saved_fov=has_saved_fov,
        effective_observer=effective_observer,
        normalized_observer=normalized_observer,
        normalized_instrument=normalized_instrument,
        implicit_response_mode=implicit_mode,
    )


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Evaluate launcher observer/FOV/response policy for a model input.")
    parser.add_argument("--model-path", type=Path, required=True)
    parser.add_argument("--instrument", type=str, default="")
    parser.add_argument("--observer", type=str, default="")
    parser.add_argument("--response-sav", type=str, default="")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    try:
        result = evaluate_script_policy(
            model_path=args.model_path,
            instrument=args.instrument,
            observer=args.observer,
            response_sav=args.response_sav,
        )
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    print(f"MODEL_OBSERVER={result.model_observer}")
    print(f"HAS_SAVED_FOV={1 if result.has_saved_fov else 0}")
    print(f"EFFECTIVE_OBSERVER={result.effective_observer}")
    print(f"NORMALIZED_OBSERVER={result.normalized_observer}")
    print(f"NORMALIZED_INSTRUMENT={result.normalized_instrument}")
    print(f"IMPLICIT_RESPONSE_MODE={result.implicit_response_mode}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
