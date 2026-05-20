from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Literal, Protocol


ObserverName = Literal["earth", "stereo-a", "stereo-b", "solar orbiter", "unknown"]


@dataclass(slots=True)
class ObserverFovRequest:
    model_path: Path
    loader: str
    cli_args: Any
    prefer_execute_center: bool = True


@dataclass(slots=True)
class ObserverFovResolution:
    model: Any
    model_dt: Any
    model_metadata: dict[str, Any]
    observer_geometry: Any
    center_source: str
    xc_auto: float
    yc_auto: float
    model_w_arcsec: float
    model_h_arcsec: float
    applied_overrides: dict[str, float]
    decisions: list[str] = field(default_factory=list)


@dataclass(slots=True)
class EUVResponseRequest:
    args: Any
    obs_time_iso: str


@dataclass(slots=True)
class EUVResponseResolution:
    response: Any
    response_dt: Any
    response_meta: Any
    mode: str
    source: str
    decisions: list[str] = field(default_factory=list)


class EUVResponseProvider(Protocol):
    name: str
    priority: int

    def supports(self, instrument: str, obs_time_iso: str) -> bool:
        ...

    def build(self, instrument: str, obs_time_iso: str, channels: list[str] | None) -> tuple[Any, Any, Any]:
        ...
