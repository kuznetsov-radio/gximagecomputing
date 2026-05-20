from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from gxrender.policy.contracts import EUVResponseProvider


_REGISTERED_PROVIDERS: list[EUVResponseProvider] = []


@dataclass(slots=True)
class _PyEUVToolsPayloadProvider:
    name: str
    priority: int
    instrument: str
    builder: Any
    builder_kwargs: dict[str, Any]

    def supports(self, instrument: str, obs_time_iso: str) -> bool:
        return str(instrument or "").strip().lower() == self.instrument

    def build(self, instrument: str, obs_time_iso: str, channels: list[str] | None):
        from gxrender.euv import EUVResponseMeta

        kwargs = dict(self.builder_kwargs)
        kwargs["obstime"] = obs_time_iso or None
        if channels is not None:
            kwargs["channels"] = [str(channel) for channel in channels]
        response, response_dt, payload_meta = self.builder(**kwargs)
        payload_meta = payload_meta if isinstance(payload_meta, dict) else {}
        meta_channels = payload_meta.get("channels", tuple(channels or ()))
        response_meta = EUVResponseMeta(
            instrument=str(payload_meta.get("instrument", str(instrument).upper())).upper(),
            channels=[str(channel) for channel in meta_channels],
            source=str(payload_meta.get("source", self.name)),
            mode="python_provider_time_dependent",
        )
        return response, response_dt, response_meta


def register_euv_response_provider(provider: EUVResponseProvider) -> None:
    """Register an optional EUV time-dependent response provider.

    Phase-1 scaffolding only: registration does not affect current workflows yet.
    """
    _REGISTERED_PROVIDERS.append(provider)


def clear_registered_euv_response_providers() -> None:
    """Test helper to clear runtime-registered providers."""
    _REGISTERED_PROVIDERS.clear()


def _discover_external_providers() -> list[EUVResponseProvider]:
    """Best-effort external provider discovery.

    Future-ready hook: if external packages expose a compatible provider list,
    we ingest it automatically without hard dependency coupling.
    """
    discovered: list[EUVResponseProvider] = []
    try:
        import pyeuvtools  # type: ignore
        from pyeuvtools import response as pyeuv_response  # type: ignore

        candidate = getattr(pyeuvtools, "get_time_dependent_response_providers", None)
        if candidate is None:
            candidate = getattr(pyeuv_response, "get_time_dependent_response_providers", None)
        if callable(candidate):
            providers = candidate()
            if isinstance(providers, list):
                for provider in providers:
                    if provider is not None:
                        discovered.append(provider)

        builder_specs = [
            ("pyeuvtools-euvi-a", "stereo-a", getattr(pyeuv_response, "build_euvi_temperature_response_gx_payload", None), {"spacecraft": "ahead"}),
            ("pyeuvtools-euvi-b", "stereo-b", getattr(pyeuv_response, "build_euvi_temperature_response_gx_payload", None), {"spacecraft": "behind"}),
            ("pyeuvtools-eui-fsi", "solo-fsi", getattr(pyeuv_response, "build_eui_temperature_response_gx_payload", None), {"detector": "fsi"}),
            ("pyeuvtools-eui-hri", "solo-hri", getattr(pyeuv_response, "build_eui_temperature_response_gx_payload", None), {"detector": "hri"}),
            ("pyeuvtools-trace", "trace", getattr(pyeuv_response, "build_trace_temperature_response_gx_payload", None), {}),
            ("pyeuvtools-sxt", "sxt", getattr(pyeuv_response, "build_sxt_temperature_response_gx_payload", None), {}),
        ]
        for index, (name, instrument, builder, builder_kwargs) in enumerate(builder_specs, start=100):
            if callable(builder):
                discovered.append(
                    _PyEUVToolsPayloadProvider(
                        name=name,
                        priority=index,
                        instrument=instrument,
                        builder=builder,
                        builder_kwargs=builder_kwargs,
                    )
                )
    except Exception:
        return []
    return discovered


def get_registered_euv_response_providers() -> list[EUVResponseProvider]:
    """Return provider list in priority order.

    The currently active behavior does not consume this list yet. It is added in
    phase-1 to avoid redesign when non-AIA time-dependent providers land.
    """
    combined: list[EUVResponseProvider] = []
    combined.extend(_REGISTERED_PROVIDERS)
    combined.extend(_discover_external_providers())

    # Keep stable deterministic ordering for future resolver usage.
    return sorted(combined, key=lambda p: int(getattr(p, "priority", 1000)))


def find_supporting_euv_response_provider(instrument: str, obs_time_iso: str) -> EUVResponseProvider | None:
    """Return the first registered provider that supports the requested instrument."""
    instrument_name = str(instrument or "").strip().lower()
    for provider in get_registered_euv_response_providers():
        try:
            if provider.supports(instrument_name, obs_time_iso):
                return provider
        except Exception:
            continue
    return None


__all__ = [
    "register_euv_response_provider",
    "clear_registered_euv_response_providers",
    "get_registered_euv_response_providers",
    "find_supporting_euv_response_provider",
]
