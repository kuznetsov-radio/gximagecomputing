from __future__ import annotations

from gxrender.policy.contracts import EUVResponseRequest, EUVResponseResolution
from gxrender.policy.response_providers import find_supporting_euv_response_provider


def apply_default_response_selection(args, *, observer_geometry):
    """Phase-1 wrapper around existing EUV instrument-default policy."""
    from gxrender.workflows import render_euv

    return render_euv._apply_default_response_selection(args, observer_geometry=observer_geometry)


def resolve_euv_response(request: EUVResponseRequest) -> EUVResponseResolution:
    """Resolve EUV response inputs, preferring any compatible provider registry entry."""
    from gxrender.workflows import render_euv

    args = request.args
    explicit_prebuilt = getattr(args, "response", None) is not None
    explicit_sav = getattr(args, "response_sav", None) is not None

    if not explicit_prebuilt and not explicit_sav:
        instrument = str(getattr(args, "instrument", "") or "").strip().lower()
        channels = getattr(args, "channels", None)
        if channels is not None:
            channels = [str(channel) for channel in channels]

        if instrument:
            provider = find_supporting_euv_response_provider(instrument, request.obs_time_iso)
            if provider is not None:
                try:
                    response, response_dt, response_meta = provider.build(instrument, request.obs_time_iso, channels)
                    source = str(getattr(response_meta, "source", getattr(provider, "name", "provider")) or "")
                    mode = str(getattr(response_meta, "mode", "python_provider_time_dependent") or "")
                    return EUVResponseResolution(
                        response=response,
                        response_dt=response_dt,
                        response_meta=response_meta,
                        mode=mode,
                        source=source,
                        decisions=[f"provider:{getattr(provider, 'name', 'unknown')}"] ,
                    )
                except Exception:
                    pass

    response, response_dt, response_meta = render_euv._resolve_response_inputs(
        args,
        obs_time_iso=request.obs_time_iso,
    )

    source = ""
    mode = ""
    if hasattr(response_meta, "source"):
        source = str(getattr(response_meta, "source", "") or "")
    if hasattr(response_meta, "mode"):
        mode = str(getattr(response_meta, "mode", "") or "")

    return EUVResponseResolution(
        response=response,
        response_dt=response_dt,
        response_meta=response_meta,
        mode=mode,
        source=source,
        decisions=[],
    )
