#!/usr/bin/env python3

from __future__ import annotations

import ctypes
from dataclasses import dataclass

import numpy as np
import scipy.io as sio

from .radio import GXRadioImageComputing


@dataclass
class EUVResponseMeta:
    instrument: str
    channels: list[str]
    source: str = ""
    mode: str = "default"


def build_default_euv_response(
    instrument: str = "aia",
    channels: list[str] | None = None,
    nt: int = 161,
    logte_min: float = 4.8,
    logte_max: float = 7.8,
) -> tuple[np.ndarray, np.dtype, EUVResponseMeta]:
    """Build a synthetic, self-contained EUV response table.

    This is a practical Python fallback for exercising the ComputeEUV pipeline
    when an instrument response table is not provided from IDL tooling.
    """
    if channels is None:
        channels = ["94", "131", "171", "193", "211", "304", "335"]

    centers = {
        "94": 6.8,
        "131": 7.0,
        "171": 5.9,
        "193": 6.2,
        "211": 6.3,
        "304": 4.9,
        "335": 6.4,
    }
    widths = {
        "94": 0.18,
        "131": 0.25,
        "171": 0.14,
        "193": 0.16,
        "211": 0.16,
        "304": 0.10,
        "335": 0.18,
    }

    logte = np.linspace(float(logte_min), float(logte_max), int(nt), dtype=np.float64)
    nchan = len(channels)
    all_resp = np.zeros((nchan, nt), dtype=np.float64)
    for i, ch in enumerate(channels):
        mu = centers.get(ch, 6.1)
        sig = widths.get(ch, 0.18)
        all_resp[i, :] = np.exp(-0.5 * ((logte - mu) / sig) ** 2)

    # IDL uses ds=0.36 for AIA-like response setup.
    ds = 0.36

    response_dt = np.dtype(
        [
            ("ds", np.float64),
            ("NT", np.int32),
            ("Nchannels", np.int32),
            # Channel-major layout: contiguous block of NT values per channel.
            ("logte", np.float64, (nt,)),
            ("all", np.float64, (nchan, nt)),
        ]
    )
    response = np.zeros(1, dtype=response_dt)
    response["ds"] = ds
    response["NT"] = nt
    response["Nchannels"] = nchan
    response["logte"] = logte
    response["all"] = all_resp
    return response, response_dt, EUVResponseMeta(instrument=instrument.upper(), channels=[str(c) for c in channels])


def load_euv_response_sav(path: str):
    """Load real IDL GX response from SAV (e.g. aia_response.sav)."""
    d = sio.readsav(path, python_dict=True, verbose=False)
    gx = None
    gx_var_name = None
    gx_field_map = None

    # Prefer the canonical IDL variable name when present, but accept any SAV
    # variable whose restored struct exposes the expected response fields.
    candidate_names = ["gxresponse", *[k for k in d.keys() if k != "gxresponse"]]
    required_fields = {"logte", "all", "channels"}
    instrument_field_aliases = ("name", "instrument")
    for name in candidate_names:
        try:
            arr = np.asarray(d[name])
        except Exception:
            continue
        if arr.size == 0 or arr.dtype.names is None:
            continue
        field_map = {str(f).lower(): f for f in arr.dtype.names}
        if not required_fields.issubset(field_map):
            continue
        if not any(alias in field_map for alias in instrument_field_aliases):
            continue
        gx = np.asarray(d[name])[0]
        gx_var_name = name
        gx_field_map = field_map
        break

    if gx is None:
        keys = ", ".join(sorted(str(k) for k in d.keys()))
        raise ValueError(
            f"Unsupported EUV response SAV (no response-like struct with LOGTE/ALL/CHANNELS and NAME|INSTRUMENT): {path}; keys=[{keys}]"
        )

    assert gx_field_map is not None
    logte = np.asarray(gx[gx_field_map["logte"]], dtype=np.float64).reshape(-1)
    all_resp = np.asarray(gx[gx_field_map["all"]], dtype=np.float64)
    if all_resp.ndim != 2:
        raise ValueError(f"gxresponse.ALL must be 2-D, got shape={all_resp.shape}")
    # Ensure channel-major [Nchannels, NT]
    if all_resp.shape[1] != logte.size and all_resp.shape[0] == logte.size:
        all_resp = all_resp.T
    if all_resp.shape[1] != logte.size:
        raise ValueError(f"gxresponse.ALL incompatible with LOGTE: ALL={all_resp.shape}, LOGTE={logte.size}")

    channels_raw = np.asarray(gx[gx_field_map["channels"]]).reshape(-1)
    channels = [x.decode("utf-8", "ignore") if isinstance(x, (bytes, np.bytes_)) else str(x) for x in channels_raw]
    instrument_field = "name" if "name" in gx_field_map else "instrument"
    name_raw = gx[gx_field_map[instrument_field]]
    instrument = name_raw.decode("utf-8", "ignore") if isinstance(name_raw, (bytes, np.bytes_)) else str(name_raw)

    nt = int(logte.size)
    nchan = int(all_resp.shape[0])
    response_dt = np.dtype(
        [
            ("ds", np.float64),
            ("NT", np.int32),
            ("Nchannels", np.int32),
            ("logte", np.float64, (nt,)),
            ("all", np.float64, (nchan, nt)),
        ]
    )
    response = np.zeros(1, dtype=response_dt)
    response["ds"] = 0.36
    response["NT"] = nt
    response["Nchannels"] = nchan
    response["logte"] = logte
    response["all"] = all_resp
    mode = "sav_gxresponse" if str(gx_var_name).lower() == "gxresponse" else f"sav_struct:{gx_var_name}"
    meta = EUVResponseMeta(instrument=instrument.upper(), channels=channels, source=path, mode=mode)
    return response, response_dt, meta


class GXEUVImageComputing(GXRadioImageComputing):
    """Python binding for EUV rendering entrypoints in RenderGRFF library."""

    def __init__(self, libname=None):
        super().__init__(libname=libname)
        self.euvfunc = self._lib.pyComputeEUV
        self.euvfunc.restype = ctypes.c_int
        self.euvfuncsh = self._lib.pyComputeEUV_SH
        self.euvfuncsh.restype = ctypes.c_int

    @staticmethod
    def _make_euv_simbox(box_nx, box_ny, box_xc, box_yc, box_dx, box_dy, parallel=False, exact=False, nthreads=0):
        projection = 0
        if parallel:
            projection |= 1
        if exact:
            projection |= 2
        if int(nthreads) > 0:
            projection |= (int(nthreads) << 16)
        dt_s = np.dtype(
            [
                ("Nx", np.int32),
                ("Ny", np.int32),
                ("xc", np.float64),
                ("yc", np.float64),
                ("dx", np.float64),
                ("dy", np.float64),
                ("projection", np.int32),
            ]
        )
        simbox = np.zeros(1, dtype=dt_s)
        simbox["Nx"] = int(box_nx)
        simbox["Ny"] = int(box_ny)
        simbox["xc"] = float(box_xc)
        simbox["yc"] = float(box_yc)
        simbox["dx"] = float(box_dx)
        simbox["dy"] = float(box_dy)
        simbox["projection"] = int(projection)
        return simbox, dt_s

    @staticmethod
    def _reserve_euv_output(box_nx: int, box_ny: int, nchan: int):
        dt_o = np.dtype(
            [
                ("flagsAll", np.int32, (6,)),
                ("flagsCorona", np.int32, (6,)),
                ("fluxCorona", np.float64, (box_nx, box_ny, nchan)),
                ("fluxTR", np.float64, (box_nx, box_ny, nchan)),
            ]
        )
        out = np.zeros(1, dtype=dt_o)
        return out, dt_o

    def synth_euv(
        self,
        model,
        model_dt,
        ebtel,
        ebtel_dt,
        response,
        response_dt,
        box_nx,
        box_ny,
        box_xc,
        box_yc,
        box_dx,
        box_dy,
        tbase,
        nbase,
        q0,
        a,
        b,
        mode=0,
        parallel=False,
        exact=False,
        nthreads=0,
        shtable=None,
    ):
        simbox, dt_s = self._make_euv_simbox(
            box_nx=box_nx,
            box_ny=box_ny,
            box_xc=box_xc,
            box_yc=box_yc,
            box_dx=box_dx,
            box_dy=box_dy,
            parallel=parallel,
            exact=exact,
            nthreads=nthreads,
        )
        dt_c = np.dtype(
            [
                ("Tbase", np.float64),
                ("nbase", np.float64),
                ("Q0", np.float64),
                ("a", np.float64),
                ("b", np.float64),
                ("mode", np.int32),
            ]
        )
        cparms = np.zeros(1, dtype=dt_c)
        cparms["Tbase"] = float(tbase)
        cparms["nbase"] = float(nbase)
        cparms["Q0"] = float(q0)
        cparms["a"] = float(a)
        cparms["b"] = float(b)
        cparms["mode"] = int(mode)

        nchan = int(response["Nchannels"][0])
        out, dt_o = self._reserve_euv_output(int(box_nx), int(box_ny), nchan)

        if shtable is not None:
            shtable = np.asarray(shtable, dtype=np.float64)
            if shtable.shape != (7, 7):
                raise ValueError("SHtable must be shape (7,7)")
            status = self.euvfuncsh(
                model.ctypes.data_as(ctypes.c_void_p),
                ebtel.ctypes.data_as(ctypes.c_void_p),
                response.ctypes.data_as(ctypes.c_void_p),
                simbox.ctypes.data_as(ctypes.c_void_p),
                cparms.ctypes.data_as(ctypes.c_void_p),
                out.ctypes.data_as(ctypes.c_void_p),
                shtable.ctypes.data_as(ctypes.c_void_p),
            )
        else:
            status = self.euvfunc(
                model.ctypes.data_as(ctypes.c_void_p),
                ebtel.ctypes.data_as(ctypes.c_void_p),
                response.ctypes.data_as(ctypes.c_void_p),
                simbox.ctypes.data_as(ctypes.c_void_p),
                cparms.ctypes.data_as(ctypes.c_void_p),
                out.ctypes.data_as(ctypes.c_void_p),
            )
        if status != 0:
            raise RuntimeError(f"ComputeEUV failed with status={status}")

        # Match the IDL memory layout convention used by the DLL (column-major).
        # This mirrors the MW path in `radio.py`.
        flux_cor = np.reshape(np.ravel(out["fluxCorona"][0]), (int(box_nx), int(box_ny), nchan), order="F").swapaxes(0, 1)
        flux_tr = np.reshape(np.ravel(out["fluxTR"][0]), (int(box_nx), int(box_ny), nchan), order="F").swapaxes(0, 1)
        return {
            "simbox": simbox,
            "simbox_dtype": dt_s,
            "outspace": out,
            "outspace_dtype": dt_o,
            "flux_corona": flux_cor,
            "flux_tr": flux_tr,
        }
