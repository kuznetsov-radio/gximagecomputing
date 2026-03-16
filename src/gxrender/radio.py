#!/usr/bin/env python3

import ctypes
import platform
import warnings
import numpy as np
from pathlib import Path
from .io.ebtel import load_ebtel as io_load_ebtel, load_ebtel_none as io_load_ebtel_none
from .io.model import (
    load_model_dict as io_load_model_dict,
    load_model_hdf as io_load_model_hdf,
    load_model_hdf_with_observer as io_load_model_hdf_with_observer,
    load_model_sav as io_load_model_sav,
    load_model_sav_with_observer as io_load_model_sav_with_observer,
)


class GXRadioImageComputing:
    def __init__(self, libname=None):
        libc_mw = None
        last_error = None
        candidates = self._resolve_library_candidates(libname)
        for candidate in candidates:
            try:
                libc_mw = ctypes.CDLL(str(candidate))
                self.libname = str(candidate)
                break
            except OSError as exc:
                last_error = exc
        if libc_mw is None:
            raise OSError(
                "Could not load any RenderGRFF library candidate. "
                f"Tried: {[str(c) for c in candidates]}"
            ) from last_error
        self._lib = libc_mw

        self.mwfunc=libc_mw.pyComputeMW
        self.mwfunc.restype=ctypes.c_int

        self.mwfuncsh=libc_mw.pyComputeMW_SH
        self.mwfuncsh.restype=ctypes.c_int

    @classmethod
    def _resolve_library_candidates(cls, libname):
        if libname is not None:
            return [Path(libname)]

        module_dir = Path(__file__).resolve().parent
        repo_root = module_dir.parent.parent
        binaries_dir = repo_root / "binaries"

        machine = platform.machine().lower()
        if machine in ("x86_64", "amd64"):
            arch = "x86_64"
        elif machine in ("arm64", "aarch64"):
            arch = "arm64"
        else:
            arch = machine

        system = platform.system()
        candidates = []

        if system == "Windows":
            candidates.extend([
                module_dir / "RenderGRFF.pyd",
                module_dir / "RenderGRFF_64.dll",
                module_dir / "RenderGRFF_32.dll",
                binaries_dir / "RenderGRFF_64.dll",
                binaries_dir / "RenderGRFF_32.dll",
            ])
        elif system == "Darwin":
            candidates.extend([
                module_dir / f"RenderGRFF_{arch}.so",
                binaries_dir / f"RenderGRFF_{arch}.so",
                module_dir / "RenderGRFF.so",
                binaries_dir / "RenderGRFF.so",
            ])
        else:
            candidates.extend([
                module_dir / "RenderGRFF.so",
                binaries_dir / "RenderGRFF.so",
            ])

        candidates.extend(sorted(module_dir.glob("RenderGRFF*.so")))
        candidates.extend(sorted(module_dir.glob("RenderGRFF*.dylib")))
        candidates.extend(sorted(module_dir.glob("RenderGRFF*.dll")))
        candidates.extend(sorted(module_dir.glob("RenderGRFF*.pyd")))
        return [candidate for candidate in candidates if candidate.exists()]

    def load_ebtel(self, ebtel_file):
        return io_load_ebtel(ebtel_file)

    @staticmethod
    def load_ebtel_none():
        return io_load_ebtel_none()
    
    def load_model_dict(self, model_dict, header):
        return io_load_model_dict(model_dict, header)

    def load_model_hdf(
        self,
        file_name,
        DSun=None,
        lonC=None,
        b0Sun=None,
        recompute_observer_ephemeris: bool = False,
        observer_name: str | None = None,
        ):
        return io_load_model_hdf(
            file_name,
            DSun=DSun,
            lonC=lonC,
            b0Sun=b0Sun,
            recompute_observer_ephemeris=recompute_observer_ephemeris,
            observer_name=observer_name,
        )

    def load_model_hdf_with_observer(
        self,
        file_name,
        DSun=None,
        lonC=None,
        b0Sun=None,
        recompute_observer_ephemeris: bool = False,
        observer_name: str | None = None,
    ):
        return io_load_model_hdf_with_observer(
            file_name,
            DSun=DSun,
            lonC=lonC,
            b0Sun=b0Sun,
            recompute_observer_ephemeris=recompute_observer_ephemeris,
            observer_name=observer_name,
        )

    def load_model_sav(
        self,
        file_name,
        DSun=None,
        lonC=None,
        b0Sun=None,
        recompute_observer_ephemeris: bool = False,
        observer_name: str | None = None,
        ):
        return io_load_model_sav(
            file_name,
            DSun=DSun,
            lonC=lonC,
            b0Sun=b0Sun,
            recompute_observer_ephemeris=recompute_observer_ephemeris,
            observer_name=observer_name,
        )

    def load_model_sav_with_observer(
        self,
        file_name,
        DSun=None,
        lonC=None,
        b0Sun=None,
        recompute_observer_ephemeris: bool = False,
        observer_name: str | None = None,
    ):
        return io_load_model_sav_with_observer(
            file_name,
            DSun=DSun,
            lonC=lonC,
            b0Sun=b0Sun,
            recompute_observer_ephemeris=recompute_observer_ephemeris,
            observer_name=observer_name,
        )

    def synth_model(
        self,
        model,
        model_dt,
        ebtel,
        ebtel_dt,
        freqlist,
        box_Nx,
        box_Ny,
        box_xc,
        box_yc,
        box_dx,
        box_dy,
        Tbase,
        nbase,
        Q0,
        a,
        b,
        SHtable=None,
        rot=0,
        projection=0,
        mode=0,
        warn_defaults=True,
    ):
        """
            mode: bit-mask matching IDL DefineCoronaParams.pro
                  0 - default
                  1 - force_isothermal
                  2 - interpolB
                  4 - analyticalNT

            projection:
                  1 - parallel
                  2 - exact
        """
        if warn_defaults:
            if mode == 0:
                warnings.warn(
                    "GXRadioImageComputing.synth_model() is using mode=0. "
                    "Pass mode explicitly for full control of the coronal-heating path.",
                    stacklevel=2,
                )
            if int(projection) == 0:
                warnings.warn(
                    "GXRadioImageComputing.synth_model() is using projection=0. "
                    "Pass projection flags explicitly for full control of the LOS geometry.",
                    stacklevel=2,
                )
            if SHtable is None:
                warnings.warn(
                    "GXRadioImageComputing.synth_model() received no SHtable, so the non-SH DLL entrypoint will be used.",
                    stacklevel=2,
                )

        dt_s=np.dtype([('Nx', np.int32),
                       ('Ny', np.int32),
                       ('Nf', np.int32),
                       ('projection', np.int32),
                       ('xc', np.float64),
                       ('yc', np.float64),
                       ('dx', np.float64),
                       ('dy', np.float64),
                       ('rot', np.float64),
                       ('freqlist', np.float64, len(freqlist))])
        simbox=np.zeros(1, dtype=dt_s)
        simbox['Nx']=box_Nx
        simbox['Ny']=box_Ny
        simbox['Nf']=len(freqlist)
        simbox['projection']=projection
        simbox['xc']=box_xc
        simbox['yc']=box_yc
        simbox['dx']=box_dx
        simbox['dy']=box_dy
        simbox['rot']=rot
        simbox['freqlist']=freqlist

        dt_c=np.dtype([('Tbase', np.float64),
                       ('nbase', np.float64),
                       ('Q0', np.float64),
                       ('a', np.float64),
                       ('b', np.float64),
                       ('mode', np.int32)])
        cparms=np.zeros(1, dtype=dt_c)
        cparms['Tbase']=Tbase
        cparms['nbase']=nbase
        cparms['Q0']=Q0
        cparms['a']=a
        cparms['b']=b
        cparms['mode']=mode

        dt_o=np.dtype([('flagsAll', np.int32, 6),
                       ('flagsCorona', np.int32, 6),
                       ('TI', np.float64, (box_Nx, box_Ny, len(freqlist))),
                       ('TV', np.float64, (box_Nx, box_Ny, len(freqlist)))])
        outspace=np.zeros(1, dtype=dt_o)

        _dt_model = np.ctypeslib.ndpointer(dtype=model_dt)
        _dt_ebtel = np.ctypeslib.ndpointer(dtype=ebtel_dt)
        _dt_s  = np.ctypeslib.ndpointer(dtype=dt_s)
        _dt_c  = np.ctypeslib.ndpointer(dtype=dt_c)
        _dt_o  = np.ctypeslib.ndpointer(dtype=dt_o)

        if SHtable is None:
            #SHtable = np.ones((7, 7), dtype=np.float64)
            func = self.mwfunc
            func.argtypes=[_dt_model, _dt_ebtel, _dt_s, _dt_c, _dt_o]
            r=func(model, ebtel, simbox, cparms, outspace)
        else:
            func = self.mwfuncsh
            dt_sh = np.dtype([("SHtable", np.float64, (7, 7))])
            shtable = np.zeros(1, dtype=dt_sh)
            shtable['SHtable'] = SHtable
            _dt_sh = np.ctypeslib.ndpointer(dtype=dt_sh)
            func.argtypes=[_dt_model, _dt_ebtel, _dt_s, _dt_c, _dt_o, _dt_sh]
            r=func(model, ebtel, simbox, cparms, outspace, shtable)

        o=outspace[0]
        TI=np.reshape(np.ravel(o['TI']), (box_Nx, box_Ny, len(freqlist)), order='F').swapaxes(0, 1)
        TV=np.reshape(np.ravel(o['TV']), (box_Nx, box_Ny, len(freqlist)), order='F').swapaxes(0, 1)
    
        return {"TI": TI, "TV": TV}
    
