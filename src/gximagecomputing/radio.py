#!/usr/bin/env python3

import ctypes
import platform
import re
import shlex
import numpy as np
import scipy.io as io
import h5py
from astropy.coordinates import SkyCoord
from sunpy.coordinates import frames, sun
import astropy.units as u
from astropy.time import Time
from pathlib import Path
from dataclasses import dataclass
from typing import Any, Dict
from .io.voxel_id import gx_box2id


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

        self.mwfunc=libc_mw.pyComputeMW
        self.mwfunc.restype=ctypes.c_int

        self.mwfuncsh=libc_mw.pyComputeMW_SH
        self.mwfuncsh.restype=ctypes.c_int

    @staticmethod
    def _decode_if_bytes(value):
        if isinstance(value, (bytes, np.bytes_)):
            return value.decode("utf-8")
        return value

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

    @staticmethod
    def _coerce_time(obs_time):
        if isinstance(obs_time, Time):
            return obs_time
        if isinstance(obs_time, (bytes, np.bytes_)):
            obs_time = obs_time.decode("utf-8")
        return Time(obs_time)

    @staticmethod
    def _sanitize_status_mask(start_idx, end_idx, chromo_mask):
        max_index = int(np.prod(chromo_mask.shape))
        start_idx = np.asarray(start_idx, dtype=np.int64).copy()
        end_idx = np.asarray(end_idx, dtype=np.int64).copy()
        # IDL-compatible behavior for linear indexing into CHROMO_MASK.
        start_idx = np.clip(start_idx, 0, max_index - 1)
        end_idx = np.clip(end_idx, 0, max_index - 1)
        return start_idx, end_idx

    @staticmethod
    def _load_h5_group(group: h5py.Group) -> Dict[str, np.ndarray]:
        loaded = {}
        for key in group.keys():
            try:
                loaded[key] = group[key][:]
            except ValueError:
                loaded[key] = group[key][()]
        return loaded

    @staticmethod
    def _normalize_cube_xyz(arr: np.ndarray, nx_hint: int | None = None, ny_hint: int | None = None) -> np.ndarray:
        """Normalize 3D cube to internal (x, y, z) layout."""
        a = np.asarray(arr)
        if a.ndim != 3:
            return a

        if nx_hint is not None and ny_hint is not None:
            if a.shape[0] == nx_hint and a.shape[1] == ny_hint:
                return a
            if a.shape[0] == ny_hint and a.shape[1] == nx_hint:
                return a.transpose((1, 0, 2))
            if a.shape[1] == ny_hint and a.shape[2] == nx_hint:
                return a.transpose((2, 1, 0))
            if a.shape[1] == nx_hint and a.shape[2] == ny_hint:
                return a.transpose((1, 2, 0))

        # Legacy fallback (pre-v0.2 files were typically y/x swapped).
        return a.transpose((1, 0, 2))

    @classmethod
    def _normalize_vector_cube_xyzc(
        cls, arr: np.ndarray, nx_hint: int | None = None, ny_hint: int | None = None
    ) -> np.ndarray:
        """Normalize vector cube to internal (x, y, z, c) layout."""
        a = np.asarray(arr)
        if a.ndim != 4:
            return a

        if a.shape[-1] == 3:
            comps = [cls._normalize_cube_xyz(a[..., i], nx_hint, ny_hint) for i in range(3)]
            return np.stack(comps, axis=-1)
        if a.shape[0] == 3:
            comps = [cls._normalize_cube_xyz(a[i, ...], nx_hint, ny_hint) for i in range(3)]
            return np.stack(comps, axis=-1)
        return a

    @classmethod
    def _components_to_vector_cube_xyzc(
        cls,
        bx: np.ndarray,
        by: np.ndarray,
        bz: np.ndarray,
        nx_hint: int | None = None,
        ny_hint: int | None = None,
    ) -> np.ndarray:
        bx_i = cls._normalize_cube_xyz(bx, nx_hint, ny_hint)
        by_i = cls._normalize_cube_xyz(by, nx_hint, ny_hint)
        bz_i = cls._normalize_cube_xyz(bz, nx_hint, ny_hint)
        if bx_i.shape != by_i.shape or bx_i.shape != bz_i.shape:
            raise ValueError(
                f"Incompatible vector component shapes: bx={bx_i.shape}, by={by_i.shape}, bz={bz_i.shape}"
            )
        return np.stack((bx_i, by_i, bz_i), axis=-1)

    @staticmethod
    def _scalar_from_any(value: Any) -> Any:
        arr = np.asarray(value)
        if arr.ndim == 0:
            return arr.item()
        if arr.size == 1:
            return arr.reshape(-1)[0].item()
        return value

    @staticmethod
    def _parse_execute_time(execute_text: str) -> str | None:
        try:
            tokens = shlex.split(execute_text)
        except ValueError:
            tokens = execute_text.split()
        for i, tok in enumerate(tokens):
            if tok == "--time" and i + 1 < len(tokens):
                return tokens[i + 1]
            if tok.startswith("--time="):
                return tok.split("=", 1)[1]
        return None

    @staticmethod
    def _decode_dataset_scalar(ds) -> Any:
        value = ds[()]
        if isinstance(value, (bytes, np.bytes_)):
            return value.decode("utf-8", "ignore")
        return value

    @classmethod
    def _fill_header_from_base_index(cls, model_f: h5py.File, header: Dict[str, Any]) -> None:
        if "base" not in model_f or "index" not in model_f["base"]:
            return
        raw = cls._decode_dataset_scalar(model_f["base"]["index"])
        text = raw if isinstance(raw, str) else str(raw)

        def _extract_fits_card(key: str) -> str | None:
            # Match FITS-header style cards, e.g.:
            #   CRVAL1  = -17.0599
            #   DATE-OBS= '2025-11-26T15:34:31.000'
            m = re.search(rf"(?m)^\s*{re.escape(key)}\s*=\s*([^\n/]+)", text)
            if not m:
                return None
            value = m.group(1).strip()
            if value.startswith("'"):
                m_str = re.match(r"'([^']*)'", value)
                if m_str:
                    return m_str.group(1).strip()
            return value.strip().strip("'")

        if "lon" not in header:
            mlon = re.search(r",\s*([+-]?\d+(?:\.\d+)?)\s*,\s*b'CRLN-CEA'", text)
            if mlon:
                header["lon"] = float(mlon.group(1))
            else:
                crval1 = _extract_fits_card("CRVAL1")
                if crval1 is not None:
                    try:
                        header["lon"] = float(crval1)
                    except ValueError:
                        pass
        if "lat" not in header:
            mlat = re.search(r",\s*([+-]?\d+(?:\.\d+)?)\s*,\s*b'CRLT-CEA'", text)
            if mlat:
                header["lat"] = float(mlat.group(1))
            else:
                crval2 = _extract_fits_card("CRVAL2")
                if crval2 is not None:
                    try:
                        header["lat"] = float(crval2)
                    except ValueError:
                        pass
        if "obs_time" not in header:
            date_obs = _extract_fits_card("DATE-OBS") or _extract_fits_card("DATE_OBS")
            if date_obs:
                header["obs_time"] = date_obs
            else:
                mtime = re.search(r"\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}(?:\.\d+)?", text)
                if mtime:
                    header["obs_time"] = mtime.group(0)
        if "dsun_obs" not in header:
            dsun_obs = _extract_fits_card("DSUN_OBS")
            if dsun_obs is not None:
                try:
                    header["dsun_obs"] = float(dsun_obs)
                except ValueError:
                    pass
            else:
                mds = re.search(
                    r"\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}(?:\.\d+)?'\s*,\s*([0-9]+(?:\.[0-9]+)?(?:[eE][+-]?[0-9]+)?)",
                    text,
                )
                if mds:
                    header["dsun_obs"] = float(mds.group(1))
        if "lonC" not in header:
            hglon = _extract_fits_card("HGLN_OBS")
            if hglon is not None:
                try:
                    header["lonC"] = float(hglon)
                except ValueError:
                    pass

    def load_ebtel(self, ebtel_file):
        ebtel_data = io.readsav(ebtel_file)
        s = ebtel_data["lrun"].shape
    
        ebtel_dtypes = [("DEM_on", np.int32),
                        ("DDM_on", np.int32),
                        ("NQ", np.int32),
                        ("NL", np.int32),
                        ("NT", np.int32),
                        ("Qrun", np.float32, ebtel_data["qrun"].shape),
                        ("Lrun", np.float32, ebtel_data["lrun"].shape),
                        ("logtdem", np.float32, ebtel_data["logtdem"].size)
                      ]
        if "dem_cor_run" in ebtel_data.keys():
            ebtel_dtypes.append(("DEM_cor_run", np.float32, ebtel_data["dem_cor_run"].shape))
        
        if "ddm_cor_run" in ebtel_data.keys():
            ebtel_dtypes.append(("DDM_cor_run", np.float32, ebtel_data["ddm_cor_run"].shape))
        
        ebtel_dt = np.dtype(ebtel_dtypes)
        ebtel_c = np.zeros(1, dtype=ebtel_dt)
        
        ebtel_c["DEM_on"] = 1 if "dem_cor_run" in ebtel_data.keys() else 0
        ebtel_c["DDM_on"] = 1 if "ddm_cor_run" in ebtel_data.keys() else 0
        # Python/C ABI in this project expects NQ/NL in this order.
        # (Swapping to match IDL helper text causes large rendering mismatch.)
        ebtel_c["NQ"] = s[1]
        ebtel_c["NL"] = s[0]
        ebtel_c["NT"] = ebtel_data["logtdem"].size
        ebtel_c["Qrun"] = np.float32(ebtel_data["qrun"])
        ebtel_c["Lrun"] = np.float32(ebtel_data["lrun"])
        ebtel_c["logtdem"] = np.float32(ebtel_data["logtdem"])
        
        if "dem_cor_run" in ebtel_data.keys():
            ebtel_c["DEM_cor_run"] = np.float32(ebtel_data["DEM_cor_run"])
        
        if "ddm_cor_run" in ebtel_data.keys():
            ebtel_c["DDM_cor_run"] = np.float32(ebtel_data["DDM_cor_run"])
        
        return ebtel_c, ebtel_dt

    @staticmethod
    def load_ebtel_none():
        """Return a minimal EBTEL struct with DEM/DDM disabled.

        This matches the legacy IDL behavior for ebtelfile='':
        heating tables are not used and the native library falls back to
        isothermal/hydrostatic coronal plasma handling.
        """
        ebtel_dt = np.dtype(
            [
                ("DEM_on", np.int32),
                ("DDM_on", np.int32),
                ("NQ", np.int32),
                ("NL", np.int32),
                ("NT", np.int32),
                ("Qrun", np.float32, (1, 1)),
                ("Lrun", np.float32, (1, 1)),
                ("logtdem", np.float32, (1,)),
            ]
        )
        ebtel_c = np.zeros(1, dtype=ebtel_dt)
        ebtel_c["DEM_on"] = 0
        ebtel_c["DDM_on"] = 0
        ebtel_c["NQ"] = 1
        ebtel_c["NL"] = 1
        ebtel_c["NT"] = 1
        return ebtel_c, ebtel_dt

    @dataclass
    class ChromoModelData:
        header: Dict[str, Any]
        model: Dict[str, np.ndarray]

    def _build_hdf_chromo_data(self, file_name):
        with h5py.File(file_name, "r") as model_f:
            if "chromo" not in model_f:
                raise KeyError("Missing required '/chromo' group in HDF5 model.")
            chromo_box = model_f["chromo"]
            model_dict = self._load_h5_group(chromo_box)
            header = {k: self._decode_if_bytes(v) for k, v in dict(chromo_box.attrs).items()}

            ny_hint = None
            nx_hint = None
            if "base" in model_f and "bx" in model_f["base"]:
                base_shape = np.asarray(model_f["base"]["bx"][:]).shape
                if len(base_shape) == 2:
                    ny_hint, nx_hint = int(base_shape[0]), int(base_shape[1])

            # pyAMPP v0.2.0 stores GEN line metadata in '/lines' rather than '/chromo'.
            if "lines" in model_f:
                lines_box = self._load_h5_group(model_f["lines"])
                for key in ("av_field", "phys_length", "voxel_status", "start_idx", "end_idx"):
                    if key not in model_dict and key in lines_box:
                        model_dict[key] = lines_box[key]

            # v0.2.0 keeps coronal geometry in '/corona'.
            if "corona" in model_f:
                corona_box = self._load_h5_group(model_f["corona"])
                if "dr" not in model_dict and "dr" in corona_box:
                    model_dict["dr"] = corona_box["dr"]
                if "corona_base" not in model_dict and "corona_base" in corona_box:
                    model_dict["corona_base"] = corona_box["corona_base"]
                if "bcube" not in model_dict:
                    if all(k in corona_box for k in ("bx", "by", "bz")):
                        model_dict["bcube"] = self._components_to_vector_cube_xyzc(
                            corona_box["bx"], corona_box["by"], corona_box["bz"], nx_hint, ny_hint
                        )
                    elif "bcube" in corona_box:
                        model_dict["bcube"] = corona_box["bcube"]

            # Fallbacks for legacy files where metadata lived under '/chromo'.
            if "dr" not in model_dict and "dr" in chromo_box:
                model_dict["dr"] = chromo_box["dr"][:]
            if "corona_base" not in model_dict and "corona_base" in chromo_box:
                model_dict["corona_base"] = chromo_box["corona_base"][()]
            if "chromo_bcube" not in model_dict:
                if all(k in model_dict for k in ("bx", "by", "bz")):
                    model_dict["chromo_bcube"] = self._components_to_vector_cube_xyzc(
                        model_dict["bx"], model_dict["by"], model_dict["bz"], nx_hint, ny_hint
                    )
                elif "chromo_bcube" in chromo_box:
                    model_dict["chromo_bcube"] = chromo_box["chromo_bcube"][:]
            if "bcube" not in model_dict and "bcube" in chromo_box:
                model_dict["bcube"] = chromo_box["bcube"][:]

            if "chromo_mask" not in model_dict and "base" in model_f and "chromo_mask" in model_f["base"]:
                model_dict["chromo_mask"] = model_f["base"]["chromo_mask"][:]
            if "base" in model_f and "wcs_header" in model_f["base"]:
                wcs_raw = self._decode_if_bytes(model_f["base"]["wcs_header"][()])
                if isinstance(wcs_raw, str):
                    m_lonc = re.search(r"CRVAL1\s*=\s*([+-]?\d+(?:\.\d+)?)", wcs_raw)
                    if m_lonc:
                        header["lonC"] = float(m_lonc.group(1))
            # Fallback metadata if chromo attrs are incomplete.
            if "corona" in model_f:
                for key, value in dict(model_f["corona"].attrs).items():
                    if key not in header:
                        header[key] = self._decode_if_bytes(value)
            if "metadata" in model_f and "execute" in model_f["metadata"] and "obs_time" not in header:
                execute_text = self._decode_dataset_scalar(model_f["metadata"]["execute"])
                if isinstance(execute_text, str):
                    t_fallback = self._parse_execute_time(execute_text)
                    if t_fallback:
                        header["obs_time"] = t_fallback
            self._fill_header_from_base_index(model_f, header)

        # Normalize H5 array orientation to the same internal convention used
        # by SAV loading/IDL (`x,y,z[,c]`).
        if "dz" in model_dict and np.asarray(model_dict["dz"]).ndim == 3:
            model_dict["dz"] = self._normalize_cube_xyz(np.asarray(model_dict["dz"]), nx_hint, ny_hint)
        if "bcube" in model_dict and np.asarray(model_dict["bcube"]).ndim == 4:
            model_dict["bcube"] = self._normalize_vector_cube_xyzc(
                np.asarray(model_dict["bcube"]), nx_hint, ny_hint
            )
        if "chromo_bcube" in model_dict and np.asarray(model_dict["chromo_bcube"]).ndim == 4:
            model_dict["chromo_bcube"] = self._normalize_vector_cube_xyzc(
                np.asarray(model_dict["chromo_bcube"]), nx_hint, ny_hint
            )
        if "chromo_layers" in model_dict:
            model_dict["chromo_layers"] = int(self._scalar_from_any(model_dict["chromo_layers"]))
        if "corona_base" in model_dict:
            model_dict["corona_base"] = int(self._scalar_from_any(model_dict["corona_base"]))

        required = (
            "dr",
            "dz",
            "bcube",
            "chromo_bcube",
            "chromo_layers",
            "corona_base",
            "chromo_idx",
            "chromo_n",
            "n_p",
            "n_hi",
            "chromo_t",
            "chromo_mask",
            "av_field",
            "phys_length",
            "voxel_status",
            "start_idx",
            "end_idx",
        )
        missing = [key for key in required if key not in model_dict]
        if missing:
            raise KeyError(
                "Missing required CHR datasets after compatibility resolution: "
                + ", ".join(missing)
            )

        required_header = ("lon", "lat", "dsun_obs", "obs_time")
        missing_header = [k for k in required_header if k not in header]
        if missing_header:
            raise KeyError(
                "Missing required observation metadata fields: " + ", ".join(missing_header)
            )
        header["obs_time"] = self._coerce_time(header["obs_time"])
        return self.ChromoModelData(header=header, model=model_dict)

    def _build_sav_chromo_data(self, file_name):
        model_data = io.readsav(file_name)
        box = model_data.box
        lon = box.index[0].CRVAL1[0]
        lat = box.index[0].CRVAL2[0]
        aptime = Time(box.index[0]["DATE_OBS"][0])
        header = {
            "lon": lon,
            "lat": lat,
            "dsun_obs": box.index[0]["DSUN_OBS"][0],
            "obs_time": aptime,
        }
        # IDL LoadGXmodel derives model lonC via WCS conversion. For CEA SAV
        # files this is provided directly in INDEX/HGLN_OBS and gives better
        # parity with IDL than recomputing from Carrington->Stonyhurst.
        if "HGLN_OBS" in box.index[0].dtype.names:
            header["lonC"] = float(box.index[0]["HGLN_OBS"][0])

        model_dict = {
            "dr": box.dr[0],
            # Match IDL LoadGXmodel.pro conventions:
            # scipy.readsav exposes arrays with reversed axis order vs IDL.
            # These permutations make downstream .T operations produce
            # render-ready arrays equivalent to IDL.
            "dz": box.dz[0].transpose((2, 1, 0)),
            "bcube": box.bcube[0].transpose((3, 2, 1, 0)),
            "chromo_bcube": box.chromo_bcube[0].transpose((3, 2, 1, 0)),
            "chromo_layers": box.chromo_layers[0],
            "corona_base": box.corona_base[0],
            "chromo_idx": box.chromo_idx[0].astype(np.int64),
            "chromo_n": box.chromo_n[0],
            "n_p": box.n_p[0],
            "n_hi": box.n_hi[0],
            "chromo_t": box.chromo_t[0],
            "chromo_mask": box["base"][0]["chromo_mask"][0],
        }

        if "AVFIELD" in box.dtype.names:
            model_dict["av_field"] = box.avfield[0].T.reshape(-1, order="F")
            model_dict["phys_length"] = box.physlength[0].T.reshape(-1, order="F")
            model_dict["voxel_status"] = box.status[0].T.reshape(-1, order="F")
            model_dict["start_idx"] = box.startidx[0].T.reshape(-1, order="F")
            model_dict["end_idx"] = box.endidx[0].T.reshape(-1, order="F")
        else:
            sc = box.bcube[0].shape
            Nx, Ny = sc[3], sc[2]
            QB = np.zeros((Nx, Ny, sc[1]), dtype=np.float64)
            QL = np.zeros((Nx, Ny, sc[1]), dtype=np.float64)
            uu = np.zeros((Nx, Ny, sc[1]), dtype=np.uint8)
            idx = np.unravel_index(box.idx[0], QB.shape, order="F")
            QB[idx] = box.bmed[0]
            QL[idx] = box.length[0]
            uu[idx] = 4
            model_dict["av_field"] = QB.reshape(-1, order="F")
            model_dict["phys_length"] = QL.reshape(-1, order="F")
            model_dict["voxel_status"] = uu.reshape(-1, order="F")
            model_dict["start_idx"] = np.zeros(QB.size, dtype=np.int64)
            model_dict["end_idx"] = np.zeros(QB.size, dtype=np.int64)

        return self.ChromoModelData(header=header, model=model_dict)
    
    def load_model_dict(self, model_dict, header):
        lon, lat, dsun_obs, obs_time = [header.get(k) for k in ("lon", "lat", "dsun_obs", "obs_time")]
        obs_time = self._coerce_time(obs_time)
        lonC_override = header.get("lonC", None)
        b0sun_override = header.get("b0Sun", None)
        dsun_override = header.get("DSun", None)
        init_coords = SkyCoord(lon*u.deg, lat*u.deg, frame=frames.HeliographicCarrington,rsun=696000*u.km,\
                       obstime=obs_time, observer="earth")
        hgs = init_coords.transform_to(frames.HeliographicStonyhurst)
        obstime = obs_time.unix - 283996800 # according to IDL specs, anytim function

        DSun = dsun_obs*1e2 if dsun_override is None else float(dsun_override)
        RSun = init_coords.rsun.to(u.cm).value
        b0Sun = sun.B0(obs_time).value if b0sun_override is None else float(b0sun_override)

        lonC = hgs.lon.to(u.deg).value if lonC_override is None else float(lonC_override)
        latC=lat

        dr = model_dict["dr"]

        dx=dr[0]*RSun
        dy=dr[1]*RSun
        dz_uniform=dr[2]*RSun
        dz=model_dict["dz"].T*RSun

        s = dz.shape
        sc = model_dict["bcube"].T.shape
        
        Nx=s[2]
        Ny=s[1]
        Nz=s[0]
    
        chromo_layers=model_dict["chromo_layers"]
        corona_layers=Nz-chromo_layers
        corona_base=model_dict["corona_base"]

        # Some pyAMPP products can have a 1-layer discrepancy between
        # chromo_layers and corona_base metadata. Align corona_base so the
        # stitched volume dimensions remain consistent with dz/Nz.
        expected_corona_base = int(sc[1] - (Nz - chromo_layers))
        if (sc[1] - int(corona_base)) != (Nz - chromo_layers):
            if 0 <= expected_corona_base < sc[1]:
                corona_base = expected_corona_base
            else:
                corona_base = int(np.clip(corona_base, 0, max(sc[1] - 1, 0)))
    
        Bx=np.zeros(s, dtype=np.float32, order="C")
        By=np.zeros(s, dtype=np.float32, order="C")
        Bz=np.zeros(s, dtype=np.float32, order="C")
    
        chromo_bcube = model_dict["chromo_bcube"].T
        bcube = model_dict["bcube"].T
    
        Bx[0:chromo_layers, :, :]=chromo_bcube[0, :, :, :]
        By[0:chromo_layers, :, :]=chromo_bcube[1, :, :, :]
        Bz[0:chromo_layers, :, :]=chromo_bcube[2, :, :, :]
    
        cor_target = Nz - chromo_layers
        cor_source = max(sc[1] - corona_base, 0)
        cor_copy = min(cor_target, cor_source)
        if cor_copy > 0:
            Bx[chromo_layers : chromo_layers + cor_copy, :, :] = bcube[0, corona_base : corona_base + cor_copy, :, :]
            By[chromo_layers : chromo_layers + cor_copy, :, :] = bcube[1, corona_base : corona_base + cor_copy, :, :]
            Bz[chromo_layers : chromo_layers + cor_copy, :, :] = bcube[2, corona_base : corona_base + cor_copy, :, :]
    
        chr_s  = (chromo_layers,  Ny, Nx)
        cor_s  = (sc[1]-corona_base, Ny, Nx)
        cor_s2 = (corona_base,    Ny, Nx)
    
        chromo_n0 =np.zeros(chr_s)
        chromo_np =np.zeros(chr_s)
        chromo_nHI=np.zeros(chr_s)
        chromo_T0 =np.zeros(chr_s)

        chromo_idx = model_dict["chromo_idx"]
        chromo_n0.flat[chromo_idx] = model_dict["chromo_n"]
        chromo_np.flat[chromo_idx] = model_dict["n_p"]
        chromo_nHI.flat[chromo_idx]= model_dict["n_hi"]
        chromo_T0.flat[chromo_idx] = model_dict["chromo_t"]

        chromo_mask = model_dict["chromo_mask"].copy()

        voxel_box = {
            "bcube": model_dict["bcube"],
            "dr": model_dict["dr"],
            "chromo_layers": chromo_layers,
            "corona_base": corona_base,
            "chromo_idx": model_dict.get("chromo_idx"),
            "chromo_t": model_dict.get("chromo_t"),
            "chromo_n": model_dict.get("chromo_n"),
        }
        if "start_idx" in model_dict:
            voxel_box["start_idx"] = model_dict["start_idx"]
        voxel_id_xyz = gx_box2id(voxel_box)
        if voxel_id_xyz is None:
            VoxelID = np.zeros(s, dtype=np.uint8)
        else:
            # gx_box2id returns (nx, ny, nz); renderer struct stores (nz, ny, nx).
            VoxelID = np.asarray(voxel_id_xyz, dtype=np.uint8).transpose((2, 1, 0))
        startidx, endidx = self._sanitize_status_mask(model_dict["start_idx"], model_dict["end_idx"], chromo_mask)

        # Match IDL linear indexing semantics used in LoadGXmodel.pro:
        # ID1=byte(chromo_mask[startidx]), ID2=byte(chromo_mask[endidx]).
        # Under this pipeline the equivalent is C-order flattened indexing.
        chromo_mask_flat = np.asarray(chromo_mask).ravel(order="C")
        ID1 = chromo_mask_flat[startidx]
        ID2 = chromo_mask_flat[endidx]

        QB = np.asarray(model_dict["av_field"], dtype=np.float64).reshape(-1)
        QL = np.asarray(model_dict["phys_length"], dtype=np.float64).reshape(-1) * RSun
        uu = np.asarray(model_dict["voxel_status"], dtype=np.uint8).reshape(-1)

        QB[ (uu & 4) != 4] = 0
        QL[ (uu & 4) != 4] = 0
        ID1[(uu & 4) != 4] = 0
        ID2[(uu & 4) != 4] = 0

        QB   = QB.reshape((Nx, Ny, sc[1]), order="F")
        QL   = QL.reshape((Nx, Ny, sc[1]), order="F")
        ID1   = ID1.reshape((Nx, Ny, sc[1]), order="F")
        ID2   = ID2.reshape((Nx, Ny, sc[1]), order="F")

        corona_Bavg         = QB[:, :, corona_base : sc[1]].T
        chromo_uniform_Bavg = QB[:, :, 0 : corona_base].T

        corona_L         = QL[:, :, corona_base : sc[1]].T
        chromo_uniform_L = QL[:, :, 0 : corona_base].T

        corona_ID1 = ID1[:, :, corona_base : sc[1]].T
        corona_ID2 = ID2[:, :, corona_base : sc[1]].T

        chromo_uniform_ID1 = ID1[:, :, 0 : corona_base].T
        chromo_uniform_ID2 = ID2[:, :, 0 : corona_base].T

        #pdb.set_trace()

        model_dt_varlist = [
                ('Nx', np.int32),
                ('Ny', np.int32),
                ('Nz', np.int32),
                ('chromo_layers', np.int32),
                ('corona_layers', np.int32),
                ('corona_base',   np.int32),
                ('DSun', np.float64),
                ('RSun', np.float64),
                ('b0Sun', np.float64),
                ('lonC', np.float64),
                ('latC', np.float64),
                ('dx',   np.float64),
                ('dy',   np.float64),
                ('dz_uniform', np.float64),
                ('obstime',    np.float64),
                ('dz', np.float32, s),
                ('Bx', np.float32, s),
                ('By', np.float32, s),
                ('Bz', np.float32, s),
                ('chromo_n0' ,  np.float32, chr_s),
                ('chromo_np' ,  np.float32, chr_s),
                ('chromo_nHI',  np.float32, chr_s),
                ('chromo_T0',   np.float32, chr_s),
                ('corona_Bavg', np.float32, cor_s),
                ('corona_L',    np.float32, cor_s),
                ('chromo_uniform_Bavg', np.float32, cor_s2),
                ('chromo_uniform_L',    np.float32, cor_s2),
                ('VoxelID',            np.byte, s),
                ('corona_ID1',         np.byte, cor_s),
                ('corona_ID2',         np.byte, cor_s),
                ('chromo_uniform_ID1', np.byte, cor_s2),
                ('chromo_uniform_ID2', np.byte, cor_s2)
        ]

        model_dt = np.dtype(model_dt_varlist)
        model = np.zeros(1, dtype=model_dt)

        for vars in model_dt_varlist:
            k, dt = vars[0], vars[1]
            val = dt(locals()[k])
            if len(val.shape) > 1:
                model[k] = val.astype(dt, order="C")
            else:
                model[k] = val
    
        return model, model_dt

    def load_model_hdf(self, file_name):
        data = self._build_hdf_chromo_data(file_name)
        return self.load_model_dict(data.model, data.header)

    def load_model_sav(self, file_name):
        data = self._build_sav_chromo_data(file_name)
        return self.load_model_dict(data.model, data.header)

    def synth_model(self, model, model_dt, ebtel, ebtel_dt, freqlist, box_Nx, box_Ny, box_xc, box_yc, box_dx, box_dy, Tbase, nbase, Q0, a, b, SHtable=None, rot=0, projection=0, mode=0):
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
    
