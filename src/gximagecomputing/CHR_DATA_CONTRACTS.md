# CHR Data Contracts

This project supports two CHR input contracts that are normalized into one internal `ChromoModel` representation before calling `RenderGRFF`.

## 1) IDL Contract (`.sav`)

Entry point: `GXRadioImageComputing.load_model_sav(path)`

Expected fields in `BOX` (current CHR format):
- `DR`, `DZ`
- `BCUBE`, `CHROMO_BCUBE`
- `CHROMO_LAYERS`, `CORONA_BASE`
- `CHROMO_IDX`, `CHROMO_N`, `N_P`, `N_HI`, `CHROMO_T`
- `BASE.CHROMO_MASK`
- `AVFIELD`, `PHYSLENGTH`, `STATUS`, `STARTIDX`, `ENDIDX`

Metadata source:
- `BOX.INDEX` (`CRVAL1`, `CRVAL2`, `DSUN_OBS`, `DATE_OBS`)

## 2) pyAMPP HDF5 Contract (`.h5`)

Entry point: `GXRadioImageComputing.load_model_hdf(path)`

Expected datasets:
- `/chromo/dr`, `/chromo/dz`
- `/chromo/bcube`, `/chromo/chromo_bcube`
- `/chromo/chromo_layers`, `/chromo/corona_base`
- `/chromo/chromo_idx`, `/chromo/chromo_n`, `/chromo/n_p`, `/chromo/n_hi`, `/chromo/chromo_t`
- `/chromo/av_field`, `/chromo/phys_length`, `/chromo/voxel_status`
- `/chromo/start_idx`, `/chromo/end_idx`
- `/chromo/chromo_mask` (fallback: `/base/chromo_mask`)

Metadata source:
- `/chromo` attributes: `lon`, `lat`, `dsun_obs`, `obs_time`

## Normalized Internal Rules

- Both contracts are adapted to shared keys before rendering.
- `voxel_status` bitmask semantics are preserved.
- `start_idx`/`end_idx` are clamped to valid chromo-mask bounds.
- `obs_time` is normalized to `astropy.time.Time`.
