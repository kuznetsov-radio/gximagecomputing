;+
; NAME:
;   LoadGXmodel
;
; PURPOSE:
;   Load a GX model file (.sav or .h5) into the DLL-ready model structure used
;   by the microwave and EUV rendering interface.
;
; NOTES:
;   - Explicit DSun/LONC/B0SUN overrides take precedence.
;   - File-first loading uses saved observer metadata when available.
;   - `observer_struct=...` returns the saved observer metadata alongside the
;     DLL-ready model without rereading the input file.
;   - `index_struct=...` returns the preserved saved index/header structure
;     alongside the DLL-ready model without rereading the input file.
;-
function LoadGXmodel, infile, noVoxelID=noVoxelID, newTime=newTime, DSUN=dsun_kw, LONC=lonc_kw, B0SUN=b0sun_kw, $
                     recompute_observer_ephemeris=recompute_observer_ephemeris, observer_name=observer_name_kw, $
                     observer_struct=observer_out, index_struct=index_out, execute_text=execute_out, box_struct=box_out
 box=LoadGXmodel__load_box(infile)
 if arg_present(box_out) then box_out=box
 if arg_present(observer_out) then begin
  if tag_exist(box, 'OBSERVER') then observer_out=box.observer else observer_out=0
 endif
 if arg_present(index_out) then begin
  if tag_exist(box, 'INDEX') then index_out=box.index else index_out=0
 endif
 if arg_present(execute_out) then begin
  if tag_exist(box, 'EXECUTE') then execute_out=box.execute else execute_out=''
 endif

 obstime=0d
 if tag_exist(box, 'INDEX') and tag_exist(box.index, 'DATE_OBS') then obstime=anytim(box.index.date_obs)

 lonC=0d
 latC=0d
 RSun=6.96d10
 saved_observer_name='earth'
 lon_ref=!values.d_nan

 if tag_exist(box, 'INDEX') then begin
  ; Prefer robust WCS conversion when full header is available.
  catch, err
  if err eq 0 then begin
   setenv, 'WCS_RSUN=6.96d8' ; same command as GX Simulator routines
   wcs=fitshead2wcs(box.index)
   wcs_convert_from_coord, wcs, wcs.crval, 'hg', lonC, latC
   RSun=wcs_rsun(unit='cm')
   catch, /cancel
  endif else begin
   catch, /cancel
  endelse

  if tag_exist(box, 'OBSERVER') then begin
   if tag_exist(box.observer, 'NAME') then begin
    obsname_tmp=strtrim(string(box.observer.name), 2)
    if obsname_tmp ne '' then saved_observer_name=strlowcase(obsname_tmp)
   endif
  endif else if tag_exist(box.index, 'OBSERVER') then saved_observer_name=strlowcase(strtrim(string(box.index.observer), 2)) $
  else if tag_exist(box.index, 'OBSERVATORY') then saved_observer_name=strlowcase(strtrim(string(box.index.observatory), 2))

  ; Fallback to direct heliographic tags when WCS parsing is unavailable.
  if (lonC eq 0d) and tag_exist(box.index, 'CRVAL1') then lonC=double(box.index.crval1)
  if tag_exist(box.index, 'CRVAL2') then latC=double(box.index.crval2)
  if tag_exist(box.index, 'CRVAL1') then lon_ref=double(box.index.crval1)
 endif

 if exist(newTime) then begin
  obstime1=anytim(newTime)
  ddays=(obstime1-obstime)/86400d
  lonC+=diff_rot(ddays, latC, /synodic)
  obstime=obstime1
 endif

 recompute_obs=keyword_set(recompute_observer_ephemeris)
 if ~recompute_obs then begin
  DSun=!values.d_nan
  b0Sun=!values.d_nan
  has_observer=0b
  has_observer_pb0r=0b
  has_observer_ephemeris=0b
  if tag_exist(box, 'OBSERVER') then begin
   has_observer=1b
   if tag_exist(box.observer, 'EPHEMERIS') then begin
    has_observer_ephemeris=1b
    if tag_exist(box.observer.ephemeris, 'DSUN_CM') then DSun=double(box.observer.ephemeris.dsun_cm)
   endif
   if tag_exist(box.observer, 'PB0R') then begin
    has_observer_pb0r=1b
    if tag_exist(box.observer.pb0r, 'B0_DEG') then b0Sun=double(box.observer.pb0r.b0_deg)
   endif
  endif
  if ~finite(DSun) and tag_exist(box, 'INDEX') and tag_exist(box.index, 'DSUN_OBS') then DSun=double(box.index.dsun_obs)
  if ~finite(b0Sun) and tag_exist(box, 'INDEX') and tag_exist(box.index, 'SOLAR_B0') then b0Sun=double(box.index.solar_b0)
  if finite(lon_ref) and tag_exist(box, 'INDEX') and tag_exist(box.index, 'CRLN_OBS') then lonC=lon_ref-double(box.index.crln_obs)
  ; pyAMPP HDF5 index stores DSUN_OBS in meters; RenderGRFF expects centimeters.
  ; Keep backward compatibility with legacy SAV values already in centimeters.
  if finite(DSun) and (DSun lt 1d12) then DSun=DSun*100d
  if ~finite(DSun) or ~finite(b0Sun) then recompute_obs=1b
 endif

 if recompute_obs then begin
  observer_name='earth'
 if n_elements(observer_name_kw) gt 0 then observer_name=strlowcase(strtrim(string(observer_name_kw), 2)) $
 else if n_elements(saved_observer_name) gt 0 then observer_name=saved_observer_name
  obstext=strtrim(anytim(obsTime, /ccsds), 2)

  case observer_name of
   'earth': angles=pb0r(obstext, /earth, l0=l0)
   'sdo': angles=pb0r(obstext, /earth, l0=l0)
   'stereo-a': angles=pb0r(obstext, stereo='A', l0=l0)
   'stereo b': angles=pb0r(obstext, stereo='B', l0=l0)
   'stereo-b': angles=pb0r(obstext, stereo='B', l0=l0)
   'solar orbiter': angles=pb0r(obstext, /solo, l0=l0)
   'solo': angles=pb0r(obstext, /solo, l0=l0)
   else: message, 'Cannot recompute observer ephemeris for observer: '+observer_name
  endcase

  b0Sun=double(angles[1])
  DSun=RSun/tan((double(angles[2])/60d)*!dtor)
  carr_earth=(tim2carr(obstext))[0]
  observer_crln=carr_earth+double(l0)
  if finite(lon_ref) then lonC=lon_ref-observer_crln else lonC=-observer_crln
 endif

 if n_elements(dsun_kw) gt 0 then DSun=double(dsun_kw)
 if n_elements(lonc_kw) gt 0 then lonC=double(lonc_kw)
 if n_elements(b0sun_kw) gt 0 then b0Sun=double(b0sun_kw)

 dx=box.dr[0]*RSun
 dy=box.dr[1]*RSun
 dz_uniform=box.dr[2]*RSun
 dz=box.dz*RSun

 s=size(box.dz, /dimensions)
 sc=size(box.bcube, /dimensions)

 Nx=s[0]
 Ny=s[1]
 Nz=s[2]

 chromo_layers=box.chromo_layers
 corona_layers=s[2]-box.chromo_layers
 corona_base=box.corona_base

 Bx=fltarr(s)
 By=fltarr(s)
 Bz=fltarr(s)

 Bx[*, *, 0 : chromo_layers-1]=box.chromo_bcube[*, *, *, 0]
 By[*, *, 0 : chromo_layers-1]=box.chromo_bcube[*, *, *, 1]
 Bz[*, *, 0 : chromo_layers-1]=box.chromo_bcube[*, *, *, 2]

 Bx[*, *, chromo_layers : s[2]-1]=box.bcube[*, *, box.corona_base : sc[2]-1, 0]
 By[*, *, chromo_layers : s[2]-1]=box.bcube[*, *, box.corona_base : sc[2]-1, 1]
 Bz[*, *, chromo_layers : s[2]-1]=box.bcube[*, *, box.corona_base : sc[2]-1, 2]

 chromo_n0=fltarr(s[0], s[1], chromo_layers)
 chromo_np=fltarr(s[0], s[1], chromo_layers)
 chromo_nHI=fltarr(s[0], s[1], chromo_layers)
 chromo_T0=fltarr(s[0], s[1], chromo_layers)

 chromo_n0[box.chromo_idx]=box.chromo_n
 chromo_np[box.chromo_idx]=box.n_p
 chromo_nHI[box.chromo_idx]=box.n_hi
 chromo_T0[box.chromo_idx]=box.chromo_T

 if tag_exist(box, 'BMED') then begin ;old version
  QB=fltarr(s[0], s[1], sc[2])
  QB[box.idx]=box.bmed
  QL=fltarr(s[0], s[1], sc[2])
  QL[box.idx]=box.length*RSun
  ID1=bytarr(s[0], s[1], s[2])
  ID1[box.idx]=box.base.chromo_mask[box.foot1]
  ID2=bytarr(s[0], s[1], s[2])
  ID2[box.idx]=box.base.chromo_mask[box.foot2]
 endif else if tag_exist(box, 'AVFIELD') then begin ;new version
  QB=box.avfield
  QL=box.physlength*RSun
  ID1=byte(box.base.chromo_mask[box.startidx])
  ID2=byte(box.base.chromo_mask[box.endidx])
  u=where((box.status and 4) ne 4, k)
  if k gt 0 then begin
   QB[u]=0
   QL[u]=0
   ID1[u]=0
   ID2[u]=0
  endif
 endif

 corona_Bavg=QB[*, *, box.corona_base : sc[2]-1]
 chromo_uniform_Bavg=QB[*, *, 0 : box.corona_base-1]

 corona_L=QL[*, *, box.corona_base : sc[2]-1]
 chromo_uniform_L=QL[*, *, 0 : box.corona_base-1]

 if ~keyword_set(NoVoxelID) then VoxelID=gx_box2id(box) $
 else VoxelID=bytarr(s)

 corona_ID1=ID1[*, *, box.corona_base : sc[2]-1]
 corona_ID2=ID2[*, *, box.corona_base : sc[2]-1]
 chromo_uniform_ID1=ID1[*, *, 0 : box.corona_base-1]
 chromo_uniform_ID2=ID2[*, *, 0 : box.corona_base-1]

 model={Nx: long(Nx), $
        Ny: long(Ny), $
        Nz: long(Nz), $
        chromo_layers: long(chromo_layers), $
        corona_layers: long(corona_layers), $
        corona_base: long(corona_base), $
        DSun: double(DSun), $
        RSun: double(Rsun), $
        b0Sun: double(b0Sun), $
        lonC: double(lonC), $
        latC: double(latC), $
        dx: double(dx), $
        dy: double(dy), $
        dz_uniform: double(dz_uniform), $
        obstime: double(obstime), $
        dz: float(dz), $
        Bx: float(Bx), $
        By: float(By), $
        Bz: float(Bz), $
        chromo_n0: float(chromo_n0), $
        chromo_np: float(chromo_np), $
        chromo_nHI: float(chromo_nHI), $
        chromo_T0: float(chromo_T0), $
        corona_Bavg: float(corona_Bavg), $
        corona_L: float(corona_L), $
        chromo_uniform_Bavg: float(chromo_uniform_Bavg), $
        chromo_uniform_L: float(chromo_uniform_L), $
        VoxelID: byte(VoxelID), $
        corona_ID1: byte(corona_ID1), $
        corona_ID2: byte(corona_ID2), $
        chromo_uniform_ID1: byte(chromo_uniform_ID1), $
        chromo_uniform_ID2: byte(chromo_uniform_ID2)}

 return, model
end
