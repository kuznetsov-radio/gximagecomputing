;+
; NAME:
;   GXObserverGeometry helpers
;
; PURPOSE:
;   Shared observer and FOV geometry helpers for gximagecomputing IDL render
;   workflows. These routines mirror the Python observer_geometry module and
;   keep DLL model-triad and simbox computations in one place.
;-

function GXObserverGeometry__normalize_name, name
 if n_elements(name) eq 0 then return, ''
 s=strlowcase(strtrim(string(name), 2))
 case s of
  'terra': return, 'earth'
  'sdo': return, 'earth'
  'stereo a': return, 'stereo-a'
  'stereoa': return, 'stereo-a'
  'stereo ahead': return, 'stereo-a'
  'stereo b': return, 'stereo-b'
  'stereob': return, 'stereo-b'
  'stereo behind': return, 'stereo-b'
  'solar orbiter': return, 'solar orbiter'
  'solar-orbiter': return, 'solar orbiter'
  'solarorbiter': return, 'solar orbiter'
  'solo': return, 'solar orbiter'
  else: return, s
 endcase
end

function GXObserverGeometry__pb0r_arcsec_to_dsun, rsun_cm, rsun_arcsec
 if n_elements(rsun_cm) eq 0 or n_elements(rsun_arcsec) eq 0 then return, !values.d_nan
 if ~finite(double(rsun_cm)) or ~finite(double(rsun_arcsec)) then return, !values.d_nan
 if double(rsun_arcsec) le 0d then return, !values.d_nan
 return, double(rsun_cm) / tan(double(rsun_arcsec) * !dtor / 3600d)
end

function GXObserverGeometry__obstime_text, obstime
 if n_elements(obstime) eq 0 then return, ''
 return, strtrim(anytim(obstime, /ccsds), 2)
end

function GXObserverGeometry__observer_crln, obstime, observer_name, l0_deg, b0_deg, dsun_cm
 name=GXObserverGeometry__normalize_name(observer_name)
 obstext=GXObserverGeometry__obstime_text(obstime)

 if name ne '' then begin
  catch, err
  if err eq 0 then begin
   spice_name=''
   case name of
    'earth': spice_name='Earth'
    'stereo-a': spice_name='A'
    'stereo-b': spice_name='B'
    'solar orbiter': spice_name='SOLO'
    else: spice_name=''
   endcase
   if (spice_name ne '') && (obstext ne '') then begin
    errmsg=''
    pos=get_sunspice_lonlat(obstext, spice_name, system='Carrington', /degrees, errmsg=errmsg)
    if (errmsg eq '') && (n_elements(pos) ge 2) then begin
     crln=double(pos[1])
     catch, /cancel
     return, crln
    endif
   endif
   catch, /cancel
  endif else begin
   catch, /cancel
  endelse
 endif

 if finite(double(l0_deg)) then begin
  catch, err
  if err eq 0 then begin
   errmsg=''
   earth_pos=get_sunspice_lonlat(obstext, 'Earth', system='Carrington', /degrees, errmsg=errmsg)
   if (errmsg eq '') && (n_elements(earth_pos) ge 2) then begin
    crln=double(earth_pos[1]) + double(l0_deg)
    catch, /cancel
    return, crln
   endif
   catch, /cancel
  endif else begin
   catch, /cancel
  endelse
 endif

 carr_earth=(tim2carr(obstext))[0]
 return, double(carr_earth) + double(l0_deg)
end

function GXObserverGeometry__get_saved_observer_state, observer_struct, rsun_cm_default
 empty={observer_name:'', l0_deg:!values.d_nan, b0_deg:!values.d_nan, dsun_cm:!values.d_nan, $
        rsun_cm:double(rsun_cm_default), rsun_arcsec:!values.d_nan, source:''}
 if n_elements(observer_struct) eq 0 then return, empty
 if size(observer_struct, /type) ne 8 then return, empty

 state=empty
 if tag_exist(observer_struct, 'NAME') then state.observer_name=GXObserverGeometry__normalize_name(observer_struct.name)
 if tag_exist(observer_struct, 'PB0R') then begin
  if tag_exist(observer_struct.pb0r, 'L0_DEG') then state.l0_deg=double(observer_struct.pb0r.l0_deg)
  if tag_exist(observer_struct.pb0r, 'B0_DEG') then state.b0_deg=double(observer_struct.pb0r.b0_deg)
  if tag_exist(observer_struct.pb0r, 'RSUN_ARCSEC') then state.rsun_arcsec=double(observer_struct.pb0r.rsun_arcsec)
 endif
 if tag_exist(observer_struct, 'EPHEMERIS') then begin
  if tag_exist(observer_struct.ephemeris, 'DSUN_CM') then state.dsun_cm=double(observer_struct.ephemeris.dsun_cm)
  if tag_exist(observer_struct.ephemeris, 'RSUN_CM') then state.rsun_cm=double(observer_struct.ephemeris.rsun_cm)
 endif
 if ~finite(state.dsun_cm) and finite(state.rsun_arcsec) then $
  state.dsun_cm=GXObserverGeometry__pb0r_arcsec_to_dsun(state.rsun_cm, state.rsun_arcsec)
 if state.observer_name eq '' then state.observer_name='custom'
 if finite(state.l0_deg) and finite(state.b0_deg) and finite(state.dsun_cm) then state.source='saved_observer_metadata'
 return, state
end

function GXObserverGeometry__resolve_named_observer, obstime, observer_name, rsun_cm_default
 name=GXObserverGeometry__normalize_name(observer_name)
 obstext=GXObserverGeometry__obstime_text(obstime)
 state={observer_name:name, l0_deg:!values.d_nan, b0_deg:!values.d_nan, dsun_cm:!values.d_nan, $
        rsun_cm:double(rsun_cm_default), rsun_arcsec:!values.d_nan, source:'named_observer'}
 if name eq '' then return, state

 case name of
  'earth': angles=pb0r(obstext, /earth, l0=l0, /arcsec)
  'stereo-a': angles=pb0r(obstext, stereo='A', l0=l0, /arcsec)
  'stereo-b': angles=pb0r(obstext, stereo='B', l0=l0, /arcsec)
  'solar orbiter': angles=pb0r(obstext, /solo, l0=l0, /arcsec)
  else: message, 'Cannot resolve observer geometry for observer: '+name
 endcase
 state.l0_deg=double(l0)
 state.b0_deg=double(angles[1])
 state.rsun_arcsec=double(angles[2])
 state.dsun_cm=GXObserverGeometry__pb0r_arcsec_to_dsun(state.rsun_cm, state.rsun_arcsec)
 return, state
end

function GXObserverGeometry__saved_fov, observer_struct
 if n_elements(observer_struct) eq 0 then return, 0
 if size(observer_struct, /type) ne 8 then return, 0
 if ~tag_exist(observer_struct, 'FOV') then return, 0
 fov=observer_struct.fov
 if ~tag_exist(fov, 'XC_ARCSEC') then return, 0
 if ~tag_exist(fov, 'YC_ARCSEC') then return, 0
 if ~tag_exist(fov, 'XSIZE_ARCSEC') then return, 0
 if ~tag_exist(fov, 'YSIZE_ARCSEC') then return, 0
 return, fov
end

function GXObserverGeometry__saved_square_fov, observer_struct
 if n_elements(observer_struct) eq 0 then return, 0b
 if size(observer_struct, /type) ne 8 then return, 0b
 if tag_exist(observer_struct, 'FOV') then begin
  if tag_exist(observer_struct.fov, 'SQUARE') then return, byte(observer_struct.fov.square)
 endif
 if tag_exist(observer_struct, 'FOV_BOX') then begin
  if tag_exist(observer_struct.fov_box, 'SQUARE') then return, byte(observer_struct.fov_box.square)
 endif
 return, 0b
end

function GXObserverGeometry__extract_execute_geometry, execute_text
 spec={valid:0b, coord_mode:'', observer_name:'earth', center_x:!values.d_nan, center_y:!values.d_nan, $
       nx:0L, ny:0L, nz:0L, dx_km:!values.d_nan}
 if n_elements(execute_text) eq 0 then return, spec
 text=strtrim(string(execute_text), 2)
 if text eq '' then return, spec
 tokens=strsplit(text, ' ', /extract)
 nt=n_elements(tokens)
 if nt le 0 then return, spec
 i=0L
 while i lt nt do begin
  tok=strlowcase(tokens[i])
  case tok of
   '--coords': if (i+2L) lt nt then begin
     spec.center_x=double(tokens[i+1L])
     spec.center_y=double(tokens[i+2L])
     i+=2L
   endif
   '--hpc': spec.coord_mode='hpc'
   '--hgc': spec.coord_mode='hgc'
   '--hgs': spec.coord_mode='hgs'
   '--box-dims': if (i+3L) lt nt then begin
     spec.nx=long(tokens[i+1L])
     spec.ny=long(tokens[i+2L])
     spec.nz=long(tokens[i+3L])
     i+=3L
   endif
  '--dx-km': if (i+1L) lt nt then begin
     spec.dx_km=double(tokens[i+1L])
     i+=1L
   endif
   '--observer-name': if (i+1L) lt nt then begin
     spec.observer_name=GXObserverGeometry__normalize_name(tokens[i+1L])
     i+=1L
   endif
   else:
  endcase
  i++
 endwhile
 if spec.coord_mode eq '' then spec.coord_mode='hpc'
 if ~finite(spec.center_x) or ~finite(spec.center_y) then return, spec
 if ~finite(spec.dx_km) then return, spec
 if (spec.nx le 0L) or (spec.ny le 0L) or (spec.nz le 0L) then return, spec
 spec.valid=1b
 return, spec
end

function GXObserverGeometry__vector_norm, v
 return, sqrt(total(double(v)^2))
end

function GXObserverGeometry__cross, a, b
 return, [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]
end

function GXObserverGeometry__unit_from_lonlat, lon_deg, lat_deg
 lon_rad=double(lon_deg)*!dtor
 lat_rad=double(lat_deg)*!dtor
 return, [cos(lat_rad)*cos(lon_rad), cos(lat_rad)*sin(lon_rad), sin(lat_rad)]
end

function GXObserverGeometry__model_box_corners_hgs, model
 lon=double(model.lonC)
 lat=double(model.latC)
 rsun=double(model.RSun)
 nx=long(model.Nx)
 ny=long(model.Ny)
 dx=double(model.dx)
 dy=double(model.dy)
 dz=double(model.dz)

 z_top=max(total(dz, 3))
 half_x=0.5d*double(nx)*dx
 half_y=0.5d*double(ny)*dy

 radial=GXObserverGeometry__unit_from_lonlat(lon, lat)
 center=radial*rsun

 earth_state=GXObserverGeometry__resolve_named_observer(model.obstime, 'earth', rsun)
 earth_unit=GXObserverGeometry__unit_from_lonlat(0d, earth_state.b0_deg)
 earth_pos=earth_unit*earth_state.dsun_cm

 los=earth_pos-center
 los_norm=GXObserverGeometry__vector_norm(los)
 if ~(finite(los_norm) && (los_norm gt 0d)) then message, 'Invalid Earth LOS vector for model box construction.'
 los=los/los_norm

 x_axis=GXObserverGeometry__cross(los, radial)
 x_norm=GXObserverGeometry__vector_norm(x_axis)
 if ~(finite(x_norm) && (x_norm gt 0d)) then begin
  north_ref=[0d, 0d, 1d]
  y_axis=north_ref-total(north_ref*radial)*radial
  y_norm=GXObserverGeometry__vector_norm(y_axis)
  if ~(finite(y_norm) && (y_norm gt 0d)) then message, 'Failed to construct tangent-plane basis.'
  y_axis=y_axis/y_norm
  x_axis=GXObserverGeometry__cross(y_axis, radial)
  x_axis=x_axis/GXObserverGeometry__vector_norm(x_axis)
 endif else begin
  x_axis=x_axis/x_norm
  y_axis=GXObserverGeometry__cross(radial, x_axis)
  y_axis=y_axis/GXObserverGeometry__vector_norm(y_axis)
 endelse

 corners=dblarr(8, 3)
 idx=0L
 for iz=0, 1 do begin
  z_cm=(iz eq 0) ? 0d : z_top
  for iy=0, 1 do begin
   y_cm=(iy eq 0) ? -half_y : half_y
   for ix=0, 1 do begin
    x_cm=(ix eq 0) ? -half_x : half_x
    corners[idx, *]=center + x_cm*x_axis + y_cm*y_axis + z_cm*radial
    idx++
   endfor
  endfor
 endfor
 return, corners
end

function GXObserverGeometry__corners_hg_hecr, corners_cm
 dims=size(corners_cm, /dimensions)
 n=dims[0]
 out=dblarr(n, 3)
 for i=0L, n-1L do begin
  x=double(corners_cm[i,0])
  y=double(corners_cm[i,1])
  z=double(corners_cm[i,2])
  r=sqrt(x*x+y*y+z*z)
  out[i,2]=r/100d
  if r gt 0d then begin
   out[i,1]=asin(z/r)/!dtor
   out[i,0]=atan(y, x)/!dtor
  endif
 endfor
 return, out
end

function GXObserverGeometry__corners_hg_to_cartesian_cm, hg
 dims=size(hg, /dimensions)
 if n_elements(dims) lt 2 then return, 0
 if dims[1] lt 3 then return, 0
 n=dims[0]
 out=dblarr(n, 3)
 for i=0L, n-1L do begin
  lon=double(hg[i,0])*!dtor
  lat=double(hg[i,1])*!dtor
  r_cm=double(hg[i,2])*100d
  out[i,0]=r_cm*cos(lat)*cos(lon)
  out[i,1]=r_cm*cos(lat)*sin(lon)
  out[i,2]=r_cm*sin(lat)
 endfor
 return, out
end

function GXObserverGeometry__execute_box_corners_hg, model, geometry, execute_text
 if n_elements(execute_text) eq 0 then return, 0
 spec=GXObserverGeometry__extract_execute_geometry(execute_text)
 if spec.valid eq 0b then return, 0

 geom_obs=GXObserverGeometry__resolve_named_observer(model.obstime, spec.observer_name, double(model.RSun))
 if ~finite(double(geom_obs.dsun_cm)) then return, 0

 hecr0=!values.d_nan
 case spec.coord_mode of
  'hpc': begin
   wcs_conv_hpc_hg, [double(spec.center_x)], [double(spec.center_y)], lon0_arr, lat0_arr, hecr_arr, $
    dsun_obs=double(geom_obs.dsun_cm)/100d, b0_angle=double(geom_obs.b0_deg), l0_angle=double(geom_obs.l0_deg), /arcseconds
   lon0=double(lon0_arr[0]) & lat0=double(lat0_arr[0]) & hecr0=double(hecr_arr[0])
  end
  'hgc': begin
   simwcs=wcs_2d_simulate(1, date_obs=anytim(model.obstime, /vms), $
    dsun_obs=double(geom_obs.dsun_cm)/100d, hgln_obs=double(geom_obs.l0_deg), hglt_obs=double(geom_obs.b0_deg))
   lon_tmp=[double(spec.center_x)]
   wcs_conv_cr_hg, lon_tmp, lon0_arr, wcs=simwcs
   lon0=double(lon0_arr[0]) & lat0=double(spec.center_y) & hecr0=double(model.RSun)/100d
  end
  else: begin
   lon0=double(spec.center_x) & lat0=double(spec.center_y) & hecr0=double(model.RSun)/100d
  end
 endcase
 if ~finite(hecr0) then hecr0=double(model.RSun)/100d

 anchor_r_cm=double(hecr0)*100d
 radial=GXObserverGeometry__unit_from_lonlat(lon0, lat0)
 anchor_cm=radial*anchor_r_cm

 north_ref=[0d, 0d, 1d]
 y_axis=north_ref-total(north_ref*radial)*radial
 y_norm=GXObserverGeometry__vector_norm(y_axis)
 if ~(finite(y_norm) && (y_norm gt 0d)) then begin
  east_ref=[0d, 1d, 0d]
  y_axis=east_ref-total(east_ref*radial)*radial
  y_norm=GXObserverGeometry__vector_norm(y_axis)
  if ~(finite(y_norm) && (y_norm gt 0d)) then return, 0
 endif
 y_axis=y_axis/y_norm
 x_axis=GXObserverGeometry__cross(y_axis, radial)
 x_axis=x_axis/GXObserverGeometry__vector_norm(x_axis)

 xdim=double(spec.nx)*double(spec.dx_km)*1d5
 ydim=double(spec.ny)*double(spec.dx_km)*1d5
 zdim=double(spec.nz)*double(spec.dx_km)*1d5
 half_x=0.5d*xdim
 half_y=0.5d*ydim

 corners_cm=dblarr(8, 3)
 idx=0L
 for iz=0,1 do begin
  z_off=(iz eq 0) ? 0d : zdim
  for iy=0,1 do begin
   y_off=(iy eq 0) ? -half_y : half_y
   for ix=0,1 do begin
    x_off=(ix eq 0) ? -half_x : half_x
    corners_cm[idx,*]=anchor_cm + x_off*x_axis + y_off*y_axis + z_off*radial
    idx++
   endfor
  endfor
 endfor
 return, GXObserverGeometry__corners_hg_hecr(corners_cm)
end

function GXObserverGeometry__index_box_corners_hg, model, box_struct, index_struct
 if n_elements(box_struct) eq 0 then return, 0
 if size(box_struct, /type) ne 8 then return, 0
 if n_elements(index_struct) eq 0 then begin
  if tag_exist(box_struct, 'INDEX') then index_struct=box_struct.index else return, 0
 endif
 if size(index_struct, /type) ne 8 then return, 0
 if ~tag_exist(index_struct, 'CRVAL1') then return, 0
 if ~tag_exist(index_struct, 'CRVAL2') then return, 0
 if ~tag_exist(box_struct, 'DR') then return, 0
 if ~tag_exist(box_struct, 'BCUBE') then return, 0

 lon_ref=double(index_struct.crval1)
 lat_ref=double(index_struct.crval2)
 dsun_obs=!values.d_nan
 if tag_exist(index_struct, 'DSUN_OBS') then dsun_obs=double(index_struct.dsun_obs)
 if ~finite(dsun_obs) then return, 0
 if dsun_obs gt 1d12 then dsun_obs=dsun_obs/100d
 hgln_obs=0d
 if tag_exist(index_struct, 'HGLN_OBS') then hgln_obs=double(index_struct.hgln_obs)
 hglt_obs=!values.d_nan
 if tag_exist(index_struct, 'HGLT_OBS') then hglt_obs=double(index_struct.hglt_obs) $
  else if tag_exist(index_struct, 'SOLAR_B0') then hglt_obs=double(index_struct.solar_b0)
 if ~finite(hglt_obs) then hglt_obs=0d

 lon0=lon_ref
 if tag_exist(index_struct, 'CTYPE1') then begin
  ctype1=strupcase(strtrim(string(index_struct.ctype1), 2))
  if strpos(ctype1, 'CRLN') ge 0 then begin
   simwcs=wcs_2d_simulate(1, date_obs=anytim(model.obstime, /vms), $
    dsun_obs=dsun_obs, hgln_obs=hgln_obs, hglt_obs=hglt_obs)
   tmp=[lon_ref]
   wcs_conv_cr_hg, tmp, lon0_arr, wcs=simwcs
   lon0=double(lon0_arr[0])
  endif
 endif

 radial=GXObserverGeometry__unit_from_lonlat(lon0, lat_ref)
 anchor_cm=radial*double(model.RSun)
 north_ref=[0d, 0d, 1d]
 y_axis=north_ref-total(north_ref*radial)*radial
 y_norm=GXObserverGeometry__vector_norm(y_axis)
 if ~(finite(y_norm) && (y_norm gt 0d)) then begin
  east_ref=[0d, 1d, 0d]
  y_axis=east_ref-total(east_ref*radial)*radial
  y_norm=GXObserverGeometry__vector_norm(y_axis)
  if ~(finite(y_norm) && (y_norm gt 0d)) then return, 0
 endif
 y_axis=y_axis/y_norm
 x_axis=GXObserverGeometry__cross(y_axis, radial)
 x_axis=x_axis/GXObserverGeometry__vector_norm(x_axis)

 dims=size(box_struct.bcube, /dimensions)
 if n_elements(dims) lt 3 then return, 0
 xdim=double(dims[0])*double(box_struct.dr[0])*double(model.RSun)
 ydim=double(dims[1])*double(box_struct.dr[1])*double(model.RSun)
 zdim=double(dims[2])*double(box_struct.dr[2])*double(model.RSun)
 half_x=0.5d*xdim
 half_y=0.5d*ydim

 corners_cm=dblarr(8, 3)
 idx=0L
 for iz=0,1 do begin
  z_off=(iz eq 0) ? 0d : zdim
  for iy=0,1 do begin
   y_off=(iy eq 0) ? -half_y : half_y
   for ix=0,1 do begin
    x_off=(ix eq 0) ? -half_x : half_x
    corners_cm[idx,*]=anchor_cm + x_off*x_axis + y_off*y_axis + z_off*radial
    idx++
   endfor
  endfor
 endfor
 return, GXObserverGeometry__corners_hg_hecr(corners_cm)
end

function GXObserverGeometry__projected_corners, model, geometry, execute_text=execute_text, box_struct=box_struct, index_struct=index_struct
 hg=GXObserverGeometry__index_box_corners_hg(model, box_struct, index_struct)
 have_exec=0b
 if size(hg, /n_dimensions) eq 2 then begin
  dims=size(hg, /dimensions)
  if (n_elements(dims) ge 2) and (dims[1] ge 3) and (dims[0] gt 0) then have_exec=1b
 endif
 if ~have_exec then begin
  hg=GXObserverGeometry__execute_box_corners_hg(model, geometry, execute_text)
  if size(hg, /n_dimensions) eq 2 then begin
   dims=size(hg, /dimensions)
   if (n_elements(dims) ge 2) and (dims[1] ge 3) and (dims[0] gt 0) then have_exec=1b
  endif
 endif
 if have_exec then corners_cm=GXObserverGeometry__corners_hg_to_cartesian_cm(hg) else begin
  corners_cm=GXObserverGeometry__model_box_corners_hgs(model)
  hg=GXObserverGeometry__corners_hg_hecr(corners_cm)
 endelse
 n=n_elements(hg[*,0])
 x_hcc=dblarr(n)
 y_hcc=dblarr(n)
 z_hcc=dblarr(n)
 hpln=dblarr(n)
 hplt=dblarr(n)
 distance=dblarr(n)
 ; Use the explicit HG -> HCC -> HPC route and disable obscuration masking so
 ; all eight model-box corners contribute to the observer-footprint rectangle,
 ; matching the Python SkyCoord.transform_to(...) behavior.
 wcs_conv_hg_hcc, hg[*,0], hg[*,1], x_hcc, y_hcc, z_hcc, hecr=hg[*,2], $
  dsun_obs=double(geometry.dsun_cm)/100d, b0_angle=double(geometry.b0_deg), $
  l0_angle=double(geometry.l0_deg)
 wcs_conv_hcc_hpc, x_hcc, y_hcc, hpln, hplt, distance, solz=z_hcc, $
  dsun_obs=double(geometry.dsun_cm)/100d, /arcseconds, nomask=1b
 return, {corners_cm: corners_cm, tx_arcsec: hpln, ty_arcsec: hplt}
end

function GXResolveObserverGeometry, model, index_struct=index_struct, observer_struct=observer_struct, $
                                   DSUN=dsun_kw, LONC=lonc_kw, B0SUN=b0sun_kw, $
                                   observer_name=observer_name_kw, $
                                   recompute_observer_ephemeris=recompute_observer_ephemeris
 rsun_cm=double(model.RSun)
 saved=GXObserverGeometry__get_saved_observer_state(observer_struct, rsun_cm)
 have_saved=(saved.source ne '')

 lon_ref=!values.d_nan
 if n_elements(index_struct) gt 0 then begin
  if size(index_struct, /type) eq 8 then begin
   if tag_exist(index_struct, 'CRVAL1') then lon_ref=double(index_struct.crval1)
  endif
 endif
 if ~finite(lon_ref) then lon_ref=double(model.lonC)

 source='default_earth'
 cli_observer_name=''
 if n_elements(observer_name_kw) gt 0 then cli_observer_name=GXObserverGeometry__normalize_name(observer_name_kw)
 if cli_observer_name ne '' then begin
  state=GXObserverGeometry__resolve_named_observer(model.obstime, cli_observer_name, rsun_cm)
  source='cli_observer'
 endif else begin
  if keyword_set(recompute_observer_ephemeris) then begin
   name=(have_saved) ? saved.observer_name : 'earth'
   if cli_observer_name ne '' then name=cli_observer_name
   state=GXObserverGeometry__resolve_named_observer(model.obstime, name, rsun_cm)
   source='recompute_observer_ephemeris'
  endif else begin
   if have_saved then begin
    state=saved
    source=saved.source
   endif else begin
    state=GXObserverGeometry__resolve_named_observer(model.obstime, 'earth', rsun_cm)
    source='default_earth'
   endelse
  endelse
 endelse

 render_lonc=double(model.lonC)
 render_dsun=double(model.DSun)
 render_b0=double(model.b0Sun)

 if source ne 'default_earth' then begin
  if finite(state.dsun_cm) then render_dsun=double(state.dsun_cm)
  if finite(state.b0_deg) then render_b0=double(state.b0_deg)
 endif

 if finite(lon_ref) && finite(state.l0_deg) && (source ne 'default_earth') then begin
  observer_crln=GXObserverGeometry__observer_crln(model.obstime, state.observer_name, state.l0_deg, state.b0_deg, state.dsun_cm)
  render_lonc=double(lon_ref)-double(observer_crln)
 endif

 if n_elements(dsun_kw) gt 0 then begin
  render_dsun=double(dsun_kw)
  source='cli_triad'
 endif
 if n_elements(lonc_kw) gt 0 then begin
  render_lonc=double(lonc_kw)
  source='cli_triad'
 endif
 if n_elements(b0sun_kw) gt 0 then begin
  render_b0=double(b0sun_kw)
  source='cli_triad'
 endif

 if state.observer_name eq '' then state.observer_name='earth'
 return, {observer_name: state.observer_name, $
          l0_deg: double(state.l0_deg), $
          b0_deg: double(state.b0_deg), $
          dsun_cm: double(state.dsun_cm), $
          rsun_cm: double(state.rsun_cm), $
          rsun_arcsec: double(state.rsun_arcsec), $
          render_lonc_deg: double(render_lonc), $
          render_b0_deg: double(render_b0), $
          render_dsun_cm: double(render_dsun), $
          source: source}
end

function GXComputeInscribingFOV, model, geometry, execute_text=execute_text, box_struct=box_struct, index_struct=index_struct, square_fov=square_fov
 proj=GXObserverGeometry__projected_corners(model, geometry, execute_text=execute_text, box_struct=box_struct, index_struct=index_struct)
 tx=double(proj.tx_arcsec)
 ty=double(proj.ty_arcsec)
 good=where(finite(tx) and finite(ty), ngood)
 if ngood le 0 then message, 'Failed to project model box corners into observer frame.'
 xmin=min(tx[good], max=xmax)
 ymin=min(ty[good], max=ymax)
 if keyword_set(square_fov) then begin
  side=(xmax-xmin) > (ymax-ymin)
  xc_mid=0.5d*(xmin+xmax)
  yc_mid=0.5d*(ymin+ymax)
  xmin=xc_mid-0.5d*side
  xmax=xc_mid+0.5d*side
  ymin=yc_mid-0.5d*side
  ymax=yc_mid+0.5d*side
 endif
 return, {xc_arcsec: 0.5d*(xmin+xmax), $
          yc_arcsec: 0.5d*(ymin+ymax), $
          xsize_arcsec: xmax-xmin, $
          ysize_arcsec: ymax-ymin, $
          xmin_arcsec: xmin, xmax_arcsec: xmax, $
          ymin_arcsec: ymin, ymax_arcsec: ymax}
end

function GXComputeInscribingFOVBox, model, geometry, execute_text=execute_text, box_struct=box_struct, index_struct=index_struct, square_fov=square_fov
 fov=GXComputeInscribingFOV(model, geometry, execute_text=execute_text, box_struct=box_struct, index_struct=index_struct, square_fov=square_fov)
 proj=GXObserverGeometry__projected_corners(model, geometry, execute_text=execute_text, box_struct=box_struct, index_struct=index_struct)
 corners=proj.corners_cm
 obs_unit=GXObserverGeometry__unit_from_lonlat(geometry.l0_deg, geometry.b0_deg)
 zmm=dblarr(8)
 for i=0L, 7L do zmm[i]=total(double(corners[i,*])*obs_unit)/1d8
 zmin=min(zmm, max=zmax)
 span=max(1d-6, zmax-zmin)
 zpad=0.10d*span
 return, create_struct(fov, 'zmin_mm', zmin-zpad, 'zmax_mm', zmax+zpad)
end

function GXResolveSimboxFromObserverAndModel, saved_fov=saved_fov, computed_fov=computed_fov, $
                                             XC=xc, YC=yc, XSIZE=xsize, YSIZE=ysize
 if n_elements(xc) gt 0 && n_elements(yc) gt 0 && n_elements(xsize) gt 0 && n_elements(ysize) gt 0 then $
  return, {source:'explicit', xc_arcsec:double(xc), yc_arcsec:double(yc), $
           xsize_arcsec:double(xsize), ysize_arcsec:double(ysize)}

 if n_elements(saved_fov) gt 0 then begin
  if size(saved_fov, /type) eq 8 then begin
   if tag_exist(saved_fov, 'XC_ARCSEC') and tag_exist(saved_fov, 'YC_ARCSEC') and $
      tag_exist(saved_fov, 'XSIZE_ARCSEC') and tag_exist(saved_fov, 'YSIZE_ARCSEC') then begin
    return, {source:'saved_observer_fov', xc_arcsec:double(saved_fov.xc_arcsec), yc_arcsec:double(saved_fov.yc_arcsec), $
             xsize_arcsec:double(saved_fov.xsize_arcsec), ysize_arcsec:double(saved_fov.ysize_arcsec)}
   endif
  endif
 endif

 if n_elements(computed_fov) gt 0 then begin
  if size(computed_fov, /type) eq 8 then begin
   return, {source:'inscribing_fov', xc_arcsec:double(computed_fov.xc_arcsec), yc_arcsec:double(computed_fov.yc_arcsec), $
            xsize_arcsec:double(computed_fov.xsize_arcsec), ysize_arcsec:double(computed_fov.ysize_arcsec)}
  endif
 endif

 return, 0
end
