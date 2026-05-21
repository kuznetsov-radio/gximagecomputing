pro RenderExampleEUV, MODelfile=modelfile, EBTELfile=ebtelfile, RESPonsefile=responsefile, LIBname=libname, $
                    INSTRument=instrument, DSUN=dsun_kw, LONC=lonc_kw, B0SUN=b0sun_kw, $
                    observer_name=observer_name_kw, recompute_observer_ephemeris=recompute_observer_ephemeris, $
                    DUMP_DLL_INPUTS=dump_dll_inputs, $
                    AUTO_FOV=auto_fov, $
                    USE_SAVED_FOV=use_saved_fov, $
                    XC=xc, YC=yc, DX=dx, DY=dy, NX=nx, NY=ny, $
                    OUTdir=outdir, OUTfile=outfile, NO_PLOT=no_plot, _extra=_extra

 exdir=file_dirname(routine_filepath('RenderExampleEUV'))
 repo_root=file_dirname(file_dirname(exdir))
 local_idlcodedir=repo_root+'/idlcode'
 local_bindir=repo_root+'/binaries'
 env_ebtel=getenv('GXIMAGECOMPUTING_EBTEL_PATH')

 if n_elements(modelfile) eq 0 then begin
  message, 'Please provide MODelfile=<path to CHR .sav or .h5 model>.', /continue
  return
 endif
 if n_elements(ebtelfile) eq 0 then begin
  if n_elements(env_ebtel) gt 0 then ebtelfile=env_ebtel else ebtelfile=''
 endif
 if n_elements(libname) eq 0 then begin
  if file_test(local_bindir+'/RenderGRFF_arm64.so') then begin
   libname=local_bindir+'/RenderGRFF_arm64.so'
  endif else if file_test(local_bindir+'/RenderGRFF_x86_64.so') then begin
   libname=local_bindir+'/RenderGRFF_x86_64.so'
  endif else if file_test(local_bindir+'/RenderGRFF.so') then begin
   libname=local_bindir+'/RenderGRFF.so'
  endif else if file_test(local_bindir+'/RenderGRFF_64.dll') then begin
   libname=local_bindir+'/RenderGRFF_64.dll'
  endif else begin
   message, 'No RenderGRFF library found in '+local_bindir+'. Please provide LIBname=<path>.', /continue
   return
  endelse
 endif
 if n_elements(instrument) eq 0 then instrument='aia'
 if n_elements(dx) eq 0 then dx=2.0
 if n_elements(dy) eq 0 then dy=2.0
 if n_elements(outdir) eq 0 then begin
  if strpos(strlowcase(!version.os_family), 'windows') ge 0 then begin
   outdir='C:/Temp/gximagecomputing_validation_groundtruth'
  endif else begin
   outdir='/tmp/gximagecomputing_validation_groundtruth'
  endelse
 endif
 if n_elements(outfile) eq 0 then outfile='EUVmaps.sav'
 spawn, 'mkdir -p "'+outdir+'"'
 if strmid(outdir, strlen(outdir)-1, 1) eq '/' then begin
  outpath=outdir+outfile
 endif else begin
  outpath=outdir+'/'+outfile
 endelse

 ; Ensure local helper routines are preferred.
 if file_test(local_idlcodedir, /directory) then begin
  if strpos(':'+!path+':', ':'+local_idlcodedir+':') lt 0 then !path=local_idlcodedir+':'+!path
 endif

 ; Analytical corona/heating defaults.
 Tbase=1d6
 nbase=1d8
 Q0=0.0217
 a=0.3
 b=2.7
 SHtable=dblarr(7, 7)
 w=[1.0d, 1.0d, 1.0d, 1.1d, 1.2d, 1.3d, 1.4d]
 for i=0, 6 do for j=0, 6 do SHtable[i, j]=w[i]*w[j]
 SHtable[6, 6]=0.1d

 print, 'Model:     ', modelfile
 if strtrim(ebtelfile, 2) eq '' then print, 'EBTEL:     none (DEM/heating disabled)'
 if strtrim(ebtelfile, 2) ne '' then print, 'EBTEL:     ', ebtelfile
 if n_elements(responsefile) gt 0 then print, 'Response:  ', responsefile
 if n_elements(dsun_kw) gt 0 then print, 'DSun override (cm): ', dsun_kw
 if n_elements(lonc_kw) gt 0 then print, 'lonC override (deg): ', lonc_kw
 if n_elements(b0sun_kw) gt 0 then print, 'b0Sun override (deg): ', b0sun_kw
 print, 'Library:   ', libname

 forward_function LoadGXmodel, LoadEBTEL, LoadEUVresponse, MakeSimulationBoxEUV, DefineCoronaParams, ReserveOutputSpaceEUV
 forward_function GXResolveObserverGeometry, GXComputeInscribingFOV, GXResolveSimboxFromObserverAndModel, $
  GXObserverGeometry__saved_fov, GXObserverGeometry__saved_square_fov
 resolve_routine, 'GXDumpDLLInputs', /either

 tm=systime(1)
 execute_text=''
 box_struct=0
 model=LoadGXmodel(modelfile, observer_struct=observer_struct, index_struct=index_struct, execute_text=execute_text, box_struct=box_struct)

 geom=GXResolveObserverGeometry(model, index_struct=index_struct, observer_struct=observer_struct, $
  DSun=dsun_kw, LONC=lonc_kw, B0SUN=b0sun_kw, observer_name=observer_name_kw, $
  recompute_observer_ephemeris=recompute_observer_ephemeris)
 model.DSun=double(geom.render_dsun_cm)
 model.lonC=double(geom.render_lonc_deg)
 model.b0Sun=double(geom.render_b0_deg)

 has_cli_observer_override=(n_elements(observer_name_kw) gt 0) or (n_elements(dsun_kw) gt 0) or $
  (n_elements(lonc_kw) gt 0) or (n_elements(b0sun_kw) gt 0)
 has_cli_view_override=(n_elements(xc) gt 0) or (n_elements(yc) gt 0) or (n_elements(nx) gt 0) or (n_elements(ny) gt 0)
 saved_fov=0
 square_fov=GXObserverGeometry__saved_square_fov(observer_struct)
 prefer_saved_fov=(~has_cli_observer_override and ~has_cli_view_override)
 if keyword_set(auto_fov) then prefer_saved_fov=0b
 if n_elements(use_saved_fov) gt 0 then prefer_saved_fov=(long(use_saved_fov) ne 0L)
 if prefer_saved_fov then saved_fov=GXObserverGeometry__saved_fov(observer_struct)
 computed_fov=GXComputeInscribingFOV(model, geom, execute_text=execute_text, box_struct=box_struct, index_struct=index_struct, square_fov=square_fov)
 simbox_default=GXResolveSimboxFromObserverAndModel(saved_fov=saved_fov, computed_fov=computed_fov)
 if n_elements(simbox_default) gt 0 then begin
  xc_auto=double(simbox_default.xc_arcsec)
  yc_auto=double(simbox_default.yc_arcsec)
  fov_x_model=double(simbox_default.xsize_arcsec)
  fov_y_model=double(simbox_default.ysize_arcsec)
  center_source=strtrim(string(simbox_default.source), 2)
 endif else begin
  dsun_km=double(model.DSun)/1d5
  km_to_arcsec=206265d0/dsun_km
  fov_x_model=double(model.Nx)*double(model.dx)/1d5*km_to_arcsec
  fov_y_model=double(model.Ny)*double(model.dy)/1d5*km_to_arcsec
  xc_auto=0d & yc_auto=0d
  center_source='fallback_zero'
 endelse

 if n_elements(xc) eq 0 then xc=xc_auto
 if n_elements(yc) eq 0 then yc=yc_auto
 if n_elements(nx) eq 0 then nx=long(ceil(fov_x_model/double(dx)))
 if n_elements(ny) eq 0 then ny=long(ceil(fov_y_model/double(dy)))
 if long(nx) lt 16L then nx=16L
 if long(ny) lt 16L then ny=16L

 print, 'Observer geometry source: ', strtrim(string(geom.source), 2)
 print, 'Observer triad: L0=', double(geom.l0_deg), ' B0=', double(geom.b0_deg), ' DSun=', double(geom.dsun_cm)
 print, 'Render triad: lonC=', double(model.lonC), ' b0Sun=', double(model.b0Sun), ' DSun=', double(model.DSun)
 print, 'Center source: ', center_source
 print, 'Window:    xc=', xc, ' yc=', yc, ' dx=', dx, ' dy=', dy, ' Nx=', nx, ' Ny=', ny

 ebtel=LoadEBTEL(ebtelfile)
 if n_elements(responsefile) gt 0 then begin
  restore, responsefile
  if n_elements(gxresponse) gt 0 then response=gxresponse
  if n_elements(response) eq 0 then begin
   message, 'Response SAV must contain variable "response" or "gxresponse".', /continue
   return
  endif
 endif else begin
  response=LoadEUVresponse(model.obstime, instrument=instrument)
 endelse
 effective_instrument=strupcase(instrument)
 if tag_exist(response, 'INSTRUMENT') then effective_instrument=strtrim(string(response.instrument), 2)
 print, 'Instrument:', effective_instrument
 simbox=MakeSimulationBoxEUV(xc, yc, dx, dy, nx, ny)
 coronaparms=DefineCoronaParams(Tbase, nbase, Q0, a, b)
 outspace=ReserveOutputSpaceEUV(simbox, response)
 if n_elements(dump_dll_inputs) gt 0 then begin
  dumpfile=''
  if size(dump_dll_inputs, /type) eq 7 then dumpfile=strtrim(string(dump_dll_inputs), 2)
  if (dumpfile eq '') or (dumpfile eq '1') then dumpfile=outpath+'.dll_input.sav'
    GXDumpDLLInputs, dumpfile, model, simbox, geom=geom, response=response, ebtel=ebtel
 endif
 print, 'Elapsed time (loading): ', systime(1)-tm, ' s'

 tm=systime(1)
 r=call_external(libname, 'ComputeEUV', model, ebtel, response, simbox, coronaparms, outspace, SHtable)
 print, 'Elapsed time (main): ', systime(1)-tm, ' s'

 ConvertToMapsEUV, outspace, simbox, model, response, mapCorona, mapTR, $
  B0=double(geom.b0_deg), L0=double(geom.l0_deg), RSun=double(geom.rsun_arcsec)
 map=obj_new('map')
 nChan=mapCorona->get(/count)
 for k=0L, nChan-1L do begin
  m=mapTR.getmap(k)
  m.id='GX (TR) '+m.id
  map.setmap, k, m
 endfor
 for k=0L, nChan-1L do begin
  m=mapCorona.getmap(k)
  m.id='GX (Corona) '+m.id
  map.setmap,nChan+k, m
 endfor
 ; Preserve the historical direct-save path because some IDL map workflows
 ; restore the per-channel EUV objects correctly while a repacked wrapper may not.
 save, map, mapCorona, mapTR, filename=outpath, /compress
 print, 'Outputs:'
 print, '- saved_sav: ', outpath
 print, '- flagsAll: ', outspace.flagsAll
 print, '- flagsCorona: ', outspace.flagsCorona

 plot_idx=0L
 plot_idx_tr=0L
 found_nonzero=0b
 found_nonzero_tr=0b
 print, 'EUV channel stats:'
 for k=0L, nChan-1L do begin
  mc=mapCorona.getmap(k)
  mt=mapTR.getmap(k)
  cmin=min(mc.data, max=cmax)
  tmin=min(mt.data, max=tmax)
  print, '  [', k, '] ', mc.id, '  corona=', cmin, ' .. ', cmax, '  TR=', tmin, ' .. ', tmax
  if (~found_nonzero) && (finite(cmax)) && (cmax gt 0d) then begin
   plot_idx=k
   found_nonzero=1b
  endif
  if (~found_nonzero_tr) && (finite(tmax)) && (tmax gt 0d) then begin
   plot_idx_tr=k
   found_nonzero_tr=1b
  endif
 endfor
 
 if keyword_set(no_plot) eq 0 then begin
  window, 1, title='EUV map (corona)'
  wset, 1
  loadct, 13, /silent
  m=mapCorona.getmap(plot_idx)
  plot_map, m, cbar=1,_extra=_extra

  window, 2, title='EUV map (TR)'
  wset, 2
  loadct, 13, /silent
  if found_nonzero_tr then begin
   m=mapTR.getmap(plot_idx_tr)
   plot_map, m, cbar=1,_extra=_extra
  endif else begin
   print, 'TR plot skipped: all TR channels are zero.'
  endelse

  window, 3, title='EUV map (corona + TR)'
  wset, 3
  loadct, 13, /silent
  m=mapCorona.getmap(plot_idx)
  if found_nonzero_tr then begin
   mtr=mapTR.getmap(plot_idx_tr)
   m.data+=mtr.data
  endif
  plot_map, m, cbar=1,_extra=_extra
 endif
end
