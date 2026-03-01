pro RenderExampleEUV, MODelfile=modelfile, EBTELfile=ebtelfile, RESPonsefile=responsefile, LIBname=libname, $
                    INSTRument=instrument, DSUN=dsun_kw, LONC=lonc_kw, B0SUN=b0sun_kw, $
                    XC=xc, YC=yc, DX=dx, DY=dy, NX=nx, NY=ny, $
                    OUTfile=outfile, NO_PLOT=no_plot, _extra=_extra

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
 if n_elements(outfile) eq 0 then outfile='EUVmaps.sav'

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
 print, 'Instrument:', strupcase(instrument)

 forward_function LoadGXmodel, LoadEBTEL, LoadEUVresponse, MakeSimulationBoxEUV, DefineCoronaParams, ReserveOutputSpaceEUV

 tm=systime(1)
 model=LoadGXmodel(modelfile, DSun=dsun_kw, LONC=lonc_kw, B0SUN=b0sun_kw)

 ; Derive default center/FOV from model when not explicitly provided (same logic as MW example).
 dsun_km=double(model.DSun)/1d5
 km_to_arcsec=206265d0/dsun_km
 fov_x_model=double(model.Nx)*double(model.dx)/1d5*km_to_arcsec
 fov_y_model=double(model.Ny)*double(model.dy)/1d5*km_to_arcsec

 xc_auto=0d & yc_auto=0d
 catch, err
 if err eq 0 then begin
  ; wcs_conv_hg_hpc expects dsun_obs in meters; model.DSun is in centimeters.
  wcs_conv_hg_hpc, double(model.lonC), double(model.latC), xc_auto, yc_auto, $
                  dsun_obs=double(model.DSun)/100d, b0_angle=double(model.b0Sun), l0_angle=0d, /arcseconds
  catch, /cancel
 endif else begin
  catch, /cancel
 endelse

 if n_elements(xc) eq 0 then xc=xc_auto
 if n_elements(yc) eq 0 then yc=yc_auto
 if n_elements(nx) eq 0 then nx=long(ceil(fov_x_model/double(dx)))
 if n_elements(ny) eq 0 then ny=long(ceil(fov_y_model/double(dy)))
 if long(nx) lt 16L then nx=16L
 if long(ny) lt 16L then ny=16L

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
  response=LoadEUVresponse(model, instrument=instrument)
 endelse
 simbox=MakeSimulationBoxEUV(xc, yc, dx, dy, nx, ny)
 coronaparms=DefineCoronaParams(Tbase, nbase, Q0, a, b)
 outspace=ReserveOutputSpaceEUV(simbox, response)
 print, 'Elapsed time (loading): ', systime(1)-tm, ' s'

 tm=systime(1)
 r=call_external(libname, 'ComputeEUV', model, ebtel, response, simbox, coronaparms, outspace, SHtable)
 print, 'Elapsed time (main): ', systime(1)-tm, ' s'

 ConvertToMapsEUV, outspace, simbox, model, response, mapCorona, mapTR
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
 save, map, filename=outfile, /compress
 print, 'Saved: ', outfile
 
 if keyword_set(no_plot) eq 0 then begin
  window, 1, title='EUV map (corona)'
  wset, 1
  loadct, 13, /silent
  m=mapCorona.getmap(2)
  plot_map, m, cbar=1,_extra=_extra

  window, 2, title='EUV map (TR)'
  wset, 2
  loadct, 13, /silent
  m=mapTR.getmap(2)
  plot_map, m, cbar=1,_extra=_extra

  window, 3, title='EUV map (corona + TR)'
  wset, 3
  loadct, 13, /silent
  m=mapCorona.getmap(2)
  mtr=mapTR.getmap(2)
  m.data+=mtr.data
  plot_map, m, cbar=1,_extra=_extra
 endif
end
