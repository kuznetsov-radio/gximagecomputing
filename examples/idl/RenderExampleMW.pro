pro RenderExampleMW, MODelfile=modelfile, EBTELfile=ebtelfile, LIBname=libname, $
                   IDLCODEdir=idlcodedir, $
                   OUTdir=outdir, OUTfile=outfile, $
                   DSUN=dsun_kw, LONC=lonc_kw, B0SUN=b0sun_kw, $
                   XC=xc, YC=yc, DX=dx, DY=dy, NX=nx, NY=ny, XRANGE=xrange, YRANGE=yrange, $
                   NO_PLOT=no_plot,_extra=_extra

 ; Temporary validation workflow:
 ; - Accept CHR SAV or H5 path directly.
 ; - Save output in same tmp folder used by Python comparisons by default.
 exdir=file_dirname(routine_filepath('RenderExampleMW'))
 repo_root=file_dirname(file_dirname(exdir))
 local_idlcodedir=repo_root+'/idlcode'
 local_bindir=repo_root+'/binaries'
 env_ebtel=getenv('GXIMAGECOMPUTING_EBTEL_PATH')

 if n_elements(modelfile) eq 0 then begin
   message, 'Please provide MODelfile=<path to CHR .sav or .h5 model>.', /continue
   return
 endif
 if n_elements(ebtelfile) eq 0 then begin
   if n_elements(env_ebtel) gt 0 then begin
     ebtelfile=env_ebtel
   endif else begin
     message, 'Please provide EBTELfile=<path> or set GXIMAGECOMPUTING_EBTEL_PATH.', /continue
     return
   endelse
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
   endif else if file_test(local_bindir+'/RenderGRFF_32.dll') then begin
     libname=local_bindir+'/RenderGRFF_32.dll'
   endif else begin
     message, 'No RenderGRFF library found in '+local_bindir+'. Please provide LIBname=<path>.', /continue
     return
   endelse
 endif
 if n_elements(idlcodedir) eq 0 then idlcodedir=local_idlcodedir
 if n_elements(outdir) eq 0 then begin
   if strpos(strlowcase(!version.os_family), 'windows') ge 0 then begin
     outdir='C:/Temp/gximagecomputing_validation_groundtruth'
   endif else begin
     outdir='/tmp/gximagecomputing_validation_groundtruth'
   endelse
 endif
 if n_elements(outfile) eq 0 then outfile=file_basename(modelfile)+'_idl_mw_maps.sav'

 ; Ensure required RenderGRFF IDL helpers are on IDL path.
 ; Put local repository idlcode first so updated routines win over SSW copies.
 if file_test(local_idlcodedir, /directory) then begin
   if strpos(':'+!path+':', ':'+local_idlcodedir+':') lt 0 then !path=local_idlcodedir+':'+!path
 endif
 if file_test(idlcodedir, /directory) then begin
   if strpos(':'+!path+':', ':'+idlcodedir+':') lt 0 then !path=!path+':'+idlcodedir
 endif else begin
   message, 'IDL helper folder not found: '+idlcodedir, /continue
 endelse

 ; Use same 16 frequencies as Python (every other from dense 5.8..12.0, step 0.2).
 freqlist=[5.8, 6.2, 6.6, 7.0, 7.4, 7.8, 8.2, 8.6, 9.0, 9.4, 9.8, 10.2, 10.6, 11.0, 11.4, 11.8]

 ; Default pixel scale.
 if n_elements(dx) eq 0 then dx=2.0d
 if n_elements(dy) eq 0 then dy=2.0d

 ; Analytical corona parameters.
 Tbase=1d6
 nbase=1d8
 Q0=0.0217d
 a=0.3d
 b=2.7d

 ; Selective heating table.
 SHtable=dblarr(7, 7)
 w=[1.0d, 1.0d, 1.0d, 1.1d, 1.2d, 1.3d, 1.4d]
 for i=0, 6 do for j=0, 6 do SHtable[i, j]=w[i]*w[j]
 SHtable[6, 6]=0.1d

 spawn, 'mkdir -p "'+outdir+'"'
 if strmid(outdir, strlen(outdir)-1, 1) eq '/' then begin
   outpath=outdir+outfile
 endif else begin
   outpath=outdir+'/'+outfile
 endelse

 print, 'Model:     ', modelfile
 print, 'EBTEL:     ', ebtelfile
 if n_elements(dsun_kw) gt 0 then print, 'DSun override (cm): ', dsun_kw
 if n_elements(lonc_kw) gt 0 then print, 'lonC override (deg): ', lonc_kw
 if n_elements(b0sun_kw) gt 0 then print, 'b0Sun override (deg): ', b0sun_kw
 print, 'Library:   ', libname
 print, 'Output:    ', outpath

 model_lower=strlowcase(string(modelfile))
 nmlen=strlen(model_lower)
 model_ext3=''
 model_ext4=''
 model_ext5=''
 if nmlen ge 3 then model_ext3=strmid(model_lower, nmlen-3, 3)
 if nmlen ge 4 then model_ext4=strmid(model_lower, nmlen-4, 4)
 if nmlen ge 5 then model_ext5=strmid(model_lower, nmlen-5, 5)
 if (model_ext4 eq '.sav') or (model_ext4 eq '.xdr') then begin
   print, 'Model load mode: direct restore (.sav/.xdr) via LoadGXmodel'
 endif else if (model_ext3 eq '.h5') or (model_ext5 eq '.hdf5') then begin
   print, 'Model load mode: HDF5 conversion (.h5/.hdf5) via ConvertToGX -> LoadGXmodel'
 endif else begin
   print, 'Model load mode: unknown extension; LoadGXmodel will decide'
 endelse

 forward_function LoadGXmodel, LoadEBTEL, MakeSimulationBox, DefineCoronaParams, ReserveOutputSpace
 resolve_routine, 'LoadGXmodel', /either
 resolve_routine, 'LoadEBTEL', /either
 resolve_routine, 'MakeSimulationBox', /either
 resolve_routine, 'DefineCoronaParams', /either
 resolve_routine, 'ReserveOutputSpace', /either
 resolve_routine, 'ConvertToMaps', /either

 tm=systime(1)
 model=LoadGXmodel(modelfile, DSun=dsun_kw, LONC=lonc_kw, B0SUN=b0sun_kw)
 ebtel=LoadEBTEL(ebtelfile)

 ; Derive default center/FOV from model when not explicitly provided.
 dsun_km=double(model.DSun)/1d5
 km_to_arcsec=206265d0/dsun_km
 fov_x_model=double(model.Nx)*double(model.dx)/1d5*km_to_arcsec
 fov_y_model=double(model.Ny)*double(model.dy)/1d5*km_to_arcsec

 xc_auto=0d & yc_auto=0d
 catch, err
 if err eq 0 then begin
  wcs_conv_hg_hpc, double(model.lonC), double(model.latC), xc_auto, yc_auto, $
                  dsun_obs=double(model.DSun)/100, b0_angle=double(model.b0Sun), l0_angle=0d, /arcseconds
  catch, /cancel
 endif else begin
  catch, /cancel
 endelse

 if n_elements(xrange) gt 0 then begin
  if n_elements(xrange) ne 2 then message, 'XRANGE must have exactly 2 elements [xmin, xmax]'
  if double(xrange[1]) le double(xrange[0]) then message, 'XRANGE must satisfy xmax > xmin'
  xc=0.5d*(double(xrange[0])+double(xrange[1]))
  nx=long(ceil((double(xrange[1])-double(xrange[0]))/double(dx)))
 endif else begin
  if n_elements(xc) eq 0 then xc=xc_auto
  if n_elements(nx) eq 0 then nx=long(ceil(fov_x_model/double(dx)))
 endelse

 if n_elements(yrange) gt 0 then begin
  if n_elements(yrange) ne 2 then message, 'YRANGE must have exactly 2 elements [ymin, ymax]'
  if double(yrange[1]) le double(yrange[0]) then message, 'YRANGE must satisfy ymax > ymin'
  yc=0.5d*(double(yrange[0])+double(yrange[1]))
  ny=long(ceil((double(yrange[1])-double(yrange[0]))/double(dy)))
 endif else begin
  if n_elements(yc) eq 0 then yc=yc_auto
  if n_elements(ny) eq 0 then ny=long(ceil(fov_y_model/double(dy)))
 endelse

 if long(nx) lt 16L then nx=16L
 if long(ny) lt 16L then ny=16L
 if n_elements(xc) ne 1 then message, 'XC could not be resolved to a scalar'
 if n_elements(yc) ne 1 then message, 'YC could not be resolved to a scalar'

 fov_x=double(nx)*double(dx)
 fov_y=double(ny)*double(dy)
 print, 'Window: xc=', xc, ' yc=', yc, ' dx=', dx, ' dy=', dy, ' Nx=', nx, ' Ny=', ny
 print, 'FOV: ', fov_x, ' x ', fov_y, ' arcsec'

 simbox=MakeSimulationBox(xc, yc, dx, dy, nx, ny, freqlist)
 coronaparms=DefineCoronaParams(Tbase, nbase, Q0, a, b)
 outspace=ReserveOutputSpace(simbox)
 print, 'Elapsed time (loading): ', systime(1)-tm, ' s'

 tm=systime(1)
 r=call_external(libname, 'ComputeMW', model, ebtel, simbox, coronaparms, outspace, SHtable)
 print, 'Elapsed time (main): ', systime(1)-tm, ' s'

 ConvertToMaps, outspace, simbox, model, mapI, mapV
 map=obj_new('map')
 for k=0L, n_elements(freqlist)-1L do map.setmap, k, mapI.getmap(k)
 for k=0L, n_elements(freqlist)-1L do map.setmap, n_elements(freqlist)+k, mapV.getmap(k)

 if keyword_set(no_plot) eq 0 then begin
   window, 1, title='IDL MW I'
   wset, 1
   loadct, 13, /silent
   m=map.getmap(0)
   plot_map, m, cbar=1,_extra=_extra

   window, 2, title='IDL MW V'
   wset, 2
   loadct, 33, /silent
   m=map.getmap(n_elements(freqlist))
   plot_map, m, cbar=1,_extra=_extra
 endif

 ; Save only output maps in a single container: first I maps, then V maps.
 save, map, filename=outpath, /compress
 nf=n_elements(freqlist)
 print, 'Outputs:'
 print, '- saved_sav: ', outpath
 print, '- outdir: ', outdir
 print, '- map_dims: Nx=', nx, ' Ny=', ny
 print, '- freqlist_n: ', nf
 if nf gt 0 then print, '- freqlist_range_GHz: ', freqlist[0], ' .. ', freqlist[nf-1]
end
