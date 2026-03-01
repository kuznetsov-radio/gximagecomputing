pro DumpComputeEUVInputs, MODelfile=modelfile, EBTELfile=ebtelfile, RESPonsefile=responsefile, $
                         DSUN=dsun_kw, LONC=lonc_kw, B0SUN=b0sun_kw, $
                         INSTRument=instrument, XC=xc, YC=yc, DX=dx, DY=dy, NX=nx, NY=ny, $
                         OUTfile=outfile
  compile_opt idl2

  exdir=file_dirname(routine_filepath('DumpComputeEUVInputs'))
  repodir=file_dirname(exdir)
  local_idlcodedir=repodir+'/idlcode'
  env_ebtel=getenv('GXIMAGECOMPUTING_EBTEL_PATH')

  if n_elements(modelfile) eq 0 then begin
    modelfile=repodir+'/test_data/test.chr.sav'
  endif
  if n_elements(ebtelfile) eq 0 then begin
    if n_elements(env_ebtel) gt 0 then ebtelfile=env_ebtel else ebtelfile=''
  endif
  if n_elements(instrument) eq 0 then instrument='aia'
  if n_elements(dx) eq 0 then dx=2.0
  if n_elements(dy) eq 0 then dy=2.0
  if n_elements(outfile) eq 0 then outfile='/tmp/gximagecomputing_validation_groundtruth/computeeuv_inputs_idl.sav'

  ; Ensure local helper routines are preferred.
  if file_test(local_idlcodedir, /directory) then begin
    if strpos(':'+!path+':', ':'+local_idlcodedir+':') lt 0 then !path=local_idlcodedir+':'+!path
  endif

  ; Analytical corona/heating defaults (must match RenderExampleEUV.pro).
  Tbase=1d6
  nbase=1d8
  Q0=0.0217
  a=0.3
  b=2.7
  SHtable=dblarr(7, 7)
  w=[1.0d, 1.0d, 1.0d, 1.1d, 1.2d, 1.3d, 1.4d]
  for i=0, 6 do for j=0, 6 do SHtable[i, j]=w[i]*w[j]
  SHtable[6, 6]=0.1d

  forward_function LoadGXmodel, LoadEBTEL, LoadEUVresponse, MakeSimulationBoxEUV, DefineCoronaParams, ReserveOutputSpaceEUV

  model=LoadGXmodel(modelfile, DSun=dsun_kw, LONC=lonc_kw, B0SUN=b0sun_kw)

  ; Derive defaults exactly as in RenderExampleEUV.pro.
  dsun_km=double(model.DSun)/1d5
  km_to_arcsec=206265d0/dsun_km
  fov_x_model=double(model.Nx)*double(model.dx)/1d5*km_to_arcsec
  fov_y_model=double(model.Ny)*double(model.dy)/1d5*km_to_arcsec

  xc_auto=0d & yc_auto=0d
  catch, err
  if err eq 0 then begin
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

  debugmeta={ $
    modelfile: string(modelfile), $
    ebtelfile: string(ebtelfile), $
    responsefile: (n_elements(responsefile) gt 0 ? string(responsefile) : ''), $
    dsun_override: (n_elements(dsun_kw) gt 0 ? double(dsun_kw) : !values.d_nan), $
    lonc_override: (n_elements(lonc_kw) gt 0 ? double(lonc_kw) : !values.d_nan), $
    b0sun_override: (n_elements(b0sun_kw) gt 0 ? double(b0sun_kw) : !values.d_nan), $
    instrument: strupcase(string(instrument)), $
    xc: double(xc), yc: double(yc), dx: double(dx), dy: double(dy), $
    nx: long(nx), ny: long(ny) $
  }

  file_mkdir, file_dirname(outfile)
  save, model, ebtel, response, simbox, coronaparms, outspace, SHtable, debugmeta, filename=outfile, /compress
  print, 'Saved ComputeEUV pre-call inputs: ', outfile
end
