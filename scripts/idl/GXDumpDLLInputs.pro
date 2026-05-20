pro GXDumpDLLInputs, dumpfile, model, simbox, geom=geom, response=response, ebtel=ebtel
 if n_elements(dumpfile) eq 0 then message, 'Please provide DUMPFILE.'

 dump_model=model
 dump_simbox=simbox
 if n_elements(geom) gt 0 then dump_geom=geom else dump_geom=0
 if n_elements(response) gt 0 then dump_response=response else dump_response=0
 if n_elements(ebtel) gt 0 then dump_ebtel=ebtel else dump_ebtel=0

 voxel_count=long(n_elements(model.VoxelID))
 voxel_bits=ulong(model.VoxelID)
 chromo_count=0L
 tr_count=0L
 corona_count=0L
 euvtr_count=0L
 if voxel_count gt 0 then begin
  chromo_count=long(total((voxel_bits and 1UL) ne 0))
  tr_count=long(total((voxel_bits and 2UL) ne 0))
  corona_count=long(total((voxel_bits and 4UL) ne 0))
  euvtr_count=long(total((voxel_bits and 8UL) ne 0))
 endif

 ebtel_dem_on=0L
 ebtel_ddm_on=0L
 ebtel_has_dem_tr=0L
 ebtel_has_ddm_tr=0L
 if n_elements(ebtel) gt 0 then begin
  if tag_exist(ebtel, 'DEM_ON') then ebtel_dem_on=long(ebtel.DEM_on)
  if tag_exist(ebtel, 'DDM_ON') then ebtel_ddm_on=long(ebtel.DDM_on)
  ebtel_has_dem_tr=long(tag_exist(ebtel, 'DEM_TR_RUN'))
  ebtel_has_ddm_tr=long(tag_exist(ebtel, 'DDM_TR_RUN'))
 endif

 response_ds=!values.d_nan
 response_nt=0L
 response_nchannels=0L
 response_instrument=''
 if n_elements(response) gt 0 then begin
  if tag_exist(response, 'DS') then response_ds=double(response.ds)
  if tag_exist(response, 'NT') then response_nt=long(response.NT)
  if tag_exist(response, 'NCHANNELS') then response_nchannels=long(response.Nchannels)
  if tag_exist(response, 'INSTRUMENT') then response_instrument=strtrim(string(response.instrument), 2)
 endif

 summary={model_nx: long(model.Nx), $
          model_ny: long(model.Ny), $
          model_nz: long(model.Nz), $
          model_lonc: double(model.lonC), $
          model_b0sun: double(model.b0Sun), $
          model_dsun: double(model.DSun), $
          model_dx_cm: double(model.dx), $
          model_dy_cm: double(model.dy), $
          model_dz_uniform_cm: double(model.dz_uniform), $
          model_dz_min_cm: double(min(model.dz, max=dz_max)), $
          model_dz_max_cm: double(dz_max), $
          simbox_nx: long(simbox.Nx), $
          simbox_ny: long(simbox.Ny), $
          simbox_xc: double(simbox.xc), $
          simbox_yc: double(simbox.yc), $
          simbox_dx: double(simbox.dx), $
          simbox_dy: double(simbox.dy), $
          simbox_projection: long(simbox.projection), $
          voxel_count: long(voxel_count), $
          voxel_chromo_count: long(chromo_count), $
          voxel_tr_count: long(tr_count), $
          voxel_corona_count: long(corona_count), $
          voxel_euvtr_count: long(euvtr_count), $
          has_response: long(n_elements(response) gt 0), $
          response_ds: double(response_ds), $
          response_nt: long(response_nt), $
          response_nchannels: long(response_nchannels), $
          ebtel_dem_on: long(ebtel_dem_on), $
          ebtel_ddm_on: long(ebtel_ddm_on), $
          ebtel_has_dem_tr: long(ebtel_has_dem_tr), $
          ebtel_has_ddm_tr: long(ebtel_has_ddm_tr)}

 if n_elements(response) gt 0 then begin
  summary.has_response=1L
 endif

 save, dump_model, dump_simbox, dump_geom, dump_response, dump_ebtel, summary, filename=dumpfile, /compress
 print, 'Saved DLL input dump: ', dumpfile
 print, '  model dims: ', summary.model_nx, ' x ', summary.model_ny, ' x ', summary.model_nz
 print, '  model render: lonC=', summary.model_lonc, ' b0Sun=', summary.model_b0sun, ' DSun=', summary.model_dsun
 print, '  model dz(cm): uniform=', summary.model_dz_uniform_cm, ' min=', summary.model_dz_min_cm, ' max=', summary.model_dz_max_cm
 print, '  simbox: xc=', summary.simbox_xc, ' yc=', summary.simbox_yc, ' dx=', summary.simbox_dx, ' dy=', summary.simbox_dy, ' Nx=', summary.simbox_nx, ' Ny=', summary.simbox_ny, ' projection=', summary.simbox_projection
 print, '  voxel bits: chromo=', summary.voxel_chromo_count, ' tr=', summary.voxel_tr_count, ' corona=', summary.voxel_corona_count, ' euvtr=', summary.voxel_euvtr_count
 print, '  EBTEL: DEM_on=', summary.ebtel_dem_on, ' DDM_on=', summary.ebtel_ddm_on, ' DEM_tr=', summary.ebtel_has_dem_tr, ' DDM_tr=', summary.ebtel_has_ddm_tr
 if n_elements(response) gt 0 then print, '  response: instrument=', response_instrument, ' channels=', summary.response_nchannels, ' NT=', summary.response_nt, ' ds=', summary.response_ds
end