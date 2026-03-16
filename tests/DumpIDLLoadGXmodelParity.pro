; Compile via:
;   @/Users/gelu/Library/CloudStorage/Dropbox/@Projects/@SUNCAST-ORG/gximagecomputing/examples/idl/compile_local_idl
;
; Then run:
;   DumpIDLLoadGXmodelParity
;
function dumpidlloadgxmodelparity__testdata_root
 compile_opt idl2
 root=getenv('GXRENDER_TEST_DATA_ROOT')
 if strtrim(root, 2) eq '' then root=getenv('GXIMAGECOMPUTING_TEST_DATA_ROOT')
 if strtrim(root, 2) ne '' then begin
  root=file_expand_path(root)
  if file_test(root+'/raw', /directory) and file_test(root+'/scripts', /directory) then return, root+'/raw'
  return, root
 endif

 exdir=file_dirname(routine_filepath('DumpIDLLoadGXmodelParity'))
 repodir=file_dirname(exdir)
 return, file_expand_path(repodir+'/../pyGXrender-test-data/raw')
end

function dumpidlloadgxmodelparity__models_dir, root
 compile_opt idl2
 model_dirs=file_search(root+'/models/*', /test_directory, count=ndirs)
 if ndirs gt 0 then return, model_dirs[0]
 return, root+'/models'
end

pro DumpIDLLoadGXmodelParity, TEST_DATA_DIR=test_data_dir, OUTPUT_DIR=output_dir, $
                              OVERRIDE_DSUN_CM=override_dsun_cm, $
                              OVERRIDE_LONC_DEG=override_lonc_deg, $
                              OVERRIDE_B0SUN_DEG=override_b0sun_deg

 compile_opt idl2
 resolve_routine, 'LoadGXmodel', /either

 if n_elements(test_data_dir) eq 0 then begin
  testdata_root=dumpidlloadgxmodelparity__testdata_root()
  test_data_dir=dumpidlloadgxmodelparity__models_dir(testdata_root)
 endif
 if n_elements(output_dir) eq 0 then output_dir='/tmp/gximagecomputing_loader_dumps/idl'
 if n_elements(override_dsun_cm) eq 0 then override_dsun_cm=1.4321098765d13
 if n_elements(override_lonc_deg) eq 0 then override_lonc_deg=23.456789d
 if n_elements(override_b0sun_deg) eq 0 then override_b0sun_deg=-6.54321d
 requested_override_dsun_cm=double(override_dsun_cm)
 requested_override_lonc_deg=double(override_lonc_deg)
 requested_override_b0sun_deg=double(override_b0sun_deg)

 test_data_dir=file_expand_path(test_data_dir)
 output_dir=file_expand_path(output_dir)
 file_mkdir, output_dir

 sav_files=file_search(test_data_dir+'/*chr*.sav', count=nsav)
 h5_files=file_search(test_data_dir+'/*chr*.h5', count=nh5)

 if (nsav + nh5) eq 0 then begin
  message, 'No *chr*.sav or *chr*.h5 files found under '+test_data_dir
 endif

 if nsav gt 0 then files=sav_files else files=h5_files
 if (nsav gt 0) and (nh5 gt 0) then files=[sav_files, h5_files]

 nan_value=!values.d_nan

 for i=0L, n_elements(files)-1L do begin
  infile=files[i]
  source_name=file_basename(infile)
  if strpos(strlowcase(source_name), '.sav') ne -1 then source_kind='sav' else source_kind='h5'

  case_name='default'
  override_applied=0B
  recompute_observer_ephemeris=0B
  observer_name=''
  override_dsun_cm=nan_value
  override_lonc_deg=nan_value
  override_b0sun_deg=nan_value
  model=LoadGXmodel(infile)
  outfile=output_dir+'/'+source_name+'.'+case_name+'.loaderdump.sav'
  source_file=infile
  save, source_file, source_name, source_kind, case_name, override_applied, recompute_observer_ephemeris, observer_name, $
        override_dsun_cm, override_lonc_deg, override_b0sun_deg, model, $
        filename=outfile

  case_name='recompute_earth'
  override_applied=0B
  recompute_observer_ephemeris=1B
  observer_name='earth'
  override_dsun_cm=nan_value
  override_lonc_deg=nan_value
  override_b0sun_deg=nan_value
  model=LoadGXmodel(infile, /recompute_observer_ephemeris, observer_name=observer_name)
  outfile=output_dir+'/'+source_name+'.'+case_name+'.loaderdump.sav'
  source_file=infile
  save, source_file, source_name, source_kind, case_name, override_applied, recompute_observer_ephemeris, observer_name, $
        override_dsun_cm, override_lonc_deg, override_b0sun_deg, model, $
        filename=outfile

  case_name='override'
  override_applied=1B
  recompute_observer_ephemeris=0B
  observer_name=''
  override_dsun_cm=requested_override_dsun_cm
  override_lonc_deg=requested_override_lonc_deg
  override_b0sun_deg=requested_override_b0sun_deg
  model=LoadGXmodel(infile, DSUN=override_dsun_cm, LONC=override_lonc_deg, B0SUN=override_b0sun_deg)
  outfile=output_dir+'/'+source_name+'.'+case_name+'.loaderdump.sav'
  source_file=infile
  save, source_file, source_name, source_kind, case_name, override_applied, recompute_observer_ephemeris, observer_name, $
        override_dsun_cm, override_lonc_deg, override_b0sun_deg, model, $
        filename=outfile
 endfor

 print, 'IDL loader dumps written to: '+output_dir
end
