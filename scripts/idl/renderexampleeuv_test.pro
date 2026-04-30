function renderexampleeuv_test__testdata_root
  compile_opt idl2
  root=getenv('GXRENDER_TEST_DATA_ROOT')
  if strtrim(root, 2) eq '' then root=getenv('GXIMAGECOMPUTING_TEST_DATA_ROOT')
  if strtrim(root, 2) ne '' then begin
    root=file_expand_path(root)
    if file_test(root+'/raw', /directory) and file_test(root+'/scripts', /directory) then return, root+'/raw'
    return, root
  endif

  exdir=file_dirname(routine_filepath('RenderExampleEUV_test'))
  repodir=file_dirname(file_dirname(exdir))
  return, file_expand_path(repodir+'/../pyGXrender-test-data/raw')
end

function renderexampleeuv_test__find_fixture, root, category, pattern
  compile_opt idl2
  files=file_search(root+'/'+category+'/*/'+pattern, count=nfiles)
  if nfiles eq 0 then files=file_search(root+'/'+category+'/'+pattern, count=nfiles)
  if nfiles eq 0 then return, ''
  return, files[0]
end

pro RenderExampleEUV_test,_extra=_extra
  testdata_root=renderexampleeuv_test__testdata_root()
  modelfile=getenv('GXIMAGECOMPUTING_IDL_MODEL_PATH')
  if strtrim(modelfile, 2) eq '' then modelfile=renderexampleeuv_test__find_fixture(testdata_root, 'models', '*.sav')
  if strtrim(modelfile, 2) eq '' then message, 'Could not locate an IDL SAV model fixture under '+testdata_root+'/models. The current pyGXrender-test-data default model bundle contains H5 models only; set GXIMAGECOMPUTING_IDL_MODEL_PATH to a SAV model for IDL examples.'

  ebtelfile=getenv('GXIMAGECOMPUTING_EBTEL_PATH')
  if strtrim(ebtelfile, 2) eq '' then ebtelfile=renderexampleeuv_test__find_fixture(testdata_root, 'ebtel', 'ebtel.sav')
  if strtrim(ebtelfile, 2) eq '' then message, 'Could not locate ebtel.sav under '+testdata_root
  instrument='aia'

  RenderExampleEUV, $
    MODelfile=modelfile, $
    EBTELfile=ebtelfile, $
    INSTRument=instrument, $
    OUTfile='/tmp/gximagecomputing_validation_groundtruth/'+file_basename(modelfile)+'_idl_euv_maps.sav',_extra=_extra
end
