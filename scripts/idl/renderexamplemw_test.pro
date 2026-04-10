function renderexamplemw_test__testdata_root
  compile_opt idl2
  root=getenv('GXRENDER_TEST_DATA_ROOT')
  if strtrim(root, 2) eq '' then root=getenv('GXIMAGECOMPUTING_TEST_DATA_ROOT')
  if strtrim(root, 2) ne '' then begin
    root=file_expand_path(root)
    if file_test(root+'/raw', /directory) and file_test(root+'/scripts', /directory) then return, root+'/raw'
    return, root
  endif

  exdir=file_dirname(routine_filepath('RenderExampleMW_test'))
  repodir=file_dirname(file_dirname(exdir))
  return, file_expand_path(repodir+'/../pyGXrender-test-data/raw')
end

function renderexamplemw_test__find_fixture, root, category, pattern
  compile_opt idl2
  files=file_search(root+'/'+category+'/*/'+pattern, count=nfiles)
  if nfiles eq 0 then files=file_search(root+'/'+category+'/'+pattern, count=nfiles)
  if nfiles eq 0 then return, ''
  return, files[0]
end

pro RenderExampleMW_test,_extra=_extra
  testdata_root=renderexamplemw_test__testdata_root()
  modelfile=renderexamplemw_test__find_fixture(testdata_root, 'models', 'test.chr.sav')
  if strtrim(modelfile, 2) eq '' then message, 'Could not locate test.chr.sav under '+testdata_root

  ebtelfile=getenv('GXIMAGECOMPUTING_EBTEL_PATH')
  if strtrim(ebtelfile, 2) eq '' then ebtelfile=renderexamplemw_test__find_fixture(testdata_root, 'ebtel', 'ebtel.sav')
  if strtrim(ebtelfile, 2) eq '' then message, 'Could not locate ebtel.sav under '+testdata_root

  RenderExampleMW, $
    MODelfile=modelfile, $
    EBTELfile=ebtelfile, $
    OUTdir='/tmp/gximagecomputing_validation_groundtruth', $
    OUTfile='test.chr.sav_idl_mw_maps.sav',_extra=_extra
end
