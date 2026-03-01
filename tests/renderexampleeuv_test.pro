pro RenderExampleEUV_test,_extra=_extra
  RenderExampleEUV, $
    MODelfile='/Users/gelu/Library/CloudStorage/Dropbox/@Projects/@SUNCAST-ORG/gximagecomputing/test_data/test.chr.sav', $
    EBTELfile=gx_findfile('ebtel.sav'), $   
    OUTfile='/tmp/gximagecomputing_validation_groundtruth/test.chr.sav_idl_euv_maps.sav',_extra=_extra
end
