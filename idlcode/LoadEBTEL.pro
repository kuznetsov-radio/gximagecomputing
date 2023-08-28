function LoadEBTEL, infile
 ebtel={DEM_on: 0L, $
        DDM_on: 0L}
        
 if infile ne '' then begin
  restore, infile
  
  s=size(Lrun, /dimensions)
  NQ=s[0]
  NL=s[1]
  NT=n_elements(logtdem)
  
  ebtel=create_struct(ebtel, $
                      'NQ', long(NQ), $
                      'NL', long(NL), $
                      'NT', long(NT), $
                      'Qrun', float(Qrun), $
                      'Lrun', float(Lrun), $
                      'logtdem', float(logtdem))
                      
  if exist(DEM_cor_run) then begin
   ebtel.DEM_on=1
   ebtel=create_struct(ebtel, $
                       'DEM_cor_run', float(DEM_cor_run)) 
  endif       
  
  if exist(DDM_cor_run) then begin
   ebtel.DDM_on=1
   ebtel=create_struct(ebtel, $
                       'DDM_cor_run', float(DDM_cor_run)) 
  endif               
 endif       
        
 return, ebtel
end
