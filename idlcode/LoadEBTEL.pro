function LoadEBTEL, infile, DEM=DEM, DDM=DDM
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
                      
  aDEM=1
  aDDM=1
  
  if exist(DEM_cor_run) && exist(DDM_cor_run) then begin
   if keyword_set(DEM) && ~keyword_set(DDM) then aDDM=0
   if keyword_set(DDM) && ~keyword_set(DEM) then aDEM=0
  endif
                      
  if exist(DEM_cor_run) && aDEM then begin
   ebtel.DEM_on=1
   ebtel=create_struct(ebtel, $
                       'DEM_cor_run', float(DEM_cor_run)) 
  endif       
  
  if exist(DDM_cor_run) && aDDM then begin
   ebtel.DDM_on=1
   ebtel=create_struct(ebtel, $
                       'DDM_cor_run', float(DDM_cor_run)) 
  endif          
  
  if exist(DEM_tr_run) && ebtel.DEM_on then ebtel=create_struct(ebtel, 'DEM_tr_run', float(DEM_tr_run)) 
  if exist(DDM_tr_run) && ebtel.DDM_on then ebtel=create_struct(ebtel, 'DDM_tr_run', float(DDM_tr_run)) 
 endif       
        
 return, ebtel
end