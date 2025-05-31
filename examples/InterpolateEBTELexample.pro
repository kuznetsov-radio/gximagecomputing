pro InterpolateEBTELexample
 ;load the EBTEL data:
 restore, 'C:\ssw\packages\gx_simulator\euv\ebtel\ebtel_ss.sav'
 
 ;specify some arrays of (Q, L) values of interest:
 Qarr=[0.05d0, 0.06d0, 0.07d0, 0.08d0, 0.09d0]
 Larr=[1d8, 1.5d8, 2d8, 2.5d8, 3d8]

 ;run the code
 InterpolateEBTEL, Qrun, Lrun, logtdem, Qarr, Larr, flag, $
                   DEM_run=DEM_cor_run, DDM_run=DDM_cor_run, $
                   DEM_arr=DEM_arr, DDM_arr=DDM_arr, $
                   n_DEM=n_DEM, T_DEM=T_DEM, n_DDM=n_DDM, T_DDM=T_DDM
        
 ;plot the results:           
 window, 1, title='DEM'
 wset, 1
 plot, logtdem, DEM_arr[*, 0], yrange=[0, max(DEM_arr)], /nodata, xstyle=1
 for i=0, n_elements(Qarr)-1 do oplot, logtdem, DEM_arr[*, i]
 
 window, 2, title='DDM'
 wset, 2
 plot, logtdem, DDM_arr[*, 0], yrange=[0, max(DDM_arr)], /nodata, xstyle=1
 for i=0, n_elements(Qarr)-1 do oplot, logtdem, DDM_arr[*, i] 
 
 ;print the distribution moments:
 print, 'n_DEM=', n_DEM
 print, 'T_DEM=', T_DEM
 print, 'n_DDM=', n_DDM
 print, 'T_DDM=', T_DDM
 
 ;now run the code to compute the distribution moments only:
 InterpolateEBTEL, Qrun, Lrun, logtdem, Qarr, Larr, flag, $
                   DEM_run=DEM_cor_run, DDM_run=DDM_cor_run, $
                   n_DEM=n_DEM, T_DEM=T_DEM, n_DDM=n_DDM, T_DDM=T_DDM, /NTonly
                   
 ;print the distribution moments:
 print, 'n_DEM=', n_DEM
 print, 'T_DEM=', T_DEM
 print, 'n_DDM=', n_DDM
 print, 'T_DDM=', T_DDM                   
end