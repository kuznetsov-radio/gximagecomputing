pro InterpolateEBTEL, Qrun, Lrun, Qarr, Larr, flag, $
                      DEM_run=DEM_run, DDM_run=DDM_run, DEM_arr=DEM_arr, DDM_arr=DDM_arr
;The code computes the local DEM or/and DDM distributions for the specified values of the heating rate Q and 
;magnetic loop length L, using the pre-computed EBTEL tables. The bilinear interpolation is used.
;
;The code calls the external library (DLL or SO), specified by the libname variable below. Set that variable to the
;appropriate library name.  
;
;Input parameters:
; Qrun and Lrun - the original EBTEL grids in Q and L, respectively, 2D arrays, float.
; DEM_run or/and DDM_run - the original EBTEL DEM or/and DDM tables, 3D arrays, float. At least one of these
;                          tables must be specified.
; Qarr and Larr - the (Q, L) values where the DEM or/and DDM should be computed, scalars or 1D arrays
;                 (with the same sizes), double.
;Note that the parameters Qrun, Lrun, DEM_run, and DDM_run must be in the single-precision (float) format, while
;the Qarr and Larr should be in the double-precision (double) format; a type mismatch will result in a crash.
;
;Output parameters:
; DEM_arr or/and DDM_arr - the resulting DEM or/and DDM distributions, 2D arrays (NT * NP elements , where NT is 
;                          the number of the temperature grid nodes, and NP is the size of the Q and L arrays), double.
; flag - 1D array (with NP elements), byte. An element of this array equals 1, if the interpolation for the 
;        corresponding (Q, L) pair was successful, and 0, if that point was beyond the EBTEL grid.
                      
 libname=gx_libpath('rendergrff')
 
 DEM_on=exist(DEM_run)
 DDM_on=exist(DDM_run) 
                      
 if DEM_on then s=size(DEM_run, /dimensions) $
 else if DDM_on then s=size(DDM_run, /dimensions) $
 else begin
  print, 'Either DEM_run or DDM_run must be specified!'
  stop
 endelse    
 
 NT=s[0]
 NQ=s[1]
 NL=s[2]
 NP=n_elements(Qarr)                 

 Lparms=long([NP, NQ, NL, NT, DEM_on, DDM_on])
 
 flag=bytarr(NP)
 
 if DEM_on then DEM_arr=dblarr(NT, NP)
 if DDM_on then DDM_arr=dblarr(NT, NP)
 
 res=call_external(libname, $
                   'InterpolateEBTEL', $
                   Lparms, $
                   Qrun, $
                   Lrun, $
                   DEM_on ? DEM_run : 0, $
                   DDM_on ? DDM_run : 0, $
                   Qarr, $
                   Larr, $
                   DEM_on ? DEM_arr : 0, $
                   DDM_on ? DDM_arr : 0, $
                   flag)
end