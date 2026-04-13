pro FillVolumeNT, ModelFileName, EBTELFileName, Q0, a, b, n_base, T_base, $
                  x, y, z, Bx, By, Bz, n_DEM, T_DEM, n_DDM, T_DDM, flag
;This code computes the thermal plasma density and temperature in an active region for a 
;specified heating model.
;
;Input parameters:
; ModelFileName - name of the .sav file that contains the GX Simulator model.
; EBTELFileName - name of the .sav file that contains the EBTEL tables.
;                 if EBTELFileName eq 'Analytical' then analytical formulae by J. Klimchuk are used
;                 to compute the plasma density and temperature.
; Q0, a, b - parameters of the heating model (same as in CHMP).
; n_base, T_base - plasma density at the base of the corona (in cm^{-3}) and constant temperature (in K)
;                  specifying the plasma parameters in the voxels where the heating model fails.
;                  The temperature there equals T_base, and the density is computed using the barometric
;                  formula, n = n_base * exp(-z / T_base / 4.73d3).
;                  
;Output parameters:
; x, y, z - 1D arrays specifying the voxel coordinates (in cm). Coordinates correspond to the centers of
;           the voxels.
; Bx, By, Bz - 3D arrays (Nx * Ny * Nz) containing the x, y, z components of the magnetic field, in G.
; n_DEM, T_DEM - 3D arrays (Nx * Ny * Nz) containing the plasma density (in cm^{-3}) and temperature (in K)
;                computed as the moments of the DEM distribution.
; n_DDM, T_DDM - 3D arrays (Nx * Ny * Nz) containing the plasma density (in cm^{-3}) and temperature (in K)
;                computed as the moments of the DDM distribution.
; flag - 3D integer array (Nx * Ny * Nz) indicating the computation status for each voxel:
;           -1: the voxel corresponds to an open field line, the plasma density and temperature are chosen
;               according to the isothermal barometric model;
;            0: the voxel corresponds to a closed field line, the plasma density and temperature are computed
;               as the moments of the DEM/DDM distributions;
;            1: the voxel corresponds to a closed field line, but the line parameters Q and/or L are beyond
;               the EBTEL DEM/DDM table; the plasma density and temperature are chosen according to the 
;               isothermal barometric model.
;        If EBTELFileName eq 'Analytical', then the analytical formulae are applied to all closed field lines
;        (flag le 0), and n_DEM eq n_DDM and T_DEM eq T_DDM.
                  
 restore, ModelFileName
 
 Bx=reform(box.bcube[*, *, *, 0])
 By=reform(box.bcube[*, *, *, 1])
 Bz=reform(box.bcube[*, *, *, 2])
 
 s=size(Bx, /dimensions)
 Nx=s[0]
 Ny=s[1]
 Nz=s[2]
 
 x=(dindgen(Nx)-0.5*Nx+0.5)*box.dr[0]*696000d5
 y=(dindgen(Ny)-0.5*Ny+0.5)*box.dr[1]*696000d5
 z=(dindgen(Nz)+0.5)*box.dr[2]*696000d5
 
 nbar=n_base*exp(-z/(T_base*4.73d3))
 n_DEM=dblarr(Nx, Ny, Nz)
 n_DDM=dblarr(Nx, Ny, Nz)
 for i=0, Nx-1 do for j=0, Ny-1 do begin
  n_DEM[i, j, *]=nbar
  n_DDM[i, j, *]=nbar
 endfor
 T_DEM=make_array(Nx, Ny, Nz, /double, value=T_base)
 T_DDM=make_array(Nx, Ny, Nz, /double, value=T_base)
 
 Bavg=fltarr(Nx, Ny, Nz)
 Lline=fltarr(Nx, Ny, Nz)
 
 if tag_exist(box, 'BMED') then begin ;old version
  Bavg[box.idx]=box.bmed
  Lline[box.idx]=box.length*696000e5
 endif else begin ;new version
  u=where((box.status and 4) eq 4, k)
  if k gt 0 then begin
   Bavg[u]=box.avfield[u]
   Lline[u]=box.physlength[u]*696000e5
  endif
 endelse 
 
 flag=intarr(Nx, Ny, Nz)
 
 u=where((Lline gt 0) and (Bavg gt 0), k, complement=nu, ncomplement=nk)
 
 if nk gt 0 then flag[nu]=-1
 
 if k gt 0 then begin
  Q=Q0*(Bavg[u]/100d0)^a/(Lline[u]/2/1d9)^b
  L=Lline[u]/2
  
  if EBTELFileName ne 'Analytical' then begin
   restore, EBTELFileName
   libname=gx_libpath('rendergrff')
   InterpolateEBTEL, Qrun, Lrun, logtdem, double(Q), double(L), ieflag, /NTonly, $
                     DEM_run=DEM_cor_run, DDM_run=DDM_cor_run, n_DEM=n_DEM1, T_DEM=T_DEM1, n_DDM=n_DDM1, T_DDM=T_DDM1
   u1=where(ieflag, k1, complement=nu1, ncomplement=nk1)
   if nk1 gt 0 then flag[u[nu1]]=1
   if k1 gt 0 then begin
    n_DEM[u[u1]]=n_DEM1[u1]
    T_DEM[u[u1]]=T_DEM1[u1]
    n_DDM[u[u1]]=n_DDM1[u1]
    T_DDM[u[u1]]=T_DDM1[u1]
   endif
  endif else begin
   KlimchukNT, Q, L, n_ana, T_ana
   n_DEM[u]=n_ana
   T_DEM[u]=T_ana
   n_DDM[u]=n_ana
   T_DDM[u]=T_ana
  endelse
 endif
end

pro KlimchukNT, Q, L, n, T
  T=(Q*L^2/4.3d-7)^(2d0/7)

  l0=dblarr(n_elements(T))
  b=dblarr(n_elements(T))

  u=where(alog10(T) le 4.97, k)
  if k gt 0 then begin
    l0[u]=1.09d-31
    b[u]=2
  endif

  u=where((alog10(T) gt 4.97) and (alog10(T) le 5.67), k)
  if k gt 0 then begin
    l0[u]=8.87d-17
    b[u]=-1
  endif

  u=where((alog10(T) gt 5.67) and (alog10(T) le 6.18), k)
  if k gt 0 then begin
    l0[u]=1.90d-22
    b[u]=0
  endif

  u=where((alog10(T) gt 6.18) and (alog10(T) le 6.55), k)
  if k gt 0 then begin
    l0[u]=3.53d-13
    b[u]=-1.5
  endif

  u=where((alog10(T) gt 6.55) and (alog10(T) le 6.90), k)
  if k gt 0 then begin
    l0[u]=3.46d-25
    b[u]=1d0/3
  endif

  u=where((alog10(T) gt 6.90) and (alog10(T) le 7.63), k)
  if k gt 0 then begin
    l0[u]=5.49d-16
    b[u]=-1
  endif

  u=where(alog10(T) gt 7.63, k)
  if k gt 0 then begin
    l0[u]=1.96d-27
    b[u]=0.5
  endif

  n=sqrt(Q/3/l0/T^b)
end