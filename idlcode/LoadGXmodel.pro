function LoadGXmodel, infile, noVoxelID=noVoxelID, newTime=newTime
 restore, infile
 
 obstime=anytim(box.index.date_obs)
 
 setenv, 'WCS_RSUN=6.96d8' ;the same command as in GX Simulator routines
 wcs=fitshead2wcs(box.index)
 wcs_convert_from_coord, wcs, wcs.crval, 'hg', lonC, latC
 RSun=wcs_rsun(unit='cm')
 
 if exist(newTime) then begin
  obstime1=anytim(newTime)
  ddays=(obstime1-obstime)/86400
  lonC+=diff_rot(ddays, latC, /synodic)
  obstime=obstime1
 endif
 
 a=get_sun(anytim(obsTime, /ex))
 DSun=a[0]*1.495978707d13 
 b0Sun=a[11]
 
 dx=box.dr[0]*RSun
 dy=box.dr[1]*RSun
 dz_uniform=box.dr[2]*RSun
 dz=box.dz*RSun
 
 s=size(box.dz, /dimensions)
 sc=size(box.bcube, /dimensions)
 
 Nx=s[0]
 Ny=s[1]
 Nz=s[2]
 
 chromo_layers=box.chromo_layers
 corona_layers=s[2]-box.chromo_layers
 corona_base=box.corona_base
 
 Bx=fltarr(s)
 By=fltarr(s)
 Bz=fltarr(s)
 
 Bx[*, *, 0 : chromo_layers-1]=box.chromo_bcube[*, *, *, 0]
 By[*, *, 0 : chromo_layers-1]=box.chromo_bcube[*, *, *, 1]
 Bz[*, *, 0 : chromo_layers-1]=box.chromo_bcube[*, *, *, 2]
 
 Bx[*, *, chromo_layers : s[2]-1]=box.bcube[*, *, box.corona_base : sc[2]-1, 0]
 By[*, *, chromo_layers : s[2]-1]=box.bcube[*, *, box.corona_base : sc[2]-1, 1]
 Bz[*, *, chromo_layers : s[2]-1]=box.bcube[*, *, box.corona_base : sc[2]-1, 2]
 
 chromo_n0=fltarr(s[0], s[1], chromo_layers)
 chromo_np=fltarr(s[0], s[1], chromo_layers)
 chromo_nHI=fltarr(s[0], s[1], chromo_layers)
 chromo_T0=fltarr(s[0], s[1], chromo_layers)
 
 chromo_n0[box.chromo_idx]=box.chromo_n
 chromo_np[box.chromo_idx]=box.n_p
 chromo_nHI[box.chromo_idx]=box.n_hi
 chromo_T0[box.chromo_idx]=box.chromo_T
 
 if tag_exist(box, 'BMED') then begin ;old version
  QB=fltarr(s[0], s[1], sc[2])
  QB[box.idx]=box.bmed
  QL=fltarr(s[0], s[1], sc[2])
  QL[box.idx]=box.length*RSun
  ID1=bytarr(s[0], s[1], s[2])
  ID1[box.idx]=box.base.chromo_mask[box.foot1]
  ID2=bytarr(s[0], s[1], s[2])
  ID2[box.idx]=box.base.chromo_mask[box.foot2]  
 endif else if tag_exist(box, 'AVFIELD') then begin ;new version
  QB=box.avfield
  QL=box.physlength*RSun
  ID1=byte(box.base.chromo_mask[box.startidx])
  ID2=byte(box.base.chromo_mask[box.endidx])
  u=where((box.status and 4) ne 4, k)
  if k gt 0 then begin
   QB[u]=0
   QL[u]=0
   ID1[u]=0
   ID2[u]=0
  endif
 endif 
 
 corona_Bavg=QB[*, *, box.corona_base : sc[2]-1]
 chromo_uniform_Bavg=QB[*, *, 0 : box.corona_base-1]
 
 corona_L=QL[*, *, box.corona_base : sc[2]-1]
 chromo_uniform_L=QL[*, *, 0 : box.corona_base-1]
 
 if ~keyword_set(NoVoxelID) then VoxelID=gx_box2id(box) $
 else VoxelID=bytarr(s)
 
 corona_ID1=ID1[*, *, box.corona_base : sc[2]-1]
 corona_ID2=ID2[*, *, box.corona_base : sc[2]-1]
 chromo_uniform_ID1=ID1[*, *, 0 : box.corona_base-1]
 chromo_uniform_ID2=ID2[*, *, 0 : box.corona_base-1]
 
 model={Nx: long(Nx), $
        Ny: long(Ny), $
        Nz: long(Nz), $
        chromo_layers: long(chromo_layers), $
        corona_layers: long(corona_layers), $
        corona_base: long(corona_base), $
        DSun: double(DSun), $
        RSun: double(Rsun), $
        b0Sun: double(b0Sun), $
        lonC: double(lonC), $
        latC: double(latC), $
        dx: double(dx), $
        dy: double(dy), $
        dz_uniform: double(dz_uniform), $
        obstime: double(obstime), $
        dz: float(dz), $
        Bx: float(Bx), $
        By: float(By), $
        Bz: float(Bz), $
        chromo_n0: float(chromo_n0), $
        chromo_np: float(chromo_np), $
        chromo_nHI: float(chromo_nHI), $
        chromo_T0: float(chromo_T0), $
        corona_Bavg: float(corona_Bavg), $
        corona_L: float(corona_L), $
        chromo_uniform_Bavg: float(chromo_uniform_Bavg), $
        chromo_uniform_L: float(chromo_uniform_L), $
        VoxelID: byte(VoxelID), $
        corona_ID1: byte(corona_ID1), $
        corona_ID2: byte(corona_ID2), $
        chromo_uniform_ID1: byte(chromo_uniform_ID1), $
        chromo_uniform_ID2: byte(chromo_uniform_ID2)}
        
 return, model
end