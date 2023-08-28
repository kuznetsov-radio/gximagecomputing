function LoadGXmodel, infile
 restore, infile
 
 obstime=anytim(box.index.date_obs)
 
 setenv, 'WCS_RSUN=6.96d8' ;the same command as in GX Simulator routines
 wcs=fitshead2wcs(box.index)
 wcs_convert_from_coord, wcs, wcs.crval, 'hcc', xC, yC, zC, length_units='cm'
 DSun=wcs.position.dsun_obs*1d2
 RSun=wcs_rsun(unit='cm')
 lonC=atan(xC, zC)/!dpi*180
 latC=asin(yC/RSun)/!dpi*180
 
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
 
 corona_Bavg=fltarr(s[0], s[1], corona_layers)
 corona_L=fltarr(s[0], s[1], corona_layers)
 chromo_uniform_Bavg=fltarr(s[0], s[1], box.corona_base)
 chromo_uniform_L=fltarr(s[0], s[1], box.corona_base)
 
 Q=fltarr(s[0], s[1], sc[2])
 Q[box.idx]=box.bmed
 corona_Bavg=Q[*, *, box.corona_base : sc[2]-1]
 chromo_uniform_Bavg=Q[*, *, 0 : box.corona_base-1]
 
 Q[*]=0
 Q[box.idx]=box.length*RSun
 corona_L=Q[*, *, box.corona_base : sc[2]-1]
 chromo_uniform_L=Q[*, *, 0 : box.corona_base-1]
 
 model={Nx: long(Nx), $
        Ny: long(Ny), $
        Nz: long(Nz), $
        chromo_layers: long(chromo_layers), $
        corona_layers: long(corona_layers), $
        corona_base: long(corona_base), $
        DSun: double(DSun), $
        RSun: double(Rsun), $
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
        chromo_uniform_L: float(chromo_uniform_L)}
        
 return, model
end
