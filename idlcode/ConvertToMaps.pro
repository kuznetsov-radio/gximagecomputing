pro ConvertToMaps, out, box, model, mapI, mapV, flux=flux
 sfu=1d-19
 kB=1.380649d-16
 c=2.99792458d10

 flux_on=keyword_set(flux)
 id=flux_on ? 'Stokes ' : 'Tb_'

 mapI=obj_new('map')
 mapV=obj_new('map')
 
 for k=0, box.Nf-1 do begin
  r=flux_on ? sfu*c^2/(2.0*kB*(box.freqlist[k]*1d9)^2)/(box.dx*box.dy*(!dpi/180/60/60)^2) : 1d0
  m=make_map(out.TI[*, *, k]/r, xc=box.xc, yc=box.yc, dx=box.dx, dy=box.dy, $
             id=id+'I '+string(box.freqlist[k], format='(F5.2)')+' GHz', time=anytim(model.obstime, /vms), $
             freq=box.freqlist[k])
  mapI->setmap, k, m
  m=make_map(out.TV[*, *, k]/r, xc=box.xc, yc=box.yc, dx=box.dx, dy=box.dy, $
             id=id+'V '+string(box.freqlist[k], format='(F5.2)')+' GHz', time=anytim(model.obstime, /vms), $
             freq=box.freqlist[k])
  mapV->setmap, k, m
 endfor
end