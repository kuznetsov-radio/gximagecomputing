pro ConvertToMaps, out, box, model, mapI, mapV
 mapI=obj_new('map')
 mapV=obj_new('map')
 
 for k=0, box.Nf-1 do begin
  m=make_map(out.TI[*, *, k], xc=box.xc, yc=box.yc, dx=box.dx, dy=box.dy, $
             id='Tb_I '+string(box.freqlist[k], format='(F5.2)')+' GHz', time=anytim(model.obstime, /vms))
  mapI->setmap, k, m
  m=make_map(out.TV[*, *, k], xc=box.xc, yc=box.yc, dx=box.dx, dy=box.dy, $
             id='Tb_V '+string(box.freqlist[k], format='(F5.2)')+' GHz', time=anytim(model.obstime, /vms))
  mapV->setmap, k, m
 endfor
end