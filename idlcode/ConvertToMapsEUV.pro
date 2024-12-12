pro ConvertToMapsEUV, out, box, model, response, mapEUV
 mapEUV=obj_new('map')
 
 for k=0, response.Nchannels-1 do begin
  m=make_map(out.flux[*, *, k], xc=box.xc, yc=box.yc, dx=box.dx, dy=box.dy, $
             id=response.instrument+' '+response.channels[k], time=anytim(model.obstime, /vms), $
             channel=response.channels[k])
  mapEUV->setmap, k, m
 endfor
end