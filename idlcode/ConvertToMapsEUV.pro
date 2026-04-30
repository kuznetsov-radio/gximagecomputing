pro ConvertToMapsEUV, out, box, model, response, mapCorona, mapTR, B0=b0, L0=l0, RSun=rsun
 mapCorona=obj_new('map')
 mapTR=obj_new('map')
 
 for k=0, response.Nchannels-1 do begin
  m=make_map(out.fluxCorona[*, *, k], xc=box.xc, yc=box.yc, dx=box.dx, dy=box.dy, $
             id=response.instrument+' '+response.channels[k], time=anytim(model.obstime, /vms), $
             channel=response.channels[k], units='DN s^-1 pix^-1', B0=b0, L0=l0, RSun=rsun)
  mapCorona->setmap, k, m
  m=make_map(out.fluxTR[*, *, k], xc=box.xc, yc=box.yc, dx=box.dx, dy=box.dy, $
             id=response.instrument+' '+response.channels[k], time=anytim(model.obstime, /vms), $
             channel=response.channels[k], units='DN s^-1 pix^-1', B0=b0, L0=l0, RSun=rsun)
  mapTR->setmap, k, m
 endfor
end
