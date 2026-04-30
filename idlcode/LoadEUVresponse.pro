function LoadEUVresponse, model, instrument=instrument, evenorm=evenorm, chiantifix=chiantifix
 if ~keyword_set(instrument) then instrument='aia'     
 instrument=strlowcase(instrument)
 ds=-1                 
      
 while ds lt 0 do case instrument of      
  'aia': begin
          if ~exist(evenorm) then evenorm=1
          if ~exist(chiantifix) then chiantifix=1                         
          r=aia_get_response(timedepend_date=anytim(model.obstime, /vms), $
                             /temperature, /dn, evenorm=evenorm, chiantifix=chiantifix)
          ds=0.36d0                   
        end  
  'aia2': begin
           r=gx_get_aia_response(anytim(model.obstime, /vms))
           ds=0.36d0                   
          end      
  'trace': begin
            restore, gx_findfile('trace_response.sav')
            r=response
            ds=1d0
           end        
  'sxt': begin
         restore, gx_findfile('sxt_response.sav')
         r=response
         ds=6.0025d0
        end
  'solo-fsi': begin
               r=gx_get_eui_response(anytim(model.obstime, /vms), /fsi)
               ds=4.4401245d0^2
              end
  'solo-hri': begin
               r=gx_get_eui_response(anytim(model.obstime, /vms), /hri)
               ds=0.49200001d0^2
              end   
  'stereo-a': begin
               r=gx_get_euvi_response(anytim(model.obstime, /vms), /a)
               ds=2.5281d0
              end                                   
  'stereo-b': begin
               r=gx_get_euvi_response(anytim(model.obstime, /vms), /b)
               ds=2.5281d0
              end
  else: instrument='aia'
 endcase             
           
 NT=n_elements(r.logte)         
 Nchannels=n_elements(r.channels)
 
 response={ds: double(ds), $ ;pixel size, arcsec^2
           NT: long(NT), $
           Nchannels: long(Nchannels), $
           logte: double(r.logte), $
           all: double(r.all), $
           channels: r.channels, $
           instrument: strupcase(instrument)}
          
 return, response         
end