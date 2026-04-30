function LoadEUVresponse__resolve_obstime, obs_input
 if n_params() eq 0 then message, 'Please provide OBS_TIME or a model structure containing OBSTIME.'

 obs_time=obs_input
 if size(obs_input, /type) eq 8 then begin
  tags=tag_names(obs_input)
  idx=where(strupcase(tags) eq 'OBSTIME', count)
  if count le 0 then message, 'MODEL input to LoadEUVresponse must contain OBSTIME.'
  obs_time=obs_input.(idx[0])
 endif

 return, anytim(obs_time, /vms)
end

function LoadEUVresponse, obs_input, instrument=instrument, evenorm=evenorm, chiantifix=chiantifix
 obs_time_vms=LoadEUVresponse__resolve_obstime(obs_input)
 if ~keyword_set(instrument) then begin
  instrument='aia'
  message, 'WARNING: No explicit instrument name was provided, AIA was assumed.', /continue
 endif
 instrument=strlowcase(instrument)
 ds=-1                 
      
 while ds lt 0 do case instrument of      
  'aia': begin
          if ~exist(evenorm) then evenorm=1
          if ~exist(chiantifix) then chiantifix=1                         
          r=aia_get_response(timedepend_date=obs_time_vms, $
                             /temperature, /dn, evenorm=evenorm, chiantifix=chiantifix)
          ds=0.36d0                   
        end  
  'aia2': begin
           r=gx_get_aia_response(obs_time_vms)
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
               r=gx_get_eui_response(obs_time_vms, /fsi)
               ds=4.4401245d0^2
              end
  'solo-hri': begin
               r=gx_get_eui_response(obs_time_vms, /hri)
               ds=0.49200001d0^2
              end   
  'stereo-a': begin
               r=gx_get_euvi_response(obs_time_vms, /a)
               ds=2.5281d0
              end                                   
  'stereo-b': begin
               r=gx_get_euvi_response(obs_time_vms, /b)
               ds=2.5281d0
              end
  else: begin
         message, 'WARNING: Unknown instrument "'+strtrim(instrument, 2)+'"; AIA was assumed.', /continue
         instrument='aia'
        end
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
