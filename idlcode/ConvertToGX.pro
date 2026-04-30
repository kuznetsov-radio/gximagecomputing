function ConvertToGX__read_dataset, file_id, dataset_path, required=required, ok=ok
 ok=0b
 catch, err
 if err ne 0 then begin
  catch, /cancel
  if keyword_set(required) then message, 'Missing required dataset: '+dataset_path
  return, 0
 endif

 dset_id=h5d_open(file_id, dataset_path)
 data=h5d_read(dset_id)
 h5d_close, dset_id
 ok=1b
 return, data
end

function ConvertToGX__read_dataset_fallback, file_id, dataset_paths, required=required, ok=ok, used_path=used_path
 ok=0b
 used_path=''
 n=n_elements(dataset_paths)
 if n eq 0 then begin
  if keyword_set(required) then message, 'Missing required dataset: <empty fallback list>'
  return, 0
 endif

 for i=0L, n-1L do begin
  data=ConvertToGX__read_dataset(file_id, dataset_paths[i], ok=oktry)
  if oktry then begin
   ok=1b
   used_path=dataset_paths[i]
   return, data
  endif
 endfor

 if keyword_set(required) then message, 'Missing required dataset (fallbacks): '+strjoin(dataset_paths, ', ')
 return, 0
end

function ConvertToGX__read_attr, obj_id, attr_name, ok=ok
 ok=0b
 catch, err
 if err ne 0 then begin
  catch, /cancel
  return, ''
 endif
 aid=h5a_open_name(obj_id, attr_name)
 v=h5a_read(aid)
 h5a_close, aid
 ok=1b
 return, v
end

function ConvertToGX__index_from_header, hdr_raw, ok=ok
 ok=0b
 cards=strsplit(string(hdr_raw), string(10b), /extract)
 if n_elements(cards) eq 0 then return, 0
 catch, err
 if err ne 0 then begin
  catch, /cancel
  return, 0
 endif
 idx=fitshead2struct(cards)
 catch, /cancel
 ok=1b
 return, idx
end

function ConvertToGX__card_int, key, value
 key8=strmid(strupcase(strtrim(key,2))+'        ', 0, 8)
 card=key8+'= '+string(long(value), format='(I20)')
 if strlen(card) lt 80 then card=card+strmid('                                                                                ', 0, 80-strlen(card))
 if strlen(card) gt 80 then card=strmid(card, 0, 80)
 return, card
end

function ConvertToGX__ensure_header_axes, hdr_cards, nx, ny
 cards=hdr_cards
 if n_elements(cards) eq 0 then cards=strarr(0)

 has_simple=0b
 has_bitpix=0b
 has_naxis=0b
 has_naxis1=0b
 has_naxis2=0b
 for i=0L, n_elements(cards)-1L do begin
  key=strupcase(strtrim(strmid(cards[i], 0, 8), 2))
  if key eq 'SIMPLE' then has_simple=1b
  if key eq 'BITPIX' then has_bitpix=1b
  if key eq 'NAXIS' then has_naxis=1b
  if key eq 'NAXIS1' then has_naxis1=1b
  if key eq 'NAXIS2' then has_naxis2=1b
 endfor

 if ~has_simple then cards=[cards, 'SIMPLE  =                    1                                                  ']
 if ~has_bitpix then cards=[cards, 'BITPIX  =                  -32                                                  ']
 if ~has_naxis then cards=[cards, ConvertToGX__card_int('NAXIS', 2L)]
 if ~has_naxis1 then cards=[cards, ConvertToGX__card_int('NAXIS1', long(nx))]
 if ~has_naxis2 then cards=[cards, ConvertToGX__card_int('NAXIS2', long(ny))]

 return, cards
end

function ConvertToGX__list_refmap_names, file_id, infile=infile
 names=''
 nfound=0L
 catch, err
 if err ne 0 then begin
  catch, /cancel
  names=''
 endif

 ; Try native IDL HDF5 enumeration first.
 catch, err1
 if err1 eq 0 then begin
  grp_id=h5g_open(file_id, '/refmaps')
  nobj=long(h5g_get_num_objs(grp_id))
  if nobj gt 0 then begin
   for i=0L, nobj-1L do begin
    nm=strtrim(string(h5g_get_objname_by_idx(grp_id, i)), 2)
    if nm eq '' then continue
    if nfound eq 0 then names=[nm] else names=[names, nm]
    nfound+=1L
   endfor
  endif
  h5g_close, grp_id
  catch, /cancel
 endif else begin
  catch, /cancel
 endelse

 ; Fallback for IDL builds without H5G_GET_OBJNAME_BY_IDX:
 ; query /refmaps keys via Python+h5py.
 if (nfound eq 0) and keyword_set(infile) then begin
  cmd='python -c "import h5py;f=h5py.File('''+string(infile)+''',''r'');g=f.get(''/refmaps'');print(''\\n''.join(list(g.keys())) if g is not None else '''')"'
  spawn, cmd, out
  if n_elements(out) gt 0 then begin
   w=where(strtrim(out,2) ne '', nw)
   if nw gt 0 then begin
    names=out[w]
    nfound=nw
   endif
  endif
 endif

 if nfound eq 0 then return, ''
 return, names
end

function ConvertToGX__build_refmaps, file_id, infile=infile
 ; Build an SSW map object container from /refmaps/*/data datasets.
 ; Returns a map object reference (possibly empty).
 resolve_routine, 'fitshead2wcs', /either
 resolve_routine, 'wcs2map', /either
 resolve_routine, 'valid_map', /either
 resolve_routine, 'make_map', /either
 refmaps=obj_new('map')
 nadded=0L
 names=ConvertToGX__list_refmap_names(file_id, infile=infile)
 if (n_elements(names) eq 0) or ((n_elements(names) eq 1) and (strtrim(string(names[0]),2) eq '')) then return, refmaps
 ; Preserve source ordering when available.
 nname=n_elements(names)
 ord=lonarr(nname) + 2147483647L
 has_order=bytarr(nname)
 for io=0L, nname-1L do begin
  catch, eo
  if eo ne 0 then begin
   catch, /cancel
   continue
  endif
  rgid=h5g_open(file_id, '/refmaps/'+names[io])
  oidx=ConvertToGX__read_attr(rgid, 'order_index', ok=okord)
  h5g_close, rgid
  catch, /cancel
  if okord then begin
   ord[io]=long(oidx)
   has_order[io]=1b
  endif
 endfor
 if total(has_order) gt 0 then begin
  order_idx=sort(ord)
  names=names[order_idx]
 endif
 for i=0L, n_elements(names)-1L do begin
  name=names[i]
  data=ConvertToGX__read_dataset(file_id, '/refmaps/'+name+'/data', ok=okdata)
  if ~okdata then continue

  ; Preferred path: build map from WCS header via SSW WCS tools.
  m=0
  idx=0
  map_id=strtrim(string(name), 2)
  t_obs='1970-01-01T00:00:00'
  okidx=0b
  fallback_reason=''
  hdr_raw=ConvertToGX__read_dataset(file_id, '/refmaps/'+name+'/wcs_header', ok=okhdr)
  if okhdr then begin
   hdr_text=string(hdr_raw)
   hdr_cards=strsplit(hdr_text, string(10b), /extract)
   sdat=size(data, /dimensions)
   if n_elements(sdat) ge 2 then hdr_cards=ConvertToGX__ensure_header_axes(hdr_cards, sdat[0], sdat[1])
   if (n_elements(hdr_cards) eq 0) then begin
    fallback_reason='empty wcs_header'
   endif else begin
    catch, err
    if err eq 0 then begin
     emsg=''
     wcs=fitshead2wcs(hdr_cards, errmsg=emsg)
     if emsg ne '' then begin
      m=0
      fallback_reason='fitshead2wcs: '+emsg
     endif else begin
      idx=ConvertToGX__index_from_header(hdr_raw, ok=okidx)
      if okidx then begin
       if tag_exist(idx, 'ID') then begin
        id_try=strtrim(string(idx.id), 2)
        if id_try ne '' then map_id=id_try
       endif
       if tag_exist(idx, 'DATE_OBS') then t_obs=string(idx.date_obs)
      endif
      emsg=''
      wcs2map, data, wcs, m, id=map_id, errmsg=emsg
      if emsg ne '' then begin
       emsg_t=''
       data_t=transpose(data)
       wcs2map, data_t, wcs, m, id=map_id, errmsg=emsg_t
       if emsg_t ne '' then begin
        m=0
        fallback_reason='wcs2map: '+emsg+' | transpose: '+emsg_t
       endif
      endif
     endelse
     catch, /cancel
    endif else begin
     catch, /cancel
     m=0
     fallback_reason='WCS conversion exception'
    endelse
   endelse
  endif else fallback_reason='missing wcs_header'

  ; Fallback: still create an SSW map structure if WCS conversion failed.
  if (datatype(m, 1) ne 'Structure') or (~valid_map(m)) then begin
   ; Build a minimal, guaranteed-valid map structure.
   sx=size(data, /dimensions)
   xc=0d & yc=0d & dx=1d & dy=1d & tobs='1970-01-01T00:00:00'
   if okidx then begin
    if tag_exist(idx, 'CRVAL1') then xc=double(idx.crval1)
    if tag_exist(idx, 'CRVAL2') then yc=double(idx.crval2)
    if tag_exist(idx, 'CDELT1') then dx=double(idx.cdelt1)
    if tag_exist(idx, 'CDELT2') then dy=double(idx.cdelt2)
    if tag_exist(idx, 'DATE_OBS') then tobs=string(idx.date_obs) else tobs=t_obs
   endif else begin
    tobs=t_obs
   endelse
   m=make_map(float(data), xc=xc, yc=yc, dx=dx, dy=dy, $
              xunits='arcsec', yunits='arcsec', id=map_id, time=tobs)
   print, 'REFMAP ', map_id, ': fallback make_map (', fallback_reason, ')'
  endif else begin
   print, 'REFMAP ', map_id, ': wcs2map OK'
  endelse
  if ~valid_map(m) then begin
   print, 'REFMAP ', map_id, ': SKIP (invalid map structure after fallback)'
   continue
  endif
  if okidx then begin
   ; Append map + index explicitly to ensure MAP object count is updated.
   refmaps->set, map=m, index=idx, /add
  endif else begin
   refmaps->set, map=m, /add
  endelse
  nadded+=1L
 endfor

 return, refmaps
end

function ConvertToGX__infer_line_depth, arr, nx, ny, fallback_depth, name
 nxy=long(nx)*long(ny)
 if nxy le 0 then return, long(fallback_depth)
 nelt=n_elements(arr)
 if (nelt mod nxy) ne 0 then begin
  message, name+' has incompatible element count: '+strtrim(string(nelt),2)+' for Nx*Ny='+strtrim(string(nxy),2)
 endif
 return, long(nelt / nxy)
end

function ConvertToGX, infile
 inlower=strlowcase(string(infile))
 nch=strlen(inlower)
 ext3=''
 ext4=''
 ext5=''
 if nch ge 3 then ext3=strmid(inlower, nch-3, 3)
 if nch ge 4 then ext4=strmid(inlower, nch-4, 4)
 if nch ge 5 then ext5=strmid(inlower, nch-5, 5)
 if (ext4 eq '.sav') or (ext4 eq '.xdr') then begin
  restore, infile
  return, box
 endif
 if (ext3 ne '.h5') and (ext5 ne '.hdf5') then message, 'Unsupported model format: '+infile

 resolve_routine, 'fitshead2struct', /either
 file_id=h5f_open(infile)
 gch_id=h5g_open(file_id, '/chromo')

 dr_raw=ConvertToGX__read_dataset_fallback(file_id, ['/corona/dr', '/chromo/dr'], /required)
 dz_raw=ConvertToGX__read_dataset(file_id, '/chromo/dz', /required)
 bcube_raw=ConvertToGX__read_dataset_fallback(file_id, ['/chromo/bcube', '/corona/bcube'], ok=okbcube)
 if ~okbcube then begin
  cor_bx=ConvertToGX__read_dataset(file_id, '/corona/bx', /required)
  cor_by=ConvertToGX__read_dataset(file_id, '/corona/by', /required)
  cor_bz=ConvertToGX__read_dataset(file_id, '/corona/bz', /required)
  sc=size(cor_bx, /dimensions)
  if n_elements(sc) ne 3 then message, 'Unsupported /corona component shape.'
  ; cor_bx is interpreted by IDL as (Nx, Ny, Nz) here.
  nxc=sc[0]
  nyc=sc[1]
  nzc=sc[2]
  ; Keep raw ordering expected by:
  ; bcube = transpose(bcube_raw, [1,2,3,0]) -> (Nx, Ny, Nz, 3)
  bcube_raw=fltarr(3, nxc, nyc, nzc)
  bcube_raw[0, *, *, *]=float(cor_bx)
  bcube_raw[1, *, *, *]=float(cor_by)
  bcube_raw[2, *, *, *]=float(cor_bz)
 endif

 chromo_bcube_raw=ConvertToGX__read_dataset(file_id, '/chromo/chromo_bcube', ok=okchrb)
 if ~okchrb then begin
  chr_bx=ConvertToGX__read_dataset(file_id, '/chromo/bx', /required)
  chr_by=ConvertToGX__read_dataset(file_id, '/chromo/by', /required)
  chr_bz=ConvertToGX__read_dataset(file_id, '/chromo/bz', /required)
  sch=size(chr_bx, /dimensions)
  if n_elements(sch) ne 3 then message, 'Unsupported /chromo component shape.'
  nxz=sch[0]
  nyz=sch[1]
  nzz=sch[2]
  ; Keep the legacy raw axis order expected by:
  ; chromo_bcube=transpose(chromo_bcube_raw, [2, 3, 1, 0]) -> (Nx, Ny, Nz, 3)
  chromo_bcube_raw=fltarr(3, nzz, nxz, nyz)
  chromo_bcube_raw[0, *, *, *]=transpose(float(chr_bx), [2, 0, 1])
  chromo_bcube_raw[1, *, *, *]=transpose(float(chr_by), [2, 0, 1])
  chromo_bcube_raw[2, *, *, *]=transpose(float(chr_bz), [2, 0, 1])
 endif

 av_field_raw=double(ConvertToGX__read_dataset_fallback(file_id, ['/lines/av_field', '/chromo/av_field'], /required))
 phys_length_raw=double(ConvertToGX__read_dataset_fallback(file_id, ['/lines/phys_length', '/chromo/phys_length'], /required))
 voxel_status_raw=byte(ConvertToGX__read_dataset_fallback(file_id, ['/lines/voxel_status', '/chromo/voxel_status'], /required))
 startidx_raw=long(ConvertToGX__read_dataset_fallback(file_id, ['/lines/start_idx', '/chromo/start_idx'], /required))
 endidx_raw=long(ConvertToGX__read_dataset_fallback(file_id, ['/lines/end_idx', '/chromo/end_idx'], /required))

 chromo_idx=long(ConvertToGX__read_dataset(file_id, '/chromo/chromo_idx', /required))
 chromo_n=float(ConvertToGX__read_dataset(file_id, '/chromo/chromo_n', /required))
 chromo_t=float(ConvertToGX__read_dataset(file_id, '/chromo/chromo_t', /required))
 n_p=float(ConvertToGX__read_dataset(file_id, '/chromo/n_p', /required))
 n_hi=float(ConvertToGX__read_dataset(file_id, '/chromo/n_hi', /required))
 n_htot=float(ConvertToGX__read_dataset(file_id, '/chromo/n_htot', /required))
 tr=long(ConvertToGX__read_dataset(file_id, '/chromo/tr', /required))
 tr_h=float(ConvertToGX__read_dataset(file_id, '/chromo/tr_h', /required))
 chromo_layers_raw=long(ConvertToGX__read_dataset(file_id, '/chromo/chromo_layers', /required))
 corona_base_raw=long(ConvertToGX__read_dataset_fallback(file_id, ['/corona/corona_base', '/chromo/corona_base'], /required))
 chromo_layers=chromo_layers_raw[0]
 corona_base=corona_base_raw[0]

 chromo_mask=ConvertToGX__read_dataset(file_id, '/chromo/chromo_mask', ok=okmask)
 if okmask eq 0b then chromo_mask=ConvertToGX__read_dataset(file_id, '/base/chromo_mask', /required)
 chromo_mask=long(chromo_mask)

 base_bx=double(ConvertToGX__read_dataset(file_id, '/base/bx', ok=okbx))
 base_by=double(ConvertToGX__read_dataset(file_id, '/base/by', ok=okby))
 base_bz=double(ConvertToGX__read_dataset(file_id, '/base/bz', ok=okbz))
 base_ic=double(ConvertToGX__read_dataset(file_id, '/base/ic', ok=okic))

 ; Normalize orientation to match IDL CHR SAV conventions.
 ; IMPORTANT: IDL HDF5 read order differs from h5py, and for these products
 ; dz_raw is already in (x, y, z) order as seen by IDL.
 dr=double(dr_raw)
 dz=double(dz_raw)
 bcube=transpose(float(bcube_raw), [1, 2, 3, 0])
 chromo_bcube=transpose(float(chromo_bcube_raw), [2, 3, 1, 0])

 s=size(dz, /dimensions)
 sc=size(bcube, /dimensions)
 Nx=s[0]
 Ny=s[1]
 Nz=s[2]
 Ncor=sc[2]

 if ~okbx then base_bx=dblarr(Nx, Ny)
 if ~okby then base_by=dblarr(Nx, Ny)
 if ~okbz then base_bz=dblarr(Nx, Ny)
 if ~okic then base_ic=dblarr(Nx, Ny)

 Nline=ConvertToGX__infer_line_depth(av_field_raw, Nx, Ny, Ncor, 'av_field')
 Nline_phys=ConvertToGX__infer_line_depth(phys_length_raw, Nx, Ny, Nline, 'phys_length')
 if Nline_phys ne Nline then message, 'phys_length depth differs from av_field'
 Nline_status=ConvertToGX__infer_line_depth(voxel_status_raw, Nx, Ny, Nline, 'voxel_status')
 if Nline_status ne Nline then message, 'voxel_status depth differs from av_field'
 Nline_start=ConvertToGX__infer_line_depth(startidx_raw, Nx, Ny, Nline, 'start_idx')
 if Nline_start ne Nline then message, 'start_idx depth differs from av_field'
 Nline_end=ConvertToGX__infer_line_depth(endidx_raw, Nx, Ny, Nline, 'end_idx')
 if Nline_end ne Nline then message, 'end_idx depth differs from av_field'

 avfield=reform(av_field_raw, Nx, Ny, Nline, /overwrite)
 physlength=reform(phys_length_raw, Nx, Ny, Nline, /overwrite)
 status=reform(voxel_status_raw, Nx, Ny, Nline, /overwrite)
 startidx=reform(startidx_raw, Nx, Ny, Nline, /overwrite)
 endidx=reform(endidx_raw, Nx, Ny, Nline, /overwrite)

 lon=ConvertToGX__read_attr(gch_id, 'lon', ok=oklon)
 lat=ConvertToGX__read_attr(gch_id, 'lat', ok=oklat)
 dsun_obs=ConvertToGX__read_attr(gch_id, 'dsun_obs', ok=okdsun)
 obs_time=ConvertToGX__read_attr(gch_id, 'obs_time', ok=oktime)
 h5g_close, gch_id

 if ~oklon then lon=0d else lon=double(lon)
 if ~oklat then lat=0d else lat=double(lat)
 if ~okdsun then dsun_obs=1.495978707d13 else dsun_obs=double(dsun_obs)
 if ~oktime then obs_time='1970-01-01T00:00:00.000' else obs_time=string(obs_time)

 ; Convert Carrington-like longitude (0..360) to Stonyhurst-style range when needed.
 if lon gt 180d then lon=lon-360d

 model_id=ConvertToGX__read_dataset(file_id, '/metadata/id', ok=okid)
 execute=ConvertToGX__read_dataset(file_id, '/metadata/execute', ok=okexec)
 if ~okid then model_id=file_basename(string(infile)) else model_id=string(model_id)
 if ~okexec then execute='h5->sav conversion' else execute=string(execute)
 index_header=ConvertToGX__read_dataset(file_id, '/base/index', ok=okindexhdr)
 if ~okindexhdr then index_header=ConvertToGX__read_dataset(file_id, '/base/index_header', ok=okindexhdr)
 if ~okindexhdr then index_header=ConvertToGX__read_dataset(file_id, '/base/wcs_header', ok=okindexhdr)

 refmaps=ConvertToGX__build_refmaps(file_id, infile=infile)
 refmaps_ptr=ptr_new(refmaps)
 h5f_close, file_id

 if okindexhdr then begin
  index=ConvertToGX__index_from_header(index_header, ok=okindex)
 endif else okindex=0b

 if ~okindex then begin
  index={SIMPLE:1B, $
         BITPIX:0L, $
         NAXIS:2L, $
         NAXIS1:long(Nx), $
         NAXIS2:long(Ny), $
         WCSNAME:'HGC-CEA', $
         CRPIX1:(double(Nx)/2.0d)+0.5d, $
         CRVAL1:double(lon), $
         CTYPE1:'CRLN-CEA', $
         CUNIT1:'deg', $
         CDELT1:double(dr[0])*180d/!dpi, $
         CRPIX2:(double(Ny)/2.0d)+0.5d, $
         CRVAL2:double(lat), $
         CTYPE2:'CRLT-CEA', $
         CUNIT2:'deg', $
         CDELT2:double(dr[1])*180d/!dpi, $
         CROTA2:0d, $
         DATE_OBS:string(obs_time), $
         DSUN_OBS:double(dsun_obs), $
         SOLAR_B0:0d, $
         HGLN_OBS:double(lon), $
         HGLT_OBS:double(lat), $
         RSUN_REF:6.96d8, $
         DATE:string(obs_time), $
         TELESCOP:'SDO', $
         INSTRUME:'HMI', $
         WAVELNTH:6173d}
 endif

 base={BX:double(base_bx), $
       BY:double(base_by), $
       BZ:double(base_bz), $
       IC:double(base_ic), $
       CHROMO_MASK:long(chromo_mask)}

 box={DR:double(dr), $
      ADD_BASE_LAYER:0, $
      INDEX:index, $
      REFMAPS:refmaps_ptr, $
      ID:string(model_id), $
      EXECUTE:string(execute), $
      STATUS:byte(status), $
      STARTIDX:long(startidx), $
      ENDIDX:long(endidx), $
      AVFIELD:double(avfield), $
      PHYSLENGTH:double(physlength), $
      BCUBE:float(bcube), $
      CHROMO_IDX:long(chromo_idx), $
      CHROMO_BCUBE:float(chromo_bcube), $
      N_HTOT:float(n_htot), $
      N_HI:float(n_hi), $
      N_P:float(n_p), $
      DZ:double(dz), $
      CHROMO_N:float(chromo_n), $
      CHROMO_T:float(chromo_t), $
      CHROMO_LAYERS:long(chromo_layers), $
      TR:long(tr), $
      TR_H:double(tr_h), $
      CORONA_BASE:long(corona_base), $
      BASE:base}

 return, box
end
