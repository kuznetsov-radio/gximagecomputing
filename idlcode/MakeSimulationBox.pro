function MakeSimulationBox, xc, yc, dx, dy, Nx, Ny, freqlist, rot=rot, $
                            parallel=parallel, exact=exact, Nthreads=Nthreads
 Nf=n_elements(freqlist)

 if ~exist(rot) then rot=0d0
 
 projection=0L
 if keyword_set(parallel) then projection=projection or 1L
 if keyword_set(exact) then projection=projection or 2L
 
 if exist(Nthreads) then begin
  Nt=long(Nthreads)
  if Nthreads gt 0 then projection=projection or ishft(Nt, 16)
 endif 

 simbox={Nx: long(Nx), $
         Ny: long(Ny), $
         Nf: long(Nf), $
         projection: long(projection), $
         xc: double(xc), $
         yc: double(yc), $
         dx: double(dx), $
         dy: double(dy), $
         rot: double(rot), $
         freqlist: double(freqlist)}
         
 return, simbox        
end