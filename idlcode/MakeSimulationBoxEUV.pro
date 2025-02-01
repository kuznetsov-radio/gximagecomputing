function MakeSimulationBoxEUV, xc, yc, dx, dy, Nx, Ny, $
                               parallel=parallel, exact=exact, Nthreads=Nthreads
 projection=0L
 if keyword_set(parallel) then projection=projection or 1L
 if keyword_set(exact) then projection=projection or 2L
 
 if exist(Nthreads) then begin
  Nt=long(Nthreads)
  if Nthreads gt 0 then projection=projection or ishft(Nt, 16)
 endif 

 simbox={Nx: long(Nx), $
         Ny: long(Ny), $
         xc: double(xc), $
         yc: double(yc), $
         dx: double(dx), $
         dy: double(dy), $
         projection: long(projection)}
         
 return, simbox        
end