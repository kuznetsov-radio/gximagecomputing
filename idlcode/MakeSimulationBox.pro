function MakeSimulationBox, xc, yc, dx, dy, Nx, Ny, freqlist
 Nf=n_elements(freqlist)

 simbox={Nx: long(Nx), $
         Ny: long(Ny), $
         Nf: long(Nf), $
         _r1: 0L, $
         xc: double(xc), $
         yc: double(yc), $
         dx: double(dx), $
         dy: double(dy), $
         freqlist: double(freqlist)}
         
 return, simbox        
end
