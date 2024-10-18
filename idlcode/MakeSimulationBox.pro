function MakeSimulationBox, xc, yc, dx, dy, Nx, Ny, freqlist, rot=rot
 Nf=n_elements(freqlist)

 if ~exist(rot) then rot=0d0

 simbox={Nx: long(Nx), $
         Ny: long(Ny), $
         Nf: long(Nf), $
         _r1: 0L, $
         xc: double(xc), $
         yc: double(yc), $
         dx: double(dx), $
         dy: double(dy), $
         rot: double(rot), $
         freqlist: double(freqlist)}
         
 return, simbox        
end