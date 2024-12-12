function MakeSimulationBoxEUV, xc, yc, dx, dy, Nx, Ny
 simbox={Nx: long(Nx), $
         Ny: long(Ny), $
         xc: double(xc), $
         yc: double(yc), $
         dx: double(dx), $
         dy: double(dy)}
         
 return, simbox        
end