pro ExportLOSparms, libname, model, ebtel, simbox, coronaparms, $
                    freqlist, Nvoxels_M, Lparms_M, Rparms_M, Parms_M, $
                    xc, yc, dx, dy, SHtable=SHtable
 cp=coronaparms
 cp.mode=cp.mode or 1024L
 outspace={NvoxMax: 0L, $
           r1: 0L, $
           Nvoxels_M: lonarr(simbox.Nx, simbox.Ny)} 
           
 if exist(SHtable) then r=call_external(libname, 'ComputeMW', model, ebtel, simbox, cp, outspace, SHtable) $
                   else r=call_external(libname, 'ComputeMW', model, ebtel, simbox, cp, outspace)
                   
 NvoxMax=long(max(outspace.Nvoxels_M))
 
 cp.mode=cp.mode or 2048L        
 r=((simbox.Nx*simbox.Ny mod 2) eq 0) ? 2 : 1  
 outspace={NvoxMax: NvoxMax, $
           r1: 0L, $
           Nvoxels_M: lonarr(simbox.Nx, simbox.Ny), $
           r2: lonarr(r), $
           Lparms_X: lonarr(5, simbox.Nx, simbox.Ny), $
           r3: lonarr(r), $
           Rparms_M: dblarr(3, simbox.Nx, simbox.Ny), $
           Parms_M: dblarr(17, NvoxMax, simbox.Nx, simbox.Ny)}                    
           
 if exist(SHtable) then r=call_external(libname, 'ComputeMW', model, ebtel, simbox, cp, outspace, SHtable) $
                   else r=call_external(libname, 'ComputeMW', model, ebtel, simbox, cp, outspace)
                   
 freqlist=simbox.freqlist
 Nvoxels_M=outspace.Nvoxels_M
 Lparms_M=[simbox.Nx*simbox.Ny, NvoxMax, simbox.Nf, 0, 0, 0]
 Rparms_M=outspace.Rparms_M
 Parms_M=outspace.Parms_M
 xc=simbox.xc
 yc=simbox.yc
 dx=simbox.dx
 dy=simbox.dy
end