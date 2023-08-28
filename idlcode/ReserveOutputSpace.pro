function ReserveOutputSpace, simbox
 out={flagsAll: lonarr(6), $
      flagsCorona: lonarr(6), $
      TI: dblarr(simbox.Nx, simbox.Ny, simbox.Nf), $
      TV: dblarr(simbox.Nx, simbox.Ny, simbox.Nf)}
      
 return, out
end
