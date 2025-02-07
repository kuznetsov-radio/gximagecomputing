function ReserveOutputSpaceEUV, simbox, response
 out={flagsAll: lonarr(6), $
      flagsCorona: lonarr(6), $
      fluxCorona: dblarr(simbox.Nx, simbox.Ny, response.Nchannels), $
      fluxTR:     dblarr(simbox.Nx, simbox.Ny, response.Nchannels)}
      
 return, out
end