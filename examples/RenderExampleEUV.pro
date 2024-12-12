pro RenderExampleEUV
 ;GX model:
 modelfile='C:\OneDrive\EUV\model20220131_0420.sav'
 
 ;file with EBTEL table (can be '' if DEM, DDM and heating model are not needed):
 ebtelfile='C:\ssw\packages\gx_simulator\euv\ebtel\ebtel_ss.sav'
  
 ;window parameters:
 xc=170.0
 yc=390.0
 dx=2.0
 dy=2.0
 Nx=150
 Ny=150
 
 ;analytical corona parameters:          
 Tbase=1d6
 nbase=1d8
 
 ;heating model parameters:
 Q0=0.0217
 a=0.3
 b=2.7
 
 ;selective heating:
 SHtable=dblarr(7, 7)
 w=[1.0, 1.0, 1.0, 1.1, 1.2, 1.3, 1.4]
 for i=0, 6 do for j=0, 6 do SHtable[i, j]=w[i]*w[j]
 SHtable[6, 6]=0.1 ;UMB2UMB
           
 forward_function LoadGXmodel, LoadEBTEL, LoadEUVresponse, MakeSimulationBoxEUV, DefineCoronaParams, ReserveOutputSpaceEUV

 tm=systime(1) 
 model=LoadGXmodel(modelfile)
 ebtel=LoadEBTEL(ebtelfile)
 response=LoadEUVresponse(model)
 simbox=MakeSimulationBoxEUV(xc, yc, dx, dy, Nx, Ny)       
 coronaparms=DefineCoronaParams(Tbase, nbase, Q0, a, b)
 outspace=ReserveOutputSpaceEUV(simbox, response)
 print, 'Elapsed time (loading): ', systime(1)-tm, ' s'
 
 tm=systime(1)
 r=call_external('C:\pCloud\Codes\RenderGR\RenderX\x64\Release\RenderGRFF_64.dll', 'ComputeEUV', $
                 model, ebtel, response, simbox, coronaparms, outspace, SHtable)
 print, 'Elapsed time (main): ', systime(1)-tm, ' s'            
 
 ConvertToMapsEUV, outspace, simbox, model, response, mapEUV
 
 window, 1, title='EUV map'
 wset, 1
 loadct, 13, /silent
 m=mapEUV.getmap(2)
 plot_map, m, cbar=1

 save, mapEUV, filename='C:\OneDrive\EUV\EUVmaps.sav', /compress
end