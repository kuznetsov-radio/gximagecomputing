pro RenderExampleMW
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
 
 ;working frequencies:
 freqlist=[5.8, 6.2, 6.6, 7.0, 7.4, 7.8, 8.2, 8.6, 9.0, 9.4, 9.8, 10.2, 10.6, 11.0, 11.4, 11.8]
           
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
           
 forward_function LoadGXmodel, LoadEBTEL, MakeSimulationBox, DefineCoronaParams, ReserveOutputSpace
 
 tm=systime(1) 
 model=LoadGXmodel(modelfile)
 ebtel=LoadEBTEL(ebtelfile)
 simbox=MakeSimulationBox(xc, yc, dx, dy, Nx, Ny, freqlist)       
 coronaparms=DefineCoronaParams(Tbase, nbase, Q0, a, b)
 outspace=ReserveOutputSpace(simbox)
 print, 'Elapsed time (loading): ', systime(1)-tm, ' s'
 
 tm=systime(1)
 r=call_external('C:\pCloud\Codes\RenderGR\RenderX\x64\Release\RenderGRFF_64.dll', 'ComputeMW', $
                 model, ebtel, simbox, coronaparms, outspace, SHtable)
 print, 'Elapsed time (main): ', systime(1)-tm, ' s'                  
 
 ConvertToMaps, outspace, simbox, model, mapI, mapV
 
 window, 1, title='I'
 wset, 1
 loadct, 13, /silent
 m=mapI.getmap(0)
 plot_map, m, cbar=1

 window, 2, title='V'
 wset, 2
 loadct, 33, /silent
 m=mapV.getmap(0)
 plot_map, m, cbar=1 
 
 save, mapI, mapV, filename='C:\OneDrive\EUV\MWmaps.sav', /compress
end