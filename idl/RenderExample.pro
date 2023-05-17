pro RenderExample
 ;GX model:
 modelfile='C:\MCloud\CoronalMW\AR-SRH\Data\Models\model20220130_0415.sav'
 
 ;file with EBTEL table (can be '' if DEM, DDM and heating model are not needed):
 ebtelfile='C:\MCloud\CoronalMW\AR-SRH\Data\ebtel\ebtelDEM.sav'
  
 ;window parameters:
 xc=-50.0
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
 Q0=4.5d-3
 a=1.5
 b=2.5
           
 forward_function LoadGXmodel, LoadEBTEL, MakeSimulationBox, DefineCoronaParams, ReserveOutputSpace
 
 model=LoadGXmodel(modelfile)
 ebtel=LoadEBTEL(ebtelfile)
 simbox=MakeSimulationBox(xc, yc, dx, dy, Nx, Ny, freqlist)       
 coronaparms=DefineCoronaParams(Tbase, nbase, Q0, a, b, /force_isothermal)
 outspace=ReserveOutputSpace(simbox)
 
 tm=systime(1)
 r=call_external('RenderX.dll', 'ComputeMW', model, ebtel, simbox, coronaparms, outspace)
 print, 'Elapsed time: ', systime(1)-tm, ' s'                  
 
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
 
 save, mapI, mapV, filename='ExampleMaps.sav', /compress
end