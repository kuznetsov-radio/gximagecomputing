pro ConvertModelFiles
 ;GX model:
 modelfile='C:\MCloud\CoronalMW\AR-SRH\Data\Models\model20220130_0415.sav'
 
 ;file with EBTEL table (can be '' if DEM, DDM and heating model are not needed):
 ebtelfile='C:\MCloud\CoronalMW\AR-SRH\Data\ebtel\ebtelDEM.sav'
           
 forward_function LoadGXmodel, LoadEBTEL
 
 model=LoadGXmodel(modelfile)
 ebtel=LoadEBTEL(ebtelfile)
 
 openw, 1, 'model.bin'
 writeu, 1, model
 close, 1
 
 openw, 1, 'ebtel.bin'
 writeu, 1, ebtel
 close, 1
end
