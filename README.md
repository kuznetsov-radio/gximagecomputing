This code computes the 2D maps of the solar gyroresonance and free-free microwave emission using the models of active regions created by the GX Simulator. The code is called from IDL (requires the SolarSoft package).

To compute the emission maps, you firstly need to create the input data blocks by calling the following functions (in any order):

1. Load the GX Simulator model:<br/>
   model=LoadGXmodel(modelfile)<br/>
   where modelfile is the name of the GX Simulator model file (the model must contain the field line information and the chromospheric part).
   
2. Load the EBTEL table:<br/>
   ebtel=LoadEBTEL(ebtelfile)<br/>
   where ebtelfile is the name of the GX Simulator file containing the EBTEL table that defines the DEM and/or DDM.<br/> 
   If ebtelfile='' then the DEM, DDM, and coronal heating model are not used, and the coronal plasma is described by a model with a constant temperature and a barometric height profile of the plasma density (see below).

3. Define the size and position of the required radio maps, as well as the emission frequencies:<br/>
   simbox=MakeSimulationBox(xc, yc, dx, dy, Nx, Ny, freqlist)<br/>
   where:<br/>
   xc and yc are the x and y coordinates of the map center (in the helioprojective coordinate system, in arcseconds);<br/>
   dx and dy are the x and y resolutions of the map (in arcseconds);<br/>
   Nx and Ny are the x and y sizes of the map (in pixels);<br/>
   freqlist is the list (1D array) of the emission frequencies (in GHz) where the maps are to be computed.
   
4. Define the parameters of the coronal plasma:<br/>
   coronaparms=DefineCoronaParams(Tbase, nbase, Q0, a, b [, /force_isothermal])<br/>
   where:<br/>
   Tbase and nbase define the "default" plasma distribution; they are respectively the plasma temperature (in K) and the base plasma density at the bottom of the simulation box (in cm^{-3}). These parameters are used to find the plasma parameters in the voxels where the heating model is not applicable, i.e., either the voxel is associated with an open field line, or the heating parameters are beyond the EBTEL table. In such voxels, the plasma temperature is set to Tbase, and the plasma density is computed using nbase, Tbase, and the barometric formula.<br/>
   Q0, a, and b define the coronal heating model (which is applied to the closed field lines).
