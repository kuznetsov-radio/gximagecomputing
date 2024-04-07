This code computes the 2D maps of the solar gyroresonance and free-free microwave emission using the models of active regions created by the GX Simulator. The code is called from IDL (requires the SolarSoft package).

Quick start: see the file /examples/RenderExample.pro (the sample GX Simulator model and EBTEL data are not included).

To compute the emission maps, you firstly need to create the input data blocks by calling the following functions:

1. Load the GX Simulator model:<br/>
   model=LoadGXmodel(modelfile)<br/>
   where modelfile is the name of the GX Simulator model file (the model must contain the field line information and the chromospheric part).
   
2. Load the EBTEL table:<br/>
   ebtel=LoadEBTEL(ebtelfile [, DEM=DEM, DDM=DDM])<br/>
   where ebtelfile is the name of the GX Simulator file containing the EBTEL table(s) that define the DEM and/or DDM.<br/> 
   If ebtelfile='' then the DEM, DDM, and coronal heating model are not used, and the coronal plasma is described by a model with a constant temperature and a barometric height profile of the plasma density (see below).<br/>
   The keywords /DEM and /DDM are only applicable if the chosen file contains both the DEM and DDM tables. In this case, if the /DEM keyword is set, the code loads the DEM table only (the DDM table is ignored). Similarly, if the /DDM keyword is set, the code loads the DDM table only (the DEM table is ignored). If both /DEM and /DDM keywords (or none of them) are set, the code loads both tables.<br/>
   If the chosen file contains only one EBTEL table (either DEM or DDM), the code loads that table; the /DEM and /DDM keywords are ignored.

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
   Tbase and nbase define the "default" plasma distribution; they are respectively the plasma temperature (in K) and the base plasma density at the bottom of the simulation box (in cm^{-3}). These parameters are used to find the plasma parameters in the voxels where the heating model is not applicable, i.e., either the voxel is associated with an open field line, or the heating parameters are beyond the boundaries of the EBTEL table. In such voxels, the plasma temperature is set to Tbase, and the plasma density is computed using nbase, Tbase, and the barometric formula.<br/>
   Q0, a, and b define the coronal heating model (which is applied to the closed field lines). The heating rate Q at each field line is computed as Q=Q0*(B/B0)^a/(L/L0)^b, where B is the average magnetic field along the line, L is the line half-length, and B0 and L0 are some pre-defined constants (the same as in GX Simulator).<br/>
   /force_isothermal - if set, the multi-thermal formulae given in the paper of Fleishman, Kuznetsov & Landi (2021) are not used, and the emission is computed using the moments of the DEM or DDM distribution (if both DEM and DDM are provided, the DDM moments are used). This option improves the computation speed greatly, although the results become less accurate.
   
5. Prepare the memory structure for the simulation results:<br/>
   outspace=ReserveOutputSpace(simbox)<br/>
   where simbox is the structure returned by the MakeSimulationBox function.
   
When the input data are ready, the computation is performed by calling the main executable module (RenderGRFF_32.dll, RenderGRFF_64.dll, or RenderGRFF.so) via the call_external function:<br/>
r=call_external(libname, 'ComputeMW', model, ebtel, simbox, coronaparms, outspace)<br/>
where libname is the name of the appropriate executable library, and model, ebtel, simbox, coronaparms, and outspace are the structures returned by the above-mentioned functions.

The output structure outspace contains the fields outspace.TI and outspace.TV, which represent the brightness temperatures corresponding respectively to the Stokes parameters I and V of the computed emission (in K). Each field is a 3D array with Nx * Ny * Nf elements, where Nx and Ny are the x and y sizes of the computed maps, and Nf is the number of the emission frequencies. These data can be processed directly, or can be converted into the SolarSoft map objects via the procedure<br/>
ConvertToMaps, outspace, simbox, model, mapI, mapV<br/>
where the input parameters outspace, simbox, and model are the structures returned by the above-mentioned functions, and the output parameters mapI and mapV are the resulting SolarSoft (multi-frequency) map objects which represent the brightness temperatures corresponding respectively to the Stokes parameters I and V of the computed emission (in K).

An example of using the code is given in the file /examples/RenderExample.pro (the sample GX Simulator model and EBTEL data are not included).
