This code computes the 2D maps of the solar microwave (gyroresonance and free-free) and EUV (spectral lines) emission using the models of active regions created by the GX Simulator. The code is called from IDL (requires the SolarSoft GX_Simulator package).

Quick start: see the files /examples/RenderExampleMW.pro and /examples/RenderExampleEUV.pro (the sample GX Simulator model and EBTEL data are not included).

To compute the microwave emission maps, you firstly need to create the input data blocks by calling the following functions:

1. Load the GX Simulator model:<br/>
   model=LoadGXmodel(modelfile)<br/>
   where modelfile is the name of the GX Simulator model file (the model must contain the field line information and the chromospheric part).<br/>
   Note: the parameter tr_mask (see below) has no effect for the microwave emission.
   
2. Load the EBTEL tables:<br/>
   ebtel=LoadEBTEL(ebtelfile [, DEM=DEM, DDM=DDM])<br/>
   where ebtelfile is the name of the GX Simulator file containing the EBTEL table(s) that define the DEM and/or DDM.<br/> 
   If ebtelfile='' then the DEM, DDM, and coronal heating model are not used, and the coronal plasma is described by a model with a constant temperature and a barometric height profile of the plasma density (see below).<br/>
   The keywords /DEM and /DDM are only applicable if the chosen file contains both the DEM and DDM tables. In this case, if the /DEM keyword is set, the code loads the DEM table only (the DDM table is ignored). Similarly, if the /DDM keyword is set, the code loads the DDM table only (the DEM table is ignored). If both /DEM and /DDM keywords (or none of them) are set, the code loads both tables.<br/>
   If the chosen file contains only one EBTEL table (either DEM or DDM), the code loads that table; the /DEM and /DDM keywords are ignored.

3. Define the size and position of the required radio maps, as well as the emission frequencies:<br/>
   simbox=MakeSimulationBox(xc, yc, dx, dy, Nx, Ny, freqlist [, rot=rot])<br/>
   where:<br/>
   xc and yc are the x and y coordinates of the map center (in the helioprojective coordinate system, in arcseconds).<br/>
   dx and dy are the x and y resolutions of the map (in arcseconds).<br/>
   Nx and Ny are the x and y sizes of the map (in pixels).<br/>
   freqlist is the list (1D array) of the emission frequencies (in GHz) where the maps are to be computed.<br/>
   rot is the optional rotation angle of the simulated radio maps (in degrees). If set, the simulation window is rotated counter-clockwise by this angle; in this case, the center of the resulting map is still at (xc, yc), and the x and y coordinates correspond to the rotated coordinate system. Default: 0.
   
5. Define the parameters of the coronal plasma:<br/>
   coronaparms=DefineCoronaParams(Tbase, nbase, Q0, a, b [, /force_isothermal])<br/>
   where:<br/>
   Tbase and nbase define the "default" plasma distribution; they are respectively the plasma temperature (in K) and the base plasma density at the bottom of the simulation box (in cm^{-3}). These parameters are used to find the plasma parameters in the voxels where the heating model is not applicable, i.e., either the voxel is associated with an open field line, or the heating parameters are beyond the boundaries of the EBTEL table. In such voxels, the plasma temperature is set to Tbase, and the plasma density is computed using nbase, Tbase, and the barometric formula.<br/>
   Q0, a, and b define the coronal heating model (which is applied to the closed field lines). The heating rate Q at each field line is computed as Q=Q0*(B/B0)^a/(L/L0)^b, where B is the average magnetic field along the line, L is the line half-length, and B0 and L0 are some pre-defined constants (the same as in GX Simulator).<br/>
   /force_isothermal - if set, the multi-thermal formulae given in the paper of Fleishman, Kuznetsov & Landi (2021) are not used, and the emission is computed using the moments of the DEM or DDM distribution (if both DEM and DDM are provided, the DDM moments are used). This option improves the computation speed greatly, although the results become less accurate.<br/>
   Note: the option /AddTR (see below) has no effect for the microwave emission.
   
6. Prepare the memory structure for the simulation results:<br/>
   outspace=ReserveOutputSpace(simbox)<br/>
   where simbox is the structure returned by the MakeSimulationBox function.

7. Optionally: prepare the table SHtable defining selective heating of the coronal magnetic field lines (see below).
   
When the input data are ready, the computation is performed by calling the main executable module (RenderGRFF_32.dll, RenderGRFF_64.dll, or RenderGRFF.so) via the call_external function:<br/>
r=call_external(libname, 'ComputeMW', model, ebtel, simbox, coronaparms, outspace [, SHtable])<br/>
where libname is the name of the appropriate executable library, and model, ebtel, simbox, coronaparms, and outspace are the structures returned by the above-mentioned functions.

The output structure outspace contains the fields outspace.TI and outspace.TV, which represent the brightness temperatures corresponding respectively to the Stokes parameters I and V of the computed emission (in K). Each field is a 3D array with Nx * Ny * Nf elements, where Nx and Ny are the x and y sizes of the computed maps, and Nf is the number of the emission frequencies. These data can be processed directly, or can be converted into the SolarSoft map objects via the procedure<br/>
ConvertToMaps, outspace, simbox, model, mapI, mapV [, /flux]<br/>
where the input parameters outspace, simbox, and model are the structures returned by the above-mentioned functions, and the output parameters mapI and mapV are the resulting SolarSoft (multi-frequency) map objects which represent the brightness temperatures corresponding respectively to the Stokes parameters I and V of the computed emission (in K or sfu/pix, if the keyword /flux is not set or set, respectively).

An example of using the code is given in the file /examples/RenderExampleMW.pro (the sample GX Simulator model and EBTEL data are not included).

Computing the EUV emission maps is similar to that for the microwave emission, with a few differences. Firstly, you need to create the input data blocks by calling the following functions:

1. Load the GX Simulator model:<br/>
   model=LoadGXmodel(modelfile [,tr_mask=tr_mask])<br/>
   where modelfile is the name of the GX Simulator model file (the model must contain the field line information and the chromospheric part).<br/>
   The optional parameter tr_mask is a 2D array specifying the 'transition region mask', i.e., defining the regions where the emission from the transition region makes a contribution to the resulting EUV flux (if the /AddTR keyword is set in the DefineCoronaParams function); see GX Simulator (panel 'Transition Region Attributes') for details. By default, the entire transition region is considered.
   
2. Load the EBTEL tables:<br/>
   ebtel=LoadEBTEL(ebtelfile)<br/>
   where ebtelfile is the name of the GX Simulator file containing the EBTEL tables that define the DEM for the corona and transition region.<br/> 
   If ebtelfile='' then the DEM and coronal heating model are not used, and the coronal plasma is described by a model with a constant temperature and a barometric height profile of the plasma density.<br/>
   The keyword /DEM (see the above description of the LoadEBTEL function for the microwave emission) can be used but will have no effect. The keyword /DDM should not be used, because the EUV emission depends on the DEM only.

3. Load the instrumental response function:<br/>
   response=LoadEUVresponse(model [, instrument, evenorm=evenorm, chiantifix=chiantifix])<br/>
   where:<br/>
   model is the structure returned by the LoadGXmodel function.<br/>
   instrument is the name of the chosen instrument. Currently, the following instruments are supported: 'AIA', 'AIA2', 'TRACE', 'SXT', 'SOLO-FSI', 'SOLO-HRI', 'STEREO-A', 'STEREO-B'. Default: 'AIA'. Note that the code computes the emission as if observed from the Earth; thus for Solar Orbiter and STEREO some additional adjustments of the EUV flux will be needed.<br/>
   evenorm and chiantifix: these parameters are applicable to the AIA instrument only (see the SolarSoft function aia_get_response.pro). Default: evenorm=1, chiantifix=1.

5. Define the size and position of the required EUV maps:<br/>
   simbox=MakeSimulationBoxEUV(xc, yc, dx, dy, Nx, Ny)<br/>
   where:<br/>
   xc and yc are the x and y coordinates of the map center (in the helioprojective coordinate system, in arcseconds).<br/>
   dx and dy are the x and y resolutions of the map (in arcseconds).<br/>
   Nx and Ny are the x and y sizes of the map (in pixels).<br/>
   Note that the particular EUV channels cannot be selected: the code computes the emission for all channels specified by the instrumental response table. Also, the code computes the emission as if observed from the Earth; thus for Solar Orbiter and STEREO the map position and pixel size should be corrected accordingly.<br/>
   
6. Define the parameters of the coronal plasma:<br/>
   coronaparms=DefineCoronaParams(Tbase, nbase, Q0, a, b [, /AddTR])<br/>
   where:<br/>
   Tbase and nbase define the "default" plasma distribution; they are respectively the plasma temperature (in K) and the base plasma density at the bottom of the simulation box (in cm^{-3}). These parameters are used to find the plasma parameters in the voxels where the heating model is not applicable, i.e., either the voxel is associated with an open field line, or the heating parameters are beyond the boundaries of the EBTEL table. In such voxels, the plasma temperature is set to Tbase, and the plasma density is computed using nbase, Tbase, and the barometric formula.<br/>
   Q0, a, and b define the coronal heating model (which is applied to the closed field lines). The heating rate Q at each field line is computed as Q=Q0*(B/B0)^a/(L/L0)^b, where B is the average magnetic field along the line, L is the line half-length, and B0 and L0 are some pre-defined constants (the same as in GX Simulator).<br/>
   /AddTR - if set, the EUV contribution from the transition region is added to the emission maps. The emissions from the corona and transition region are computed using the respective DEM tables.<br/>
   Note: the option /force_isothermal (see above) has no effect for the EUV emission.
   
7. Prepare the memory structure for the simulation results:<br/>
   outspace=ReserveOutputSpaceEUV(simbox, response)<br/>
   where simbox is the structure returned by the MakeSimulationBox function, and response is the structure returned by the LoadEUVresponse function.

8. Optionally: prepare the table SHtable defining selective heating of the coronal magnetic field lines (see below).
   
When the input data are ready, the computation is performed by calling the main executable module (RenderGRFF_32.dll, RenderGRFF_64.dll, or RenderGRFF.so) via the call_external function:<br/>
r=call_external(libname, 'ComputeEUV', model, ebtel, response, simbox, coronaparms, outspace [, SHtable])<br/>
where libname is the name of the appropriate executable library, and model, ebtel, response, simbox, coronaparms, and outspace are the structures returned by the above-mentioned functions.

The output structure outspace contains the field flux, which represents the computed EUV flux (in DN s^{-1} pix^{-1}). It is a 3D array with Nx * Ny * Nchannels elements, where Nx and Ny are the x and y sizes of the computed maps, and Nchannels is the number of the EUV channels (defined by the selected instrumental response table). These data can be processed directly, or can be converted into the SolarSoft map object via the procedure<br/>
ConvertToMapsEUV, outspace, simbox, model, response, mapEUV<br/>
where the input parameters outspace, simbox, model, and response are the structures returned by the above-mentioned functions, and the output parameter mapEUV is the resulting SolarSoft (multi-channel) map object which represents the computed emission.

An example of using the code is given in the file /examples/RenderExampleEUV.pro (the sample GX Simulator model and EBTEL data are not included).

Optionally, for both the microwave and EUV emissions, selective heating of the coronal magnetic field lines can be defined via the SHtable parameter. This parameter is a 2D array (double precision floating point) with 7 * 7 elements. Each element of that table represents the factor applied to the heating rate Q for the field lines connecting specific regions at the photosphere; see the 'Selective Heating Mask' panel in GX Simulator. The SHtable table is supposed to be symmetric, i.e., SHtable[j, i]=SHtable[i, j]; asymmetric tables are accepted but the result will likely have no sense. By default (if the parameter SHtable is not supplied in calls to ComputeMW or ComputeEUV), all elements of that table are assumed to equal 1.
