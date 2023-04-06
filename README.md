This code computes the 2D maps of the solar gyroresonance and free-free microwave emission using the models of active regions created by the GX Simulator. The code is called from IDL (requires the SolarSoft package).

To compute the emission maps, you firstly need to create the input data blocks by calling the following functions:

1. Load the GX Simulator model:<br/>
   model=LoadGXmodel(modelfile)

   where modelfile is 
