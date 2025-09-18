#include "Rmain.h"
#include "EUVmain.h"

#ifdef WINDOWS
extern "C" __declspec(dllexport) int pyComputeMW(void *model, void *ebtel, void *simbox, void *cparms, void *out)
#else
extern "C" int pyComputeMW(void* model, void* ebtel, void* simbox, void* cparms, void* out)
#endif
{
 void *ARGV[5];

 ARGV[0]=model;
 ARGV[1]=ebtel;
 ARGV[2]=simbox;
 ARGV[3]=cparms;
 ARGV[4]=out;

 return ComputeMW(5, ARGV);
}

#ifdef WINDOWS
extern "C" __declspec(dllexport) int pyComputeMW_SH(void *model, void *ebtel, void *simbox, void *cparms, void *out, 
	                                                void *shtable)
#else
extern "C" int pyComputeMW_SH(void* model, void* ebtel, void* simbox, void* cparms, void* out, 
	                          void* shtable)
#endif
{
 void *ARGV[6];

 ARGV[0]=model;
 ARGV[1]=ebtel;
 ARGV[2]=simbox;
 ARGV[3]=cparms;
 ARGV[4]=out;
 ARGV[5]=shtable;

 return ComputeMW(6, ARGV);
}

#ifdef WINDOWS
extern "C" __declspec(dllexport) int pyComputeEUV(void *model, void *ebtel, void *response, void *simbox, void *cparms, void *out)
#else
extern "C" int pyComputeEUV(void* model, void* ebtel, void* response, void* simbox, void* cparms, void* out)
#endif
{
 void *ARGV[6];

 ARGV[0]=model;
 ARGV[1]=ebtel;
 ARGV[2]=response;
 ARGV[3]=simbox;
 ARGV[4]=cparms;
 ARGV[5]=out;

 return ComputeEUV(6, ARGV);
}

#ifdef WINDOWS
extern "C" __declspec(dllexport) int pyComputeEUV_SH(void *model, void *ebtel, void *response, void *simbox, void *cparms, void *out,
	                                                 void *shtable)
#else
extern "C" int pyComputeEUV_SH(void* model, void* ebtel, void* response, void* simbox, void* cparms, void* out, 
	                           void* shtable)
#endif
{
 void *ARGV[7];

 ARGV[0]=model;
 ARGV[1]=ebtel;
 ARGV[2]=response;
 ARGV[3]=simbox;
 ARGV[4]=cparms;
 ARGV[5]=out;
 ARGV[6]=shtable;

 return ComputeEUV(7, ARGV);
}