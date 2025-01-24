#include <stdlib.h>
#include <math.h>
#include "ExtMath.h"
#include "IDLinterface.h"

#ifndef LINUX
#include <ppl.h>
#else
#include <omp.h>
#define __int32 int32_t
#endif

#ifndef LINUX
extern "C" __declspec(dllexport) int InterpolateEBTEL(int argc, void **argv)
#else
extern "C" int InterpolateEBTEL(int argc, void **argv)
#endif
{
 int *Lparms=(int*)argv[0];
 int Npoints=Lparms[0];
 int NQ=Lparms[1];
 int NL=Lparms[2];
 int NT=Lparms[3];
 int DEM_on=Lparms[4];
 int DDM_on=Lparms[5];

 float *Qrun=(float*)argv[1];
 float *Lrun=(float*)argv[2];
 float *logtdem=(float*)argv[3];
 float *DEM_run=(float*)argv[4];
 float *DDM_run=(float*)argv[5];

 double *Qarr=(double*)argv[6];
 double *Larr=(double*)argv[7];

 double *DEM_arr=(double*)argv[8];
 double *DDM_arr=(double*)argv[9];

 double *n_DEM=(double*)argv[10];
 double *T_DEM=(double*)argv[11];
 double *n_DDM=(double*)argv[12];
 double *T_DDM=(double*)argv[13];

 char *flag=(char*)argv[14];

 double *LgridL=(double*)malloc(NL*sizeof(double));
 for (int j=0; j<NL; j++) LgridL[j]=log((double)Lrun[D2(NQ, 0, j)]);

 double *QgridL=(double*)malloc(NQ*NL*sizeof(double));
 for (int i=0; i<NQ*NL; i++) QgridL[i]=log((double)Qrun[i]);

 double *Tarr=(double*)malloc(NT*sizeof(double));
 double *lnTarr=(double*)malloc(NT*sizeof(double));
 double ln10=log(10.0);
 for (int l=0; l<NT; l++) 
 {
  lnTarr[l]=logtdem[l]*ln10;
  Tarr[l]=exp(lnTarr[l]);
 }

 #ifndef LINUX
 concurrency::parallel_for(0, Npoints, [&](int k)
 {
 #else
 #pragma omp parallel for
 for (int k=0; k<Npoints; k++)
 {
 #endif

  flag[k]=0;

  double Qlog=log(Qarr[k]);
  double Llog=log(Larr[k]);

  double *y1=(double*)malloc(NT*sizeof(double));
  double *y2=(double*)malloc(NT*sizeof(double));

  int Lind=value_locate(LgridL, NL, Llog);

  if (Lind>=0 && Lind<(NL-1))
  {
   int Qind1=value_locate(QgridL+Lind*NQ, NQ, Qlog);
   int Qind2=value_locate(QgridL+(Lind+1)*NQ, NQ, Qlog);

   if (Qind1>=0 && Qind1<(NQ-1) && Qind2>=0 && Qind2<(NQ-1))
   {
	flag[k]=1;

	double dL=(Llog-LgridL[Lind])/(LgridL[Lind+1]-LgridL[Lind]);
    double dQ1=(Qlog-QgridL[D2(NQ, Qind1, Lind)])/(QgridL[D2(NQ, Qind1+1, Lind)]-QgridL[D2(NQ, Qind1, Lind)]);
    double dQ2=(Qlog-QgridL[D2(NQ, Qind2, Lind+1)])/(QgridL[D2(NQ, Qind2+1, Lind+1)]-QgridL[D2(NQ, Qind2, Lind+1)]);

	if (DEM_on)
    {
	 for (int l=0; l<NT; l++)
     {
	  double DEM=DEM_arr[D2(NT, l, k)]=DEM_run[D3(NT, NQ, l, Qind1, Lind)]*(1.0-dL)*(1.0-dQ1)+
                                       DEM_run[D3(NT, NQ, l, Qind1+1, Lind)]*(1.0-dL)*dQ1+
                                       DEM_run[D3(NT, NQ, l, Qind2, Lind+1)]*dL*(1.0-dQ2)+
                                       DEM_run[D3(NT, NQ, l, Qind2+1, Lind+1)]*dL*dQ2;
      y1[l]=DEM*Tarr[l];
      y2[l]=y1[l]*Tarr[l];
     }

     double n2_avg=IntTabulated(lnTarr, y1, NT);
     n_DEM[k]=sqrt(n2_avg);
     T_DEM[k]=(n2_avg>0) ? IntTabulated(lnTarr, y2, NT)/n2_avg : 0.0;
    }

	if (DDM_on)
    {
 	 for (int l=0; l<NT; l++)
     {
	  double DDM=DDM_arr[D2(NT, l, k)]=DDM_run[D3(NT, NQ, l, Qind1, Lind)]*(1.0-dL)*(1.0-dQ1)+
                                       DDM_run[D3(NT, NQ, l, Qind1+1, Lind)]*(1.0-dL)*dQ1+
                                       DDM_run[D3(NT, NQ, l, Qind2, Lind+1)]*dL*(1.0-dQ2)+
                                       DDM_run[D3(NT, NQ, l, Qind2+1, Lind+1)]*dL*dQ2; 
      y1[l]=DDM*Tarr[l];
      y2[l]=y1[l]*Tarr[l];
     }

     double n_avg=IntTabulated(lnTarr, y1, NT);
     n_DDM[k]=n_avg;
     T_DDM[k]=(n_avg>0) ? IntTabulated(lnTarr, y2, NT)/n_avg : 0.0;
    }
   }
  }

  free(y2);
  free(y1);

 #ifndef LINUX
 });
 #else
 }
 #endif

 free(lnTarr);
 free(Tarr);
 free(LgridL);
 free(QgridL);

 return 0;
}