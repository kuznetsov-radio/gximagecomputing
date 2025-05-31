#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include "IDLinterface.h"
#include "Plasma.h"
#include "ExtMath.h"
#include "RenderIrregular.h"
#include "MWtransfer.h"
#include "GXdefs.h"

#ifdef WINDOWS
#include <ppl.h>
#include <concrtrm.h>
#else
#include <omp.h>
#define __int32 int32_t
#endif

#define InSize 15
#define OutSize 7

#ifdef WINDOWS
extern "C" __declspec(dllexport) int ComputeMW_fragment(int argc, void **argv)
#else
extern "C" double ComputeMW_fragment(int argc, void **argv)
#endif
{
 __int32 *m32=(__int32*)argv[0];
 int m_Nx=*(m32++);
 int m_Ny=*(m32++);
 int m_Nz=*(m32++); 
 int m_chromo_layers=*(m32++);
 int m_corona_layers=*(m32++); 
 int m_corona_base=*(m32++); 

 double *m64=(double*)m32;
 double m_DSun=*(m64++); 
 double m_RSun=*(m64++); 
 double m_b0Sun=*(m64++);
 double m_lonC=*(m64++);
 double m_latC=*(m64++); 
 double m_dx=*(m64++); 
 double m_dy=*(m64++); 
 double m_dz_uniform=*(m64++); 
 m64++; //skip obstime

 float *m_dz=(float*)(m64);
 float *m_Bx=m_dz+m_Nx*m_Ny*m_Nz;
 float *m_By=m_Bx+m_Nx*m_Ny*m_Nz;
 float *m_Bz=m_By+m_Nx*m_Ny*m_Nz;
 float *m_chromo_n0=m_Bz+m_Nx*m_Ny*m_Nz;
 float *m_chromo_np=m_chromo_n0+m_Nx*m_Ny*m_chromo_layers;
 float *m_chromo_nHI=m_chromo_np+m_Nx*m_Ny*m_chromo_layers;
 float *m_chromo_T0=m_chromo_nHI+m_Nx*m_Ny*m_chromo_layers;
 float *m_corona_Bavg=m_chromo_T0+m_Nx*m_Ny*m_chromo_layers;
 float *m_corona_L=m_corona_Bavg+m_Nx*m_Ny*m_corona_layers;
 float *m_chromo_uniform_Bavg=m_corona_L+m_Nx*m_Ny*m_corona_layers;
 float *m_chromo_uniform_L=m_chromo_uniform_Bavg+m_Nx*m_Ny*m_corona_base;

 char *m_VoxelID=(char*)(m_chromo_uniform_L+m_Nx*m_Ny*m_corona_base);

 char *m_corona_ID1=m_VoxelID+m_Nx*m_Ny*m_Nz;
 char *m_corona_ID2=m_corona_ID1+m_Nx*m_Ny*m_corona_layers;
 char *m_chromo_uniform_ID1=m_corona_ID2+m_Nx*m_Ny*m_corona_layers;
 char *m_chromo_uniform_ID2=m_chromo_uniform_ID1+m_Nx*m_Ny*m_corona_base;

 //-------------------------------------

 __int32 *e32=(__int32*)argv[1];
 int DEM_on=*(e32++);
 int DDM_on=*(e32++);

 int e_NQ, e_NL, e_NT;
 float *e_Qrun, *e_Lrun, *e_logtdem, *e_DEM_cor_run, *e_DDM_cor_run;
 e_NQ=e_NL=e_NT=0;
 e_Qrun=e_Lrun=e_logtdem=e_DEM_cor_run=e_DDM_cor_run=0;

 if (DEM_on || DDM_on)
 {
  e_NQ=*(e32++);
  e_NL=*(e32++);
  e_NT=*(e32++);

  e_Qrun=(float*)e32;
  e_Lrun=e_Qrun+e_NQ*e_NL;
  e_logtdem=e_Lrun+e_NQ*e_NL;

  float *e_cor=e_logtdem+e_NT;

  if (DEM_on)
  {
   e_DEM_cor_run=e_cor;
   e_cor+=e_NQ*e_NL*e_NT;
  }

  if (DDM_on) e_DDM_cor_run=e_cor;
 }

 //-------------------------------------

 __int32 *b32=(__int32*)argv[2]; 
 int b_Nx=*(b32++); 
 int b_Ny=*(b32++); 
 int b_Nf=*(b32++); 
 int b_projection=*(b32++);

 double *b64=(double*)b32; 
 double b_xc=*(b64++);
 double b_yc=*(b64++);
 double b_dx=*(b64++);
 double b_dy=*(b64++);
 double b_rot=*(b64++);

 double *b_freqlist=b64;

 int ProjectionParallel=(b_projection & 1)!=0;
 int ProjectionExact=(b_projection & 2)!=0;

 //-------------------------------------

 double *cp64=(double*)argv[3];
 double cp_Tbase=*(cp64++);
 double cp_nbase=*(cp64++);
 double cp_Q0=*(cp64++);
 double cp_a=*(cp64++);
 double cp_b=*(cp64++);

 __int32 *cp32=(__int32*)cp64;
 int cp_mode=*cp32;
 int force_isothermal=(cp_mode & 1)!=0;
 int interpolB=(cp_mode & 2)!=0;
 int aNT=(cp_mode & 4)!=0;

 //-------------------------------------

 __int32 *o_flagsAll=(__int32*)argv[4];
 __int32 *o_flagsCorona=o_flagsAll+6;

 double *o_TI=(double*)(o_flagsCorona+6);
 double *o_TV=o_TI+b_Nx*b_Ny*b_Nf;

 //-------------------------------------

 double *SHtable=(double*)argv[5];

 //-------------------------------------

 int *fragment=(int*)argv[6];
 int i_start=fragment[0];
 int i_end=fragment[1];
 int j_start=fragment[2];
 int j_end=fragment[3];

 char *flagsGlobal=(char*)argv[7];

 #ifdef WINDOWS
 concurrency::critical_section *cs=(concurrency::critical_section*)argv[8];
 #endif

 //-------------------------------------

 double srot=sin(b_rot/180*M_PI);
 double crot=cos(b_rot/180*M_PI);
 double *wx=(double*)malloc(b_Nx*b_Ny*sizeof(double));
 double *wy=(double*)malloc(b_Nx*b_Ny*sizeof(double));

 for (int i=0; i<b_Nx; i++) for (int j=0; j<b_Ny; j++)
 {
  double x1=-b_dx*(b_Nx-1)/2+i*b_dx;
  double y1=-b_dy*(b_Ny-1)/2+j*b_dy;
  wx[D2(b_Nx, i, j)]=b_xc+x1*crot-y1*srot;
  wy[D2(b_Nx, i, j)]=b_yc+x1*srot+y1*crot;
 }

 double *I_L=(double*)malloc(b_Nx*b_Ny*b_Nf*sizeof(double));
 double *I_R=(double*)malloc(b_Nx*b_Ny*b_Nf*sizeof(double));
 memset(I_L, 0, b_Nx*b_Ny*b_Nf*sizeof(double));
 memset(I_R, 0, b_Nx*b_Ny*b_Nf*sizeof(double));

 double *dz=(double*)malloc(m_Nx*m_Ny*m_Nz*sizeof(double));
 double *h=(double*)malloc(m_Nx*m_Ny*m_Nz*sizeof(double));

 for (int i=0; i<m_Nx*m_Ny*m_Nz; i++) dz[i]=(double)m_dz[i];

 for (int i=0; i<m_Nx; i++) for (int j=0; j<m_Ny; j++)
 {
  h[D3(m_Nx, m_Ny, i, j, 0)]=dz[D3(m_Nx, m_Ny, i, j, 0)];
  for (int k=1; k<m_Nz; k++) h[D3(m_Nx, m_Ny, i, j, k)]=h[D3(m_Nx, m_Ny, i, j, k-1)]+dz[D3(m_Nx, m_Ny, i, j, k)]; //cumulative sum

  for (int k=0; k<m_Nz; k++) h[D3(m_Nx, m_Ny, i, j, k)]-=dz[D3(m_Nx, m_Ny, i, j, k)]/2;
 }

 double z1=-m_RSun;
 double z2=m_RSun*2;
 double D_exact=sqrt(sqr(m_RSun*sin(m_latC*M_PI/180))+
	                 sqr(m_RSun*cos(m_latC*M_PI/180)*sin(m_lonC*M_PI/180))+
	                 sqr(m_RSun*cos(m_latC*M_PI/180)*cos(m_lonC*M_PI/180)-m_DSun));

 double rb[3]; //the bottom front left box boundary
 rb[0]=-m_dx*(m_Nx-1)/2;
 rb[1]=-m_dy*(m_Ny-1)/2;
 rb[2]=m_RSun;

 double *LgridL, *QgridL, *Tgrid;
 LgridL=QgridL=Tgrid=0;

 if (DEM_on || DDM_on)
 {
  LgridL=(double*)malloc(e_NL*sizeof(double));
  for (int j=0; j<e_NL; j++) LgridL[j]=log((double)e_Lrun[D2(e_NQ, 0, j)]);

  QgridL=(double*)malloc(e_NQ*e_NL*sizeof(double));
  for (int i=0; i<e_NQ*e_NL; i++) QgridL[i]=log((double)e_Qrun[i]);

  Tgrid=(double*)malloc(e_NT*sizeof(double));
  for (int i=0; i<e_NT; i++) Tgrid[i]=pow(10.0, (double)e_logtdem[i]);
 }

 double H_corona=corona_Hscale*cp_Tbase;

 char *flags=(char*)malloc(m_Nx*m_Ny*m_Nz*sizeof(char));
 memset(flags, 0, m_Nx*m_Ny*m_Nz*sizeof(char));

 int Nvoxels=0;
 int *VoxList=(int*)malloc(arrN*sizeof(int));
 double *ds=(double*)malloc(arrN*sizeof(double));
 double *xmid=(double*)malloc(arrN*sizeof(double));
 double *ymid=(double*)malloc(arrN*sizeof(double));
 double *zmid=(double*)malloc(arrN*sizeof(double));

 int Rdim[3];
 Rdim[0]=m_Nx;
 Rdim[1]=m_Ny;
 Rdim[2]=m_Nz;

 double Rdxdy[2];
 Rdxdy[0]=m_dx;
 Rdxdy[1]=m_dy;

 int Lparms[5]={0, 0, 0, 0, 0};
 Lparms[1]=b_Nf;
 Lparms[2]=(DEM_on || DDM_on) ? e_NT : 0;

 double Rparms[3]={0, 0, 0};
 Rparms[0]=b_dx*b_dy;

 int NvoxMax=300;
 double *Parms=(double*)malloc(InSize*NvoxMax*sizeof(double));
 double *DEM_arr=(DEM_on) ? (double*)malloc(e_NT*NvoxMax*sizeof(double)) : 0;
 double *DDM_arr=(DDM_on) ? (double*)malloc(e_NT*NvoxMax*sizeof(double)) : 0;

 double *RL=(double*)malloc(OutSize*b_Nf*sizeof(double));
 
 void *ARGV[11];

 for (int i=i_start; i<=i_end; i++) for (int j=j_start; j<=j_end; j++)
 {
  double r1[3], r2[3], LOS[3], norm_x[3], norm_y[3];

  double spx=sin(M_PI/648000*wx[D2(b_Nx, i, j)]);
  double spy=sin(M_PI/648000*wy[D2(b_Nx, i, j)]);
  double q=sqrt(1.0-sqr(spx)-sqr(spy));
  double xD=spx/q;
  double yD=spy/q;

  if (!ProjectionParallel)
  {
   r1[0]=(m_DSun-z1)*xD;
   r1[1]=(m_DSun-z1)*yD;
   r1[2]=z1;

   r2[0]=(m_DSun-z2)*xD;
   r2[1]=(m_DSun-z2)*yD;
   r2[2]=z2;

   LOS[0]=(z1-z2)*xD;
   LOS[1]=(z1-z2)*yD;
   LOS[2]=z2-z1;
  }
  else
  {
   double D=ProjectionExact ? D_exact : m_DSun;

   r1[0]=r2[0]=D*wx[D2(b_Nx, i, j)]*M_PI/648000;
   r1[1]=r2[1]=D*wy[D2(b_Nx, i, j)]*M_PI/648000;
   r1[2]=z1;
   r2[2]=z2;

   LOS[0]=LOS[1]=0;
   LOS[2]=z2-z1;
  }

  double aLOS=sqrt(sqr(LOS[0])+sqr(LOS[1])+sqr(LOS[2]));
  for (int l=0; l<3; l++) LOS[l]/=aLOS;

  double jvec[3]={0, 1, 0};
  CrossP(jvec, LOS, norm_x);
  CrossP(LOS, norm_x, norm_y);

  rotC(r1, m_latC, m_lonC, m_b0Sun);
  rotC(r2, m_latC, m_lonC, m_b0Sun);
  rotC(LOS, m_latC, m_lonC, m_b0Sun);
  rotC(norm_x, m_latC, m_lonC, m_b0Sun);
  rotC(norm_y, m_latC, m_lonC, m_b0Sun); 

  for (int l=0; l<3; l++)
  {
   r1[l]-=rb[l];
   r2[l]-=rb[l];
  }

  ARGV[0]=(void*)Rdim;
  ARGV[1]=(void*)Rdxdy;
  ARGV[2]=(void*)dz;
  ARGV[3]=(void*)r1;
  ARGV[4]=(void*)r2;
  ARGV[5]=(void*)(&Nvoxels);
  ARGV[6]=(void*)VoxList;
  ARGV[7]=(void*)ds;
  ARGV[8]=(void*)xmid;
  ARGV[9]=(void*)ymid;
  ARGV[10]=(void*)zmid;

  int res=RENDER(11, ARGV);

  if (Nvoxels>0)
  {
   Lparms[0]=Nvoxels;

   if (Nvoxels>NvoxMax)
   {
	NvoxMax=Nvoxels+100;
	Parms=(double*)realloc(Parms, InSize*NvoxMax*sizeof(double));
	if (DEM_on) DEM_arr=(double*)realloc(DEM_arr, e_NT*NvoxMax*sizeof(double));
	if (DDM_on) DDM_arr=(double*)realloc(DDM_arr, e_NT*NvoxMax*sizeof(double));
   }

   memset(Parms, 0, InSize*Nvoxels*sizeof(double));

   for (int k=0; k<Nvoxels; k++)
   {
	flags[VoxList[k]]|=1; //voxels crossed by LOSs

	Parms[D2(InSize, 0, k)]=ds[k]; //voxel depths
	Parms[D2(InSize, 1, k)]=cp_Tbase; //default temperature
	Parms[D2(InSize, 2, k)]=cp_nbase*exp(-h[VoxList[k]]/H_corona); //default density
	Parms[D2(InSize, 6, k)]=(force_isothermal) ? 8 : 0; //turn isothermal approximation on/off
	Parms[D2(InSize, 7, k)]=10; //max harmonic number
	Parms[D2(InSize, 11, k)]=1; //DEM off
	Parms[D2(InSize, 12, k)]=1; //DDM off

	double Bavg, Lline;
	Bavg=Lline=0;

    int ID1, ID2;
    ID1=ID2=0;

	int idx_i=VoxList[k] % m_Nx;
	int idx_j=(VoxList[k]/m_Nx) % m_Ny;
	int idx_k=(VoxList[k]/m_Nx)/m_Ny;

	if (idx_k<m_chromo_layers)
	{
	 if (m_chromo_T0[D3(m_Nx, m_Ny, idx_i, idx_j, idx_k)]>0)
	 {
	  Parms[D2(InSize, 1, k)]=m_chromo_T0[D3(m_Nx, m_Ny, idx_i, idx_j, idx_k)]; //temperature
	  Parms[D2(InSize, 2, k)]=m_chromo_n0[D3(m_Nx, m_Ny, idx_i, idx_j, idx_k)]; //total density
	  Parms[D2(InSize, 8, k)]=m_chromo_np[D3(m_Nx, m_Ny, idx_i, idx_j, idx_k)]; //proton density
	  Parms[D2(InSize, 9, k)]=m_chromo_nHI[D3(m_Nx, m_Ny, idx_i, idx_j, idx_k)]; //neutral hydrogen density
	  Parms[D2(InSize, 13, k)]=2; //photospheric abundance by Scott
	  flags[D3(m_Nx, m_Ny, idx_i, idx_j, idx_k)]|=2; //chromosphere voxels
	 }
	 else
	 {
	  int l=int(h[VoxList[k]]/m_dz_uniform);
	  Bavg=m_chromo_uniform_Bavg[D3(m_Nx, m_Ny, idx_i, idx_j, l)];
	  Lline=m_chromo_uniform_L[D3(m_Nx, m_Ny, idx_i, idx_j, l)];
	  ID1=m_chromo_uniform_ID1[D3(m_Nx, m_Ny, idx_i, idx_j, l)];
	  ID2=m_chromo_uniform_ID2[D3(m_Nx, m_Ny, idx_i, idx_j, l)];
	 }
	}
	else
	{
	 Bavg=m_corona_Bavg[D3(m_Nx, m_Ny, idx_i, idx_j, idx_k-m_chromo_layers)];
	 Lline=m_corona_L[D3(m_Nx, m_Ny, idx_i, idx_j, idx_k-m_chromo_layers)];
	 ID1=m_corona_ID1[D3(m_Nx, m_Ny, idx_i, idx_j, idx_k-m_chromo_layers)];
	 ID2=m_corona_ID2[D3(m_Nx, m_Ny, idx_i, idx_j, idx_k-m_chromo_layers)];
	}

	if (Bavg>0 && (DEM_on || DDM_on || aNT))
	{
	 flags[VoxList[k]]|=4; //voxels with B and L known

	 Lline/=2; //switch from line length to line half-length, for consistency with GX Simulator

	 double Q=cp_Q0*pow(Bavg/Bavg0, cp_a)/pow(Lline/Lline0, cp_b);
	 if (SHtable) Q*=SHtable[D2(SHsize, ID1-1, ID2-1)];
	 
	 if (DEM_on || DDM_on)
	 {
	  double Qlog=log(Q);
	  double Llog=log(Lline);

	  int EBTEL_hit=0;

	  int Lind=value_locate(LgridL, e_NL, Llog);

	  if (Lind>=0 && Lind<(e_NL-1))
	  {
	   int Qind1=value_locate(QgridL+Lind*e_NQ, e_NQ, Qlog);
	   int Qind2=value_locate(QgridL+(Lind+1)*e_NQ, e_NQ, Qlog);

	   if (Qind1>=0 && Qind1<(e_NQ-1) && Qind2>=0 && Qind2<(e_NQ-1))
	   {
	    flags[VoxList[k]]|=8; //EBTEL table hit

	    EBTEL_hit=1;

	    double dL=(Llog-LgridL[Lind])/(LgridL[Lind+1]-LgridL[Lind]);
        double dQ1=(Qlog-QgridL[D2(e_NQ, Qind1, Lind)])/(QgridL[D2(e_NQ, Qind1+1, Lind)]-QgridL[D2(e_NQ, Qind1, Lind)]);
        double dQ2=(Qlog-QgridL[D2(e_NQ, Qind2, Lind+1)])/(QgridL[D2(e_NQ, Qind2+1, Lind+1)]-QgridL[D2(e_NQ, Qind2, Lind+1)]);

	    if (DEM_on)
	    {
		 for (int l=0; l<e_NT; l++)
		  DEM_arr[D2(e_NT, l, k)]=e_DEM_cor_run[D3(e_NT, e_NQ, l, Qind1, Lind)]*(1.0-dL)*(1.0-dQ1)+
                                  e_DEM_cor_run[D3(e_NT, e_NQ, l, Qind1+1, Lind)]*(1.0-dL)*dQ1+
                                  e_DEM_cor_run[D3(e_NT, e_NQ, l, Qind2, Lind+1)]*dL*(1.0-dQ2)+
                                  e_DEM_cor_run[D3(e_NT, e_NQ, l, Qind2+1, Lind+1)]*dL*dQ2;
		 Parms[D2(InSize, 11, k)]=0; //DEM on
	    }

	    if (DDM_on)
	    {
		 for (int l=0; l<e_NT; l++)
		  DDM_arr[D2(e_NT, l, k)]=e_DDM_cor_run[D3(e_NT, e_NQ, l, Qind1, Lind)]*(1.0-dL)*(1.0-dQ1)+
                                  e_DDM_cor_run[D3(e_NT, e_NQ, l, Qind1+1, Lind)]*(1.0-dL)*dQ1+
                                  e_DDM_cor_run[D3(e_NT, e_NQ, l, Qind2, Lind+1)]*dL*(1.0-dQ2)+
                                  e_DDM_cor_run[D3(e_NT, e_NQ, l, Qind2+1, Lind+1)]*dL*dQ2;
		 Parms[D2(InSize, 12, k)]=0; //DDM on
	    }
	   }
	   else flags[VoxList[k]]|=32; //EBTEL table miss (Q) 
	  }
	  else flags[VoxList[k]]|=16; //EBTEL table miss (L)

	  if (!EBTEL_hit && aNT)
	  {
	   double n0a, T0a;
	   FindAnalyticalNT(Q, Lline, &n0a, &T0a);
	   Parms[D2(InSize, 1, k)]=T0a; 
	   Parms[D2(InSize, 2, k)]=n0a; 
	  }
	 }
	 else if (aNT)
	 {
	  double n0a, T0a;
	  FindAnalyticalNT(Q, Lline, &n0a, &T0a);
	  Parms[D2(InSize, 1, k)]=T0a; 
	  Parms[D2(InSize, 2, k)]=n0a; 
	 }
	}

	double Bx, By, Bz;
	if (interpolB) InterpolateTrilinear(m_Nx, m_Ny, m_Nz, m_dx, m_dy, h, 
                                        m_Bx, m_By, m_Bz, xmid[k], ymid[k], zmid[k],
		                                &Bx, &By, &Bz);
	else
	{
	 Bx=m_Bx[VoxList[k]];
	 By=m_By[VoxList[k]];
	 Bz=m_Bz[VoxList[k]];
	}

	Parms[D2(InSize, 3, k)]=sqrt(sqr(Bx)+sqr(By)+sqr(Bz)); //magnetic field
	Parms[D2(InSize, 4, k)]=(Parms[D2(InSize, 3, k)]>0) ?
		                    acos((Bx*LOS[0]+By*LOS[1]+Bz*LOS[2])/Parms[D2(InSize, 3, k)])/M_PI*180 : 0; //viewing angle
    Parms[D2(InSize, 5, k)]=(Parms[D2(InSize, 3, k)]>0) ?
		                    atan2(Bx*norm_y[0]+By*norm_y[1]+Bz*norm_y[2], 
                                  Bx*norm_x[0]+By*norm_x[1]+Bz*norm_x[2])/M_PI*180 : 0; //azimuthal angle
   }

   memset(RL, 0, OutSize*b_Nf*sizeof(double));
   for (int l=0; l<b_Nf; l++) RL[D2(OutSize, 0, l)]=b_freqlist[l];

   ARGV[0]=(void*)Lparms;
   ARGV[1]=(void*)Rparms;
   ARGV[2]=(void*)Parms;
   ARGV[3]=(DEM_on || DDM_on) ? (void*)Tgrid : 0;
   ARGV[4]=(DEM_on) ? (void*)DEM_arr : 0;
   ARGV[5]=(DDM_on) ? (void*)DDM_arr : 0;
   ARGV[6]=(void*)RL;

   res=GET_MW(7, ARGV);

   for (int l=0; l<b_Nf; l++)
   {
	I_L[D3(b_Nx, b_Ny, i, j, l)]=RL[D2(OutSize, 5, l)];
	I_R[D3(b_Nx, b_Ny, i, j, l)]=RL[D2(OutSize, 6, l)];
   }
  }
 }

 for (int k=0; k<b_Nf; k++)
 {
  double r=sfu*c*c/(2.0*kB*sqr(b_freqlist[k]*1e9))/Rparms[0]*AU*AU;

  for (int i=i_start; i<=i_end; i++) for (int j=j_start; j<=j_end; j++)
  {
   o_TI[D3(b_Nx, b_Ny, i, j, k)]=(I_R[D3(b_Nx, b_Ny, i, j, k)]+I_L[D3(b_Nx, b_Ny, i, j, k)])*r;
   o_TV[D3(b_Nx, b_Ny, i, j, k)]=(I_R[D3(b_Nx, b_Ny, i, j, k)]-I_L[D3(b_Nx, b_Ny, i, j, k)])*r;
  }
 }

 #ifdef WINDOWS
 if (cs) cs->lock();
 for (int i=0; i<m_Nx*m_Ny*m_Nz; i++) flagsGlobal[i]|=flags[i];
 if (cs) cs->unlock();
 #else
 #pragma omp critical (UpdateFlags)
 for (int i=0; i<m_Nx*m_Ny*m_Nz; i++) flagsGlobal[i]|=flags[i];
 #endif

 free(RL);
 if (DDM_arr) free(DDM_arr);
 if (DEM_arr) free(DEM_arr);
 free(Parms);
 free(zmid);
 free(ymid);
 free(xmid);
 free(ds);
 free(VoxList);
 free(flags);
 if (Tgrid) free(Tgrid);
 if (QgridL) free(QgridL);
 if (LgridL) free(LgridL);
 free(h);
 free(dz);
 free(I_R);
 free(I_L);
 free(wy);
 free(wx);

 return 0;
}

#ifdef WINDOWS
extern "C" __declspec(dllexport) int ComputeMW(int argc, void **argv)
#else
extern "C" int ComputeMW(int argc, void **argv)
#endif
{
 __int32 *m32=(__int32*)argv[0];
 int m_Nx=*(m32++);
 int m_Ny=*(m32++);
 int m_Nz=*(m32++); 
 int m_chromo_layers=*(m32++);

 __int32 *b32=(__int32*)argv[2]; 
 int b_Nx=*(b32++); 
 int b_Ny=*(b32++); 
 b32++;
 int projection=*(b32++);
 int Nthreads=projection>>16;
                     
 __int32 *o_flagsAll=(__int32*)argv[4];
 __int32 *o_flagsCorona=o_flagsAll+6;

 char *flags=(char*)malloc(m_Nx*m_Ny*m_Nz*sizeof(char));
 memset(flags, 0, m_Nx*m_Ny*m_Nz*sizeof(char));

 void *SHtable=(argc>5) ? argv[5] : 0;

 #ifdef WINDOWS
 concurrency::critical_section cs;

 int NtMax=concurrency::GetProcessorCount();
 if (Nthreads>0) NtMax=min(Nthreads, NtMax);
 int K=int(b_Ny/NtMax);
 int Np=b_Ny % NtMax;
 int Nm=NtMax-Np;
 
 concurrency::parallel_for(0, NtMax, [&](int j)
 {
  void *ARGV[9];
  memcpy(ARGV, argv, sizeof(void*)*5);
  ARGV[5]=SHtable;

  int fragment[4];
  ARGV[6]=(void*)fragment;

  ARGV[7]=(void*)flags;
  ARGV[8]=(void*)&cs;

  fragment[0]=0;
  fragment[1]=b_Nx-1;

  fragment[2]=(j<Nm) ? K*j : K*Nm+(K+1)*(j-Nm);
  fragment[3]=fragment[2]+((j<Nm) ? K : K+1)-1;

  ComputeMW_fragment(9, ARGV);
 });
 #else
 int NtMax=omp_get_max_threads();
 if (Nthreads>NtMax) Nthreads=NtMax;
 if (Nthreads>0) omp_set_num_threads(Nthreads);

 #pragma omp parallel for
 for (int j=0; j<b_Ny; j++)
 {
  void *ARGV[9];
  memcpy(ARGV, argv, sizeof(void*)*5);
  ARGV[5]=SHtable;

  int fragment[4];
  ARGV[6]=(void*)fragment;

  ARGV[7]=(void*)flags;

  fragment[0]=0;
  fragment[1]=b_Nx-1;
  fragment[2]=fragment[3]=j;

  ComputeMW_fragment(9, ARGV);
 }
 #endif

 for (int i=0; i<m_Nx; i++) for (int j=0; j<m_Ny; j++) for (int k=0; k<m_Nz; k++) for (int m=0; m<6; m++)
 {
  if ((flags[D3(m_Nx, m_Ny, i, j, k)] & (char(1)<<m))!=0)
  {
   o_flagsAll[m]++;
   if (k>=m_chromo_layers) o_flagsCorona[m]++;
  }
 }

 free(flags);

 return 0;
}

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