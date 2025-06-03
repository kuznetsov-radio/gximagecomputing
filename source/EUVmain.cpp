#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include "IDLinterface.h"
#include "ExtMath.h"
#include "RenderIrregular.h"
#include "Plasma.h"
#include "GXdefs.h"

#ifdef WINDOWS
#include <ppl.h>
#include <concrtrm.h>
#else
#include <omp.h>
#define __int32 int32_t
#endif

#ifdef WINDOWS
extern "C" __declspec(dllexport) int ComputeEUV_fragment(int argc, void **argv)
#else
extern "C" int ComputeEUV_fragment(int argc, void **argv)
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
 float *e_Qrun, *e_Lrun, *e_logtdem, *e_DEM_cor_run, *e_DDM_cor_run, *e_DEM_tr_run, *e_DDM_tr_run;
 e_NQ=e_NL=e_NT=0;
 e_Qrun=e_Lrun=e_logtdem=e_DEM_cor_run=e_DDM_cor_run=e_DEM_tr_run=e_DDM_tr_run=0;

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

  if (DDM_on) 
  {
   e_DDM_cor_run=e_cor;
   e_cor+=e_NQ*e_NL*e_NT;
  }

  if (DEM_on)
  {
   e_DEM_tr_run=e_cor;
   e_cor+=e_NQ*e_NL*e_NT;
  }

  if (DDM_on) e_DDM_tr_run=e_cor;
 }

 //-------------------------------------

 double *rs64=(double*)argv[2];

 double rs_ds=*(rs64++);

 __int32 *rs32=(__int32*)rs64;

 int rs_NT=*(rs32++);
 int rs_Nch=*(rs32++);

 double *rs_logte=(double*)rs32;
 double *rs_all=rs_logte+rs_NT;
        
 //-------------------------------------

 __int32 *b32=(__int32*)argv[3]; 
 int b_Nx=*(b32++); 
 int b_Ny=*(b32++); 

 double *b64=(double*)b32; 
 double b_xc=*(b64++);
 double b_yc=*(b64++);
 double b_dx=*(b64++);
 double b_dy=*(b64++);

 b32=(__int32*)b64;
 int b_projection=*(b32++);
 int ProjectionParallel=(b_projection & 1)!=0;
 int ProjectionExact=(b_projection & 2)!=0;
         
 //-------------------------------------

 double *cp64=(double*)argv[4];
 double cp_Tbase=*(cp64++);
 double cp_nbase=*(cp64++);
 double cp_Q0=*(cp64++);
 double cp_a=*(cp64++);
 double cp_b=*(cp64++);

 __int32 *cp32=(__int32*)cp64;
 int cp_mode=*cp32;
 int force_isothermal=(cp_mode & 1)!=0;
 int aNT=(cp_mode & 4)!=0;

 //-------------------------------------

 __int32 *o_flagsAll=(__int32*)argv[5];
 __int32 *o_flagsCorona=o_flagsAll+6;

 double *o_fluxCorona=(double*)(o_flagsCorona+6);
 double *o_fluxTR=o_fluxCorona+b_Nx*b_Ny*rs_Nch;

 //-------------------------------------

 double *SHtable=(double*)argv[6];

 //-------------------------------------

 int *fragment=(int*)argv[7];
 int i_start=fragment[0];
 int i_end=fragment[1];
 int j_start=fragment[2];
 int j_end=fragment[3];

 char *flagsGlobal=(char*)argv[8];

 #ifdef WINDOWS
 concurrency::critical_section *cs=(concurrency::critical_section*)argv[9];
 #endif

 //-------------------------------------

 double *wx=(double*)malloc(b_Nx*sizeof(double));
 wx[0]=b_xc-b_dx*(b_Nx-1)/2;
 for (int i=1; i<b_Nx; i++) wx[i]=wx[i-1]+b_dx;

 double *wy=(double*)malloc(b_Ny*sizeof(double));
 wy[0]=b_yc-b_dy*(b_Ny-1)/2;
 for (int j=1; j<b_Ny; j++) wy[j]=wy[j-1]+b_dy;

 double *fluxCorona=(double*)malloc(b_Nx*b_Ny*rs_Nch*sizeof(double));
 double *fluxTR=(double*)malloc(b_Nx*b_Ny*rs_Nch*sizeof(double));
 memset(fluxCorona, 0, b_Nx*b_Ny*rs_Nch*sizeof(double));
 memset(fluxTR, 0, b_Nx*b_Ny*rs_Nch*sizeof(double));

 double *dz, *h;

 if (argc>=12)
 {
  dz=(double*)argv[10];
  h=(double*)argv[11];
 }
 else
 {
  dz=(double*)malloc(m_Nx*m_Ny*m_Nz*sizeof(double));
  h=(double*)malloc(m_Nx*m_Ny*m_Nz*sizeof(double));

  for (int i=0; i<m_Nx*m_Ny*m_Nz; i++) dz[i]=(double)m_dz[i];

  for (int i=0; i<m_Nx; i++) for (int j=0; j<m_Ny; j++)
  {
   h[D3(m_Nx, m_Ny, i, j, 0)]=dz[D3(m_Nx, m_Ny, i, j, 0)];
   for (int k=1; k<m_Nz; k++) h[D3(m_Nx, m_Ny, i, j, k)]=h[D3(m_Nx, m_Ny, i, j, k-1)]+dz[D3(m_Nx, m_Ny, i, j, k)]; //cumulative sum

   for (int k=0; k<m_Nz; k++) h[D3(m_Nx, m_Ny, i, j, k)]-=dz[D3(m_Nx, m_Ny, i, j, k)]/2;
  }
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

 int rs_NT1=0;
 int DEM_idx1=0;
 double *LgridL, *QgridL, *Tgrid, *DEM_local_corona, *DEM_local_TR, *rs_logte1, *rs_te1, *rs_all1, *EUV_integrand;
 LgridL=QgridL=Tgrid=DEM_local_corona=DEM_local_TR=rs_logte1=rs_te1=rs_all1=EUV_integrand=0;

 Spline **rs_spl_arr=(Spline**)malloc(rs_Nch*sizeof(Spline*));
 for (int j=0; j<rs_Nch; j++) rs_spl_arr[j]=new Spline(rs_NT, rs_logte, rs_all+j*rs_NT);

 if (DEM_on)
 {
  LgridL=(double*)malloc(e_NL*sizeof(double));
  for (int j=0; j<e_NL; j++) LgridL[j]=log((double)e_Lrun[D2(e_NQ, 0, j)]);

  QgridL=(double*)malloc(e_NQ*e_NL*sizeof(double));
  for (int i=0; i<e_NQ*e_NL; i++) QgridL[i]=log((double)e_Qrun[i]);

  Tgrid=(double*)malloc(e_NT*sizeof(double));
  for (int i=0; i<e_NT; i++) Tgrid[i]=pow(10.0, (double)e_logtdem[i]);

  DEM_local_corona=(double*)malloc(e_NT*sizeof(double));
  DEM_local_TR=(double*)malloc(e_NT*sizeof(double));

  double Tmin=max(rs_logte[0], e_logtdem[0]);
  double Tmax=min(rs_logte[rs_NT-1], e_logtdem[e_NT-1]);

  int i1=-1;
  int i2=0;
  for (int i=0; i<e_NT; i++)
  {
   if ((e_logtdem[i]>=Tmin) && (i1<0)) i1=i;
   if (e_logtdem[i]<=Tmax) i2=i;
  }
  rs_NT1=i2-i1+1;
  DEM_idx1=i1;
  rs_logte1=(double*)malloc(rs_NT1*sizeof(double));
  rs_te1=(double*)malloc(rs_NT1*sizeof(double));
  for (int i=0; i<rs_NT1; i++) 
  {
   rs_logte1[i]=e_logtdem[i+i1];
   rs_te1[i]=pow(10.0, rs_logte1[i]);
  }
     
  rs_all1=(double*)malloc(rs_NT1*rs_Nch*sizeof(double));
  for (int j=0; j<rs_Nch; j++) for (int i=0; i<rs_NT1; i++) rs_spl_arr[j]->Interpolate(rs_logte1[i], rs_all1+i+j*rs_NT1, 0);

  EUV_integrand=(double*)malloc(rs_NT1*sizeof(double));
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

 void *ARGV[11];

 int Rdim[3];
 Rdim[0]=m_Nx;
 Rdim[1]=m_Ny;
 Rdim[2]=m_Nz;

 double Rdxdy[2];
 Rdxdy[0]=m_dx;
 Rdxdy[1]=m_dy;

 for (int i=i_start; i<=i_end; i++) for (int j=j_start; j<=j_end; j++)
 {
  double r1[3], r2[3], LOS[3];

  if (!ProjectionParallel)
  {
   double spx=sin(M_PI/648000*wx[i]);
   double spy=sin(M_PI/648000*wy[j]);
   double q=sqrt(1.0-sqr(spx)-sqr(spy));
   double xD=spx/q;
   double yD=spy/q;

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

   r1[0]=r2[0]=D*wx[i]*M_PI/648000;
   r1[1]=r2[1]=D*wy[j]*M_PI/648000;
   r1[2]=z1;
   r2[2]=z2;

   LOS[0]=LOS[1]=0;
   LOS[2]=z2-z1;
  }

  double aLOS=sqrt(sqr(LOS[0])+sqr(LOS[1])+sqr(LOS[2]));
  for (int l=0; l<3; l++) LOS[l]/=aLOS;

  rotC(r1, m_latC, m_lonC, m_b0Sun);
  rotC(r2, m_latC, m_lonC, m_b0Sun);
  rotC(LOS, m_latC, m_lonC, m_b0Sun);

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
       
  int done=0;
  int TR_on=0;
  if (Nvoxels>0) for (int k=Nvoxels-1; k>=0; k--) if (!done)
  {               
   flags[VoxList[k]]|=1; //voxels crossed by LOSs

   double T_iso=cp_Tbase; //default temperature
   double n_iso=cp_nbase*exp(-h[VoxList[k]]/H_corona); //default density
   int useDEM=0; //default: isothermal
                     
   double Bavg, Lline;
   Bavg=Lline=0;
   
   int ID1, ID2;
   ID1=ID2=0;

   int idx_i=VoxList[k] % m_Nx;
   int idx_j=(VoxList[k]/m_Nx) % m_Ny;
   int idx_k=(VoxList[k]/m_Nx)/m_Ny;
           
   if (idx_k<m_chromo_layers)
   {
	if (m_chromo_nHI[D3(m_Nx, m_Ny, idx_i, idx_j, idx_k)]>0)
	{
	 T_iso=n_iso=0;
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

   if (Bavg>0 && (DEM_on || aNT))
   {
	flags[VoxList[k]]|=4; //voxels with B and L known

	Lline/=2; //switch from line length to line half-length, for consistency with GX Simulator

	double Q=cp_Q0*pow(Bavg/Bavg0, cp_a)/pow(Lline/Lline0, cp_b);
	if (SHtable) Q*=SHtable[D2(SHsize, ID1-1, ID2-1)];
	 
    if (DEM_on)
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
	    useDEM=1;

	    if ((m_VoxelID[VoxList[k]] & 4)!=0) //corona
	    {
	     for (int l=0; l<e_NT; l++) 
		  DEM_local_corona[l]=e_DEM_cor_run[D3(e_NT, e_NQ, l, Qind1, Lind)]*(1.0-dL)*(1.0-dQ1)+
                              e_DEM_cor_run[D3(e_NT, e_NQ, l, Qind1+1, Lind)]*(1.0-dL)*dQ1+
                              e_DEM_cor_run[D3(e_NT, e_NQ, l, Qind2, Lind+1)]*dL*(1.0-dQ2)+
                              e_DEM_cor_run[D3(e_NT, e_NQ, l, Qind2+1, Lind+1)]*dL*dQ2;
	    }

	    if ((m_VoxelID[VoxList[k]] & 2)!=0) //TR
	    {
		 TR_on=1;

		 double costheta=abs(m_Bz[VoxList[k]])/sqrt(sqr(m_Bx[VoxList[k]])+sqr(m_By[VoxList[k]])+sqr(m_Bz[VoxList[k]]));
		 double cosphi=abs(LOS[2]);
		 double cosphi0=sin(0.5*acos(1.0-dz[VoxList[k]]/m_RSun));
		 double TRfactor=costheta/max(cosphi, cosphi0);

		 for (int l=0; l<e_NT; l++) 
		 {
		  DEM_local_TR[l]=e_DEM_tr_run[D3(e_NT, e_NQ, l, Qind1, Lind)]*(1.0-dL)*(1.0-dQ1)+
                          e_DEM_tr_run[D3(e_NT, e_NQ, l, Qind1+1, Lind)]*(1.0-dL)*dQ1+
                          e_DEM_tr_run[D3(e_NT, e_NQ, l, Qind2, Lind+1)]*dL*(1.0-dQ2)+
                          e_DEM_tr_run[D3(e_NT, e_NQ, l, Qind2+1, Lind+1)]*dL*dQ2;
		  DEM_local_TR[l]*=TRfactor;
		 }
	    } 
	   }
	  }
	  else flags[VoxList[k]]|=32; //EBTEL table miss (Q) 
	 }
	 else flags[VoxList[k]]|=16; //EBTEL table miss (L)

	 if (!EBTEL_hit && aNT) FindAnalyticalNT(Q, Lline, &n_iso, &T_iso);
	}
    else if (aNT) FindAnalyticalNT(Q, Lline, &n_iso, &T_iso);
   }

   if (useDEM)
   {
    for (int l=0; l<rs_Nch; l++)
	{
	 for (int m=0; m<rs_NT1; m++) EUV_integrand[m]=rs_all1[D2(rs_NT1, m, l)]*DEM_local_corona[DEM_idx1+m]*rs_te1[m];
	 fluxCorona[D3(b_Nx, b_Ny, i, j, l)]+=(IntTabulated(rs_logte1, EUV_integrand, rs_NT1)*log(10.0)*ds[k]);
	}
   }
   else if (n_iso>0)
   {
	double logT_iso=log10(T_iso);
	for (int l=0; l<rs_Nch; l++)
	{
	 double rs_local;
	 rs_spl_arr[l]->Interpolate(logT_iso, &rs_local, 0);
	 fluxCorona[D3(b_Nx, b_Ny, i, j, l)]+=(sqr(n_iso)*ds[k]*rs_local);
	}
   }

   if ((m_VoxelID[VoxList[k]] & 2)!=0) done=1; //arrived to TR
  }

  if (TR_on)
  {
   for (int l=0; l<rs_Nch; l++)
   {
	for (int m=0; m<rs_NT1; m++) EUV_integrand[m]=rs_all1[D2(rs_NT1, m, l)]*DEM_local_TR[DEM_idx1+m]*rs_te1[m];
	fluxTR[D3(b_Nx, b_Ny, i, j, l)]=IntTabulated(rs_logte1, EUV_integrand, rs_NT1)*log(10.0);
   }
  }
 }

 double norm=b_dx*b_dy/rs_ds;
 for (int i=i_start; i<=i_end; i++) for (int j=j_start; j<=j_end; j++) for (int l=0; l<rs_Nch; l++)
 {
  o_fluxCorona[D3(b_Nx, b_Ny, i, j, l)]=fluxCorona[D3(b_Nx, b_Ny, i, j, l)]*norm;
  o_fluxTR[D3(b_Nx, b_Ny, i, j, l)]=fluxTR[D3(b_Nx, b_Ny, i, j, l)]*norm;
 }

 #ifdef WINDOWS
 if (cs) cs->lock();
 for (int i=0; i<m_Nx*m_Ny*m_Nz; i++) flagsGlobal[i]|=flags[i];
 if (cs) cs->unlock();
 #else
 #pragma omp critical (UpdateFlags)
 for (int i=0; i<m_Nx*m_Ny*m_Nz; i++) flagsGlobal[i]|=flags[i];
 #endif

 free(zmid);
 free(ymid);
 free(xmid);
 free(ds);
 free(VoxList);
 free(flags);
 for (int j=0; j<rs_Nch; j++) delete rs_spl_arr[j];
 free(rs_spl_arr);
 if (DEM_on)
 {
  free(EUV_integrand);
  free(rs_all1);
  free(rs_te1);
  free(rs_logte1);
  free(DEM_local_TR);
  free(DEM_local_corona);
  free(Tgrid);
  free(QgridL);
  free(LgridL);
 }
 if (argc<12)
 {
  free(dz);
  free(h);
 }
 free(fluxTR);
 free(fluxCorona);
 free(wy);
 free(wx);

 return 0;
}

#ifdef WINDOWS
extern "C" __declspec(dllexport) int ComputeEUV(int argc, void **argv)
#else
extern "C" int ComputeEUV(int argc, void **argv)
#endif
{
 __int32 *m32=(__int32*)argv[0];
 int m_Nx=*(m32++);
 int m_Ny=*(m32++);
 int m_Nz=*(m32++); 
 int m_chromo_layers=*(m32++);

 float *m_dz=(float*)(m32+20);

 __int32 *b32=(__int32*)argv[3]; 
 int b_Nx=*(b32++); 
 int b_Ny=*(b32++); 
 b32++;
 int projection=*(b32++);
 int Nthreads=projection>>16;

 __int32 *o_flagsAll=(__int32*)argv[5];
 __int32 *o_flagsCorona=o_flagsAll+6;

 char *flags=(char*)malloc(m_Nx*m_Ny*m_Nz*sizeof(char));
 memset(flags, 0, m_Nx*m_Ny*m_Nz*sizeof(char));

 void *SHtable=(argc>6) ? argv[6] : 0;

 double *dz=(double*)malloc(m_Nx*m_Ny*m_Nz*sizeof(double));
 double *h=(double*)malloc(m_Nx*m_Ny*m_Nz*sizeof(double));

 for (int i=0; i<m_Nx*m_Ny*m_Nz; i++) dz[i]=(double)m_dz[i];

 for (int i=0; i<m_Nx; i++) for (int j=0; j<m_Ny; j++)
 {
  h[D3(m_Nx, m_Ny, i, j, 0)]=dz[D3(m_Nx, m_Ny, i, j, 0)];
  for (int k=1; k<m_Nz; k++) h[D3(m_Nx, m_Ny, i, j, k)]=h[D3(m_Nx, m_Ny, i, j, k-1)]+dz[D3(m_Nx, m_Ny, i, j, k)]; //cumulative sum

  for (int k=0; k<m_Nz; k++) h[D3(m_Nx, m_Ny, i, j, k)]-=dz[D3(m_Nx, m_Ny, i, j, k)]/2;
 }

 #ifdef WINDOWS
 concurrency::critical_section cs;

 int NtMax=concurrency::GetProcessorCount();
 if (Nthreads>0) NtMax=min(Nthreads, NtMax);
 int K=int(b_Ny/NtMax);
 int Np=b_Ny % NtMax;
 int Nm=NtMax-Np;
 
 concurrency::parallel_for(0, NtMax, [&](int j)
 {
  void *ARGV[12];
  memcpy(ARGV, argv, sizeof(void*)*6);
  ARGV[6]=SHtable;

  int fragment[4];
  ARGV[7]=(void*)fragment;

  ARGV[8]=(void*)flags;
  ARGV[9]=(void*)&cs;
  ARGV[10]=(void*)dz;
  ARGV[11]=(void*)h;

  fragment[0]=0;
  fragment[1]=b_Nx-1;

  fragment[2]=(j<Nm) ? K*j : K*Nm+(K+1)*(j-Nm);
  fragment[3]=fragment[2]+((j<Nm) ? K : K+1)-1;

  ComputeEUV_fragment(10, ARGV);
 });
 #else
 int NtMax=omp_get_max_threads();
 if (Nthreads>NtMax) Nthreads=NtMax;
 if (Nthreads>0) omp_set_num_threads(Nthreads);

 #pragma omp parallel for
 for (int j=0; j<b_Ny; j++)
 {
  void *ARGV[12];
  memcpy(ARGV, argv, sizeof(void*)*6);
  ARGV[6]=SHtable;

  int fragment[4];
  ARGV[7]=(void*)fragment;

  ARGV[8]=(void*)flags;
  ARGV[10]=(void*)dz;
  ARGV[11]=(void*)h;

  fragment[0]=0;
  fragment[1]=b_Nx-1;
  fragment[2]=fragment[3]=j;

  ComputeEUV_fragment(12, ARGV);
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

 free(h);
 free(dz);
 free(flags);

 return 0;
}