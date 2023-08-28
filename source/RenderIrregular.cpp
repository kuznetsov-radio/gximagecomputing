#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "RenderIrregular.h"
#include "ExtMath.h"

int CheckCrossPixel(double xl, double xr, double yb, double yt, double x1, double x2, double y1, double y2)
{
 if ((x1>xl) && (x1<=xr) && (y1>yb) && (y1<=yt)) return 1; //1st point inside
 
 if ((x2>xl) && (x2<=xr) && (y2>yb) && (y2<=yt)) return 1; //2nd point inside
 
 double bx=x2-x1;
 if (bx==0.0) if ((x1>xl) && (x1<=xr) && (min(y1, y2)<=yb) && (max(y1, y2)>=yt)) return 1; //vertical crossing
                                                                            else return 0;
 
 double by=y2-y1;
 if (by==0.0) if ((y1>yb) && (y1<=yt) && (min(x1, x2)<=xl) && (max(x1, x2)>=xr)) return 1; //horizontal crossing
                                                                            else return 0;
 
 double yc=y1+(xl-x1)*by/bx;
 if (finite(yc)) if ((yc>=min(y1, y2)) && (yc<=max(y1, y2)) && (yc>yb) && (yc<=yt)) return 1; //crosses the left boundary
 
 yc=y1+(xr-x1)*by/bx;
 if (finite(yc)) if ((yc>=min(y1, y2)) && (yc<=max(y1, y2)) && (yc>yb) && (yc<=yt)) return 1; //crosses the right boundary
 
 double xc=x1+(yb-y1)*bx/by;
 if (finite(xc)) if ((xc>=min(x1, x2)) && (xc<=max(x1, x2)) && (xc>xl) && (xc<=xr)) return 1; //crosses the bottom boundary
 
 xc=x1+(yt-y1)*bx/by;
 if (finite(xc)) if ((xc>=min(x1, x2)) && (xc<=max(x1, x2)) && (xc>xl) && (xc<=xr)) return 1; //crosses the top boundary
 
 return 0;
}

#ifndef LINUX
extern "C" __declspec(dllexport) int RENDER(int argc, void **argv)
#else
extern "C" int RENDER(int argc, void **argv)
#endif
{
 int *Lparms=(int*)argv[0];
 int Nx=Lparms[0];
 int Ny=Lparms[1];
 int Nz=Lparms[2];

 double *dxdy=(double*)argv[1];
 double dx=dxdy[0];
 double dy=dxdy[1];

 double *dz=(double*)argv[2];

 double *LOS_r1=(double*)argv[3];
 double x1=LOS_r1[0];
 double y1=LOS_r1[1];
 double z1=LOS_r1[2];

 double *LOS_r2=(double*)argv[4];
 double x2=LOS_r2[0];
 double y2=LOS_r2[1];
 double z2=LOS_r2[2];

 int *Nvoxels=(int*)argv[5];
 int *VoxList=(int*)argv[6]; 
 double *ds=(double*)argv[7];

 double *xmid=(double*)argv[8];
 double *ymid=(double*)argv[9];
 double *zmid=(double*)argv[10];

 //------------------------------------------------------

 double *tmidarr=(double*)malloc(sizeof(double)*arrN);
 int *iarr=(int*)malloc(sizeof(int)*Nx*Ny);
 int *jarr=(int*)malloc(sizeof(int)*Nx*Ny);
 double *zbarr=(double*)malloc(sizeof(double)*(Nz+1));
 int *Ncrossings=(int*)malloc(sizeof(int)*Nz);
 double *tarr=(double*)malloc(sizeof(double)*Nz*2);

 int res=0;

 double L=sqrt(sqr(x2-x1)+sqr(y2-y1)+sqr(z2-z1));
 
 double bx=x2-x1;
 double by=y2-y1;
 double bz=z2-z1;
 
 double tcmin= 1e100;
 double tcmax=-1e100;
 
 *Nvoxels=0;
 
 int imin=max(int(floor(min(x1, x2)/dx)), 0);
 int imax=min(int(ceil(max(x1, x2)/dx)), Nx-1);
 int jmin=max(int(floor(min(y1, y2)/dy)), 0);
 int jmax=min(int(ceil(max(y1, y2)/dy)), Ny-1);
 
 if (bx==0.0) imin=max(imin-1, 0);
 if (by==0.0) jmin=max(jmin-1, 0);
 
 if ((imax>=imin) && (jmax>=jmin))
 {
  int Nxy=0;
 
  for (int i=imin; i<=imax; i++)
  {
   double xl=dx*i;
   double xr=dx*(i+1);
  
   for (int j=jmin; j<=jmax; j++)
   {
    double yb=dy*j;
    double yt=dy*(j+1);
   
    if (CheckCrossPixel(xl, xr, yb, yt, x1, x2, y1, y2))
	{
     iarr[Nxy]=i;
     jarr[Nxy]=j; 
     Nxy+=1;
	}
   }
  }
 
  if (Nxy>0)
  {
   for (int m=0; m<Nxy; m++)
   {
    int i=iarr[m];
    int j=jarr[m];

    double xl=dx*i;
    double xr=dx*(i+1);
    double yb=dy*j;
    double yt=dy*(j+1);
    
    zbarr[0]=0.0;
    for (int k=1; k<=Nz; k++) zbarr[k]=zbarr[k-1]+dz[i+(j+(k-1)*Ny)*Nx];
    
    for (int k=0; k<Nz; k++) Ncrossings[k]=0;
    
    if (bx!=0.0) 
	{
     double tc=(xl-x1)/bx;
     double yc=y1+tc*by;
     double zc=z1+tc*bz;
     
     if ((yc>yb) && (yc<=yt) && (zc>=zbarr[0]) && (zc<=zbarr[Nz]))
	 {
      int r=value_locate(zbarr, Nz+1, zc);
      if ((r>=0) && (r<Nz))
	  {
       if (Ncrossings[r]<2) tarr[r+Ncrossings[r]*Nz]=tc;
       Ncrossings[r]++;
       tcmin=min(tcmin, tc);
       tcmax=max(tcmax, tc); 
	  }
	 }
     
     tc=(xr-x1)/bx;
     yc=y1+tc*by;
     zc=z1+tc*bz;
     
     if ((yc>yb) && (yc<=yt) && (zc>=zbarr[0]) && (zc<=zbarr[Nz]))
	 {
      int r=value_locate(zbarr, Nz+1, zc);
      if ((r>=0) && (r<Nz))
	  {
       if (Ncrossings[r]<2) tarr[r+Ncrossings[r]*Nz]=tc;
       Ncrossings[r]++;
       tcmin=min(tcmin, tc);
       tcmax=max(tcmax, tc); 
	  }
	 }
	}
    
    if (by!=0.0)
	{
     double tc=(yb-y1)/by;
     double xc=x1+tc*bx;
     double zc=z1+tc*bz;
     
     if ((xc>xl) && (xc<=xr) && (zc>=zbarr[0]) && (zc<=zbarr[Nz]))
	 {
      int r=value_locate(zbarr, Nz+1, zc);
      if ((r>=0) && (r<Nz))
	  {
       if (Ncrossings[r]<2) tarr[r+Ncrossings[r]*Nz]=tc;
       Ncrossings[r]++;
       tcmin=min(tcmin, tc);
       tcmax=max(tcmax, tc); 
	  }
	 }
     
     tc=(yt-y1)/by;
     xc=x1+tc*bx;
     zc=z1+tc*bz;
     
     if ((xc>xl) && (xc<=xr) && (zc>=zbarr[0]) && (zc<=zbarr[Nz]))
	 {
      int r=value_locate(zbarr, Nz+1, zc);
      if ((r>=0) && (r<Nz))
	  {
       if (Ncrossings[r]<2) tarr[r+Ncrossings[r]*Nz]=tc;
       Ncrossings[r]++;
       tcmin=min(tcmin, tc);
       tcmax=max(tcmax, tc); 
	  }
	 }
	}
    
    if (bz!=0.0) for (int k=0; k<=Nz; k++)
	{
     double tc=(zbarr[k]-z1)/bz;
     double xc=x1+tc*bx;
     double yc=y1+tc*by ;
     
     if ((xc>xl) && (xc<=xr) && (yc>yb) && (yc<=yt)) for (int r=max(k-1, 0); r<=min(k, Nz-1); r++)
	 {
      if (Ncrossings[r]<2) tarr[r+Ncrossings[r]*Nz]=tc;
      Ncrossings[r]++;
      tcmin=min(tcmin, tc);
      tcmax=max(tcmax, tc); 
	 }
	}
    
    for (int r=0; r<Nz; r++) if ((*Nvoxels<arrN) && (Ncrossings[r]==2))
	{
     VoxList[*Nvoxels]=r*Nx*Ny+j*Nx+i;
     ds[*Nvoxels]=fabs(tarr[r]-tarr[r+Nz])*L;
     tmidarr[*Nvoxels]=(tarr[r]+tarr[r+Nz])/2;

     xmid[*Nvoxels]=x1+tmidarr[*Nvoxels]*bx;
     ymid[*Nvoxels]=y1+tmidarr[*Nvoxels]*by;
     zmid[*Nvoxels]=z1+tmidarr[*Nvoxels]*bz;
     
     (*Nvoxels)++; 
     if (*Nvoxels>=arrN) res=1;
	}
   }
  }
 }
 
 if ((*Nvoxels>0) && (*Nvoxels<arrN))
 {
  for (int p=0; p<(*Nvoxels-1); p++) for (int q=p+1; q<(*Nvoxels); q++) if (tmidarr[p]>tmidarr[q])
  {
   arrswap(tmidarr, p, q);
   arrswap(VoxList, p, q);
   arrswap(ds, p, q);
   arrswap(xmid, p, q);
   arrswap(ymid, p, q);
   arrswap(zmid, p, q);
  }
 }

 if (*Nvoxels>=arrN) *Nvoxels=0;

 free(tmidarr);
 free(iarr);
 free(jarr);
 free(zbarr);
 free(Ncrossings);
 free(tarr);

 return res;
}