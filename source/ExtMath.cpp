#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "ExtMath.h"
#include "IDLinterface.h"

void CrossP(double *a, double *b, double *axb) //cross product, c = a \times b
{
 axb[0]=a[1]*b[2]-a[2]*b[1];
 axb[1]=a[2]*b[0]-a[0]*b[2];
 axb[2]=a[0]*b[1]-a[1]*b[0];
}

void rotC(double *r, double lat, double lon, double b0)
{
 double slat=sin(M_PI/180*lat);
 double clat=cos(M_PI/180*lat);
 double slon=sin(M_PI/180*lon);
 double clon=cos(M_PI/180*lon);
 double sb0=sin(M_PI/180*b0);
 double cb0=cos(M_PI/180*b0);

 double r1[3], r2[3];

 r1[0]= r[0]; 
 r1[1]= r[2]*sb0+r[1]*cb0;
 r1[2]= r[2]*cb0-r[1]*sb0;

 r2[0]=r1[0]*clon-r1[2]*slon;
 r2[1]=r1[1]; 
 r2[2]=r1[0]*slon+r1[2]*clon;

 r[0]= r2[0]; 
 r[1]=-r2[2]*slat+r2[1]*clat;
 r[2]= r2[2]*clat+r2[1]*slat;
}

int value_locate(double *a, int n, double x)
{
 int asc=a[n-1]>a[0];

 if (asc ? x<a[0] : x>a[0]) return -1;
 if (asc ? x>=a[n-1] : x<=a[n-1]) return n-1;

 int j, j1, l;
 j=0; 
 j1=n-1; 
 while (j1-j>1) 
 { 
  l=(j+j1)>>1; 
  if (asc ? a[l]>x : a[l]<x) j1=l; 
  else j=l; 
 } 
 return j;
} 

double gammln(double xx) 
{ 
 int j;
 double x, tmp, y, ser;
 static const double cof[14]={57.1562356658629235, -59.5979603554754912, 14.1360979747417471, -0.491913816097620199, 
                              0.339946499848118887e-4, 0.465236289270485756e-4, -0.983744753048795646e-4, 0.158088703224912494e-3, 
                              -0.210264441724104883e-3, 0.217439618115212643e-3, -0.164318106536763890e-3, 0.844182239838527433e-4, 
                              -0.261908384015814087e-4, 0.368991826595316234e-5};
 y=x=xx;
 tmp=x+5.24218750000000000;
 tmp=(x+0.5)*log(tmp)-tmp;
 ser=0.999999999999997092;
 for (j=0; j<14; j++) ser+=cof[j]/++y;
 return tmp+log(2.5066282746310005*ser/x);
}

double LogFactorial(int n)
{
 #define Nmax 21

 static double a[Nmax];
 static int init=1;

 if (init)
 {
  double f=1.0;

  for (int i=0; i<Nmax; i++)
  {
   f*=((i>0) ? i : 1); 
   a[i]=log(f);
  }

  init=0;
 }

 return (n<Nmax) ? a[n] : gammln(n+1.0);
}

double IntTabulated(double *x, double *y, int N)
{
 double s=0;
 for (int i=1; i<N; i++) s+=0.5*(y[i-1]+y[i])*(x[i]-x[i-1]);
 return s;
}

double InterpolBilinear(double *arr, double *x1arr, double *x2arr, double x1, double x2, int N1, int N2)
/* Interpolation on an arbitrary grid (like the IDL interpol function). 
   Performs extrapolation, if the point is outside the range. */   
{
 int j, j1, k, k1, l;

 if (x1<x1arr[0])
 {
  j=0;
  j1=1;
 }
 else if (x1>x1arr[N1-1])
 {
  j=N1-2;
  j1=N1-1;
 }
 else
 {
  j=0;
  j1=N1-1;
  while ((j1-j)>1)
  {
   l=(j1+j) >> 1;
   if (x1arr[l]>x1) j1=l;
   else j=l;
  }
 }
 double dx1=x1arr[j1]-x1arr[j];
 double t=(x1-x1arr[j])/dx1;

 if (x2<x2arr[0])
 {
  k=0;
  k1=1;
 }
 else if (x2>x2arr[N2-1])
 {
  k=N2-2;
  k1=N2-1;
 }
 else
 {
  k=0;
  k1=N2-1;
  while ((k1-k)>1)
  {
   l=(k1+k) >> 1;
   if (x2arr[l]>x2) k1=l;
   else k=l;
  }
 }
 double dx2=x2arr[k1]-x2arr[k];
 double u=(x2-x2arr[k])/dx2;
                                                           
 double y1=arr[N2*j+k];
 double y2=arr[N2*j1+k];
 double y3=arr[N2*j1+k1];
 double y4=arr[N2*j+k1];

 return (1.0-t)*(1.0-u)*y1+t*(1.0-u)*y2+t*u*y3+(1.0-t)*u*y4;
}

double InterpolateBilinear(double *arr, double i1, double i2, int N1, int N2, double missing)
/* Interpolation on an equidistant grid (like the IDL interpolate function),
   i1 and i2 - fractional indices of the required point. */
{
 if (i1<0 || i1>(N1-1) || i2<0 || i2>(N2-1)) return missing;

 int j=int(i1);
 int k=int(i2);
 double t=i1-j;
 double u=i2-k;

 double y1=arr[N2*j+k];
 double y2=arr[N2*(j+1)+k];
 double y3=arr[N2*(j+1)+k+1];
 double y4=arr[N2*j+k+1];

 return (1.0-t)*(1.0-u)*y1+t*(1.0-u)*y2+t*u*y3+(1.0-t)*u*y4;
}

void InterpolateTrilinear(int Nx, int Ny, int Nz, double Dx, double Dy, double *zc_arr, 
                          float *Bx_arr, float *By_arr, float *Bz_arr, double x, double y, double z,
	                      double *Bx, double *By, double *Bz)
{
 double x_ind=x/Dx-0.5;
 double y_ind=y/Dy-0.5;

 int i1=min(max(int(x_ind), 0), Nx-2);
 int i2=i1+1;
 int j1=min(max(int(y_ind), 0), Ny-2);
 int j2=j1+1;

 double dx=x_ind-i1;
 double dy=y_ind-j1;

 double *zc_local=(double*)malloc(sizeof(double)*Nz);

 for (int k=0; k<Nz; k++) zc_local[k]=zc_arr[D3(Nx, Ny, i1, j1, k)];
 int k111=min(max(value_locate(zc_local, Nz, z), 0), Nz-2);
 int k112=k111+1;
 double dz11=(z-zc_local[k111])/(zc_local[k112]-zc_local[k111]);

 for (int k=0; k<Nz; k++) zc_local[k]=zc_arr[D3(Nx, Ny, i1, j2, k)];
 int k121=min(max(value_locate(zc_local, Nz, z), 0), Nz-2);
 int k122=k121+1;
 double dz12=(z-zc_local[k121])/(zc_local[k122]-zc_local[k121]);

 for (int k=0; k<Nz; k++) zc_local[k]=zc_arr[D3(Nx, Ny, i2, j1, k)];
 int k211=min(max(value_locate(zc_local, Nz, z), 0), Nz-2);
 int k212=k211+1;
 double dz21=(z-zc_local[k211])/(zc_local[k212]-zc_local[k211]);

 for (int k=0; k<Nz; k++) zc_local[k]=zc_arr[D3(Nx, Ny, i2, j2, k)];
 int k221=min(max(value_locate(zc_local, Nz, z), 0), Nz-2);
 int k222=k221+1;
 double dz22=(z-zc_local[k221])/(zc_local[k222]-zc_local[k221]);

 free(zc_local);

 double u111=(1.0-dx)*(1.0-dy)*(1.0-dz11);
 double u211=dx*(1.0-dy)*(1.0-dz21);
 double u121=(1.0-dx)*dy*(1.0-dz12);
 double u112=(1.0-dx)*(1.0-dy)*dz11;
 double u212=dx*(1.0-dy)*dz21;
 double u122=(1.0-dx)*dy*dz12;
 double u221=dx*dy*(1.0-dz22);
 double u222=dx*dy*dz22;

 *Bx=Bx_arr[D3(Nx, Ny, i1, j1, k111)]*u111+
     Bx_arr[D3(Nx, Ny, i2, j1, k211)]*u211+
     Bx_arr[D3(Nx, Ny, i1, j2, k121)]*u121+
     Bx_arr[D3(Nx, Ny, i1, j1, k112)]*u112+
     Bx_arr[D3(Nx, Ny, i2, j1, k212)]*u212+
     Bx_arr[D3(Nx, Ny, i1, j2, k122)]*u122+
     Bx_arr[D3(Nx, Ny, i2, j2, k221)]*u221+
     Bx_arr[D3(Nx, Ny, i2, j2, k222)]*u222;

 *By=By_arr[D3(Nx, Ny, i1, j1, k111)]*u111+
     By_arr[D3(Nx, Ny, i2, j1, k211)]*u211+
     By_arr[D3(Nx, Ny, i1, j2, k121)]*u121+
     By_arr[D3(Nx, Ny, i1, j1, k112)]*u112+
     By_arr[D3(Nx, Ny, i2, j1, k212)]*u212+
     By_arr[D3(Nx, Ny, i1, j2, k122)]*u122+
     By_arr[D3(Nx, Ny, i2, j2, k221)]*u221+
     By_arr[D3(Nx, Ny, i2, j2, k222)]*u222;

 *Bz=Bz_arr[D3(Nx, Ny, i1, j1, k111)]*u111+
     Bz_arr[D3(Nx, Ny, i2, j1, k211)]*u211+
     Bz_arr[D3(Nx, Ny, i1, j2, k121)]*u121+
     Bz_arr[D3(Nx, Ny, i1, j1, k112)]*u112+
     Bz_arr[D3(Nx, Ny, i2, j1, k212)]*u212+
     Bz_arr[D3(Nx, Ny, i1, j2, k122)]*u122+
     Bz_arr[D3(Nx, Ny, i2, j2, k221)]*u221+
     Bz_arr[D3(Nx, Ny, i2, j2, k222)]*u222;
}

void spline_init(double *x, double *y, int n, double yp1, double ypn, double *y2)
{
 int i, k;
 double p, qn, sig, un, *u;
 u=(double*)malloc(sizeof(double)*n);
 if (!finite(yp1)) y2[0]=u[0]=0.0; 
 else 
 { 
  y2[0]=-0.5;
  u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
 }
 for (i=1; i<n-1; i++) 
 { 
  sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
  p=sig*y2[i-1]+2.0;
  y2[i]=(sig-1.0)/p;
  u[i]=(y[i+1]-y[i])/(x[i+1]-x[i])-(y[i]-y[i-1])/(x[i]-x[i-1]);
  u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
 }
 if (!finite(ypn)) qn=un=0.0;
 else 
 { 
  qn=0.5;
  un=(3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
 }
 y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);
 for (k=n-2; k>=0; k--) y2[k]=y2[k]*y2[k+1]+u[k]; 
 free(u);
}

void spline_short(double *x, double *y, int n, double *y2)
{
 double dxl, dxr;
 double K1[3];

 dxl=x[1]-x[0];
 dxr=x[2]-x[1];
 K1[0]=-(2.0*dxl+dxr)/dxl/(dxl+dxr);
 K1[1]=(dxr+dxl)/(dxr*dxl);
 K1[2]=-dxl/dxr/(dxr+dxl);
 double y1l=K1[0]*y[0]+K1[1]*y[1]+K1[2]*y[2];

 dxl=x[n-2]-x[n-3];
 dxr=x[n-1]-x[n-2];
 K1[0]=dxr/dxl/(dxr+dxl);
 K1[1]=-(dxr+dxl)/(dxr*dxl);
 K1[2]=(2.0*dxr+dxl)/dxr/(dxr+dxl);
 double y1r=K1[0]*y[n-3]+K1[1]*y[n-2]+K1[2]*y[n-1];

 spline_init(x, y, n, y1l, y1r, y2);
}

Spline :: Spline(int _N, double *x, double *y)
{
 N=_N;

 x_arr=(double*)malloc(sizeof(double)*N);
 y_arr=(double*)malloc(sizeof(double)*N);
 y2_arr=(double*)malloc(sizeof(double)*N);

 for (int i=0; i<N; i++)
 {
  x_arr[i]=x[i];
  y_arr[i]=y[i];
 }

 spline_short(x_arr, y_arr, N, y2_arr);
}

Spline :: ~Spline()
{
 free(x_arr);
 free(y_arr);
 free(y2_arr);
}

void spline_interp(double *xa, double *ya, double *y2a, int n, double x, double *y, double *y1)
{
 int klo, khi, k;
 double h, b, a;

 if ((x<xa[0]) || (x>xa[n-1]))
 {
  if (y) *y=0;
  if (y1) *y1=0;
 }
 else
 {
  klo=0; 
  khi=n-1;
  while (khi-klo>1) 
  {
   k=(khi+klo)>>1;
   if (xa[k]>x) khi=k;
   else klo=k;
  } 
 
  h=xa[khi]-xa[klo];
  a=(xa[khi]-x)/h; 
  b=(x-xa[klo])/h; 

  if (y) *y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
  if (y1) *y1=(ya[khi]-ya[klo])/h+((1.0-3.0*a*a)*y2a[klo]+(3.0*b*b-1)*y2a[khi])*h/6.0;
 }
}

void Spline :: Interpolate(double x, double *y, double *y1)
{
 spline_interp(x_arr, y_arr, y2_arr, N, x, y, y1);
}