#pragma once

#ifdef LINUX
#define finite isfinite
#else
#define finite _finite
#endif

#define dNaN (double(HUGE_VAL))

inline double sqr(double x)
{
 return x*x;
}

inline int min(int a, int b)
{
 return (a<b) ? a : b;
}

inline int max(int a, int b)
{
 return (a>b) ? a : b;
}

inline double max(double a, double b)
{
 return (a>b) ? a : b;
}

inline double min(double a, double b)
{
 return (a<b) ? a : b;
}

inline double sign(double x)
{
 return (x<0.0) ? -1.0 : 1.0;
}

inline void arrswap(double *a, int i, int j)
{
 double tmp=a[i];
 a[i]=a[j];
 a[j]=tmp;
}

inline void arrswap(int *a, int i, int j)
{
 int tmp=a[i];
 a[i]=a[j];
 a[j]=tmp;
}

void CrossP(double *a, double *b, double *axb);
void rotC(double *r, double lat, double lon, double b0);
int value_locate(double *a, int n, double x);

double LogFactorial(int n);
double IntTabulated(double *x, double *y, int N);
double InterpolBilinear(double *arr, double *x1arr, double *x2arr, double x1, double x2, int N1, int N2);
double InterpolateBilinear(double *arr, double i1, double i2, int N1, int N2, double missing);
void InterpolateTrilinear(int Nx, int Ny, int Nz, double Dx, double Dy, double *zc_arr, 
                          float *Bx_arr, float *By_arr, float *Bz_arr, double x, double y, double z,
	                      double *Bx, double *By, double *Bz);

class Spline
{
 int N;
 double *x_arr, *y_arr, *y2_arr;
 public:
 Spline(int _N, double *x, double *y);
 ~Spline();
 void Interpolate(double x, double *y, double *y1);
};