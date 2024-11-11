#ifndef _BPM_PROTOTYPE_H_
#define _BPM_PROTOTYPE_H_

#ifdef __CUDACC__
#include <cuComplex.h>
#endif


#ifdef __NVCC__
  #include <thrust/complex.h>
  typedef thrust::complex<float> floatcomplex;
  #define I thrust::complex<float>{0,1}
  #define CEXPF(x) (thrust::exp(x))
  #define CREALF(x) (x.real())
  #define CIMAGF(x) (x.imag())
  #define MAX(x,y) (max(x,y))
  #define MIN(x,y) (min(x,y))
  #define FLOORF(x) (floor(x))
  #include <nvml.h>
  #define TILE_DIM 32
#else
  #ifdef __GNUC__ // This is defined for GCC and CLANG but not for Microsoft Visual C++ compiler
    //#define MAX(a,b) ({__typeof__ (a) _a = (a); __typeof__ (b) _b = (b); _a > _b? _a: _b;})
    //#define MIN(a,b) ({__typeof__ (a) _a = (a); __typeof__ (b) _b = (b); _a > _b? _b: _a;})
    // Prevent conflicting definition of d_complex_t
    #include <stdbool.h>
    #include <complex.h>

    //typedef float complex floatcomplex;
    typedef float _Complex floatcomplex;
    //#define I (1.0f * _Complex_I)
    // 手動で虚数単位を定義
    #define I ((float _Complex) 0.0f + 1.0f * (float _Complex) 1.0f)
    #define CEXPF(x) (cexpf(x))
    #define CREALF(x) (crealf(x))
    #define CIMAGF(x) (cimagf(x))
    //#define I std::complex<float>{0,1}
    //#define CREALF(x) (x.r)
    //#define CIMAGF(x) (x.i)
    #define FLOORF(x) (floorf(x))
  #else
    #include <algorithm>
    #include <complex>
    typedef std::complex<float> floatcomplex;
    #define I std::complex<float>{0,1}
    #define CEXPF(x) (std::exp(x))
    #define CREALF(x) (x.r)
    #define CIMAGF(x) (x.i)
    #define MAX(x,y) (std::max(x,y))
    #define MIN(x,y) (std::min(x,y))
    #define FLOORF(x) (std::floor(x))
  #endif
#endif


struct parameters {
  long Nx;
  long Ny;
  float dx;
  float dy;
  float dz;
  long iz_start;
  long iz_end;
  unsigned char xSymmetry;
  unsigned char ySymmetry;
  float taperPerStep;
  float twistPerStep;
  float d;
  float n_0;
  floatcomplex *n_in;
  long  Nx_n;
  long  Ny_n;
  long  Nz_n;
  float dz_n;
  floatcomplex *Efinal;
  floatcomplex *E1;
  floatcomplex *E2;
  floatcomplex *Eyx;
  floatcomplex *n_out;
  floatcomplex *b;
  float *multiplier;
  floatcomplex ax;
  floatcomplex ay;
  float rho_e;
  float RoC;
  float sinBendDirection;
  float cosBendDirection;
  double precisePower;
  float precisePowerDiff;
  float EfieldPower;
};

#ifdef __cplusplus
extern "C" {
#endif

// C
void substep1a(struct parameters *P_global);
void substep1b(struct parameters *P_global);
void substep2a(struct parameters *P_global);
void substep2b(struct parameters *P_global);
void applyMultiplier(struct parameters *P_global, long iz, struct debug *D);
void swapEPointers(struct parameters *P, long iz);
void updatePrecisePower(struct parameters *P);
#ifdef __cplusplus
}
#endif


// CUDA
#ifdef __CUDACC__
#endif

#endif

