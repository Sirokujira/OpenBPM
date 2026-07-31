#ifndef _BPM_PROTOTYPE_H_
#define _BPM_PROTOTYPE_H_

#ifdef __CUDACC__
#include <cuComplex.h>
#endif


#ifdef __NVCC__
  #include <thrust/complex.h>
  typedef thrust::complex<float> floatcomplex;
  #undef I
  #define I thrust::complex<float>{0,1}
  #define CEXPF(x) (thrust::exp(x))
  #define CREALF(x) (x.real())
  #define CIMAGF(x) (x.imag())
  // obpm.h が定義済みの場合はそちらを使用する (ホスト/デバイス両対応)
  #ifndef MAX
    #define MAX(x,y) (max(x,y))
  #endif
  #ifndef MIN
    #define MIN(x,y) (min(x,y))
  #endif
  #define FLOORF(x) (floor(x))
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
    // 虚数単位 (GCC 拡張の虚数定数。キャストでは虚数単位にならないことに注意)
    #undef I
    #define I (__extension__ 1.0fi)
    // macOS の complex.h は C++ では C99 API (crealf 等) を宣言しないため、
    // gcc/clang 共通のビルトイン・拡張構文で実装する
    #define CEXPF(x) (__builtin_cexpf(x))
    #define CREALF(x) (__real__ (x))
    #define CIMAGF(x) (__imag__ (x))
    //#define I std::complex<float>{0,1}
    //#define CREALF(x) (x.r)
    //#define CIMAGF(x) (x.i)
    #define FLOORF(x) (floorf(x))
  #else
    // MSVC: C99 _Complex が無いため std::complex を使用 (C++ としてコンパイルする)
    #include <algorithm>
    #include <complex>
    typedef std::complex<float> floatcomplex;
    #define I std::complex<float>{0,1}
    #define CEXPF(x) (std::exp(x))
    #define CREALF(x) ((x).real())
    #define CIMAGF(x) ((x).imag())
    #define MAX(x,y) (std::max(x,y))
    #define MIN(x,y) (std::min(x,y))
    #define FLOORF(x) (std::floor(x))
  #endif
#endif

struct debug {
  double             dbls[3];
  unsigned long long ulls[3];
};

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

#ifdef __CUDACC__

// CUDA カーネル (bpm/FDBPMpropagator.cu)
__global__ void substep1a(struct parameters *);
__global__ void substep1b(struct parameters *);
__global__ void substep2a(struct parameters *);
__global__ void substep2b(struct parameters *);
__global__ void applyMultiplier(struct parameters *, long, struct debug*);
__global__ void swapEPointers(struct parameters *, long);
__global__ void updatePrecisePower(struct parameters *);
__global__ void applyTPA(struct parameters *, const float *, double, double *);
__global__ void scalePrecisePowerByTPA(struct parameters *, const double *);

#else

#ifdef __cplusplus
extern "C" {
#endif

extern void substep1a(struct parameters *);
extern void substep1b(struct parameters *);
extern void substep2a(struct parameters *);
extern void substep2b(struct parameters *);
extern void applyMultiplier(struct parameters *, long, struct debug*);
extern void swapEPointers(struct parameters *, long);
extern void updatePrecisePower(struct parameters *);

#ifdef __cplusplus
}
#endif

#endif

//extern void applyMultiplier(struct parameters *, long, struct debug*);

// CUDA
#ifdef __CUDACC__
extern void gpuAssert(cudaError_t code, const char *file, int line);
extern void createDeviceStructs(struct parameters *P, struct parameters **P_devptr, struct debug *D, struct debug **D_devptr);
extern void retrieveAndFreeDeviceStructs(struct parameters *P, struct parameters *P_dev, struct debug *D, struct debug *D_dev);
#endif

#endif

