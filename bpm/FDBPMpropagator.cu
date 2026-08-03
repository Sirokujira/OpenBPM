/********************************************
 * FDBPMpropagator.cu, in the C programming language, written for MATLAB MEX function generation
 * 
 ** Compiling on Windows
 * Can be compiled with GCC using
 * "mex COPTIMFLAGS='$COPTIMFLAGS -Ofast -fopenmp -std=c11 -Wall' LDOPTIMFLAGS='$LDOPTIMFLAGS -Ofast -fopenmp -std=c11 -Wall' -outdir +BPMmatlab\@model\private .\src\FDBPMpropagator.c ".\src\libut.lib" -R2018a"
 * ... or the Microsoft Visual C++ compiler (MSVC) with
 * "copyfile ./src/FDBPMpropagator.c ./src/FDBPMpropagator.cpp; mex COMPFLAGS='/Zp8 /GR /EHs /nologo /MD /openmp /W4 /WX /wd4204 /wd4100' -outdir +BPMmatlab\@model\private .\src\FDBPMpropagator.cpp ".\src\libut.lib" -R2018a"
 * 
 * The source code in this file is written is such a way that it is
 * compilable by either C or C++ compilers, either with GCC, MSVC or
 * the Nvidia CUDA compiler called NVCC, which is based on MSVC. To
 * compile with CUDA GPU acceleration support, you must have MSVC
 * installed. As of January 2020, mexcuda does not work with MSVC 2019,
 * so I'd recommend MSVC 2017. You also need the Parallel Computing
 * Toolbox, which you will find in the MATLAB addon manager. To compile, run:
 * "copyfile ./src/FDBPMpropagator.c ./src/FDBPMpropagator_CUDA.cu; mexcuda -llibut COMPFLAGS='-use_fast_math -res-usage $COMPFLAGS' -outdir +BPMmatlab\@model\private .\src\FDBPMpropagator_CUDA.cu -R2018a"
 *
 ** Compiling on macOS
 * As of March 2021, the macOS compiler doesn't support libut (for ctrl+c 
 * breaking) or openmp (for multithreading).
 * "mex COPTIMFLAGS='$COPTIMFLAGS -Ofast -std=c11 -Wall' LDOPTIMFLAGS='$LDOPTIMFLAGS -Ofast -std=c11 -Wall' -outdir +BPMmatlab/@model/private ./src/FDBPMpropagator.c -R2018a"
 *
 * To get the MATLAB C compiler to work, try this:
 * 1. Install XCode from the App Store
 * 2. Type "mex -setup" in the MATLAB command window
 *
 ** Compiling on Linux
 * "mex COPTIMFLAGS='$COPTIMFLAGS -Ofast -fopenmp -std=c11 -Wall' LDOPTIMFLAGS='$LDOPTIMFLAGS -Ofast -fopenmp -std=c11 -Wall' -outdir +BPMmatlab/@model/private ./src/FDBPMpropagator.c -R2018a -lut"
 *
 * To get the MATLAB C compiler to work, try this:
 * 1. Use a package manager like apt to install GCC (on Ubuntu, part of the build-essential package)
 * 2. Type "mex -setup" in the MATLAB command window
 ********************************************/
// printf("Reached line %d...\n",__LINE__);mexEvalString("drawnow; pause(.001);");mexEvalString("drawnow; pause(.001);");mexEvalString("drawnow; pause(.001);"); // For inserting into code for debugging purposes

#include <math.h>
#include <stdint.h>
#include <stdbool.h>
#include <float.h>
//#include "mex.h"
#include "obpm.h"
//#define FLT_EPSILON 
//#define PI acosf(-1.0f)
#ifdef _OPENMP
  #include "omp.h"
#endif
#ifndef __clang__
  #ifdef __cplusplus
  extern "C"
  #endif
  //extern bool utIsInterruptPending(); // Allows catching ctrl+c while executing the mex function
#endif

#include "bpm/bpm_prototype.h"

__host__ __device__
float sqrf(float x) {return x*x;}

// #if defined(__CUDA_ARCH__) && __CUDA_ARCH__ < 600
// __device__ double atomicAdd(double* address, double val) {
//   unsigned long long int* address_as_ull = (unsigned long long int*)address;
//   unsigned long long int old = *address_as_ull, assumed;
// 
//   do {
//     assumed = old;
//     old = atomicCAS(address_as_ull, assumed,__double_as_longlong(val + __longlong_as_double(assumed)));
//   } while (assumed != old);
//   return __longlong_as_double(old);
// }
// #endif

__global__
void substep1a(struct parameters *P_global) {
  // Explicit part of substep 1 out of 2

  __shared__ char Pdummy[sizeof(struct parameters)];
  struct parameters *P = (struct parameters *)Pdummy;
  if(!threadIdx.x && !threadIdx.y) *P = *P_global; // Only let one thread per block do the copying. 
  __syncthreads(); // All threads in the block wait for the copy to have finished

  bool xAntiSymm = P->xSymmetry == 2;
  bool yAntiSymm = P->ySymmetry == 2;

  __shared__ double tiledummy[TILE_DIM][TILE_DIM+1]; // We declare with double because a double is the same size as a float complex. +1 is to avoid memory bank conflicts
  floatcomplex *tile = (floatcomplex *)tiledummy;
  
  long xTiles = (P->Nx + TILE_DIM - 1)/TILE_DIM;
  long yTiles = (P->Ny + TILE_DIM - 1)/TILE_DIM;
  for(long tileNum=blockIdx.x; tileNum<xTiles*yTiles; tileNum += gridDim.x) {
    long tilexoffset = TILE_DIM*(tileNum%xTiles);
    long tileyoffset = TILE_DIM*(tileNum/xTiles);
    long ix = tilexoffset + threadIdx.x;
    long iy = tileyoffset + threadIdx.y;

    if(ix<P->Nx && iy<P->Ny) {
      long i = ix + iy*P->Nx;
      tile[threadIdx.x + threadIdx.y*(TILE_DIM+1)] = P->E1[i];
      if(ix != 0                                 ) tile[threadIdx.x + threadIdx.y*(TILE_DIM+1)] += (P->E1[i-1]     - P->E1[i])*P->ax;
      if(ix != P->Nx-1 && (!yAntiSymm || ix != 0)) tile[threadIdx.x + threadIdx.y*(TILE_DIM+1)] += (P->E1[i+1]     - P->E1[i])*P->ax;
      if(iy != 0                                 ) tile[threadIdx.x + threadIdx.y*(TILE_DIM+1)] += (P->E1[i-P->Nx] - P->E1[i])*P->ay*2;
      if(iy != P->Ny-1 && (!xAntiSymm || iy != 0)) tile[threadIdx.x + threadIdx.y*(TILE_DIM+1)] += (P->E1[i+P->Nx] - P->E1[i])*P->ay*2;
    }
    __syncthreads();
    // Save transposed xy -> yx
    ix = tilexoffset + threadIdx.y;
    iy = tileyoffset + threadIdx.x;
    if(ix<P->Nx && iy<P->Ny) P->Eyx[iy + ix*P->Ny] = tile[threadIdx.y + threadIdx.x*(TILE_DIM+1)];
    __syncthreads();
  }
}

__global__
void substep1b(struct parameters *P_global) {
  // Implicit part of substep 1 out of 2
  long threadNum = threadIdx.x + threadIdx.y*blockDim.x + blockIdx.x*blockDim.x*blockDim.y;
  __shared__ char Pdummy[sizeof(struct parameters)];
  struct parameters *P = (struct parameters *)Pdummy;
  if(!threadIdx.x && !threadIdx.y) *P = *P_global; // Only let one thread per block do the copying. 
  __syncthreads(); // All threads in the block wait for the copy to have finished
  bool yAntiSymm = P->ySymmetry == 2;
  for(long iy=threadNum;iy<P->Ny;iy+=gridDim.x*blockDim.x*blockDim.y){
    for(long ix=0; ix<P->Nx; ix++) {
      long i = iy + ix*P->Ny;
      if     (ix == 0 && yAntiSymm) P->b[i] = 1          ;
      else if(ix == 0             ) P->b[i] = 1 +   P->ax;
      else if(ix < P->Nx-1        ) P->b[i] = 1 + 2*P->ax;
      else                          P->b[i] = 1 +   P->ax;

      if(ix > 0) {
        floatcomplex w   = -P->ax/P->b[i-P->Ny];
        P->b[i]         += w*(ix == 1 && yAntiSymm? 0: P->ax);
        P->Eyx[i]       -= w*P->Eyx[i-P->Ny];
      }
    }

    for(long ix=P->Nx-1; ix>=0 + yAntiSymm; ix--) {
      long i = iy + ix*P->Ny;
      P->Eyx[i] = (P->Eyx[i] + (ix == P->Nx-1? 0: P->ax*P->Eyx[i+P->Ny]))/P->b[i];
    }
  }
}

__global__
void substep2a(struct parameters *P_global) {
  // Explicit part of substep 2 out of 2
  __shared__ char Pdummy[sizeof(struct parameters)];
  struct parameters *P = (struct parameters *)Pdummy;
  if(!threadIdx.x && !threadIdx.y) *P = *P_global; // Only let one thread per block do the copying. 
  __syncthreads(); // All threads in the block wait for the copy to have finished
  __shared__ double tiledummy[TILE_DIM][TILE_DIM+1]; // We declare with double because a double is the same size as a float complex. +1 is to avoid memory bank conflicts
  floatcomplex *tile = (floatcomplex *)tiledummy;

  bool xAntiSymm = P->xSymmetry == 2;

  long xTiles = (P->Nx + TILE_DIM - 1)/TILE_DIM;
  long yTiles = (P->Ny + TILE_DIM - 1)/TILE_DIM;
  for(long tileNum=blockIdx.x; tileNum<xTiles*yTiles; tileNum += gridDim.x) {
    long tilexoffset = TILE_DIM*(tileNum%xTiles);
    long tileyoffset = TILE_DIM*(tileNum/xTiles);
    long ix = tilexoffset + threadIdx.y;
    long iy = tileyoffset + threadIdx.x;

    __syncthreads(); // All threads in the block wait for any previous tile usage to have completed
    if(ix<P->Nx && iy<P->Ny) tile[threadIdx.y + threadIdx.x*(TILE_DIM+1)] = P->Eyx[ix*P->Ny + iy]; // load yx data and store in shared memory in tile, which is xy
    __syncthreads(); // All threads in the block wait for the copy to have finished

    ix = tilexoffset + threadIdx.x;
    iy = tileyoffset + threadIdx.y;
    if(ix<P->Nx && iy<P->Ny) {
      long i = ix + iy*P->Nx;

      floatcomplex deltaE = 0;
      if(iy != 0                                 ) deltaE -= (P->E1[i-P->Nx] - P->E1[i])*P->ay;
      if(iy != P->Ny-1 && (!xAntiSymm || iy != 0)) deltaE -= (P->E1[i+P->Nx] - P->E1[i])*P->ay;
      P->E2[i] = tile[threadIdx.x + threadIdx.y*(TILE_DIM+1)] + deltaE;
    }
  }
}

__global__
void substep2b(struct parameters *P_global) {
  // Implicit part of substep 2 out of 2
  float EfieldPowerThread = 0.0f;
  long threadNum = threadIdx.x + threadIdx.y*blockDim.x + blockIdx.x*blockDim.x*blockDim.y;
  __shared__ char Pdummy[sizeof(struct parameters)];
  struct parameters *P = (struct parameters *)Pdummy;
  if(!threadIdx.x && !threadIdx.y) *P = *P_global; // Only let one thread per block do the copying. 
  __syncthreads(); // All threads in the block wait for the copy to have finished
  bool xAntiSymm = P->xSymmetry == 2;
  for(long ix=threadNum;ix<P->Nx;ix+=gridDim.x*blockDim.x*blockDim.y) {
    for(long iy=0; iy<P->Ny; iy++) {
      long i = ix + iy*P->Nx;
      if     (iy == 0 && xAntiSymm) P->b[i] = 1          ;
      else if(iy == 0             ) P->b[i] = 1 +   P->ay;
      else if(iy < P->Ny-1        ) P->b[i] = 1 + 2*P->ay;
      else                          P->b[i] = 1 +   P->ay;

      if(iy > 0) {
        floatcomplex w   = -P->ay/P->b[i-P->Nx];
        P->b[i]         += w*(iy == 1 && xAntiSymm? 0: P->ay);
        P->E2[i]        -= w*P->E2[i-P->Nx];
      }
    }

    for(long iy=P->Ny-1; iy>=0 + xAntiSymm; iy--) {
      long i = ix + iy*P->Nx;
      P->E2[i] = (P->E2[i] + (iy == P->Ny-1? 0: P->ay*P->E2[i+P->Nx]))/P->b[i];
      EfieldPowerThread += sqrf(CREALF(P->E2[i])) + sqrf(CIMAGF(P->E2[i]));
    }
  }

  atomicAdd(&P_global->EfieldPower,EfieldPowerThread);
}

__global__
void applyMultiplier(struct parameters *P_global, long iz, struct debug *D) {
  float precisePowerDiffThread = 0.0f;
  long threadNum = threadIdx.x + threadIdx.y*blockDim.x + blockIdx.x*blockDim.x*blockDim.y;
  __shared__ char Pdummy[sizeof(struct parameters)];
  struct parameters *P = (struct parameters *)Pdummy;
  if(!threadIdx.x && !threadIdx.y) *P = *P_global; // Only let one thread per block do the copying
  __syncthreads(); // All threads in the block wait for the copy to have finished
  float fieldCorrection = sqrtf((float)P->precisePower/P->EfieldPower);
  float cosvalue = cosf(-P->twistPerStep*iz); // Minus is because we go from the rotated frame to the source frame
  float sinvalue = sinf(-P->twistPerStep*iz);
  float scaling = 1/(1 - P->taperPerStep*iz); // Take reciprocal because we go from scaled frame to unscaled frame
  for(long i=threadNum;i<P->Nx*P->Ny;i+=gridDim.x*blockDim.x*blockDim.y) {
    long ix = i%P->Nx;
    float x = P->dx*(ix - (P->Nx-1)/2.0f*(P->ySymmetry == 0));
    long iy = i/P->Nx;
    float y = P->dy*(iy - (P->Ny-1)/2.0f*(P->xSymmetry == 0));
    floatcomplex n = 0;
    if(P->taperPerStep || P->twistPerStep) { // Rotate, scale, interpolate. If we are tapering or twisting, we know that the RIP is 2D
      float x_src = scaling*(cosvalue*x - sinvalue*y);
      float y_src = scaling*(sinvalue*x + cosvalue*y);
      float ix_src = MIN(MAX(0.0f,x_src/P->dx + (P->Nx - 1)/2.0f*(P->ySymmetry == 0)),(P->Nx - 1)*(1-FLT_EPSILON)); // Fractional index, coerced to be within the source window
      float iy_src = MIN(MAX(0.0f,y_src/P->dy + (P->Ny - 1)/2.0f*(P->xSymmetry == 0)),(P->Ny - 1)*(1-FLT_EPSILON));
      long ix_low = (long)FLOORF(ix_src);
      long iy_low = (long)FLOORF(iy_src);
      float ix_frac = ix_src - FLOORF(ix_src);
      float iy_frac = iy_src - FLOORF(iy_src);
      n = P->n_in[ix_low     + P->Nx*(iy_low    )]*(1 - ix_frac)*(1 - iy_frac) +
          P->n_in[ix_low + 1 + P->Nx*(iy_low    )]*(    ix_frac)*(1 - iy_frac) +
          P->n_in[ix_low     + P->Nx*(iy_low + 1)]*(1 - ix_frac)*(    iy_frac) +
          P->n_in[ix_low + 1 + P->Nx*(iy_low + 1)]*(    ix_frac)*(    iy_frac); // Bilinear interpolation
    } else if(P->Nz_n == 1) { // 2D RIP
      n = P->n_in[i];
    } else { // 3D RIP
      float z = iz*P->dz;
      long ix_n = MIN(MAX(0L,ix - (P->Nx - P->Nx_n)/2),P->Nx_n-1);
      long iy_n = MIN(MAX(0L,iy - (P->Ny - P->Ny_n)/2),P->Ny_n-1);
      float iz_n = MIN(MAX(0.0f,z/P->dz_n),(P->Nz_n - 1)*(1-FLT_EPSILON)); // Fractional index, coerced to be within the n window
      long iz_n_low = (long)FLOORF(iz_n);
      float iz_n_frac = iz_n - FLOORF(iz_n);
      n = P->n_in[ix_n + P->Nx_n*iy_n + P->Ny_n*P->Nx_n*(iz_n_low    )]*(1 - iz_n_frac) +
          P->n_in[ix_n + P->Nx_n*iy_n + P->Ny_n*P->Nx_n*(iz_n_low + 1)]*(    iz_n_frac); // Linear interpolation in z
    }
    if(iz == P->iz_end-1) P->n_out[i] = n;
    float n_bend = CREALF(n)*(1-(sqrf(CREALF(n))*(x*P->cosBendDirection+y*P->sinBendDirection)/2/P->RoC*P->rho_e))*exp((x*P->cosBendDirection+y*P->sinBendDirection)/P->RoC);
    floatcomplex a = P->multiplier[i]*CEXPF(P->d*(CIMAGF(n) + (sqrf(n_bend) - sqrf(P->n_0))*I/(2*P->n_0))); // Multiplier includes only the edge absorber
    P->E2[i] *= fieldCorrection*a;
    float anormsqr = sqrf(CREALF(a)) + sqrf(CIMAGF(a));
    if(anormsqr > 1 - 10*FLT_EPSILON && anormsqr < 1 + 10*FLT_EPSILON) anormsqr = 1; // To avoid accumulating power discrepancies due to rounding errors
    precisePowerDiffThread += (sqrf(CREALF(P->E2[i])) + sqrf(CIMAGF(P->E2[i])))*(1 - 1/anormsqr);
  }

  atomicAdd(&P_global->precisePowerDiff,precisePowerDiffThread);
}

__global__
void updatePrecisePower(struct parameters *P) {
P->precisePower += P->precisePowerDiff;
P->precisePowerDiff = 0;
}

// 二光子吸収 (TPA) : E2 *= exp(-(beta/2)*I*dz), I = |E2|^2 [W/m^2]
// (CPU 版 sol/solve_bpm.cpp と同一。物理スケーリング済みの界を仮定し、
//  強度の減衰率 alpha = beta*I に対し界には alpha/2 を適用する)
// sums[0] に適用前、sums[1] に適用後の電力和を積算する (電力簿記の補正用)。
__global__
void applyTPA(struct parameters *P, const float *beta, double dz, double *sums) {
  double pb = 0, pa = 0;
  for(long i = blockIdx.x*(long)blockDim.x + threadIdx.x; i < P->Nx*P->Ny; i += gridDim.x*(long)blockDim.x) {
    const floatcomplex e = P->E2[i];
    const double i2 = ((double)CREALF(e)*CREALF(e)) + ((double)CIMAGF(e)*CIMAGF(e));
    const double b = beta[i];
    pb += i2;
    if(b > 0) {
      // CPU 版 (sol/solve_bpm.cpp) と同一 : 倍精度で exp を評価してから float 化する
      const float g = (float)exp(-0.5*b*i2*dz);
      P->E2[i] = e*g;
      pa += i2*(double)g*g;
    } else {
      pa += i2;
    }
  }
  atomicAdd(&sums[0], pb);
  atomicAdd(&sums[1], pa);
}

// TPA による電力簿記の補正 : precisePower *= (適用後 / 適用前)
// (次ステップの fieldCorrection が TPA 減衰を打ち消さないようにする)
__global__
void scalePrecisePowerByTPA(struct parameters *P, const double *sums) {
  if(sums[0] > 0) P->precisePower *= sums[1]/sums[0];
}


// z 断面の統計量 (GUI 表示用の /trace) をデバイス上で集計する
// acc[0..4] = sum|E|^2, sum x|E|^2, sum y|E|^2, sum x^2|E|^2, sum y^2|E|^2
// acc[5]    = |E|^2 の最大値 (非負 double は IEEE ビット列の大小と順序が一致するため
//             unsigned long long の atomicMax で求められる)
__global__
void fieldTrace(struct parameters *P, const double *xc, const double *yc, double *acc) {
  double s = 0, sx = 0, sy = 0, sxx = 0, syy = 0, pk = 0;
  for(long i = blockIdx.x*(long)blockDim.x + threadIdx.x; i < P->Nx*P->Ny; i += gridDim.x*(long)blockDim.x) {
    const long ix = i % P->Nx;
    const long iy = i / P->Nx;
    const floatcomplex e = P->E2[i];
    const double iv = ((double)CREALF(e)*CREALF(e)) + ((double)CIMAGF(e)*CIMAGF(e));
    s   += iv;
    sx  += xc[ix]*iv;
    sy  += yc[iy]*iv;
    sxx += xc[ix]*xc[ix]*iv;
    syy += yc[iy]*yc[iy]*iv;
    if(iv > pk) pk = iv;
  }
  atomicAdd(&acc[0], s);
  atomicAdd(&acc[1], sx);
  atomicAdd(&acc[2], sy);
  atomicAdd(&acc[3], sxx);
  atomicAdd(&acc[4], syy);
  atomicMax((unsigned long long *)&acc[5], (unsigned long long)__double_as_longlong(pk));
}

__global__
void swapEPointers(struct parameters *P, long iz) {
  P->EfieldPower = 0;
  floatcomplex *temp = P->E1;
  P->E1 = P->E2;
  P->E2 = temp;
}

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line) {
  if (code != cudaSuccess) {
    printf("GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
    //mexEvalString("drawnow; pause(.001);");
    while(true) {;}
  }
}

void createDeviceStructs(struct parameters *P, struct parameters **P_devptr,
                         struct debug *D, struct debug **D_devptr) {
  long N = P->Nx*P->Ny;
  long N_n = P->Nx_n*P->Ny_n*P->Nz_n;
  struct parameters P_tempvar = *P;

  gpuErrchk(cudaMalloc(&P_tempvar.E1,N*sizeof(floatcomplex)));
  gpuErrchk(cudaMemcpy(P_tempvar.E1,P->E1,N*sizeof(floatcomplex),cudaMemcpyHostToDevice));
  gpuErrchk(cudaMalloc(&P_tempvar.E2,N*sizeof(floatcomplex)));
  gpuErrchk(cudaMalloc(&P_tempvar.Eyx,N*sizeof(floatcomplex)));

  gpuErrchk(cudaMalloc(&P_tempvar.n_out,N*sizeof(floatcomplex)));

  gpuErrchk(cudaMalloc(&P_tempvar.b,N*sizeof(floatcomplex)));
  gpuErrchk(cudaMalloc(&P_tempvar.multiplier,N*sizeof(float)));
  gpuErrchk(cudaMemcpy(P_tempvar.multiplier,P->multiplier,N*sizeof(float),cudaMemcpyHostToDevice));
  gpuErrchk(cudaMalloc(&P_tempvar.n_in,N_n*sizeof(floatcomplex)));
  gpuErrchk(cudaMemcpy(P_tempvar.n_in,P->n_in,N_n*sizeof(floatcomplex),cudaMemcpyHostToDevice));

  gpuErrchk(cudaMalloc(P_devptr, sizeof(struct parameters)));
  gpuErrchk(cudaMemcpy(*P_devptr,&P_tempvar,sizeof(struct parameters),cudaMemcpyHostToDevice));

  // Allocate and copy debug struct
  struct debug D_tempvar = *D;
  gpuErrchk(cudaMalloc(D_devptr, sizeof(struct debug)));
  gpuErrchk(cudaMemcpy(*D_devptr,&D_tempvar,sizeof(struct debug),cudaMemcpyHostToDevice));
}

void retrieveAndFreeDeviceStructs(struct parameters *P, struct parameters *P_dev,
                                  struct debug *D, struct debug *D_dev) {
  long N = P->Nx*P->Ny;
  struct parameters P_temp; gpuErrchk(cudaMemcpy(&P_temp,P_dev,sizeof(struct parameters),cudaMemcpyDeviceToHost));
  gpuErrchk(cudaMemcpy(P->Efinal,P_temp.E2,N*sizeof(floatcomplex),cudaMemcpyDeviceToHost));
  gpuErrchk(cudaMemcpy(P->n_out,P_temp.n_out,N*sizeof(floatcomplex),cudaMemcpyDeviceToHost));
  gpuErrchk(cudaFree(P_temp.n_out));
  P->precisePower = P_temp.precisePower;

  gpuErrchk(cudaFree(P_temp.E1));
  gpuErrchk(cudaFree(P_temp.E2));
  gpuErrchk(cudaFree(P_temp.Eyx));
  gpuErrchk(cudaFree(P_temp.b));
  gpuErrchk(cudaFree(P_temp.multiplier));
  gpuErrchk(cudaFree(P_temp.n_in));
  gpuErrchk(cudaFree(P_dev));


  gpuErrchk(cudaMemcpy(D, D_dev, sizeof(struct debug),cudaMemcpyDeviceToHost));
  gpuErrchk(cudaFree(D_dev));
}

/*
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[]) {
  struct parameters P_var;
  struct parameters *P = &P_var;
  P->Nx = (long)mxGetM(prhs[0]);
  P->Ny = (long)mxGetN(prhs[0]);
  P->dx = *(float *)mxGetData(mxGetField(prhs[1],0,"dx"));
  P->dy = *(float *)mxGetData(mxGetField(prhs[1],0,"dy"));
  P->dz = *(float *)mxGetData(mxGetField(prhs[1],0,"dz"));
  P->iz_start = *(long *)mxGetData(mxGetField(prhs[1],0,"iz_start"));
  P->iz_end = *(long *)mxGetData(mxGetField(prhs[1],0,"iz_end"));
  P->taperPerStep = *(float *)mxGetData(mxGetField(prhs[1],0,"taperPerStep"));
  P->twistPerStep = *(float *)mxGetData(mxGetField(prhs[1],0,"twistPerStep"));
  P->xSymmetry = *(unsigned char *)mxGetData(mxGetField(prhs[1],0,"xSymmetry"));
  P->ySymmetry = *(unsigned char *)mxGetData(mxGetField(prhs[1],0,"ySymmetry"));
  P->d = *(float *)mxGetData(mxGetField(prhs[1],0,"d"));
  P->n_0 = *(float *)mxGetData(mxGetField(prhs[1],0,"n_0"));
  P->n_in = (floatcomplex *)mxGetData(mxGetField(prhs[1],0,"n_mat"));
  mwSize nDims = mxGetNumberOfDimensions(mxGetField(prhs[1],0,"n_mat"));
  mwSize const *dimPtr = mxGetDimensions(mxGetField(prhs[1],0,"n_mat"));
  P->Nx_n = (long)dimPtr[0];
  P->Ny_n = (long)dimPtr[1];
  P->Nz_n = nDims > 2? (long)dimPtr[2]: 1;
  P->dz_n = *(float *)mxGetData(mxGetField(prhs[1],0,"dz_n"));
  P->rho_e = *(float *)mxGetData(mxGetField(prhs[1],0,"rho_e"));
  P->RoC = *(float *)mxGetData(mxGetField(prhs[1],0,"RoC"));
  P->sinBendDirection = sin(*(float *)mxGetData(mxGetField(prhs[1],0,"bendDirection"))/180*PI);
  P->cosBendDirection = cos(*(float *)mxGetData(mxGetField(prhs[1],0,"bendDirection"))/180*PI);
  P->E1 = (floatcomplex *)mxGetData(prhs[0]); // Input E field
  dimPtr = mxGetDimensions(prhs[0]);
  P->Efinal = (floatcomplex *)mxGetData(plhs[0] = mxCreateNumericArray(2,dimPtr,mxSINGLE_CLASS,mxCOMPLEX)); // Output E field
  P->n_out = (floatcomplex *)mxGetData(plhs[1] = mxCreateNumericArray(2,dimPtr,mxSINGLE_CLASS,mxCOMPLEX)); // Output refractive index
  P->precisePower = (float)mxGetScalar(mxGetField(prhs[1],0,"inputPrecisePower"));
  #ifndef __NVCC__
  P->E2 = (floatcomplex *)((P->iz_end - P->iz_start)%2? P->Efinal: malloc(P->Nx*P->Ny*sizeof(floatcomplex)));
  #endif
  P->multiplier = (float *)mxGetData(mxGetField(prhs[1],0,"multiplier")); // Array of multiplier values to apply to the E field after each step, due to the edge absorber outside the main simulation window
  P->ax = *(floatcomplex *)mxGetData(mxGetField(prhs[1],0,"ax"));
  P->ay = *(floatcomplex *)mxGetData(mxGetField(prhs[1],0,"ay"));
  
  bool ctrlc_caught = false;      // Has a ctrl+c been passed from MATLAB?
  P->EfieldPower = 0;
  P->precisePowerDiff = 0;
  #ifdef __NVCC__
  int temp, nBlocks; gpuErrchk(cudaOccupancyMaxPotentialBlockSize(&nBlocks,&temp,&substep1a,0,0));
  dim3 blockDims(TILE_DIM,TILE_DIM,1);

  struct parameters *P_dev;
  struct debug D_var = {{0.0,0.0,0.0},{0,0,0}};
  struct debug *D = &D_var;
  struct debug *D_dev;
  createDeviceStructs(P,&P_dev,D,&D_dev);
  #else
  #ifdef _OPENMP
  bool useAllCPUs = mxIsLogicalScalarTrue(mxGetField(prhs[1],0,"useAllCPUs"));
  long numThreads = useAllCPUs || omp_get_num_procs() == 1? omp_get_num_procs(): omp_get_num_procs()-1;
  #else
  long numThreads = 1;
  #endif
  P->b = (floatcomplex *)malloc(numThreads*MAX(P->Nx,P->Ny)*sizeof(floatcomplex));
  #ifdef _OPENMP
  #pragma omp parallel num_threads(useAllCPUs || omp_get_num_procs() == 1? omp_get_num_procs(): omp_get_num_procs()-1)
  #endif
  #endif
  {
    for(long iz=P->iz_start; iz<P->iz_end; iz++) {
      if(ctrlc_caught) break;
      
      #ifdef __NVCC__
      substep1a<<<nBlocks, blockDims>>>(P_dev); // xy -> yx
      substep1b<<<nBlocks, blockDims>>>(P_dev); // yx -> yx
      substep2a<<<nBlocks, blockDims>>>(P_dev); // yx -> xy
      substep2b<<<nBlocks, blockDims>>>(P_dev); // xy -> xy
      applyMultiplier<<<nBlocks, blockDims>>>(P_dev,iz,D_dev); // xy -> xy
      #else
      substep1a(P);
      substep1b(P);
      substep2a(P);
      substep2b(P);
      applyMultiplier(P,iz,NULL);
      #endif

      #ifdef _OPENMP
      #pragma omp master
      #endif
      {
        #ifdef __NVCC__
        if(iz+1 < P->iz_end) swapEPointers<<<1,1>>>(P_dev,iz);
        updatePrecisePower<<<1,1>>>(P_dev);
        gpuErrchk(cudaDeviceSynchronize()); // Wait until all kernels have finished
        #else
        if(iz+1 < P->iz_end) swapEPointers(P,iz);
        updatePrecisePower(P);
        #endif
        
        #ifndef __clang__
        if(utIsInterruptPending()) {
          ctrlc_caught = true;
          printf("\nCtrl+C detected, stopping.\n");
        }
        #endif
      }
      #ifdef _OPENMP
      #pragma omp barrier
      #endif
    }
  }
  #ifdef __NVCC__
  gpuErrchk(cudaDeviceSynchronize()); // Wait until all kernels have finished
  retrieveAndFreeDeviceStructs(P,P_dev,D,D_dev);
//   printf("\nDebug: %.18e %.18e %.18e %llu %llu %llu\n          ",D->dbls[0],D->dbls[1],D->dbls[2],D->ulls[0],D->ulls[1],D->ulls[2]);
  #else
  if(P->E1 != mxGetData(prhs[0]) && P->E1 != P->Efinal) free(P->E1); // Part of the reason for checking this is to properly handle ctrl-c cases
  free(P->b);
  #endif
  double *outputPrecisePowerPtr = (double *)mxGetData(plhs[2] = mxCreateDoubleMatrix(1,1,mxREAL));
  *outputPrecisePowerPtr = P->precisePower;
  return;
}
*/



