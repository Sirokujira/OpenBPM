/*
wabpm.cu (CUDA)

一般化 BPM 伝搬エンジン (広角 Pade(1,1) / 半ベクトル対応) の CUDA 版
bpm/wabpm.cpp と同一のアルゴリズム:

  (1 + d+ * P) E^{n+1} = (1 + d- * P) E^n
  P  = Px + Py,  Px = Lx + (1/2) k0^2 (n^2 - n0^2),  Py = Ly + (1/2) ...
  近軸      : d± = ±i dz/(4 k0 n0)
  Pade(1,1) : d± = (1 ± i k0 n0 dz)/(4 k0^2 n0^2)

陽的スイープは要素並列、陰的スイープは行 (列) 毎に 1 スレッドを割り当てた
Thomas 法 (bpm/FDBPMpropagator.cu の substep1b/2b と同じ流儀)。
半ベクトルは Stern の差分 (重みは n^2 の実部)。境界は Dirichlet。
*/

#include <cstdio>
#include <cstdlib>
#include <cuda_runtime.h>
#include <thrust/complex.h>

#include "bpm/wabpm.h"

typedef thrust::complex<double> cplxd;

static void wabpm_cuda_check(cudaError_t code, const char *file, int line)
{
	if (code != cudaSuccess) {
		fprintf(stderr, "*** CUDA error : %s (%s:%d)\n", cudaGetErrorString(code), file, line);
		exit(1);
	}
}
#define WABPM_CHECK(ans) wabpm_cuda_check((ans), __FILE__, __LINE__)

// カーネルへ値渡しするデバイス用パラメータ
struct wabpm_dev {
	int   Nx, Ny;
	double idx2, idy2;
	double k0, n02;
	cplxd dp, dm;
	int   polx, poly;
	// 対称境界 : ix=0 / iy=0 の手前 (半セル外) に鏡像面。E[-1] = sgn*E[0]
	int    hasSx, hasSy;
	double sgnx, sgny;
};

struct wabpm_gpu {
	struct wabpm_dev D;
	cplxd *d_E;     // 電界
	cplxd *d_T;     // 作業バッファ
	cplxd *d_n2;    // 複素比誘電率 (スライス)
	cplxd *d_cw;    // Thomas 法の上対角バッファ (Nx*Ny)
	float *d_mult;  // 端部吸収体
	float *d_beta;  // TPA 係数 (スライス) [m/W]
	size_t N;
};

// 偏波方向のステンシル重み (pol = 1 なら Stern, 0 なら標準ラプラシアン)
__device__ __forceinline__
static void dstencil(int pol, const cplxd *n2, long i, long stride, long idx, long num,
                     double *a, double *b, double *c)
{
	if (pol) {
		const double nc = n2[i].real();
		const double nm = (idx > 0)       ? n2[i - stride].real() : nc;
		const double np = (idx < num - 1) ? n2[i + stride].real() : nc;
		*a = 2 * nm / (nc + nm);
		*c = 2 * np / (nc + np);
		*b = -((2 * nc / (nc + nm)) + (2 * nc / (nc + np)));
	}
	else {
		*a = 1;
		*b = -2;
		*c = 1;
	}
}

// 陽的 y : T = (1 + dm*Py) E
__global__
static void k_explicit_y(struct wabpm_dev D, const cplxd *E, cplxd *T, const cplxd *n2)
{
	const long N = (long)D.Nx * D.Ny;
	for (long i = blockIdx.x * (long)blockDim.x + threadIdx.x; i < N; i += gridDim.x * (long)blockDim.x) {
		const long iy = i / D.Nx;
		double a, b, c;
		dstencil(D.poly, n2, i, D.Nx, iy, D.Ny, &a, &b, &c);
		cplxd lap = b * E[i];
		if (iy > 0)          lap += a * E[i - D.Nx];
		else if (D.hasSy)    lap += (a * D.sgny) * E[i];   // 鏡像セル
		if (iy < D.Ny - 1)   lap += c * E[i + D.Nx];
		const cplxd V2 = 0.5 * D.k0 * D.k0 * (n2[i] - D.n02);
		T[i] = E[i] + D.dm * (lap * D.idy2 + V2 * E[i]);
	}
}

// 陽的 x : E = (1 + dm*Px) T
__global__
static void k_explicit_x(struct wabpm_dev D, cplxd *E, const cplxd *T, const cplxd *n2)
{
	const long N = (long)D.Nx * D.Ny;
	for (long i = blockIdx.x * (long)blockDim.x + threadIdx.x; i < N; i += gridDim.x * (long)blockDim.x) {
		const long ix = i % D.Nx;
		double a, b, c;
		dstencil(D.polx, n2, i, 1, ix, D.Nx, &a, &b, &c);
		cplxd lap = b * T[i];
		if (ix > 0)          lap += a * T[i - 1];
		else if (D.hasSx)    lap += (a * D.sgnx) * T[i];   // 鏡像セル
		if (ix < D.Nx - 1)   lap += c * T[i + 1];
		const cplxd V2 = 0.5 * D.k0 * D.k0 * (n2[i] - D.n02);
		E[i] = T[i] + D.dm * (lap * D.idx2 + V2 * T[i]);
	}
}

// 陰的 x : (1 + dp*Px) E' = E を行毎に Thomas 法で解く (1 行 = 1 スレッド)
__global__
static void k_implicit_x(struct wabpm_dev D, cplxd *E, const cplxd *n2, cplxd *cw)
{
	for (long iy = blockIdx.x * (long)blockDim.x + threadIdx.x; iy < D.Ny; iy += gridDim.x * (long)blockDim.x) {
		cplxd *e = &E[iy * (long)D.Nx];
		cplxd *w = &cw[iy * (long)D.Nx];
		for (long ix = 0; ix < D.Nx; ix++) {
			const long i = ix + iy * (long)D.Nx;
			double a, b, c;
			dstencil(D.polx, n2, i, 1, ix, D.Nx, &a, &b, &c);
			const cplxd V2 = 0.5 * D.k0 * D.k0 * (n2[i] - D.n02);
			const cplxd sub  = D.dp * (a * D.idx2);
			cplxd       diag = 1.0 + D.dp * (b * D.idx2 + V2);
			const cplxd sup  = D.dp * (c * D.idx2);
			if (ix > 0) {
				// w / e の前要素は対角で正規化済み (c', r')
				diag  -= sub * w[ix - 1];
				e[ix] -= sub * e[ix - 1];
			}
			else if (D.hasSx) {
				diag += sub * D.sgnx;   // E[-1] = sgn*E[0] を対角へ畳み込む
			}
			w[ix] = sup / diag;
			e[ix] /= diag;
		}
		for (long ix = D.Nx - 2; ix >= 0; ix--) {
			e[ix] -= w[ix] * e[ix + 1];
		}
	}
}

// 陰的 y : (1 + dp*Py) E' = E を列毎に Thomas 法で解く (1 列 = 1 スレッド)
__global__
static void k_implicit_y(struct wabpm_dev D, cplxd *E, const cplxd *n2, cplxd *cw)
{
	for (long ix = blockIdx.x * (long)blockDim.x + threadIdx.x; ix < D.Nx; ix += gridDim.x * (long)blockDim.x) {
		cplxd *w = &cw[ix * (long)D.Ny];
		for (long iy = 0; iy < D.Ny; iy++) {
			const long i = ix + iy * (long)D.Nx;
			double a, b, c;
			dstencil(D.poly, n2, i, D.Nx, iy, D.Ny, &a, &b, &c);
			const cplxd V2 = 0.5 * D.k0 * D.k0 * (n2[i] - D.n02);
			const cplxd sub  = D.dp * (a * D.idy2);
			cplxd       diag = 1.0 + D.dp * (b * D.idy2 + V2);
			const cplxd sup  = D.dp * (c * D.idy2);
			if (iy > 0) {
				diag -= sub * w[iy - 1];
				E[i] -= sub * E[i - D.Nx];
			}
			else if (D.hasSy) {
				diag += sub * D.sgny;   // E[-1] = sgn*E[0] を対角へ畳み込む
			}
			w[iy] = sup / diag;
			E[i] /= diag;
		}
		for (long iy = D.Ny - 2; iy >= 0; iy--) {
			E[ix + iy * (long)D.Nx] -= w[iy] * E[ix + (iy + 1) * (long)D.Nx];
		}
	}
}

// 端部吸収体
__global__
static void k_multiplier(struct wabpm_dev D, cplxd *E, const float *mult)
{
	const long N = (long)D.Nx * D.Ny;
	for (long i = blockIdx.x * (long)blockDim.x + threadIdx.x; i < N; i += gridDim.x * (long)blockDim.x) {
		E[i] *= (double)mult[i];
	}
}

// 二光子吸収 (TPA) : E *= exp(-(beta/2)*I*dz), I = |E|^2 [W/m^2]
// (物理スケーリング済みの界を仮定。強度の減衰率 alpha = beta*I に対し界には alpha/2)
__global__
static void k_tpa(struct wabpm_dev D, cplxd *E, const float *beta, double dz)
{
	const long N = (long)D.Nx * D.Ny;
	for (long i = blockIdx.x * (long)blockDim.x + threadIdx.x; i < N; i += gridDim.x * (long)blockDim.x) {
		const double b = beta[i];
		if (b > 0) {
			E[i] *= exp(-0.5 * b * thrust::norm(E[i]) * dz);
		}
	}
}


// z 断面の統計量 (GUI 表示用の /trace) : 近軸版 fieldTrace と同一の集計
__global__
static void k_trace(struct wabpm_dev D, const cplxd *E, const double *xc, const double *yc, double *acc)
{
	const long N = (long)D.Nx * D.Ny;
	double s = 0, sx = 0, sy = 0, sxx = 0, syy = 0, pk = 0;
	for (long i = blockIdx.x * (long)blockDim.x + threadIdx.x; i < N; i += gridDim.x * (long)blockDim.x) {
		const long ix = i % D.Nx;
		const long iy = i / D.Nx;
		const double iv = thrust::norm(E[i]);
		s   += iv;
		sx  += xc[ix] * iv;
		sy  += yc[iy] * iv;
		sxx += xc[ix] * xc[ix] * iv;
		syy += yc[iy] * yc[iy] * iv;
		if (iv > pk) pk = iv;
	}
	atomicAdd(&acc[0], s);
	atomicAdd(&acc[1], sx);
	atomicAdd(&acc[2], sy);
	atomicAdd(&acc[3], sxx);
	atomicAdd(&acc[4], syy);
	atomicMax((unsigned long long *)&acc[5], (unsigned long long)__double_as_longlong(pk));
}

struct wabpm_gpu *wabpm_gpu_create(const wabpm_params *W,
                                   const wabpm_cplx *E_host, const float *mult_host)
{
	struct wabpm_gpu *G = (struct wabpm_gpu *)malloc(sizeof(struct wabpm_gpu));

	G->D.Nx = W->Nx;
	G->D.Ny = W->Ny;
	G->D.idx2 = 1.0 / (W->dx * W->dx);
	G->D.idy2 = 1.0 / (W->dy * W->dy);
	G->D.k0 = W->k0;
	G->D.n02 = W->n0 * W->n0;
	G->D.polx = (W->pol == 1);
	G->D.poly = (W->pol == 2);
	G->D.hasSx = (W->symx != 0);
	G->D.hasSy = (W->symy != 0);
	G->D.sgnx = (W->symx == 2) ? -1.0 : 1.0;
	G->D.sgny = (W->symy == 2) ? -1.0 : 1.0;
	if (W->wideangle) {
		G->D.dp = cplxd(1.0,  W->k0 * W->n0 * W->dz) / (4 * W->k0 * W->k0 * W->n0 * W->n0);
		G->D.dm = cplxd(1.0, -W->k0 * W->n0 * W->dz) / (4 * W->k0 * W->k0 * W->n0 * W->n0);
	}
	else {
		G->D.dp = cplxd(0.0,  W->dz / (4 * W->k0 * W->n0));
		G->D.dm = cplxd(0.0, -W->dz / (4 * W->k0 * W->n0));
	}

	G->N = (size_t)W->Nx * W->Ny;
	WABPM_CHECK(cudaMalloc(&G->d_E,    G->N * sizeof(cplxd)));
	WABPM_CHECK(cudaMalloc(&G->d_T,    G->N * sizeof(cplxd)));
	WABPM_CHECK(cudaMalloc(&G->d_n2,   G->N * sizeof(cplxd)));
	WABPM_CHECK(cudaMalloc(&G->d_cw,   G->N * sizeof(cplxd)));
	WABPM_CHECK(cudaMalloc(&G->d_mult, G->N * sizeof(float)));
	WABPM_CHECK(cudaMalloc(&G->d_beta, G->N * sizeof(float)));
	WABPM_CHECK(cudaMemcpy(G->d_E,    E_host,    G->N * sizeof(cplxd), cudaMemcpyHostToDevice));
	WABPM_CHECK(cudaMemcpy(G->d_mult, mult_host, G->N * sizeof(float), cudaMemcpyHostToDevice));

	return G;
}

// z 断面の統計量を集計してホストへ返す (acc_host は 6 要素)
// xc_dev / yc_dev は呼び出し側が確保・転送したデバイス側の座標配列
void wabpm_gpu_trace(struct wabpm_gpu *G, const double *xc_dev, const double *yc_dev,
                     double *acc_dev, double *acc_host)
{
	WABPM_CHECK(cudaMemset(acc_dev, 0, 6 * sizeof(double)));
	const int tpb = 256;
	const int nbe = (int)((G->N + tpb - 1) / tpb);
	k_trace <<<nbe, tpb>>>(G->D, G->d_E, xc_dev, yc_dev, acc_dev);
	WABPM_CHECK(cudaGetLastError());
	WABPM_CHECK(cudaMemcpy(acc_host, acc_dev, 6 * sizeof(double), cudaMemcpyDeviceToHost));
}

// TPA を 1 ステップ分適用する (beta_host : このスライスの TPA 係数 [m/W], Nx*Ny)
void wabpm_gpu_tpa(struct wabpm_gpu *G, const float *beta_host, double dz)
{
	WABPM_CHECK(cudaMemcpy(G->d_beta, beta_host, G->N * sizeof(float), cudaMemcpyHostToDevice));
	const int tpb = 256;
	const int nbe = (int)((G->N + tpb - 1) / tpb);
	k_tpa <<<nbe, tpb>>>(G->D, G->d_E, G->d_beta, dz);
	WABPM_CHECK(cudaGetLastError());
	WABPM_CHECK(cudaDeviceSynchronize());
}

void wabpm_gpu_step(struct wabpm_gpu *G, const wabpm_cplx *n2_host)
{
	WABPM_CHECK(cudaMemcpy(G->d_n2, n2_host, G->N * sizeof(cplxd), cudaMemcpyHostToDevice));

	const int tpb = 256;
	const int nbe = (int)((G->N + tpb - 1) / tpb);
	const int tpl = 64;
	const int nbx = (G->D.Ny + tpl - 1) / tpl;
	const int nby = (G->D.Nx + tpl - 1) / tpl;

	k_explicit_y <<<nbe, tpb>>>(G->D, G->d_E, G->d_T, G->d_n2);
	k_explicit_x <<<nbe, tpb>>>(G->D, G->d_E, G->d_T, G->d_n2);
	k_implicit_x <<<nbx, tpl>>>(G->D, G->d_E, G->d_n2, G->d_cw);
	k_implicit_y <<<nby, tpl>>>(G->D, G->d_E, G->d_n2, G->d_cw);
	k_multiplier <<<nbe, tpb>>>(G->D, G->d_E, G->d_mult);
	WABPM_CHECK(cudaGetLastError());
	WABPM_CHECK(cudaDeviceSynchronize());
}

void wabpm_gpu_get_row(struct wabpm_gpu *G, int iy, wabpm_cplx *row_host)
{
	WABPM_CHECK(cudaMemcpy(row_host, G->d_E + (size_t)iy * G->D.Nx,
	                       G->D.Nx * sizeof(cplxd), cudaMemcpyDeviceToHost));
}

void wabpm_gpu_get_field(struct wabpm_gpu *G, wabpm_cplx *E_host)
{
	WABPM_CHECK(cudaMemcpy(E_host, G->d_E, G->N * sizeof(cplxd), cudaMemcpyDeviceToHost));
}

void wabpm_gpu_destroy(struct wabpm_gpu *G)
{
	cudaFree(G->d_E);
	cudaFree(G->d_T);
	cudaFree(G->d_n2);
	cudaFree(G->d_cw);
	cudaFree(G->d_mult);
	cudaFree(G->d_beta);
	free(G);
}
