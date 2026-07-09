#ifndef _WABPM_H_
#define _WABPM_H_

// 一般化 BPM 伝搬エンジン (広角 Pade(1,1) / 半ベクトル対応)
//
// (1 + d+ * P) E^{n+1} = (1 + d- * P) E^n,  P = Lx + Ly + k0^2(n^2 - n0^2)
//   近軸       : d± = ±i dz/(4 k0 n0)
//   Pade(1,1)  : d± = (1 ± i k0 n0 dz)/(4 k0^2 n0^2)
// を Peaceman-Rachford 型 ADI (x/y 分離, 屈折率項は半分ずつ) で解く。
// 偏波方向の x (または y) 微分は Stern の半ベクトル差分
// d/dx[(1/n^2) d(n^2 E)/dx] を用いる (界面での n^2 E 連続を反映)。

#ifdef __cplusplus

#include <complex>

struct wabpm_params {
	int    Nx, Ny;       // 格子数
	double dx, dy, dz;   // 格子幅 [m]
	double k0;           // 真空波数 [1/m]
	double n0;           // 参照屈折率
	int    wideangle;    // 0 = 近軸, 1 = Pade(1,1)
	int    pol;          // 0 = スカラー, 1 = x 偏波 (半ベクトル), 2 = y 偏波
};

// 1 ステップ伝搬 (E, n2 は Nx*Ny, 行優先 [iy*Nx+ix])
// n2 は複素比誘電率 (= 複素屈折率の2乗, 損失は負の虚部)
void wabpm_step(const wabpm_params *W, std::complex<double> *E, const std::complex<double> *n2);

// 虚軸伝搬 1 ステップ (モード解析用)
// (1 - a P) E' = (1 + a P) E, a = dz/(4 k0 n0) を解く。P の固有値が
// 大きい (実効屈折率が高い) モードほど利得が大きく、反復でモードが析出する。
void wabpm_imagdist_step(const wabpm_params *W, std::complex<double> *E, const std::complex<double> *n2);

// ヘルムホルツ演算子 P = Lx + Ly + k0^2(n^2 - n0^2) の適用 (Rayleigh 商用)
void wabpm_apply_P(const wabpm_params *W, const std::complex<double> *E,
                   const std::complex<double> *n2, std::complex<double> *out);

// モード解析 (虚軸伝搬 + Gram-Schmidt 直交化のべき乗法, bpm/modes.cpp)
//   nModes 個の導波モードを実効屈折率の降順で求める。
//   modes : nModes * (Nx*Ny) [出力, L2 正規化済み]
//   neff  : nModes           [出力, Rayleigh 商から算出]
//   W->dz は虚軸ステップ幅 (収束パラメータ)。戻り値は収束したモード数。
int wabpm_find_modes(const wabpm_params *W, const std::complex<double> *n2,
                     int nModes, int maxIter, double tol,
                     std::complex<double> *modes, double *neff);

#endif

// ---------------------------------------------------------------
// CUDA 版 (bpm/wabpm.cu) : デバイス常駐の伝搬コンテキスト
// ---------------------------------------------------------------
#ifdef __CUDACC__

#include <thrust/complex.h>
typedef thrust::complex<double> wabpm_cplx;

struct wabpm_gpu;  // opaque

// コンテキスト生成 (E_host : 初期電界, mult_host : 端部吸収体, いずれも Nx*Ny)
extern struct wabpm_gpu *wabpm_gpu_create(const wabpm_params *W,
                                          const wabpm_cplx *E_host, const float *mult_host);
// 1 ステップ伝搬 (n2_host : このスライスの複素比誘電率) + 吸収体適用
extern void wabpm_gpu_step(struct wabpm_gpu *G, const wabpm_cplx *n2_host);
// 行 iy の電界をホストへ取得 (Ixz 記録用)
extern void wabpm_gpu_get_row(struct wabpm_gpu *G, int iy, wabpm_cplx *row_host);
// 全電界をホストへ取得
extern void wabpm_gpu_get_field(struct wabpm_gpu *G, wabpm_cplx *E_host);
// 解放
extern void wabpm_gpu_destroy(struct wabpm_gpu *G);

#endif

#endif
