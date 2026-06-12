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

#endif

#endif
