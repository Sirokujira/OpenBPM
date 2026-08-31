/*
wabpm.cpp

一般化 BPM 伝搬エンジン (広角 Pade(1,1) / 半ベクトル対応, CPU, 倍精度)

  (1 + d+ * P) E^{n+1} = (1 + d- * P) E^n
  P  = Px + Py
  Px = Lx + (1/2) k0^2 (n^2 - n0^2)
  Py = Ly + (1/2) k0^2 (n^2 - n0^2)

  近軸      : d± = ±i dz/(4 k0 n0)
  Pade(1,1) : d± = (1 ± i k0 n0 dz)/(4 k0^2 n0^2)
              (∂E/∂z = -i k0 n0 [(X/2)/(1+X/4)] E, X = P/(k0 n0)^2
               の Crank-Nicolson 離散化)

ADI 分離:
  T = (1 + d- Py) E      (陽的 y)
  E = (1 + d- Px) T      (陽的 x)
  (1 + d+ Px) E = E      (陰的 x : 行毎の三重対角)
  (1 + d+ Py) E = E      (陰的 y : 列毎の三重対角)

半ベクトル (Stern) 差分 : 偏波方向の 2 階微分を
  d/dx[(1/n^2) d(n^2 E)/dx]_i
  = [a_i E_{i-1} + b_i E_i + c_i E_{i+1}] / dx^2
  a_i = 2 n^2_{i-1}/(n^2_i + n^2_{i-1})
  c_i = 2 n^2_{i+1}/(n^2_i + n^2_{i+1})
  b_i = -[2 n^2_i/(n^2_i + n^2_{i-1}) + 2 n^2_i/(n^2_i + n^2_{i+1})]
で置き換える (重みには n^2 の実部を使用)。均一媒質では標準の
ラプラシアン (1, -2, 1) に一致する。境界は Dirichlet (領域外 E = 0)。
*/

#include <complex>
#include <cstdlib>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "bpm/wabpm.h"

typedef std::complex<double> cplx;

// 偏波方向のステンシル重み (pol_dir = 1 なら Stern, 0 なら標準)
//
// 対角項 b は「-側の面の寄与 bm」と「+側の面の寄与 bp」に分けて持つ。
// PML (複素座標伸長) 使用時は面ごとに異なる係数 gm/gp を掛けるため、
// b をまとめてしまうと PML を適用できない。
//   標準  : a = 1, bm = bp = -1, c = 1       (合計 b = -2)
//   Stern : a = 2nm/(nc+nm), bm = -2nc/(nc+nm), c = 2np/(nc+np), bp = -2nc/(nc+np)
// gm/gp が NULL のときは従来と同じ実数係数になる。
static inline void stencil(int pol_dir, const cplx *n2, long i, long stride, long idx, long num,
                           const cplx *gm, const cplx *gp,
                           cplx *a, cplx *b, cplx *c)
{
	double wa, wbm, wbp, wc;
	if (pol_dir) {
		const double nc = n2[i].real();
		const double nm = (idx > 0)       ? n2[i - stride].real() : nc;
		const double np = (idx < num - 1) ? n2[i + stride].real() : nc;
		wa  = 2 * nm / (nc + nm);
		wc  = 2 * np / (nc + np);
		wbm = -(2 * nc / (nc + nm));
		wbp = -(2 * nc / (nc + np));
	}
	else {
		wa  = 1;
		wc  = 1;
		wbm = -1;
		wbp = -1;
	}
	if (gm != NULL) {
		*a = wa * gm[idx];
		*b = (wbm * gm[idx]) + (wbp * gp[idx]);
		*c = wc * gp[idx];
	}
	else {
		*a = wa;
		*b = wbm + wbp;
		*c = wc;
	}
}

// ADI 1 ステップ (共通カーネル) : (1 + dp*P) E' = (1 + dm*P) E を解く
static void adi_step(const wabpm_params *W, cplx *E, const cplx *n2, cplx dp, cplx dm)
{
	const long Nx = W->Nx;
	const long Ny = W->Ny;
	const long N = Nx * Ny;
	const double idx2 = 1.0 / (W->dx * W->dx);
	const double idy2 = 1.0 / (W->dy * W->dy);
	const double k0 = W->k0;
	const double n02 = W->n0 * W->n0;
	const int polx = (W->pol == 1);
	const int poly = (W->pol == 2);
	// 対称境界 : ix=0 / iy=0 の手前 (半セル外) に鏡像面。E[-1] = sgn*E[0]
	const int    hasSx = (W->symx != 0), hasSy = (W->symy != 0);
	const double sgnx = (W->symx == 2) ? -1.0 : 1.0;
	const double sgny = (W->symy == 2) ? -1.0 : 1.0;

	cplx *T = new cplx[N];

	// 陽的 y : T = (1 + dm*Py) E
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (long iy = 0; iy < Ny; iy++) {
		for (long ix = 0; ix < Nx; ix++) {
			const long i = ix + iy * Nx;
			cplx a, b, c;
			stencil(poly, n2, i, Nx, iy, Ny, W->gym, W->gyp, &a, &b, &c);
			cplx lap = b * E[i];
			if (iy > 0)      lap += a * E[i - Nx];
			else if (hasSy)  lap += (a * sgny) * E[i];   // 鏡像セル
			if (iy < Ny - 1) lap += c * E[i + Nx];
			const cplx V2 = 0.5 * k0 * k0 * (n2[i] - n02);
			T[i] = E[i] + dm * (lap * idy2 + V2 * E[i]);
		}
	}

	// 陽的 x : E = (1 + dm*Px) T
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (long iy = 0; iy < Ny; iy++) {
		for (long ix = 0; ix < Nx; ix++) {
			const long i = ix + iy * Nx;
			cplx a, b, c;
			stencil(polx, n2, i, 1, ix, Nx, W->gxm, W->gxp, &a, &b, &c);
			cplx lap = b * T[i];
			if (ix > 0)      lap += a * T[i - 1];
			else if (hasSx)  lap += (a * sgnx) * T[i];   // 鏡像セル
			if (ix < Nx - 1) lap += c * T[i + 1];
			const cplx V2 = 0.5 * k0 * k0 * (n2[i] - n02);
			E[i] = T[i] + dm * (lap * idx2 + V2 * T[i]);
		}
	}

	// 陰的 x : (1 + dp*Px) E' = E を行毎に Thomas 法で解く
#ifdef _OPENMP
#pragma omp parallel
#endif
	{
		cplx *cw = new cplx[Nx > Ny ? Nx : Ny];  // 前進消去した上対角
#ifdef _OPENMP
#pragma omp for
#endif
		for (long iy = 0; iy < Ny; iy++) {
			cplx *e = &E[iy * Nx];
			for (long ix = 0; ix < Nx; ix++) {
				const long i = ix + iy * Nx;
				cplx a, b, c;
				stencil(polx, n2, i, 1, ix, Nx, W->gxm, W->gxp, &a, &b, &c);
				const cplx V2 = 0.5 * k0 * k0 * (n2[i] - n02);
				const cplx sub  = dp * (a * idx2);
				cplx       diag = 1.0 + dp * (b * idx2 + V2);
				const cplx sup  = dp * (c * idx2);
				if (ix > 0) {
					// cw / e の前要素は対角で正規化済み (c', r')
					diag  -= sub * cw[ix - 1];
					e[ix] -= sub * e[ix - 1];
				}
				else if (hasSx) {
					diag += sub * sgnx;   // E[-1] = sgn*E[0] を対角へ畳み込む
				}
				cw[ix] = sup / diag;
				e[ix] /= diag;
			}
			for (long ix = Nx - 2; ix >= 0; ix--) {
				e[ix] -= cw[ix] * e[ix + 1];
			}
		}

		// 陰的 y : (1 + dp*Py) E' = E を列毎に Thomas 法で解く
#ifdef _OPENMP
#pragma omp for
#endif
		for (long ix = 0; ix < Nx; ix++) {
			for (long iy = 0; iy < Ny; iy++) {
				const long i = ix + iy * Nx;
				cplx a, b, c;
				stencil(poly, n2, i, Nx, iy, Ny, W->gym, W->gyp, &a, &b, &c);
				const cplx V2 = 0.5 * k0 * k0 * (n2[i] - n02);
				const cplx sub  = dp * (a * idy2);
				cplx       diag = 1.0 + dp * (b * idy2 + V2);
				const cplx sup  = dp * (c * idy2);
				if (iy > 0) {
					// cw / E の前要素は対角で正規化済み (c', r')
					diag -= sub * cw[iy - 1];
					E[i] -= sub * E[i - Nx];
				}
				else if (hasSy) {
					diag += sub * sgny;   // E[-1] = sgn*E[0] を対角へ畳み込む
				}
				cw[iy] = sup / diag;
				E[i] /= diag;
			}
			for (long iy = Ny - 2; iy >= 0; iy--) {
				E[ix + iy * Nx] -= cw[iy] * E[ix + (iy + 1) * Nx];
			}
		}

		delete[] cw;
	}

	delete[] T;
}

void wabpm_step(const wabpm_params *W, cplx *E, const cplx *n2)
{
	const double k0 = W->k0;
	const double n0 = W->n0;

	cplx dp, dm;
	if (W->wideangle) {
		dp = cplx(1.0,  k0 * n0 * W->dz) / (4 * k0 * k0 * n0 * n0);
		dm = cplx(1.0, -k0 * n0 * W->dz) / (4 * k0 * k0 * n0 * n0);
	}
	else {
		dp = cplx(0.0,  W->dz / (4 * k0 * n0));
		dm = cplx(0.0, -W->dz / (4 * k0 * n0));
	}

	adi_step(W, E, n2, dp, dm);
}

void wabpm_imagdist_step(const wabpm_params *W, cplx *E, const cplx *n2)
{
	// 虚軸伝搬 (z -> -i z) : d± が実数になり、伝搬演算子
	// (1 - a P)^-1 (1 + a P) は P の固有値 (実効屈折率) が大きい
	// モードほど大きな利得を持つ -> べき乗法で基本モードへ収束する
	const double a = W->dz / (4 * W->k0 * W->n0);
	adi_step(W, E, n2, cplx(-a, 0.0), cplx(a, 0.0));
}

void wabpm_apply_P(const wabpm_params *W, const cplx *E, const cplx *n2, cplx *out)
{
	const long Nx = W->Nx;
	const long Ny = W->Ny;
	const double idx2 = 1.0 / (W->dx * W->dx);
	const double idy2 = 1.0 / (W->dy * W->dy);
	const double k0 = W->k0;
	const double n02 = W->n0 * W->n0;
	const int polx = (W->pol == 1);
	const int poly = (W->pol == 2);
	const int    hasSx = (W->symx != 0), hasSy = (W->symy != 0);
	const double sgnx = (W->symx == 2) ? -1.0 : 1.0;
	const double sgny = (W->symy == 2) ? -1.0 : 1.0;

#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (long iy = 0; iy < Ny; iy++) {
		for (long ix = 0; ix < Nx; ix++) {
			const long i = ix + iy * Nx;
			cplx a, b, c;
			stencil(polx, n2, i, 1, ix, Nx, W->gxm, W->gxp, &a, &b, &c);
			cplx lapx = b * E[i];
			if (ix > 0)      lapx += a * E[i - 1];
			else if (hasSx)  lapx += (a * sgnx) * E[i];
			if (ix < Nx - 1) lapx += c * E[i + 1];
			stencil(poly, n2, i, Nx, iy, Ny, W->gym, W->gyp, &a, &b, &c);
			cplx lapy = b * E[i];
			if (iy > 0)      lapy += a * E[i - Nx];
			else if (hasSy)  lapy += (a * sgny) * E[i];
			if (iy < Ny - 1) lapy += c * E[i + Nx];
			out[i] = lapx * idx2 + lapy * idy2 + k0 * k0 * (n2[i] - n02) * E[i];
		}
	}
}
