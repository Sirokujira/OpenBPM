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

// ADI の陽的部分 : E = (1 + dm*Px)(1 + dm*Py) E (T は Nx*Ny の作業領域)
//   polx / poly : x / y 方向の 2 階微分に Stern 差分を使うか
static void adi_explicit(const wabpm_params *W, int polx, int poly,
                         cplx *E, const cplx *n2, cplx dm, cplx *T)
{
	const long Nx = W->Nx;
	const long Ny = W->Ny;
	const double idx2 = 1.0 / (W->dx * W->dx);
	const double idy2 = 1.0 / (W->dy * W->dy);
	const double k0 = W->k0;
	const double n02 = W->n0 * W->n0;
	// 対称境界 : ix=0 / iy=0 の手前 (半セル外) に鏡像面。E[-1] = sgn*E[0]
	const int    hasSx = (W->symx != 0), hasSy = (W->symy != 0);
	const double sgnx = (W->symx == 2) ? -1.0 : 1.0;
	const double sgny = (W->symy == 2) ? -1.0 : 1.0;

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
}

// ADI の陰的部分 : (1 + dp*Py)(1 + dp*Px) E' = E を解く (E を上書き)
static void adi_implicit(const wabpm_params *W, int polx, int poly,
                         cplx *E, const cplx *n2, cplx dp)
{
	const long Nx = W->Nx;
	const long Ny = W->Ny;
	const double idx2 = 1.0 / (W->dx * W->dx);
	const double idy2 = 1.0 / (W->dy * W->dy);
	const double k0 = W->k0;
	const double n02 = W->n0 * W->n0;
	const int    hasSx = (W->symx != 0), hasSy = (W->symy != 0);
	const double sgnx = (W->symx == 2) ? -1.0 : 1.0;
	const double sgny = (W->symy == 2) ? -1.0 : 1.0;

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
}

// ADI 1 ステップ (共通カーネル) : (1 + dp*P) E' = (1 + dm*P) E を解く
static void adi_step(const wabpm_params *W, cplx *E, const cplx *n2, cplx dp, cplx dm)
{
	const int polx = (W->pol == 1);
	const int poly = (W->pol == 2);
	cplx *T = new cplx[(size_t)W->Nx * W->Ny];
	adi_explicit(W, polx, poly, E, n2, dm, T);
	adi_implicit(W, polx, poly, E, n2, dp);
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

void wabpm_imagdist_step_be(const wabpm_params *W, cplx *E, const cplx *n2)
{
	// 後退オイラー : (1 - a Px)(1 - a Py) E' = E (陽的部分なし = dm = 0)
	const double a = W->dz / (4 * W->k0 * W->n0);
	adi_step(W, E, n2, cplx(-a, 0.0), cplx(0.0, 0.0));
}

// 対角ブロックの演算子 P = Lx + Ly + k0^2(n^2 - n0^2) の適用 (polx/poly は Stern の有無)
static void apply_diag(const wabpm_params *W, int polx, int poly,
                       const cplx *E, const cplx *n2, cplx *out)
{
	const long Nx = W->Nx;
	const long Ny = W->Ny;
	const double idx2 = 1.0 / (W->dx * W->dx);
	const double idy2 = 1.0 / (W->dy * W->dy);
	const double k0 = W->k0;
	const double n02 = W->n0 * W->n0;
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

void wabpm_apply_P(const wabpm_params *W, const cplx *E, const cplx *n2, cplx *out)
{
	apply_diag(W, (W->pol == 1), (W->pol == 2), E, n2, out);
}

// ===================================================================
// 完全ベクトル BPM (pol = WABPM_POL_FULL)
//
//   d/dz [Ex; Ey] = (-i/(2 k0 n0)) [[Pxx, Pxy], [Pyx, Pyy]] [Ex; Ey]
//   Pxx Ex = d/dx[(1/n^2) d(n^2 Ex)/dx] + d2Ex/dy2 + k0^2(n^2 - n0^2) Ex
//   Pyy Ey = d2Ey/dx2 + d/dy[(1/n^2) d(n^2 Ey)/dy] + k0^2(n^2 - n0^2) Ey
//   Pxy Ey = d/dx[(1/n^2) d(n^2 Ey)/dy] - d2Ey/dxdy
//   Pyx Ex = d/dy[(1/n^2) d(n^2 Ex)/dx] - d2Ex/dydx
// (横方向のベクトル波動方程式 ∇²E + k²n²E + ∇(E·∇ln n²) = 0 の x, y 成分)
//
// 対角ブロック Pxx / Pyy は半ベクトル (Stern) と同じ演算子で、従来の ADI で解く。
// 交差項 C = [[0, Pxy], [Pyx, 0]] は中心差分 (4 隅の点) で離散化する:
//   Pxy Ey|ij = Σ_{sx,sy=±1} sx sy (n²_{i+sx,j+sy}/n²_{i+sx,j} - 1) Ey_{i+sx,j+sy} / (4 dx dy)
//   Pyx Ex|ij = Σ_{sx,sy=±1} sx sy (n²_{i+sx,j+sy}/n²_{i,j+sy} - 1) Ex_{i+sx,j+sy} / (4 dx dy)
// 一様媒質では係数が厳密に 0 になり、Ex / Ey は独立なスカラー伝搬に帰着する。
//
// 時間積分は (1 + dp(D + C)) E' = (1 + dm(D + C)) E を
//   (1 + dp D) E'^{k+1} = (1 + dm D) E + dm C E - dp C E'^k
// の不動点反復で解く (左辺の (1 + dp D) は従来の ADI 分解)。収束すれば
// 交差項を含めた Crank-Nicolson (広角では Pade(1,1)) と同じ離散式になる。
// 交差項は屈折率が x, y 両方向に変化する所 (コアの角・円周) にしか無く小さいため
// 反復は数回で収束する。
//
// 制約 (呼び出し側で弾く): 対称境界 (symx/symy) は未対応 (Ex と Ey で鏡像の
// パリティが逆になるため)。PML の伸長係数は対角ブロックにのみ掛ける
// (交差項は n が一様な PML 層内では 0)。重みは n^2 の実部を使う (Stern と同じ)。
// ===================================================================

// 交差項 cx = Pxy Ey, cy = Pyx Ex
static void apply_cross(const wabpm_params *W, const cplx *n2,
                        const cplx *Ex, const cplx *Ey, cplx *cx, cplx *cy)
{
	const long Nx = W->Nx;
	const long Ny = W->Ny;
	const double s = 1.0 / (4 * W->dx * W->dy);

#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (long iy = 0; iy < Ny; iy++) {
		for (long ix = 0; ix < Nx; ix++) {
			const long i = ix + iy * Nx;
			cplx ax = 0, ay = 0;
			for (int sx = -1; sx <= 1; sx += 2) {
				const long jx = ix + sx;
				if ((jx < 0) || (jx >= Nx)) continue;
				for (int sy = -1; sy <= 1; sy += 2) {
					const long jy = iy + sy;
					if ((jy < 0) || (jy >= Ny)) continue;
					const long   k   = jx + jy * Nx;           // 隅の点
					const double nk  = n2[k].real();
					const double wxy = nk / n2[jx + iy * Nx].real() - 1.0;  // /n²_{i+sx,j}
					const double wyx = nk / n2[ix + jy * Nx].real() - 1.0;  // /n²_{i,j+sy}
					const double sg  = (double)(sx * sy);
					ax += (sg * wxy) * Ey[k];
					ay += (sg * wyx) * Ex[k];
				}
			}
			cx[i] = ax * s;
			cy[i] = ay * s;
		}
	}
}

// 交差項が存在するか。係数 n²_k/n²_l - 1 は隣接点で n² が異なると 0 でなくなる
// (x のみに変化するスラブでも、界が y に変化すれば Pyx Ex = (∂y Ex)(∂x ln n²) が残る)
// ため、Re(n²) が全点で等しい一様媒質のときだけ 0 になる。
static int has_cross(const wabpm_params *W, const cplx *n2)
{
	const long N = (long)W->Nx * W->Ny;
	const double r0 = n2[0].real();
	for (long i = 1; i < N; i++) {
		if (n2[i].real() != r0) return 1;
	}
	return 0;
}

static int fv_adi(const wabpm_params *W, cplx *Ex, cplx *Ey, const cplx *n2, cplx dp, cplx dm)
{
	const size_t N = (size_t)W->Nx * W->Ny;
	cplx *T  = new cplx[N];
	cplx *Rx = new cplx[N];
	cplx *Ry = new cplx[N];

	// 右辺の対角部分 (1 + dm D) E
	for (size_t i = 0; i < N; i++) { Rx[i] = Ex[i]; Ry[i] = Ey[i]; }
	if (dm != cplx(0.0, 0.0)) {
		adi_explicit(W, 1, 0, Rx, n2, dm, T);
		adi_explicit(W, 0, 1, Ry, n2, dm, T);
	}

	if (!has_cross(W, n2)) {
		// 交差項なし : Ex は Pxx, Ey は Pyy の半ベクトル伝搬そのもの
		adi_implicit(W, 1, 0, Rx, n2, dp);
		adi_implicit(W, 0, 1, Ry, n2, dp);
		for (size_t i = 0; i < N; i++) { Ex[i] = Rx[i]; Ey[i] = Ry[i]; }
		delete[] T; delete[] Rx; delete[] Ry;
		return 1;
	}

	cplx *cx = new cplx[N];
	cplx *cy = new cplx[N];
	cplx *Nxf = new cplx[N];
	cplx *Nyf = new cplx[N];

	// 右辺の交差項 dm C E
	if (dm != cplx(0.0, 0.0)) {
		apply_cross(W, n2, Ex, Ey, cx, cy);
		for (size_t i = 0; i < N; i++) { Rx[i] += dm * cx[i]; Ry[i] += dm * cy[i]; }
	}

	// 不動点反復 (初期推定 E'^0 = E)
	const int    maxit = 50;
	const double tol2  = 1e-24;   // 相対変化 1e-12 (倍精度の丸め近傍)
	int it;
	for (it = 1; it <= maxit; it++) {
		apply_cross(W, n2, Ex, Ey, cx, cy);
		for (size_t i = 0; i < N; i++) {
			Nxf[i] = Rx[i] - dp * cx[i];
			Nyf[i] = Ry[i] - dp * cy[i];
		}
		adi_implicit(W, 1, 0, Nxf, n2, dp);
		adi_implicit(W, 0, 1, Nyf, n2, dp);
		double dif = 0, nrm = 0;
		for (size_t i = 0; i < N; i++) {
			dif += std::norm(Nxf[i] - Ex[i]) + std::norm(Nyf[i] - Ey[i]);
			nrm += std::norm(Nxf[i]) + std::norm(Nyf[i]);
			Ex[i] = Nxf[i];
			Ey[i] = Nyf[i];
		}
		if (dif <= tol2 * nrm) break;
	}

	delete[] T; delete[] Rx; delete[] Ry;
	delete[] cx; delete[] cy; delete[] Nxf; delete[] Nyf;
	return (it <= maxit) ? it : -maxit;
}

int wabpm_fv_step(const wabpm_params *W, cplx *Ex, cplx *Ey, const cplx *n2)
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
	return fv_adi(W, Ex, Ey, n2, dp, dm);
}

int wabpm_fv_imagdist_step(const wabpm_params *W, cplx *Ex, cplx *Ey, const cplx *n2, int backward_euler)
{
	const double a = W->dz / (4 * W->k0 * W->n0);
	return fv_adi(W, Ex, Ey, n2, cplx(-a, 0.0), backward_euler ? cplx(0.0, 0.0) : cplx(a, 0.0));
}

void wabpm_fv_apply_P(const wabpm_params *W, const cplx *Ex, const cplx *Ey,
                      const cplx *n2, cplx *outx, cplx *outy)
{
	const size_t N = (size_t)W->Nx * W->Ny;
	cplx *cx = new cplx[N];
	cplx *cy = new cplx[N];
	apply_diag(W, 1, 0, Ex, n2, outx);
	apply_diag(W, 0, 1, Ey, n2, outy);
	apply_cross(W, n2, Ex, Ey, cx, cy);
	for (size_t i = 0; i < N; i++) {
		outx[i] += cx[i];
		outy[i] += cy[i];
	}
	delete[] cx;
	delete[] cy;
}
