/*
test_fdbpm.c

スカラー近軸 BPM 伝搬カーネル (bpm/FDBPMpropagator.c, Douglas-Gunn ADI) の
単体テスト。既定の伝搬経路でありながら単体検証が無かったため追加した。
外部テストフレームワークに依存しない自己完結ハーネス。

カーネルは C (float _Complex) で書かれているため、型・リンケージを合わせて
このテストも C で記述する。

検証内容 (いずれも解析解との比較):
  1. 自由空間ガウシアンビームの回折 : w(z) = w0 sqrt(1 + (z/zR)^2)
                                      (zR = pi w0^2 n / lambda)
  2. 損失なし自由空間でのエネルギー保存 : sum|E|^2 = const
  3. 一様吸収媒質での強度減衰       : P(z)/P(0) = exp(-2 k0 n'' z)
  4. ビーム中心の保持 (対称な初期界は中心を移動しない)
  5. 曲げ (等価屈折率法) の偏向            : <x>(z) = z^2/(2 RoC)
                                      (applyMultiplier 内の n_bend 変換を通す)
  6. PML (複素座標伸長) の境界反射        : 壁の無い広領域の参照計算と内部窓の
                                      残留電力を比べる。PML の吸収を
                                      fieldCorrection (電力の再正規化) が
                                      打ち消さないこと (電力簿記の同期) も確認する

呼び出し順は sol/solve_bpm.cpp と同一 (substep1a/1b/2a/2b →
applyMultiplier → swapEPointers → updatePrecisePower)。
失敗が 1 件でもあれば終了コード 1 を返す (CTest 連携用)。
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bpm/bpm_prototype.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static int g_failures = 0;

static void check_close(const char *name, double got, double expected, double rtol)
{
	const double denom = (fabs(expected) > 1e-300) ? fabs(expected) : 1.0;
	const double relerr = fabs(got - expected) / denom;
	const int ok = (relerr <= rtol);
	printf("[%s] %-44s got=%.6g expected=%.6g relerr=%.3g (tol=%.3g)\n",
	       ok ? " OK " : "FAIL", name, got, expected, relerr, rtol);
	if (!ok) g_failures++;
}

/* 伝搬に必要な一式 (sol/solve_bpm.cpp と同じ組み立て) */
struct rig {
	struct parameters P;
	floatcomplex *E1_orig, *E2_orig;
	double dx, dy, dz, k0, n0;
	int Nx, Ny, Nz;
};

/* 一様媒質 (屈折率 nr + i*ni、ni > 0 が吸収) を伝搬する準備 */
static void rig_init(struct rig *R, int Nx, int Ny, int Nz,
                     double dx, double dz, double lambda,
                     double nr, double ni, double w0)
{
	struct parameters *P = &R->P;
	const long N = (long)Nx * Ny;

	R->Nx = Nx; R->Ny = Ny; R->Nz = Nz;
	R->dx = R->dy = dx; R->dz = dz;
	R->k0 = 2 * M_PI / lambda;
	R->n0 = nr;

	memset(P, 0, sizeof(*P));
	P->Nx = Nx;
	P->Ny = Ny;
	P->dx = (float)dx;
	P->dy = (float)dx;
	P->dz = (float)dz;
	P->iz_start = 0;
	P->iz_end = Nz;
	P->xSymmetry = 0;
	P->ySymmetry = 0;
	P->n_0 = (float)nr;
	P->d = (float)(-dz * R->k0);
	P->Nx_n = Nx;
	P->Ny_n = Ny;
	P->Nz_n = 1;
	P->dz_n = (float)dz;
	P->taperPerStep = 0;
	P->twistPerStep = 0;
	P->rho_e = 0;
	P->RoC = (float)INFINITY;   /* 直線 (曲げなし) */
	P->sinBendDirection = 0;
	P->cosBendDirection = 0;
	/* ADI 係数 : ax = dz/(4i dx^2 k0 n0) */
	P->ax = -I * (float)(dz / (4 * dx * dx * R->k0 * nr));
	P->ay = -I * (float)(dz / (4 * dx * dx * R->k0 * nr));

	P->n_in      = (floatcomplex *)malloc(N * sizeof(floatcomplex));
	P->E1        = (floatcomplex *)malloc(N * sizeof(floatcomplex));
	P->E2        = (floatcomplex *)malloc(N * sizeof(floatcomplex));
	P->Efinal    = (floatcomplex *)malloc(N * sizeof(floatcomplex));
	P->n_out     = (floatcomplex *)malloc(N * sizeof(floatcomplex));
	P->multiplier = (float *)malloc(N * sizeof(float));
	P->b         = (floatcomplex *)malloc((Nx > Ny ? Nx : Ny) * sizeof(floatcomplex));
	R->E1_orig = P->E1;
	R->E2_orig = P->E2;

	/* 一様屈折率 (損失は +imag で保持する規約) と初期ガウシアン */
	double power = 0;
	for (int iy = 0; iy < Ny; iy++) {
		for (int ix = 0; ix < Nx; ix++) {
			const long i = ix + (long)iy * Nx;
			const double x = dx * (ix - (Nx - 1) / 2.0);
			const double y = dx * (iy - (Ny - 1) / 2.0);
			const float amp = (float)exp(-(x * x + y * y) / (w0 * w0));
			floatcomplex n = (float)nr + (float)ni * I;
			floatcomplex e = amp + 0.0f * I;
			P->n_in[i] = n;
			P->E1[i] = e;
			P->multiplier[i] = 1.0f;   /* 端部吸収体なし (解析解と比較するため) */
			power += (double)amp * amp;
		}
	}
	P->precisePower = power;
	P->EfieldPower = 0;
	P->precisePowerDiff = 0;
}

/* Nz ステップ伝搬する (呼び出し順は sol/solve_bpm.cpp と同一) */
static void rig_run(struct rig *R)
{
	struct parameters *P = &R->P;
	for (long iz = P->iz_start; iz < P->iz_end; iz++) {
		substep1a(P);
		substep1b(P);
		substep2a(P);
		substep2b(P);
		applyMultiplier(P, iz, NULL);
		if (iz + 1 < P->iz_end) swapEPointers(P, iz);
		updatePrecisePower(P);
	}
	/* 最終界は E2 にある */
	if (P->E2 != P->Efinal) {
		memcpy(P->Efinal, P->E2, (size_t)P->Nx * P->Ny * sizeof(floatcomplex));
	}
}

static void rig_free(struct rig *R)
{
	struct parameters *P = &R->P;
	/* swapEPointers でポインタが入れ替わるため重複を除いて解放する */
	floatcomplex *bufs[4];
	int nbuf = 0;
	bufs[nbuf++] = R->E1_orig;
	bufs[nbuf++] = R->E2_orig;
	bufs[nbuf++] = P->E1;
	bufs[nbuf++] = P->E2;
	for (int i = 0; i < nbuf; i++) {
		int dup = 0;
		for (int j = 0; j < i; j++) if (bufs[j] == bufs[i]) dup = 1;
		if (!dup && bufs[i] != NULL && bufs[i] != P->Efinal) free(bufs[i]);
	}
	free(P->Efinal);
	free(P->n_in);
	free(P->n_out);
	free(P->multiplier);
	free(P->b);
}

/* 最終界の 2 次モーメント幅 w = 2*sigma [m] と全電力・重心を返す */
static void field_stats(const struct rig *R, double *w, double *power, double *xbar)
{
	const struct parameters *P = &R->P;
	double p = 0, sx = 0, sxx = 0;
	for (int iy = 0; iy < P->Ny; iy++) {
		for (int ix = 0; ix < P->Nx; ix++) {
			const long i = ix + (long)iy * P->Nx;
			const floatcomplex e = P->Efinal[i];
			const double in = (double)CREALF(e) * CREALF(e) + (double)CIMAGF(e) * CIMAGF(e);
			const double x = R->dx * (ix - (P->Nx - 1) / 2.0);
			p += in;
			sx += in * x;
			sxx += in * x * x;
		}
	}
	const double xm = (p > 0) ? sx / p : 0;
	const double var = (p > 0) ? (sxx / p - xm * xm) : 0;
	*w = 2 * sqrt(var > 0 ? var : 0);
	*power = p;
	*xbar = xm;
}

/* 1. 自由空間回折 + 2. エネルギー保存 + 4. 中心保持 */
static void test_free_space_diffraction(void)
{
	const double lambda = 1.55e-6, n = 1.444, w0 = 5e-6;
	const int Nx = 128, Ny = 128, Nz = 60;
	const double dx = 40e-6 / Nx, dz = 1e-6;

	struct rig R;
	rig_init(&R, Nx, Ny, Nz, dx, dz, lambda, n, 0.0, w0);
	const double p0 = R.P.precisePower;
	rig_run(&R);

	double w, p, xbar;
	field_stats(&R, &w, &p, &xbar);

	const double z = dz * Nz;
	const double zR = M_PI * w0 * w0 * n / lambda;
	const double w_ana = w0 * sqrt(1 + (z / zR) * (z / zR));

	check_close("fdbpm: free-space w(z) (2nd moment)", w, w_ana, 0.03);
	check_close("fdbpm: energy conservation", p / p0, 1.0, 1e-3);
	check_close("fdbpm: beam center preserved", xbar * 1e6 + 100.0, 100.0, 1e-3);

	rig_free(&R);
}

/* 3. 一様吸収媒質での減衰 */
static void test_uniform_absorption(void)
{
	const double lambda = 1.55e-6, n = 1.444, ni = 1e-4, w0 = 8e-6;
	const int Nx = 96, Ny = 96, Nz = 50;
	const double dx = 40e-6 / Nx, dz = 2e-6;

	struct rig R;
	rig_init(&R, Nx, Ny, Nz, dx, dz, lambda, n, ni, w0);
	const double p0 = R.P.precisePower;
	rig_run(&R);

	double w, p, xbar;
	field_stats(&R, &w, &p, &xbar);

	const double z = dz * Nz;
	const double k0 = 2 * M_PI / lambda;
	const double ratio_ana = exp(-2 * k0 * ni * z);

	check_close("fdbpm: absorption P(z)/P(0)", p / p0, ratio_ana, 0.02);
	/* カーネル内部の電力簿記 (precisePower) も同じ減衰を追うこと */
	check_close("fdbpm: precisePower bookkeeping", R.P.precisePower / p0, ratio_ana, 0.02);

	rig_free(&R);
}

/* 5. 曲げ (等価屈折率法) : 一様媒質 + n_bend = n0*exp(x/RoC) (rho_e = 0) では
 *    ビーム重心が円弧を描く : <x>(z) = z^2/(2*RoC) (Ehrenfest の定理)。
 *    test_wabpm の bend_deflection と同じ条件で、既定経路 (近軸カーネル) の
 *    applyMultiplier 内の曲げ変換を検証する。許容 2% は test_wabpm と同じ
 *    (実測 -0.75%) */
static void test_bend_deflection(void)
{
	const double lambda = 1.55e-6, n = 1.444, w0 = 8e-6;
	const int Nx = 192, Ny = 192, Nz = 400;
	/* dx = 0.75um (w0 の約 1/10)。誤差は横方向の離散化が支配し dx^2 で縮む :
	 *   dx = 1.5um で -2.2%、0.75um で -0.75% (dz を半分にしてもほぼ不変)。
	 *   領域幅 (144um / 288um) は結果に影響しない (偏向 4um に対し十分広い) */
	const double dx = 144e-6 / Nx, dz = 1e-6;
	const double RoC = 20e-3;

	struct rig R;
	rig_init(&R, Nx, Ny, Nz, dx, dz, lambda, n, 0.0, w0);
	R.P.RoC = (float)RoC;
	R.P.rho_e = 0.0f;
	R.P.cosBendDirection = 1.0f;   /* +x 方向へ曲げる */
	R.P.sinBendDirection = 0.0f;
	rig_run(&R);

	double w, p, xbar;
	field_stats(&R, &w, &p, &xbar);
	const double z = dz * Nz;
	check_close("fdbpm: bend deflection z^2/(2 RoC)", xbar, z * z / (2 * RoC), 0.02);
	rig_free(&R);
}

/* PML の面係数 (sol/solve_bpm.cpp の bpmPmlAxis と同じ規約)
 *   s(u) = 1 - i*sigma(u), sigma(d) = smax*(d/W)^3
 *   gm[i] = 1/(s(cent_i) s(face_{i-1/2})), gp[i] = 1/(s(cent_i) s(face_{i+1/2})) */
static double pml_sigma(double u, double lo, double hi, double W, double smax)
{
	double d = 0;
	if (u < lo)      d = lo - u;
	else if (u > hi) d = u - hi;
	if (d <= 0) return 0;
	if (d > W) d = W;
	return smax * pow(d / W, 3.0);
}

static void pml_axis(int N, double d, double W, double smax, floatcomplex *gm, floatcomplex *gp)
{
	const double half = 0.5 * N * d;
	for (int i = 0; i < N; i++) {
		const double uc = d * (i - (N - 1) / 2.0);
		const double sc = pml_sigma(uc,           -half + W, half - W, W, smax);
		const double sm = pml_sigma(uc - 0.5 * d, -half + W, half - W, W, smax);
		const double sp = pml_sigma(uc + 0.5 * d, -half + W, half - W, W, smax);
		/* 1/((1 - i a)(1 - i b)) = 1/((1 - ab) - i(a + b)) */
		double re = 1 - sc * sm, im = -(sc + sm), den = re * re + im * im;
		gm[i] = (float)(re / den) + (float)(-im / den) * I;
		re = 1 - sc * sp; im = -(sc + sp); den = re * re + im * im;
		gp[i] = (float)(re / den) + (float)(-im / den) * I;
	}
}

/* 15 度傾けたガウシアンを伝搬させ、内部窓 |x|,|y| < win の残留電力を返す */
static double tilt_run(int Nx, int Ny, double pmlW, double win, double *pin)
{
	const double lambda = 1.55e-6, n = 1.444, w0 = 4e-6, tilt = 15.0;
	const double dx = 0.25e-6, dz = 0.5e-6;
	const int Nz = 400;
	struct rig R;
	rig_init(&R, Nx, Ny, Nz, dx, dz, lambda, n, 0.0, w0);
	/* 傾きの位相ランプ (solve_bpm と同じ符号 : ph = -kt*x) */
	const double kt = R.k0 * n * sin(tilt * M_PI / 180.0);
	for (int iy = 0; iy < Ny; iy++) {
		for (int ix = 0; ix < Nx; ix++) {
			const long i = ix + (long)iy * Nx;
			const double x = dx * (ix - (Nx - 1) / 2.0);
			const float a = CREALF(R.P.E1[i]);
			R.P.E1[i] = (float)(a * cos(-kt * x)) + (float)(a * sin(-kt * x)) * I;
		}
	}
	floatcomplex *g[4] = {NULL, NULL, NULL, NULL};
	if (pmlW > 0) {
		const double smax = -(3.0 + 1.0) * log(1e-20) / (2 * R.k0 * n * pmlW);
		for (int k = 0; k < 4; k++) g[k] = (floatcomplex *)malloc((size_t)((k < 2) ? Nx : Ny) * sizeof(floatcomplex));
		pml_axis(Nx, dx, pmlW, smax, g[0], g[1]);
		pml_axis(Ny, dx, pmlW, smax, g[2], g[3]);
		R.P.pmlxm = g[0]; R.P.pmlxp = g[1]; R.P.pmlym = g[2]; R.P.pmlyp = g[3];
	}
	*pin = R.P.precisePower;
	rig_run(&R);
	double p = 0;
	for (int iy = 0; iy < Ny; iy++) {
		for (int ix = 0; ix < Nx; ix++) {
			const double x = dx * (ix - (Nx - 1) / 2.0);
			const double y = dx * (iy - (Ny - 1) / 2.0);
			const floatcomplex e = R.P.Efinal[ix + (long)iy * Nx];
			if (fabs(x) < win && fabs(y) < win) p += (double)CREALF(e) * CREALF(e) + (double)CIMAGF(e) * CIMAGF(e);
		}
	}
	for (int k = 0; k < 4; k++) free(g[k]);
	rig_free(&R);
	return p;
}

/* 6. PML : 小領域 (Nx = 160, ±20um) の内部窓 |x|,|y| < 16um の残留電力を、
 *    x を 4 倍に広げた参照 (壁に届かない) と比べて反射電力を測る。
 *    solve_bpm 経由の実測は 9.5e-8 (PML) / 0.73 (従来吸収体)。判定は 1e-5 / 0.1 */
static void test_pml_reflection(void)
{
	const double win = 16e-6;
	double pin = 0, pin_ref = 0, pin_h = 0;
	const double ref  = tilt_run(640, 160, 0.0,  win, &pin_ref);
	const double pml  = tilt_run(160, 160, 4e-6, win, &pin);
	const double hard = tilt_run(160, 160, 0.0,  win, &pin_h);
	const double r_pml = fabs(pml - ref) / pin, r_hard = fabs(hard - ref) / pin_h;
	const int ok = (r_pml < 1e-5) && (r_hard > 0.1);
	printf("[%s] %-44s R_pml = %.3g (< 1e-5), R_dirichlet = %.3g (> 0.1)\n",
	       ok ? " OK " : "FAIL", "fdbpm: PML reflection vs wall-free ref", r_pml, r_hard);
	if (!ok) g_failures++;
}

int main(void)
{
	printf("=== FDBPMpropagator (scalar paraxial ADI) unit tests ===\n");
	test_free_space_diffraction();
	test_uniform_absorption();
	test_bend_deflection();
	test_pml_reflection();
	printf("=== %s (%d failure(s)) ===\n",
	       g_failures ? "FAILED" : "PASSED", g_failures);
	return g_failures ? 1 : 0;
}
