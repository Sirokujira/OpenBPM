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

int main(void)
{
	printf("=== FDBPMpropagator (scalar paraxial ADI) unit tests ===\n");
	test_free_space_diffraction();
	test_uniform_absorption();
	printf("=== %s (%d failure(s)) ===\n",
	       g_failures ? "FAILED" : "PASSED", g_failures);
	return g_failures ? 1 : 0;
}
