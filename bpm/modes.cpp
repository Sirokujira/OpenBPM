/*
modes.cpp

導波モードソルバ (虚軸伝搬法 + Gram-Schmidt 直交化)

虚軸伝搬 1 ステップ (wabpm_imagdist_step) を伝搬演算子とするべき乗法で、
実効屈折率の大きい順に導波モードを求める。m 番目のモード探索では、
既に得られたモード 0..m-1 への射影を毎反復で除去 (deflation) する。

実効屈折率は Rayleigh 商
  mu = <E, P E> / <E, E>,  neff = sqrt(n0^2 + mu/k0^2)
から算出し、反復間の neff 変化が tol を下回ったら収束とみなす。

スカラー Helmholtz 演算子は実対称のため固有モードは L2 直交であり、
標準内積の Gram-Schmidt で正しく deflation できる (スカラー使用を想定。
半ベクトル Stern 差分は非対称のため直交性は近似となる)。
*/

#include <complex>
#include <cmath>
#include <cstring>

#include "bpm/wabpm.h"

typedef std::complex<double> cplx;

// 内積 <a, b> = sum conj(a)*b
static cplx dot(const cplx *a, const cplx *b, long N)
{
	cplx s = 0.0;
	for (long i = 0; i < N; i++) s += std::conj(a[i]) * b[i];
	return s;
}

static double normalize(cplx *E, long N)
{
	const double nrm = std::sqrt(dot(E, E, N).real());
	if (nrm > 0) {
		for (long i = 0; i < N; i++) E[i] /= nrm;
	}
	return nrm;
}

// 既出モードへの射影を除去 (Gram-Schmidt deflation)
static void deflate(cplx *E, const cplx *modes, int m, long N)
{
	for (int j = 0; j < m; j++) {
		const cplx *phi = &modes[(size_t)j * N];
		const cplx c = dot(phi, E, N);
		for (long i = 0; i < N; i++) E[i] -= c * phi[i];
	}
}

// Rayleigh 商から実効屈折率を算出
static double rayleigh_neff(const wabpm_params *W, const cplx *E, const cplx *n2, cplx *work)
{
	const long N = (long)W->Nx * W->Ny;
	wabpm_apply_P(W, E, n2, work);
	const double mu = dot(E, work, N).real() / dot(E, E, N).real();
	const double v = W->n0 * W->n0 + mu / (W->k0 * W->k0);
	return (v > 0) ? std::sqrt(v) : 0.0;
}

int wabpm_find_modes(const wabpm_params *W, const cplx *n2,
                     int nModes, int maxIter, double tol,
                     cplx *modes, double *neff)
{
	const long Nx = W->Nx;
	const long Ny = W->Ny;
	const long N = Nx * Ny;

	cplx *E = new cplx[N];
	cplx *work = new cplx[N];

	// 再現性のある擬似乱数 (LCG) : 全モード成分を含むシードを作る
	unsigned long long seed = 88172645463325252ULL;
	auto frand = [&seed]() {
		seed = seed * 6364136223846793005ULL + 1442695040888963407ULL;
		return (double)(seed >> 33) / (double)(1ULL << 31) - 1.0;
	};

	int found = 0;
	for (int m = 0; m < nModes; m++) {
		// 初期値 : ガウシアン包絡 × 乱数 (全対称性のモードを含む)
		const double wx = 0.25 * Nx * W->dx;
		const double wy = 0.25 * Ny * W->dy;
		for (long iy = 0; iy < Ny; iy++) {
			for (long ix = 0; ix < Nx; ix++) {
				const double x = W->dx * (ix - (Nx - 1) / 2.0);
				const double y = W->dy * (iy - (Ny - 1) / 2.0);
				E[ix + iy * Nx] = frand() * std::exp(-(x * x) / (wx * wx) - (y * y) / (wy * wy));
			}
		}
		deflate(E, modes, m, N);
		normalize(E, N);

		double neff_prev = 0.0;
		int converged = 0;
		for (int it = 0; it < maxIter; it++) {
			wabpm_imagdist_step(W, E, n2);
			deflate(E, modes, m, N);
			if (normalize(E, N) == 0.0) break;   // 部分空間が尽きた

			const double ne = rayleigh_neff(W, E, n2, work);
			if (it > 0 && std::fabs(ne - neff_prev) < tol) {
				neff_prev = ne;
				converged = 1;
				break;
			}
			neff_prev = ne;
		}
		if (!converged) break;

		std::memcpy(&modes[(size_t)m * N], E, (size_t)N * sizeof(cplx));
		neff[m] = neff_prev;
		found++;
	}

	delete[] E;
	delete[] work;
	return found;
}
