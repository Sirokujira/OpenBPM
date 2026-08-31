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

// Rayleigh 商から neff^2 に相当する量 v = n0^2 + mu/k0^2 を返す
// (クランプしない : 非導波状態では負にもなる。収束判定はこの生の値で行う。
//  クランプした値で判定すると、連続して同じ値に丸められた非導波状態を
//  「収束」と誤判定してモード探索が即座に打ち切られてしまう)
static double rayleigh_v(const wabpm_params *W, const cplx *E, const cplx *n2, cplx *work)
{
	const long N = (long)W->Nx * W->Ny;
	wabpm_apply_P(W, E, n2, work);
	const double mu = dot(E, work, N).real() / dot(E, E, N).real();
	return (W->n0 * W->n0) + (mu / (W->k0 * W->k0));
}

// 導波判定のしきい値 (= クラッド相当の neff^2) を n^2 スライスの境界から求める。
//
// 導波モードは開いた境界へ向かって減衰する必要があるので、境界セルの屈折率と
// 比較する。参照屈折率 W->n0 を使わないのは、n0 は位相基準の選択にすぎず、
// GRIN のように「n0 = 軸上最大値」を選ぶと全モードが誤って棄却されるため。
//
// **辺ごとに最大値を取り、その最小値をしきい値とする**。境界リング全体の最大値に
// すると、y 不変のスラブ導波路 (コアが y の上下端に達する) のようにコアが境界へ
// 接する構造で「しきい値 = コア屈折率」となり、正しいモードまで棄却されてしまう
// (data/slab_polarization.ofd がこれに該当する)。辺ごとの最大値の最小値なら、
// ファイバのようにクラッドが 4 辺を囲む構造では従来と同じ値になり、
// スラブでは開いた辺 (クラッド) の屈折率が選ばれる。
// 対称境界 (symx/symy) の鏡像面 (ix = 0 / iy = 0) は物理境界ではないため除外する。
double wabpm_guided_threshold(const wabpm_params *W, const cplx *n2)
{
	const long Nx = W->Nx;
	const long Ny = W->Ny;
	double side[4] = {-1.0, -1.0, -1.0, -1.0};   // x-, x+, y-, y+ の各辺の最大値
	for (long iy = 0; iy < Ny; iy++) {
		for (long ix = 0; ix < Nx; ix++) {
			const int on_xlo = (ix == 0) && !W->symx;
			const int on_xhi = (ix == Nx - 1);
			const int on_ylo = (iy == 0) && !W->symy;
			const int on_yhi = (iy == Ny - 1);
			if (!on_xlo && !on_xhi && !on_ylo && !on_yhi) continue;
			const double v = n2[ix + iy * Nx].real();
			if (on_xlo && (v > side[0])) side[0] = v;
			if (on_xhi && (v > side[1])) side[1] = v;
			if (on_ylo && (v > side[2])) side[2] = v;
			if (on_yhi && (v > side[3])) side[3] = v;
		}
	}
	double v0 = -1.0;
	for (int k = 0; k < 4; k++) {
		if (side[k] < 0) continue;                    // 対称境界で存在しない辺
		if ((v0 < 0) || (side[k] < v0)) v0 = side[k];
	}
	return (v0 < 0) ? 0.0 : v0;
}

int wabpm_find_modes(const wabpm_params *W, const cplx *n2,
                     int nModes, int maxIter, double tol,
                     cplx *modes, double *neff)
{
	const long Nx = W->Nx;
	const long Ny = W->Ny;
	const long N = Nx * Ny;

	// 反復に使う参照屈折率は「境界 (クラッド) の屈折率」に取り直す。
	//
	// 虚軸伝搬の増幅率は Cayley 写像 (1 + a*mu)/(1 - a*mu) (mu = P の固有値) で、
	// これは |a*mu| >> 1 の高周波格子モードを減衰させない (絶対値が 1 に漸近する)。
	// 一方 P には屈折率項 k0^2 (n^2 - n0^2) が入るため、n0 をコアとクラッドの
	// 中間に取ると、コア領域で正の屈折率項を持つ「x 方向 Nyquist モード」の増幅率が
	// 導波モードを上回り、そちらへ収束してしまう (実測: n1=2.0/n2=1.0 のスラブに
	// refindex = 1.5 を与えると neff^2 = -372 のスプリアス状態に収束した)。
	// n0 を境界の屈折率に取ると導波モードの mu = k0^2(neff^2 - n_clad^2) が
	// 最大固有値になり、この病理が起きない (neff 自体は n0 に依らない量)。
	//
	// 虚軸ステップ幅は呼び出し側が a = dz/(4 k0 n0) を意図して決めているので、
	// a を保つように dz を n0 の変更に合わせて取り直す。
	const double v0_thr = wabpm_guided_threshold(W, n2);
	wabpm_params Wm = *W;
	if (v0_thr > 0) {
		const double n0_new = std::sqrt(v0_thr);
		const double a_keep = W->dz / (4 * W->k0 * W->n0);
		Wm.n0 = n0_new;
		Wm.dz = a_keep * 4 * W->k0 * n0_new;
	}
	W = &Wm;

	// 虚軸ステップ幅 : 段階 1 (選択) は呼び出し側の指定、
	// 段階 2 (仕上げ) は格子スケール a = 1/(4/dx^2 + 4/dy^2)
	const wabpm_params Wsel = Wm;
	const double a_sel = Wsel.dz / (4 * Wsel.k0 * Wsel.n0);
	const double a_pol = 1.0 / ((4.0 / (Wsel.dx * Wsel.dx)) + (4.0 / (Wsel.dy * Wsel.dy)));
	wabpm_params Wpol = Wsel;
	Wpol.dz = a_pol * 4 * Wsel.k0 * Wsel.n0;
	const int usePolish = (a_pol < a_sel);

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

		// 導波条件のしきい値 (辺ごとの最大値の最小値)
		const double v0 = v0_thr;
		double v_prev = 0.0;
		int converged = 0;

		// 2 段階で反復する。
		//  段階 1 (選択)  : 呼び出し側が指定した虚軸ステップ幅。1 ステップの利得が
		//                   大きく速いが、|a*mu| >> 1 の高周波格子モードを減衰させない
		//                   (Cayley 写像の絶対値が 1 に漸近する) ため、高コントラスト
		//                   構造ではスプリアス状態に落ちることがある。
		//  段階 2 (仕上げ): 格子スケール a = 1/(4/dx^2 + 4/dy^2) で反復する。高周波成分は
		//                   確実に減衰し、ADI 分離誤差 O(a^2 Px Py) も小さくなるため、
		//                   得られる neff が格子細分化で解析解に収束するようになる
		//                   (段階 1 のみだと dx を半分にすると誤差が増えていた)。
		// 段階 1 が導波条件に届かなくても段階 2 を試す (半ベクトルの高コントラスト
		// 構造は段階 1 でスプリアスに落ちるが、段階 2 で正しいモードに復帰する)。
		for (int stage = 0; stage < 2; stage++) {
			const wabpm_params *Ws = (stage == 0) ? &Wsel : &Wpol;
			if ((stage == 1) && !usePolish) break;
			// 段階 2 は前段の結果を初期値として引き継ぐ
			converged = 0;
			int stall = 0;
			for (int it = 0; it < maxIter; it++) {
				wabpm_imagdist_step(Ws, E, n2);
				deflate(E, modes, m, N);
				if (normalize(E, N) == 0.0) { stall = 1; break; }

				const double v = rayleigh_v(Ws, E, n2, work);
				// 十分反復しても導波条件に達しないなら放射状態と判断して打ち切る
				// (段階 2 はスプリアスからの復帰に時間がかかるので余裕を持たせる)
				const int giveup = (stage == 0) ? 2000 : 4000;
				if ((it >= giveup) && (v <= v0)) {
					v_prev = v;
					stall = 1;
					break;
				}
				if ((it > 0) && (std::fabs(v - v_prev) < tol)) {
					v_prev = v;
					converged = 1;
					break;
				}
				v_prev = v;
			}
			if (stall && !usePolish) break;
			if (stall && (stage == 1)) break;
		}
		// 導波条件を満たすモードのみ採用する
		if (!converged || (v_prev <= v0)) break;

		std::memcpy(&modes[(size_t)m * N], E, (size_t)N * sizeof(cplx));
		neff[m] = std::sqrt(v_prev);
		found++;
	}

	delete[] E;
	delete[] work;
	return found;
}
