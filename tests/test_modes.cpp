/*
test_modes.cpp

モードソルバ (bpm/modes.cpp の wabpm_find_modes) の単体テスト。

ステップインデックス光ファイバの LP モードは弱導波近似で分散方程式
  u J_{l+1}(u) / J_l(u) = w K_{l+1}(w) / K_l(w),  u^2 + w^2 = V^2
の解として厳密に求まる (u = a*sqrt(k0^2 ncore^2 - beta^2),
w = a*sqrt(beta^2 - k0^2 nclad^2))。これを二分法で解いた neff と、
虚軸伝搬法によるモードソルバの結果を比較する。

検証内容:
  1. シングルモードファイバ (V=2.26) の LP01 実効屈折率
  2. 数モードファイバ (V=5) の LP01/LP11 実効屈折率と順序
  3. モードの L2 直交性
*/

#include "bpm/wabpm.h"
#include "bessel_ref.h"   // std::cyl_bessel_j/k は libc++ 未実装のため自前実装を使用

#include <complex>
#include <vector>
#include <cmath>
#include <cstdio>
#include <string>

typedef std::complex<double> cplx;

static int g_failures = 0;

static void check_close(const std::string &name, double got, double expected,
                        double atol)
{
    const double err = std::fabs(got - expected);
    const bool ok = (err <= atol);
    std::printf("[%s] %-40s got=%.8f expected=%.8f err=%.2e (tol=%.1e)\n",
                ok ? " OK " : "FAIL", name.c_str(), got, expected, err, atol);
    if (!ok) g_failures++;
}

static void check_true(const std::string &name, bool cond, const std::string &detail)
{
    std::printf("[%s] %-40s %s\n", cond ? " OK " : "FAIL", name.c_str(), detail.c_str());
    if (!cond) g_failures++;
}

// LP_l モードの分散関数 f(u) = u J_{l+1}(u)/J_l(u) - w K_{l+1}(w)/K_l(w)
static double lp_dispersion(int l, double u, double V)
{
    const double w = std::sqrt(V * V - u * u);
    return u * bessel_ref_j(l + 1, u) / bessel_ref_j(l, u)
         - w * bessel_ref_k(l + 1, w) / bessel_ref_k(l, w);
}

// LP_l1 (最低次ラジアルモード) の u を二分法で解く
// 探索区間 : (lower, upper) = (J_l の最初の零点の手前まで, V 以下)
static double solve_lp_u(int l, double V)
{
    // J_l の最初の零点 (l=0: 2.405, l=1: 3.832)
    const double jzero = (l == 0) ? 2.404825557695773 : 3.831705970207512;
    double lo = (l == 0) ? 1e-9 : 1e-9;
    double hi = std::min(V, jzero) - 1e-9;
    // f(lo) < 0 (w K/K 項が勝つ), f(hi) > 0 (J_l -> 0 で発散) を仮定して二分法
    double flo = lp_dispersion(l, lo, V);
    if (flo > 0) return -1.0;  // カットオフ以下 (モードなし)
    for (int it = 0; it < 200; it++) {
        const double mid = 0.5 * (lo + hi);
        if (lp_dispersion(l, mid, V) > 0) hi = mid; else lo = mid;
    }
    return 0.5 * (lo + hi);
}

// 厳密 LP neff (弱導波近似の分散方程式から)
static double exact_lp_neff(int l, double a, double ncore, double nclad, double lambda)
{
    const double k0 = 2.0 * M_PI / lambda;
    const double V = k0 * a * std::sqrt(ncore * ncore - nclad * nclad);
    const double u = solve_lp_u(l, V);
    if (u < 0) return -1.0;
    const double beta2 = k0 * k0 * ncore * ncore - (u / a) * (u / a);
    return std::sqrt(beta2) / k0;
}

// ステップインデックスファイバの n2 分布を作成
static std::vector<cplx> fiber_n2(int Nx, int Ny, double dx, double dy,
                                  double a, double ncore, double nclad)
{
    std::vector<cplx> n2((size_t)Nx * Ny);
    for (int iy = 0; iy < Ny; iy++) {
        for (int ix = 0; ix < Nx; ix++) {
            const double x = dx * (ix - (Nx - 1) / 2.0);
            const double y = dy * (iy - (Ny - 1) / 2.0);
            const double n = (x * x + y * y <= a * a) ? ncore : nclad;
            n2[ix + (long)iy * Nx] = cplx(n * n, 0.0);
        }
    }
    return n2;
}

// 中心化座標 (グリッド中心を原点)
static inline double coordm(long i, long N, double d)
{
    return (static_cast<double>(i) - (N - 1) / 2.0) * d;
}

static wabpm_params mode_params(int Nx, int Ny, double dx, double dy,
                                double lambda, double n0, double nmax)
{
    wabpm_params W;
    W.Nx = Nx; W.Ny = Ny;
    W.dx = dx; W.dy = dy;
    W.k0 = 2.0 * M_PI / lambda;
    W.n0 = n0;
    W.wideangle = 0; W.pol = 0;
    W.symx = 0; W.symy = 0;
    // a*mu_max ~ 0.5 となる虚軸ステップ幅
    const double mu_max = W.k0 * W.k0 * (nmax * nmax - n0 * n0);
    W.dz = 2.0 * W.k0 * W.n0 / mu_max;
    return W;
}

// シングルモードファイバ LP01 (SMF-28 相当, V = 2.26)
static void test_smf_lp01()
{
    const double lambda = 1.55e-6;
    const double a = 4.1e-6, ncore = 1.4504, nclad = 1.4440;
    const int Nx = 160, Ny = 160;
    const double dx = 40.0e-6 / Nx, dy = dx;

    const wabpm_params W = mode_params(Nx, Ny, dx, dy, lambda, nclad, ncore);
    const std::vector<cplx> n2 = fiber_n2(Nx, Ny, dx, dy, a, ncore, nclad);

    std::vector<cplx> modes((size_t)Nx * Ny);
    double neff = 0;
    const int nFound = wabpm_find_modes(&W, n2.data(), 1, 20000, 1e-12, modes.data(), &neff);

    const double neff_ana = exact_lp_neff(0, a, ncore, nclad, lambda);
    check_true("smf: found LP01", nFound == 1, "nFound = " + std::to_string(nFound));
    // 離散化誤差を考慮し neff の絶対誤差 2e-4 (コア/クラッド差 6.4e-3 の 3%)
    check_close("smf: LP01 neff", neff, neff_ana, 2e-4);
    check_true("smf: nclad < neff < ncore", neff > nclad && neff < ncore,
               "neff = " + std::to_string(neff));
}

// 数モードファイバ (V = 5) : LP01 / LP11 の neff・順序・直交性
static void test_fmf_lp01_lp11()
{
    const double lambda = 1.55e-6;
    const double a = 4.1e-6, nclad = 1.4440;
    // V = 5 となる ncore
    const double k0 = 2.0 * M_PI / lambda;
    const double NA = 5.0 / (k0 * a);
    const double ncore = std::sqrt(nclad * nclad + NA * NA);
    const int Nx = 160, Ny = 160;
    const double dx = 40.0e-6 / Nx, dy = dx;
    const long N = (long)Nx * Ny;

    const wabpm_params W = mode_params(Nx, Ny, dx, dy, lambda, nclad, ncore);
    const std::vector<cplx> n2 = fiber_n2(Nx, Ny, dx, dy, a, ncore, nclad);

    // LP01 x1 + LP11 縮退ペア x2 = 3 モード
    const int nModes = 3;
    std::vector<cplx> modes((size_t)nModes * N);
    std::vector<double> neff(nModes);
    const int nFound = wabpm_find_modes(&W, n2.data(), nModes, 20000, 1e-12,
                                        modes.data(), neff.data());

    const double neff01 = exact_lp_neff(0, a, ncore, nclad, lambda);
    const double neff11 = exact_lp_neff(1, a, ncore, nclad, lambda);

    check_true("fmf: found 3 modes", nFound == 3, "nFound = " + std::to_string(nFound));
    check_close("fmf: LP01 neff", neff[0], neff01, 3e-4);
    check_close("fmf: LP11a neff", neff[1], neff11, 3e-4);
    check_close("fmf: LP11b neff", neff[2], neff11, 3e-4);
    check_true("fmf: neff order", neff[0] > neff[1] && neff[1] >= neff[2] - 1e-6,
               "neff = " + std::to_string(neff[0]) + ", " + std::to_string(neff[1]) +
               ", " + std::to_string(neff[2]));

    // L2 直交性
    double maxovl = 0.0;
    for (int p = 0; p < nFound; p++) {
        for (int q = p + 1; q < nFound; q++) {
            cplx s = 0.0;
            for (long i = 0; i < N; i++) {
                s += std::conj(modes[(size_t)p * N + i]) * modes[(size_t)q * N + i];
            }
            maxovl = std::max(maxovl, std::abs(s));
        }
    }
    check_true("fmf: orthogonality", maxovl < 1e-6,
               "max |<m_p, m_q>| = " + std::to_string(maxovl));
}

// 高コントラスト矩形コアで複数モードが確実に見つかること (回帰テスト)
//
// 収束判定を Rayleigh 商から算出した neff で行い、非導波状態を 0 にクランプすると、
// 2 回連続で 0 に丸められた時点で「収束」と誤判定し、モード 2 以降の探索が
// 即座に打ち切られる。この退行を検出するため、モード数と neff の降順性を確認する。
// (判定はクランプしない生の v = n0^2 + mu/k0^2 で行うこと)
static void test_highcontrast_multimode()
{
    const double lambda = 1.55e-6;
    const double ncore = 2.0, nclad = 1.0;
    const int Nx = 120, Ny = 120;
    const double dx = 6.0e-6 / Nx, dy = dx;
    const double ax = 1.2e-6, ay = 0.9e-6;   // 矩形コアの半幅
    const long N = (long)Nx * Ny;

    wabpm_params W = mode_params(Nx, Ny, dx, dy, lambda, nclad, ncore);
    std::vector<cplx> n2((size_t)N);
    for (int iy = 0; iy < Ny; iy++) {
        for (int ix = 0; ix < Nx; ix++) {
            const double x = coordm(ix, Nx, dx);
            const double y = coordm(iy, Ny, dy);
            const double n = (std::fabs(x) <= ax && std::fabs(y) <= ay) ? ncore : nclad;
            n2[ix + (long)iy * Nx] = cplx(n * n, 0.0);
        }
    }

    const int nModes = 4;
    std::vector<cplx> modes((size_t)nModes * N);
    std::vector<double> neff(nModes);
    const int nFound = wabpm_find_modes(&W, n2.data(), nModes, 60000, 1e-12,
                                        modes.data(), neff.data());

    check_true("multimode: found 4 modes", nFound == 4,
               "nFound = " + std::to_string(nFound));
    if (nFound < 2) return;
    // 全モードが導波条件を満たし neff は降順
    bool ok = true;
    for (int m = 0; m < nFound; m++) {
        if (!(neff[m] > nclad && neff[m] < ncore)) ok = false;
        if (m > 0 && neff[m] > neff[m - 1] + 1e-9) ok = false;
    }
    check_true("multimode: neff guided and descending", ok,
               "neff[0] = " + std::to_string(neff[0]) +
               ", neff[" + std::to_string(nFound - 1) + "] = " + std::to_string(neff[nFound - 1]));
    // 直交性
    double maxovl = 0.0;
    for (int p = 0; p < nFound; p++) {
        for (int q = p + 1; q < nFound; q++) {
            cplx s = 0.0;
            for (long i = 0; i < N; i++) {
                s += std::conj(modes[(size_t)p * N + i]) * modes[(size_t)q * N + i];
            }
            maxovl = std::max(maxovl, std::abs(s));
        }
    }
    check_true("multimode: orthogonality", maxovl < 1e-6,
               "max |<m_p, m_q>| = " + std::to_string(maxovl));
}

// スラブ導波路の TE / TM 実効屈折率 (半ベクトル差分 = Stern の解析解検証)
//
// 対称スラブ (コア厚 2a, 芯 n1 / クラッド n2) の基本モードは
//   u = kappa*a,  V = k0*a*sqrt(n1^2-n2^2),  (gamma*a)^2 = V^2 - u^2
//   TE : u tan(u) = gamma*a
//   TM : u tan(u) = (n1/n2)^2 * gamma*a       <- 半ベクトル (Stern) 差分が効く
// で決まる。TE と TM で neff が大きく分かれる (この設定で 1.674 vs 1.358) ため、
// Stern の重みが正しいかどうかを鋭敏に検出できる。
//
// 2D ソルバーで y 不変のスラブを解くと y 方向の閉じ込め (Dirichlet) の分だけ
// neff が下がる。y は分離可能なので、離散ラプラシアンの固有値
//   lam_y = -(4/dy^2) sin^2(pi/(2(Ny+1)))
// を厳密に差し引けば 1D スラブの neff を復元できる。
//
// 許容誤差は離散化誤差から決める : dx = 25nm での実測は TE +9.8e-4 / TM +1.6e-3 で、
// dx を半分にすると 1/4 (O(dx^2)) になることを確認済み (dx=6.25nm で 6e-5 / 1e-4)。
static double slab_neff_analytic(double k0, double n1, double n2, double a, int tm)
{
    const double V = k0 * a * std::sqrt(n1 * n1 - n2 * n2);
    auto f = [&](double u) {
        const double w = std::sqrt(std::max(0.0, V * V - u * u));
        const double c = tm ? ((n1 * n1) / (n2 * n2)) : 1.0;
        return u * std::tan(u) - c * w;
    };
    double lo = 1e-12, hi = std::min(0.5 * M_PI - 1e-12, V - 1e-12);
    for (int it = 0; it < 200; it++) {
        const double mid = 0.5 * (lo + hi);
        if (f(mid) > 0) hi = mid; else lo = mid;
    }
    const double kap = 0.5 * (lo + hi) / a;
    return std::sqrt(n1 * n1 - (kap / k0) * (kap / k0));
}

static void test_slab_te_tm()
{
    const double lambda = 1.55e-6;
    const double k0 = 2 * M_PI / lambda;
    const double n1 = 2.0, n2 = 1.0;      // 高コントラスト (Si 系相当)
    const double a = 0.2e-6;              // 半幅 (コア厚 0.4um)
    const int Nx = 320, Ny = 40;
    const double dx = 8.0e-6 / Nx, dy = 10.0e-6 / Ny;
    const long N = (long)Nx * Ny;

    std::vector<cplx> n2m((size_t)N);
    for (int iy = 0; iy < Ny; iy++) {
        for (int ix = 0; ix < Nx; ix++) {
            const double x = coordm(ix, Nx, dx);
            const double n = (std::fabs(x) < a) ? n1 : n2;
            n2m[ix + (long)iy * Nx] = cplx(n * n, 0.0);
        }
    }
    // y 方向 (Dirichlet) の離散固有値 : 分離可能なので厳密に差し引ける
    const double sy = std::sin(M_PI / (2 * (Ny + 1)));
    const double lam_y = -(4.0 / (dy * dy)) * sy * sy;

    double neff_te = 0, neff_tm = 0;
    for (int pol = 1; pol <= 2; pol++) {          // 1 = x 偏波 (TM), 2 = y 偏波 (TE)
        // 参照屈折率は中間値 (1.5) にする : 以前はこの指定でスプリアス格子モードへ
        // 収束していた (find_modes が内部で境界屈折率に取り直すことの回帰検出)
        wabpm_params W = mode_params(Nx, Ny, dx, dy, lambda, 1.5, n1);
        W.pol = pol;
        std::vector<cplx> mode((size_t)N);
        double v = 0;
        const int nf = wabpm_find_modes(&W, n2m.data(), 1, 40000, 1e-12, mode.data(), &v);
        const double neff = (nf == 1) ? std::sqrt(v * v - lam_y / (k0 * k0)) : 0.0;
        if (pol == 1) neff_tm = neff; else neff_te = neff;
        check_true(std::string("slab: ") + ((pol == 1) ? "TM" : "TE") + " mode found",
                   nf == 1, "nFound = " + std::to_string(nf));
    }

    const double te_ana = slab_neff_analytic(k0, n1, n2, a, 0);
    const double tm_ana = slab_neff_analytic(k0, n1, n2, a, 1);
    check_close("slab: TE neff (dispersion eq.)", neff_te, te_ana, 3e-3);
    check_close("slab: TM neff (dispersion eq., Stern)", neff_tm, tm_ana, 3e-3);
    // 偏波分離そのもの (TE - TM = 0.3164) が再現できていること
    check_close("slab: TE-TM splitting", neff_te - neff_tm, te_ana - tm_ana, 2e-2);
}

// 方向性結合器 (同一 2 コア) の超モード (偶 / 奇)
//
// 弱く導波する高次モード (奇超モード) は mu = k0^2 (neff^2 - nclad^2) が小さいため、
// 虚軸伝搬の選択段階を Cayley 写像 (1+a mu)/(1-a mu) で行うと高周波格子モードに
// 負けてスプリアス状態に落ち、「1 / 2 converged」となっていた (data/fiber_coupler.ofd)。
// 選択段階を後退オイラー 1/(1 - a mu) にした修正の回帰検出を兼ねる。
//
// 検証 (いずれも解析的な基準):
//   1. 偶 / 奇の 2 本が見つかり、x について偶関数 / 奇関数であること (厳密な対称性)
//   2. 超モードの平均 neff が孤立コア LP01 の厳密解と一致する (分裂は 1 次で対称)。
//      許容 1e-4 : 孤立コアの離散化誤差 (smf テストで ~1e-5) + 分裂の 2 次項
//   3. 結合長 Lc = lambda / (2 (n_even - n_odd)) が結合モード理論 (Snyder & Love)
//        kappa = sqrt(2 Delta)/a * U^2/V^3 * K0(W D/a)/K1(W)^2,  Lc = pi/(2 kappa)
//      と一致する。CMT は弱結合近似 (この構造で約 3% の差) なので許容 10%。
//      BPM 伝搬から直接フィットした Lc (833.0um) とは 0.07% で一致することを確認済み
static void test_coupler_supermodes()
{
    const double lambda = 1.55e-6;
    const double k0 = 2 * M_PI / lambda;
    const double ncore = 1.4504, nclad = 1.4440;
    const double a = 3.0e-6, D = 8.0e-6;      // コア半径 / コア中心間隔
    const int Nx = 120, Ny = 120;
    const double dx = 40.0e-6 / Nx, dy = dx;
    const long N = (long)Nx * Ny;

    std::vector<cplx> n2((size_t)N);
    for (int iy = 0; iy < Ny; iy++) {
        for (int ix = 0; ix < Nx; ix++) {
            const double x = coordm(ix, Nx, dx);
            const double y = coordm(iy, Ny, dy);
            const bool in = (std::hypot(x + 0.5 * D, y) < a) || (std::hypot(x - 0.5 * D, y) < a);
            const double n = in ? ncore : nclad;
            n2[ix + (long)iy * Nx] = cplx(n * n, 0.0);
        }
    }
    const wabpm_params W = mode_params(Nx, Ny, dx, dy, lambda, nclad, ncore);
    std::vector<cplx> modes((size_t)2 * N);
    double neff[2] = {0, 0};
    const int nFound = wabpm_find_modes(&W, n2.data(), 2, 20000, 1e-12, modes.data(), neff);
    check_true("coupler: found even + odd supermodes", nFound == 2,
               "nFound = " + std::to_string(nFound));
    if (nFound < 2) return;

    // x 方向の対称性 : E(-x) = +E(x) (偶) / -E(x) (奇)
    double sym[2] = {0, 0};
    for (int m = 0; m < 2; m++) {
        const cplx *e = &modes[(size_t)m * N];
        double sp = 0, sm = 0;
        for (int iy = 0; iy < Ny; iy++) {
            for (int ix = 0; ix < Nx; ix++) {
                const cplx v = e[ix + (long)iy * Nx];
                const cplx w = e[(Nx - 1 - ix) + (long)iy * Nx];
                sp += std::norm(v + w);
                sm += std::norm(v - w);
            }
        }
        sym[m] = (m == 0) ? (sm / (sp + sm)) : (sp / (sp + sm));   // 期待と逆の成分の割合
    }
    check_true("coupler: mode 1 even / mode 2 odd in x", (sym[0] < 1e-6) && (sym[1] < 1e-6),
               "残留 = " + std::to_string(sym[0]) + ", " + std::to_string(sym[1]));

    const double n_iso = exact_lp_neff(0, a, ncore, nclad, lambda);
    check_close("coupler: mean supermode neff = isolated LP01", 0.5 * (neff[0] + neff[1]), n_iso, 1e-4);

    const double V = k0 * a * std::sqrt(ncore * ncore - nclad * nclad);
    const double U = solve_lp_u(0, V);
    const double Wv = std::sqrt(V * V - U * U);
    const double Delta = (ncore * ncore - nclad * nclad) / (2 * ncore * ncore);
    const double kb = bessel_ref_k(1, Wv);
    const double kappa = std::sqrt(2 * Delta) / a * (U * U) / (V * V * V)
                       * bessel_ref_k(0, Wv * D / a) / (kb * kb);
    const double Lc_cmt = M_PI / (2 * kappa);
    const double Lc = lambda / (2 * (neff[0] - neff[1]));
    // check_close は絶対誤差なので比 Lc/Lc_cmt で比較する (相対 10%)
    check_close("coupler: Lc / Lc(coupled-mode theory)", Lc / Lc_cmt, 1.0, 0.10);
}

int main()
{
    std::printf("=== wabpm_find_modes unit tests ===\n");
    test_smf_lp01();
    test_fmf_lp01_lp11();
    test_highcontrast_multimode();
    test_slab_te_tm();
    test_coupler_supermodes();

    std::printf("\n%d failure(s).\n", g_failures);
    return g_failures == 0 ? 0 : 1;
}
