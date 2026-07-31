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

static wabpm_params mode_params(int Nx, int Ny, double dx, double dy,
                                double lambda, double n0, double nmax)
{
    wabpm_params W;
    W.Nx = Nx; W.Ny = Ny;
    W.dx = dx; W.dy = dy;
    W.k0 = 2.0 * M_PI / lambda;
    W.n0 = n0;
    W.wideangle = 0; W.pol = 0;
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

int main()
{
    std::printf("=== wabpm_find_modes unit tests ===\n");
    test_smf_lp01();
    test_fmf_lp01_lp11();

    std::printf("\n%d failure(s).\n", g_failures);
    return g_failures == 0 ? 0 : 1;
}
