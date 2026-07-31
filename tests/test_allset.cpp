/*
test_allset.cpp

include/bpm/allset.hpp (BPM-MATLAB 互換の P_Struct API) の単体テスト。

検証内容:
  1. findModes : ステップインデックスファイバの LP01 実効屈折率が
     分散方程式の厳密解と一致し、モード形状が中心対称であること
  2. modeSuperposition / getLabeledModeIndex : モードの重ね合わせと参照
  3. offsetField : 界の重心が指定量だけ平行移動すること
  4. tiltField : 位相ランプが付与され振幅分布は不変であること
*/

#include "bpm/allset.hpp"
#include "bessel_ref.h"   // std::cyl_bessel_j/k は libc++ 未実装のため自前実装を使用

#include <cmath>
#include <cstdio>

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

// LP01 の厳密 neff (test_modes.cpp と同じ分散方程式, l = 0)
static double exact_lp01_neff(double a, double ncore, double nclad, double lambda)
{
    const double k0 = 2.0 * M_PI / lambda;
    const double V = k0 * a * std::sqrt(ncore * ncore - nclad * nclad);
    double lo = 1e-9, hi = std::min(V, 2.404825557695773) - 1e-9;
    for (int it = 0; it < 200; it++) {
        const double u = 0.5 * (lo + hi);
        const double w = std::sqrt(V * V - u * u);
        const double f = u * bessel_ref_j(1, u) / bessel_ref_j(0, u)
                       - w * bessel_ref_k(1, w) / bessel_ref_k(0, w);
        if (f > 0) hi = u; else lo = u;
    }
    const double u = 0.5 * (lo + hi);
    return std::sqrt(ncore * ncore - (u / (k0 * a)) * (u / (k0 * a)));
}

// SMF-28 相当のステップインデックスファイバを P_Struct に構築
static P_Struct make_fiber(int Nx, int Ny, double L, double a,
                           double ncore, double nclad, double lambda)
{
    P_Struct P;
    const double dx = L / Nx, dy = L / Ny;
    P.dx = (float)dx;
    P.dy = (float)dy;
    P.lambda = (float)lambda;
    P.n_0 = (float)nclad;
    P.n_background = (float)nclad;
    P.xSymmetry = "none";
    P.ySymmetry = "none";
    for (int i = 0; i < Nx; i++) P.x.push_back((float)(dx * (i - (Nx - 1) / 2.0)));
    for (int j = 0; j < Ny; j++) P.y.push_back((float)(dy * (j - (Ny - 1) / 2.0)));
    P.n.Lx = (float)L;
    P.n.Ly = (float)L;
    P.n.n = MatrixXf(Nx, Ny);
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            const double r2 = (double)P.x[i] * P.x[i] + (double)P.y[j] * P.y[j];
            P.n.n(i, j) = (r2 <= a * a) ? (float)ncore : (float)nclad;
        }
    }
    return P;
}

// 界の重心 (x 方向)
static double centroid_x(const P_Struct &P)
{
    double sumI = 0, sumX = 0;
    for (int i = 0; i < (int)P.x.size(); i++) {
        for (int j = 0; j < (int)P.y.size(); j++) {
            const double I = std::norm(P.E.field(i, j));
            sumI += I;
            sumX += P.x[i] * I;
        }
    }
    return sumX / sumI;
}

int main()
{
    std::printf("=== allset.hpp (P_Struct API) unit tests ===\n");

    const double lambda = 1.55e-6;
    const double a = 4.1e-6, ncore = 1.4504, nclad = 1.4440;
    const int Nx = 160, Ny = 160;
    P_Struct P = make_fiber(Nx, Ny, 40.0e-6, a, ncore, nclad, lambda);

    // --- findModes ---
    P = findModes(P, 1);
    check_true("findModes: 1 mode found", P.modes.size() == 1,
               "modes.size() = " + to_string(P.modes.size()));
    const double neff_ana = exact_lp01_neff(a, ncore, nclad, lambda);
    check_close("findModes: LP01 neff", P.modes[0].neff.real(), neff_ana, 2e-4);

    // モード形状 : 中心にピーク・境界で減衰
    const int ic = Nx / 2, jc = Ny / 2;
    const float peak = std::abs(P.modes[0].field(ic, jc));
    const float edge = std::abs(P.modes[0].field(0, jc));
    check_true("findModes: mode confined", peak > 100 * edge,
               "peak/edge = " + to_string(edge > 0 ? peak / edge : 1e9f));

    // --- modeSuperposition / getLabeledModeIndex ---
    const int idx = getLabeledModeIndex(P, "Mode 1");
    check_true("getLabeledModeIndex", idx == 0, "idx = " + to_string(idx));
    P = modeSuperposition(P, {0}, {2.0f});
    const float ratio = std::abs(P.E.field(ic, jc)) / peak;
    check_close("modeSuperposition: 2x amplitude", ratio, 2.0, 1e-4);

    // --- offsetField : +x へ 3um 移動 ---
    const double x0 = centroid_x(P);
    P = offsetField(P, 0.0f, 3.0e-6f);
    const double x1 = centroid_x(P);
    check_close("offsetField: centroid shift [um]", (x1 - x0) * 1e6, 3.0, 0.05);

    // --- tiltField : 振幅分布は不変・位相ランプ付与 ---
    const double amp_before = std::abs(P.E.field(ic, jc));
    P = tiltField(P, 0.0f, 5.0f);   // +x へ 5 度
    const double amp_after = std::abs(P.E.field(ic, jc));
    check_close("tiltField: amplitude preserved", amp_after / amp_before, 1.0, 1e-5);
    // 隣接格子点間の位相差 = k0 n0 sin?(小角では angle) * dx
    const double k_t = P.n_0 * 2 * M_PI / P.lambda * (5.0 * M_PI / 180);
    const double dphi = std::arg(P.E.field(ic + 1, jc) / P.E.field(ic, jc));
    check_close("tiltField: phase ramp", -dphi, k_t * P.dx, 1e-3 * k_t * P.dx);

    std::printf("\n%d failure(s).\n", g_failures);
    return g_failures == 0 ? 0 : 1;
}
