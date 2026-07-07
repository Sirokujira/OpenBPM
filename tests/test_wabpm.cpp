/*
test_wabpm.cpp

一般化 BPM 伝搬カーネル (bpm/wabpm.cpp の wabpm_step) の単体テスト。
外部テストフレームワークに依存しない自己完結ハーネス。

検証内容 (いずれも解析解との比較):
  1. 自由空間ガウシアンビームの回折 (近軸)  : w(z) = w0 sqrt(1 + (z/zR)^2)
  2. 損失なし自由空間でのエネルギー保存 (近軸): sum|E|^2 = const
  3. 一様吸収媒質での強度減衰              : P(z)/P(0) = exp(-2 k0 n'' z)
  4. 広角 (Pade(1,1)) パスでの回折         : 近軸領域では w(z) に一致
  5. 半ベクトル差分の一様媒質での整合       : pol=x/y はスカラーと一致 (Stern が標準ラプラシアンに帰着)

失敗が 1 件でもあれば終了コード 1 を返す (CTest 連携用)。
*/

#include "bpm/wabpm.h"

#include <complex>
#include <vector>
#include <cmath>
#include <cstdio>
#include <string>

typedef std::complex<double> cplx;

static int g_failures = 0;

// 相対誤差での合否判定
static void check_close(const std::string &name, double got, double expected,
                        double rtol)
{
    const double denom = (std::fabs(expected) > 1e-300) ? std::fabs(expected) : 1.0;
    const double relerr = std::fabs(got - expected) / denom;
    const bool ok = (relerr <= rtol);
    std::printf("[%s] %-44s got=%.6g expected=%.6g relerr=%.3g (tol=%.3g)\n",
                ok ? " OK " : "FAIL", name.c_str(), got, expected, relerr, rtol);
    if (!ok) g_failures++;
}

static void check_true(const std::string &name, bool cond, const std::string &detail)
{
    std::printf("[%s] %-44s %s\n", cond ? " OK " : "FAIL", name.c_str(), detail.c_str());
    if (!cond) g_failures++;
}

// 中心化座標 (グリッド中心を原点)
static inline double coord(long i, long N, double d)
{
    return (static_cast<double>(i) - (N - 1) / 2.0) * d;
}

// 強度の 2 次モーメントからビーム幅 w を推定 (I ∝ exp(-2 r^2/w^2) -> <x^2> = w^2/4)
static double beam_width(const std::vector<cplx> &E, int Nx, int Ny, double dx)
{
    double sumI = 0.0, sumX2 = 0.0;
    for (int iy = 0; iy < Ny; iy++) {
        for (int ix = 0; ix < Nx; ix++) {
            const double I = std::norm(E[ix + iy * Nx]);
            const double x = coord(ix, Nx, dx);
            sumI += I;
            sumX2 += x * x * I;
        }
    }
    return 2.0 * std::sqrt(sumX2 / sumI);
}

static double total_power(const std::vector<cplx> &E)
{
    double p = 0.0;
    for (const cplx &e : E) p += std::norm(e);
    return p;
}

// 一様な複素比誘電率 n2 = (nr - i*ni)^2 を生成
static std::vector<cplx> uniform_n2(int Nx, int Ny, double nr, double ni)
{
    const cplx n(nr, -ni);            // 損失は負の虚部 (wabpm.h の規約)
    return std::vector<cplx>((size_t)Nx * Ny, n * n);
}

// ガウシアンビーム初期値 (field 1/e 半径 = w0)
static std::vector<cplx> gaussian(int Nx, int Ny, double dx, double dy, double w0)
{
    std::vector<cplx> E((size_t)Nx * Ny);
    for (int iy = 0; iy < Ny; iy++) {
        for (int ix = 0; ix < Nx; ix++) {
            const double x = coord(ix, Nx, dx);
            const double y = coord(iy, Ny, dy);
            E[ix + iy * Nx] = std::exp(-(x * x + y * y) / (w0 * w0));
        }
    }
    return E;
}

static wabpm_params make_params(int Nx, int Ny, double dx, double dy, double dz,
                                double k0, double n0, int wideangle, int pol)
{
    wabpm_params W;
    W.Nx = Nx; W.Ny = Ny;
    W.dx = dx; W.dy = dy; W.dz = dz;
    W.k0 = k0; W.n0 = n0;
    W.wideangle = wideangle; W.pol = pol;
    return W;
}

// 自由空間ガウシアン回折テスト (wideangle = 0/1 で共用)
static void test_diffraction(int wideangle)
{
    const std::string tag = wideangle ? "diffraction(wideangle)" : "diffraction(paraxial)";
    const int Nx = 128, Ny = 128;
    const double lambda = 1.0e-6;
    const double n0 = 1.0;
    const double k0 = 2.0 * M_PI / lambda;
    const double w0 = 10.0e-6;
    const double dx = 200.0e-6 / Nx;     // 領域 200um (ビーム幅に対し十分広い)
    const double dy = dx;
    const double dz = 1.0e-6;

    const double zR = 0.5 * k0 * n0 * w0 * w0;   // レイリー長
    const int nstep = static_cast<int>(zR / dz); // z = zR まで伝搬
    const double z = nstep * dz;

    wabpm_params W = make_params(Nx, Ny, dx, dy, dz, k0, n0, wideangle, 0);
    std::vector<cplx> E = gaussian(Nx, Ny, dx, dy, w0);
    const std::vector<cplx> n2 = uniform_n2(Nx, Ny, n0, 0.0);

    for (int s = 0; s < nstep; s++) wabpm_step(&W, E.data(), n2.data());

    const double w_num = beam_width(E, Nx, Ny, dx);
    const double w_ana = w0 * std::sqrt(1.0 + (z / zR) * (z / zR));
    check_close(tag, w_num, w_ana, 0.05);
}

// 損失なし自由空間でのエネルギー保存 (近軸 ADI のユニタリ性)
static void test_energy_conservation()
{
    const int Nx = 128, Ny = 128;
    const double lambda = 1.0e-6;
    const double n0 = 1.45;
    const double k0 = 2.0 * M_PI / lambda;
    const double w0 = 12.0e-6;
    const double dx = 240.0e-6 / Nx;
    const double dy = dx;
    const double dz = 1.0e-6;
    const int nstep = 100;

    wabpm_params W = make_params(Nx, Ny, dx, dy, dz, k0, n0, 0, 0);
    std::vector<cplx> E = gaussian(Nx, Ny, dx, dy, w0);
    const std::vector<cplx> n2 = uniform_n2(Nx, Ny, n0, 0.0);

    const double p0 = total_power(E);
    for (int s = 0; s < nstep; s++) wabpm_step(&W, E.data(), n2.data());
    const double p1 = total_power(E);

    check_close("energy_conservation", p1 / p0, 1.0, 0.01);
}

// 一様吸収媒質での強度減衰 (回折はビーム全体の電力を保存するため、減衰のみが効く)
static void test_absorption()
{
    const int Nx = 96, Ny = 96;
    const double lambda = 1.55e-6;
    const double n0 = 1.45;
    const double k0 = 2.0 * M_PI / lambda;
    const double w0 = 12.0e-6;
    const double dx = 240.0e-6 / Nx;
    const double dy = dx;
    const double dz = 1.0e-6;
    const double kappa = 2.0e-4;     // n'' (吸収)
    const int nstep = 200;
    const double z = nstep * dz;

    wabpm_params W = make_params(Nx, Ny, dx, dy, dz, k0, n0, 0, 0);
    std::vector<cplx> E = gaussian(Nx, Ny, dx, dy, w0);
    const std::vector<cplx> n2 = uniform_n2(Nx, Ny, n0, kappa);

    const double p0 = total_power(E);
    for (int s = 0; s < nstep; s++) wabpm_step(&W, E.data(), n2.data());
    const double p1 = total_power(E);

    const double ratio_ana = std::exp(-2.0 * k0 * kappa * z);
    check_close("absorption_decay", p1 / p0, ratio_ana, 0.02);
}

// 曲げ (等価屈折率法) : 一様媒質 + n(x) = n0*exp(x/RoC) では
// ビーム重心が円弧 (レイ) を描く : <x>(z) = z^2/(2*RoC) (Ehrenfest の定理)
static void test_bend_deflection()
{
    const int Nx = 192, Ny = 64;
    const double lambda = 1.55e-6;
    const double n0 = 1.444;
    const double k0 = 2.0 * M_PI / lambda;
    const double w0 = 8.0e-6;
    const double dx = 288.0e-6 / Nx;    // x は偏向方向なので広め
    const double dy = 96.0e-6 / Ny;
    const double dz = 1.0e-6;
    const double RoC = 20.0e-3;
    const int nstep = 400;
    const double z = nstep * dz;

    wabpm_params W = make_params(Nx, Ny, dx, dy, dz, k0, n0, 0, 0);
    std::vector<cplx> E((size_t)Nx * Ny);
    for (int iy = 0; iy < Ny; iy++) {
        for (int ix = 0; ix < Nx; ix++) {
            const double x = coord(ix, Nx, dx);
            const double y = coord(iy, Ny, dy);
            E[ix + (long)iy * Nx] = std::exp(-(x * x + y * y) / (w0 * w0));
        }
    }
    // 等価屈折率 (solve_bpm の曲げ変換と同じ, rho_e = 0)
    std::vector<cplx> n2((size_t)Nx * Ny);
    for (int iy = 0; iy < Ny; iy++) {
        for (int ix = 0; ix < Nx; ix++) {
            const double x = coord(ix, Nx, dx);
            const double n = n0 * std::exp(x / RoC);
            n2[ix + (long)iy * Nx] = cplx(n * n, 0.0);
        }
    }

    for (int s = 0; s < nstep; s++) wabpm_step(&W, E.data(), n2.data());

    // x 重心
    double sumI = 0, sumX = 0;
    for (int iy = 0; iy < Ny; iy++) {
        for (int ix = 0; ix < Nx; ix++) {
            const double I = std::norm(E[ix + (long)iy * Nx]);
            sumI += I;
            sumX += coord(ix, Nx, dx) * I;
        }
    }
    const double xc_num = sumX / sumI;
    const double xc_ana = z * z / (2.0 * RoC);
    check_close("bend_deflection", xc_num, xc_ana, 0.02);
}

// 半ベクトル差分 (Stern) は一様媒質では標準ラプラシアンに帰着する
static void test_semivectorial_uniform()
{
    const int Nx = 64, Ny = 64;
    const double lambda = 1.55e-6;
    const double n0 = 1.45;
    const double k0 = 2.0 * M_PI / lambda;
    const double w0 = 10.0e-6;
    const double dx = 160.0e-6 / Nx;
    const double dy = dx;
    const double dz = 1.0e-6;
    const int nstep = 50;

    const std::vector<cplx> n2 = uniform_n2(Nx, Ny, n0, 0.0);
    const std::vector<cplx> E0 = gaussian(Nx, Ny, dx, dy, w0);

    auto run = [&](int pol) {
        wabpm_params W = make_params(Nx, Ny, dx, dy, dz, k0, n0, 0, pol);
        std::vector<cplx> E = E0;
        for (int s = 0; s < nstep; s++) wabpm_step(&W, E.data(), n2.data());
        return E;
    };

    const std::vector<cplx> Es = run(0);   // scalar
    const std::vector<cplx> Ex = run(1);   // x 偏波 (半ベクトル)
    const std::vector<cplx> Ey = run(2);   // y 偏波 (半ベクトル)

    double maxdiff_x = 0.0, maxdiff_y = 0.0;
    for (size_t i = 0; i < Es.size(); i++) {
        maxdiff_x = std::max(maxdiff_x, std::abs(Ex[i] - Es[i]));
        maxdiff_y = std::max(maxdiff_y, std::abs(Ey[i] - Es[i]));
    }
    check_true("semivectorial==scalar (uniform, x)", maxdiff_x < 1e-9,
               "max|Ex-Es| = " + std::to_string(maxdiff_x));
    check_true("semivectorial==scalar (uniform, y)", maxdiff_y < 1e-9,
               "max|Ey-Es| = " + std::to_string(maxdiff_y));
}

int main()
{
    std::printf("=== wabpm_step unit tests ===\n");
    test_diffraction(0);
    test_diffraction(1);
    test_energy_conservation();
    test_absorption();
    test_bend_deflection();
    test_semivectorial_uniform();

    std::printf("\n%d failure(s).\n", g_failures);
    return g_failures == 0 ? 0 : 1;
}
