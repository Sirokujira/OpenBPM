/*
 * test_wabpm.cpp
 *
 * wabpm_step() (一般化 BPM 伝搬エンジン) の単体テスト
 *
 * テスト一覧:
 *   1. test_power_conservation_paraxial   – 損失なし均一媒質、近軸、100 ステップでパワー保存
 *   2. test_power_conservation_wideangle  – 同、広角 Pade(1,1) で 100 ステップ
 *   3. test_symmetry_preserved            – x 対称な初期電界は伝搬後も x 対称を保つ (近軸)
 *   4. test_freespace_diffraction         – 自由空間ガウシアンビーム幅が解析値 w(L) の 3% 以内
 *   5. test_tilt_wideangle_vs_paraxial    – チルトビーム: 広角は z*tan(θ), 近軸は z*sin(θ)
 *   6. test_power_conservation_halfvector – 均一媒質、半ベクトル x 偏波でパワー保存
 *
 * 依存: bpm/wabpm.cpp, include/bpm/wabpm.h のみ。外部フレームワーク不要。
 *
 * 設計上の注意点 (Ny=1 のとき dy を大きく取る理由):
 *   wabpm_step は Dirichlet 境界条件を使うため、Ny=1 のとき iy=0 の節点は
 *   「隣接節点が存在しない Dirichlet 境界」として扱われ、y-ラプラシアン寄与が
 *   |dm * 2/dy^2| となる。dy=1um では 0.25/step と大きく伝搬を歪めるが、
 *   dy=1mm (= 1e-3 m) では 2.5e-7/step と無視できるレベルになる。
 *   これにより 1D (x のみ) の伝搬を 2D コードで安価に近似できる。
 */

#include <cmath>
#include <cstdio>
#include <complex>
#include <vector>

#include "bpm/wabpm.h"

typedef std::complex<double> cplx;

// -----------------------------------------------------------------------
// 軽量テストマクロ
// -----------------------------------------------------------------------

// 条件チェック: 失敗時にファイル/行と msg を表示して 1 を返す
#define CHECK(cond, msg)                                                      \
    do {                                                                      \
        if (!(cond)) {                                                        \
            fprintf(stderr, "FAIL [%s:%d] %s\n", __FILE__, __LINE__, msg);   \
            return 1;                                                         \
        }                                                                     \
    } while (0)

// 数値近似チェック: |a - b| > tol なら失敗
#define CHECK_NEAR(a, b, tol, msg)                                             \
    do {                                                                       \
        double _diff = std::abs((double)(a) - (double)(b));                    \
        if (_diff > (tol)) {                                                   \
            fprintf(stderr,                                                    \
                    "FAIL [%s:%d] %s: |%g - %g| = %g > %g\n",                 \
                    __FILE__, __LINE__, msg,                                   \
                    (double)(a), (double)(b), _diff, (double)(tol));           \
            return 1;                                                          \
        }                                                                      \
    } while (0)

// -----------------------------------------------------------------------
// ヘルパー関数
// -----------------------------------------------------------------------

/** 強度の総和 (= 出力パワーの指標) */
static double total_power(const cplx *E, long N)
{
    double p = 0;
    for (long i = 0; i < N; i++)
        p += std::norm(E[i]);
    return p;
}

/** x 方向の強度重心 [物理座標 m] */
static double centroid_x(const cplx *E, long Nx, long Ny, double dx)
{
    double sx = 0, s = 0;
    for (long iy = 0; iy < Ny; iy++) {
        for (long ix = 0; ix < Nx; ix++) {
            double I = std::norm(E[ix + iy * Nx]);
            sx += I * (ix * dx);
            s  += I;
        }
    }
    return s > 0 ? sx / s : 0;
}

/**
 * 2 次元ガウシアンビーム E(x,y) = exp[-(x-x0)^2/w0^2 - (y-y0)^2/w0^2]
 *
 * 中心は (x0, y0) で指定する。x 対称性テストでは
 * x0 = (Nx-1)/2.0 * dx (格子点の中点) を渡すこと。
 * x0 = Nx/2.0 * dx (格子点上) にすると中心が 0.5 格子ずれて
 * 初期電界が非対称になる。
 */
static void make_gaussian(cplx *E, long Nx, long Ny,
                          double dx, double dy,
                          double x0, double y0, double w0)
{
    for (long iy = 0; iy < Ny; iy++) {
        for (long ix = 0; ix < Nx; ix++) {
            double x = ix * dx - x0;
            double y = iy * dy - y0;
            E[ix + iy * Nx] = std::exp(-(x * x + y * y) / (w0 * w0));
        }
    }
}

// -----------------------------------------------------------------------
// Test 1: 損失のない均一媒質 (スカラー近軸) でパワー保存
//
// 均一・無損失媒質での近軸 ADI は酉演算子になるため、
// パワーは理論上厳密に保存される。100 ステップでの誤差 < 1e-6 を確認。
// -----------------------------------------------------------------------
static int test_power_conservation_paraxial()
{
    const long Nx = 64, Ny = 64;
    const double lambda = 1.55e-6;
    const double k0     = 2 * M_PI / lambda;
    const double n0     = 1.5;
    const double dx     = 0.5e-6, dy = 0.5e-6, dz = 1.0e-6;

    wabpm_params W = {};
    W.Nx = Nx; W.Ny = Ny;
    W.dx = dx; W.dy = dy; W.dz = dz;
    W.k0 = k0; W.n0 = n0;
    W.wideangle = 0; W.pol = 0;  // 近軸・スカラー

    std::vector<cplx> E(Nx * Ny), n2(Nx * Ny, cplx(n0 * n0, 0));
    make_gaussian(E.data(), Nx, Ny, dx, dy,
                  (Nx - 1) / 2.0 * dx, (Ny - 1) / 2.0 * dy, 2.0e-6);

    double p0 = total_power(E.data(), Nx * Ny);
    for (int s = 0; s < 100; s++)
        wabpm_step(&W, E.data(), n2.data());
    double p1 = total_power(E.data(), Nx * Ny);

    CHECK_NEAR(p1 / p0, 1.0, 1e-6,
               "パワー保存 (近軸, 100 step): |p1/p0 - 1| < 1e-6");
    printf("PASS test_power_conservation_paraxial   p1/p0 = %.9f\n", p1 / p0);
    return 0;
}

// -----------------------------------------------------------------------
// Test 2: 損失のない均一媒質 (広角 Pade(1,1)) でパワー保存
//
// 広角 ADI も酉演算子となるため、パワー保存は近軸と同様に成立する。
// -----------------------------------------------------------------------
static int test_power_conservation_wideangle()
{
    const long Nx = 64, Ny = 64;
    const double lambda = 1.55e-6;
    const double k0     = 2 * M_PI / lambda;
    const double n0     = 1.5;
    const double dx     = 0.5e-6, dy = 0.5e-6, dz = 1.0e-6;

    wabpm_params W = {};
    W.Nx = Nx; W.Ny = Ny;
    W.dx = dx; W.dy = dy; W.dz = dz;
    W.k0 = k0; W.n0 = n0;
    W.wideangle = 1; W.pol = 0;  // 広角・スカラー

    std::vector<cplx> E(Nx * Ny), n2(Nx * Ny, cplx(n0 * n0, 0));
    make_gaussian(E.data(), Nx, Ny, dx, dy,
                  (Nx - 1) / 2.0 * dx, (Ny - 1) / 2.0 * dy, 2.0e-6);

    double p0 = total_power(E.data(), Nx * Ny);
    for (int s = 0; s < 100; s++)
        wabpm_step(&W, E.data(), n2.data());
    double p1 = total_power(E.data(), Nx * Ny);

    CHECK_NEAR(p1 / p0, 1.0, 1e-6,
               "パワー保存 (広角, 100 step): |p1/p0 - 1| < 1e-6");
    printf("PASS test_power_conservation_wideangle  p1/p0 = %.9f\n", p1 / p0);
    return 0;
}

// -----------------------------------------------------------------------
// Test 3: x 対称な初期電界は伝搬後も x 対称を保つ (近軸・スカラー)
//
// ガウシアン中心を格子点の「中点」(= (Nx-1)/2 * dx) に置くと
// 初期電界が x 対称になる。Peaceman-Rachford ADI はこの対称性を保存する。
//
// max|E(ix) - E(Nx-1-ix)| < 1e-8 を確認 (浮動小数点の有限精度による
// Thomas 法での加算順序の違いが ~1e-12 程度の非対称性を生じさせうる)。
// -----------------------------------------------------------------------
static int test_symmetry_preserved()
{
    const long Nx = 64, Ny = 32;
    const double lambda = 1.55e-6;
    const double k0     = 2 * M_PI / lambda;
    const double n0     = 1.44;
    const double dx     = 0.5e-6, dy = 0.5e-6, dz = 2.0e-6;

    wabpm_params W = {};
    W.Nx = Nx; W.Ny = Ny;
    W.dx = dx; W.dy = dy; W.dz = dz;
    W.k0 = k0; W.n0 = n0;
    W.wideangle = 0; W.pol = 0;

    std::vector<cplx> E(Nx * Ny), n2(Nx * Ny, cplx(n0 * n0, 0));
    // 中心を格子点間 (= (Nx-1)/2 * dx) に置き x 対称なガウシアンを作る
    make_gaussian(E.data(), Nx, Ny, dx, dy,
                  (Nx - 1) / 2.0 * dx, (Ny - 1) / 2.0 * dy, 3.0e-6);

    for (int s = 0; s < 50; s++)
        wabpm_step(&W, E.data(), n2.data());

    // |E(ix)| と |E(Nx-1-ix)| が一致することを確認
    double max_asym = 0;
    for (long iy = 0; iy < Ny; iy++) {
        for (long ix = 0; ix < Nx / 2; ix++) {
            double diff = std::abs(E[ix + iy * Nx]) - std::abs(E[(Nx - 1 - ix) + iy * Nx]);
            if (std::abs(diff) > max_asym)
                max_asym = std::abs(diff);
        }
    }
    CHECK_NEAR(max_asym, 0.0, 1e-8,
               "x 対称性保存 (50 step): max|E(ix) - E(Nx-1-ix)| < 1e-8");
    printf("PASS test_symmetry_preserved            max_asym = %.2e\n", max_asym);
    return 0;
}

// -----------------------------------------------------------------------
// Test 4: 均一媒質でのガウシアンビーム回折 (1D 近似, 近軸)
//
// 解析解: w(z) = w0 * sqrt(1 + (z/zR)^2), zR = pi*w0^2*n0/lambda
// 1D 近似として Ny=1, dy=1mm を使用する。dy が大きいと y-ラプラシアン
// 寄与 (|dm*2/dy^2| ∝ 1/dy^2) が無視できるレベルになるため、
// 実質的に 1D (x のみ) の伝搬を計算できる。
//
// 2 次モーメントから測定幅を求め、解析値と 3% 以内で一致することを確認。
// パラメータは end-to-end テストの freespace.ofd に近い値を採用。
// -----------------------------------------------------------------------
static int test_freespace_diffraction()
{
    // Ny=1, dy=1mm で y-ラプラシアン寄与を無視できるようにする
    const long Nx = 128, Ny = 1;
    const double lambda = 1.55e-6;
    const double k0     = 2 * M_PI / lambda;
    const double n0     = 1.444;   // freespace.ofd の参照屈折率
    const double dx     = 0.5e-6;
    const double dy     = 1.0e-3;  // 大きな dy → y-ラプラシアン ∝ 1/dy^2 が微小
    const double dz     = 1.0e-6;
    const double w0     = 5.0e-6;  // freespace.ofd に合わせた入射ビームウェスト

    wabpm_params W = {};
    W.Nx = Nx; W.Ny = Ny;
    W.dx = dx; W.dy = dy; W.dz = dz;
    W.k0 = k0; W.n0 = n0;
    W.wideangle = 0; W.pol = 0;

    std::vector<cplx> E(Nx), n2(Nx, cplx(n0 * n0, 0));
    // x 対称なガウシアン (中心を格子間点に置く)
    const double x0 = (Nx - 1) / 2.0 * dx;
    for (long ix = 0; ix < Nx; ix++)
        E[ix] = std::exp(-std::pow((ix * dx - x0) / w0, 2));

    const int steps = 60;
    for (int s = 0; s < steps; s++)
        wabpm_step(&W, E.data(), n2.data());

    // 解析値: レイリー長 zR = pi*w0^2*n0/lambda, 出力幅 w(L)
    const double L       = steps * dz;
    const double zR      = M_PI * w0 * w0 * n0 / lambda;
    const double w_theory = w0 * std::sqrt(1.0 + (L / zR) * (L / zR));

    // 2 次モーメントから 1/e^2 径を計算 (重心は x0 に固定されている前提)
    double s_tot = 0, sx2 = 0;
    for (long ix = 0; ix < Nx; ix++) {
        double I = std::norm(E[ix]);
        double x = ix * dx - x0;
        s_tot += I;
        sx2   += I * x * x;
    }
    const double w_meas = 2.0 * std::sqrt(sx2 / s_tot);

    const double rel = std::abs(w_meas - w_theory) / w_theory;
    CHECK_NEAR(rel, 0.0, 0.03, "ガウシアン回折 w(L) vs 解析値 (3%)");
    printf("PASS test_freespace_diffraction         "
           "w_meas=%.3f um  w_theory=%.3f um  rel=%.2f%%\n",
           w_meas * 1e6, w_theory * 1e6, rel * 100);
    return 0;
}

// -----------------------------------------------------------------------
// Test 5: チルトビーム — 広角と近軸で重心変位が異なることを確認 (1D 近似)
//
// 横方向波数 kx = k0*n0*sin(θ) を位相ランプとして与えたガウシアンを伝搬し、
// 重心変位を測定する。
//
//   広角 (Pade(1,1)): 変位 ≈ z * tan(θ)  (厳密値)
//   近軸             : 変位 ≈ z * sin(θ)
//
// θ = 15° では tan と sin の比が 1.035 (3.5% 差) あり、これを識別する。
// Ny=1, dy=1mm で y-ラプラシアン寄与を無視 (|dm*2/dy^2| ≈ 2.5e-7/step)。
// -----------------------------------------------------------------------
static int test_tilt_wideangle_vs_paraxial()
{
    // Ny=1, dy=1mm で y-ラプラシアン寄与を無視できるようにする
    const long Nx = 320, Ny = 1;
    const double lambda = 1.55e-6;
    const double k0     = 2 * M_PI / lambda;
    const double n0     = 1.0;
    const double dx     = 0.2e-6;
    const double dy     = 1.0e-3;  // 大きな dy → y-ラプラシアン ∝ 1/dy^2 が微小
    const double dz     = 1.0e-6;
    const double w0     = 3.0e-6;
    const double theta  = 15.0 * M_PI / 180.0;
    // wabpm は exp(+iωt) 規約を使うため kx > 0 は -x 方向へドリフトさせる。
    // kx < 0 にすることでビームが +x 方向へ変位し、
    // 変位 = +z*sin(θ) (近軸) / +z*tan(θ) (広角) となる。
    const double kx     = -k0 * n0 * std::sin(theta);
    const double x0     = (Nx - 1) / 2.0 * dx;        // グリッド中心

    std::vector<cplx> n2(Nx, cplx(n0 * n0, 0));
    const int steps = 80;
    const double L = steps * dz;

    // チルト + ガウシアンの初期電界を生成
    auto make_tilt = [&](std::vector<cplx> &E) {
        for (long ix = 0; ix < Nx; ix++) {
            double x = ix * dx - x0;
            // ガウシアン × 横方向位相ランプ (= チルト入射)
            E[ix] = std::exp(-x * x / (w0 * w0)) * std::exp(cplx(0, kx * x));
        }
    };

    wabpm_params W = {};
    W.Nx = Nx; W.Ny = Ny;
    W.dx = dx; W.dy = dy; W.dz = dz;
    W.k0 = k0; W.n0 = n0;
    W.pol = 0;

    // --- 近軸 ---
    W.wideangle = 0;
    std::vector<cplx> E_paraxial(Nx);
    make_tilt(E_paraxial);
    for (int s = 0; s < steps; s++)
        wabpm_step(&W, E_paraxial.data(), n2.data());
    const double cx_paraxial = centroid_x(E_paraxial.data(), Nx, Ny, dx) - x0;

    // --- 広角 ---
    W.wideangle = 1;
    std::vector<cplx> E_wide(Nx);
    make_tilt(E_wide);
    for (int s = 0; s < steps; s++)
        wabpm_step(&W, E_wide.data(), n2.data());
    const double cx_wide = centroid_x(E_wide.data(), Nx, Ny, dx) - x0;

    const double disp_exact = L * std::tan(theta);   // z*tan(θ): 厳密値
    const double disp_parax = L * std::sin(theta);   // z*sin(θ): 近軸値

    // 広角は近軸より大きい変位 (tan(θ) > sin(θ) for 0 < θ < 90°)
    CHECK(cx_wide > cx_paraxial,
          "広角の変位 > 近軸の変位 (θ=15°)");
    // 広角は z*tan(θ) の 5% 以内
    CHECK_NEAR(cx_wide / disp_exact, 1.0, 0.05,
               "広角変位 ≈ z*tan(θ) (5%)");
    // 近軸は z*sin(θ) の 3% 以内
    CHECK_NEAR(cx_paraxial / disp_parax, 1.0, 0.03,
               "近軸変位 ≈ z*sin(θ) (3%)");

    printf("PASS test_tilt_wideangle_vs_paraxial    "
           "z*tan=%.3f um  wide=%.3f um  z*sin=%.3f um  parax=%.3f um\n",
           disp_exact * 1e6, cx_wide * 1e6, disp_parax * 1e6, cx_paraxial * 1e6);
    return 0;
}

// -----------------------------------------------------------------------
// Test 6: 半ベクトル x 偏波 (pol=1) でのパワー保存
//
// 均一媒質では半ベクトルステンシルはスカラーステンシルと等価になるため、
// パワー保存が同様に成立するはず (Stern 差分の均一媒質での退化を確認)。
// -----------------------------------------------------------------------
static int test_power_conservation_halfvector()
{
    const long Nx = 48, Ny = 48;
    const double lambda = 1.55e-6;
    const double k0     = 2 * M_PI / lambda;
    const double n0     = 1.5;
    const double dx     = 0.5e-6, dy = 0.5e-6, dz = 1.0e-6;

    wabpm_params W = {};
    W.Nx = Nx; W.Ny = Ny;
    W.dx = dx; W.dy = dy; W.dz = dz;
    W.k0 = k0; W.n0 = n0;
    W.wideangle = 0; W.pol = 1;  // 近軸・x 偏波 (半ベクトル)

    std::vector<cplx> E(Nx * Ny), n2(Nx * Ny, cplx(n0 * n0, 0));
    make_gaussian(E.data(), Nx, Ny, dx, dy,
                  (Nx - 1) / 2.0 * dx, (Ny - 1) / 2.0 * dy, 2.0e-6);

    double p0 = total_power(E.data(), Nx * Ny);
    for (int s = 0; s < 50; s++)
        wabpm_step(&W, E.data(), n2.data());
    double p1 = total_power(E.data(), Nx * Ny);

    CHECK_NEAR(p1 / p0, 1.0, 1e-6,
               "パワー保存 (半ベクトル x 偏波, 50 step): |p1/p0 - 1| < 1e-6");
    printf("PASS test_power_conservation_halfvector p1/p0 = %.9f\n", p1 / p0);
    return 0;
}

// -----------------------------------------------------------------------
// エントリポイント
// -----------------------------------------------------------------------
int main()
{
    int fail = 0;
    fail += test_power_conservation_paraxial();
    fail += test_power_conservation_wideangle();
    fail += test_symmetry_preserved();
    fail += test_freespace_diffraction();
    fail += test_tilt_wideangle_vs_paraxial();
    fail += test_power_conservation_halfvector();

    if (fail == 0)
        printf("\nAll 6 tests passed.\n");
    else
        fprintf(stderr, "\n%d test(s) FAILED.\n", fail);
    return fail;
}
