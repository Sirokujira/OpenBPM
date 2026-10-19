/*
test_pml.cpp

PML (複素座標伸長) の単体テスト。伝搬エンジンは bpm/wabpm.cpp (倍精度)。
外部テストフレームワークに依存しない自己完結ハーネス。

検証内容:
  1. 後方互換   : PML 係数を全て 1 (伸長なし) にした場合、g 配列を渡さない
                  従来経路と**ビット一致**すること
  2. 境界反射   : 傾いたガウシアンを側壁へ入射させ、壁の無い広領域の参照計算と
                  内部窓の残留電力を比べる。PML では反射が入射電力の 1e-5 以下
                  (実測 ~1e-8)、PML 無し (Dirichlet 境界) では 0.1 以上になる
  3. 理論との一致: 目標反射率 R0 を弱く (1e-3) 取ったとき、往復電力反射率が
                  連続体の理論値 R = R0^(2 sin(theta)) と桁で一致すること
                  (強吸収側は離散化由来の反射床が支配するため弱吸収側で判定する)

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

static void check_true(const std::string &name, bool cond, const std::string &detail)
{
    std::printf("[%s] %-46s %s\n", cond ? " OK " : "FAIL", name.c_str(), detail.c_str());
    if (!cond) g_failures++;
}

// 中心化座標 (グリッド中心を原点)
static inline double coord(long i, long N, double d)
{
    return (static_cast<double>(i) - (N - 1) / 2.0) * d;
}

// PML の伸長係数 (sol/solve_bpm.cpp の bpmPmlAxis と同じ規約)
//   s(u) = 1 - i*sigma(u), sigma(d) = smax*(d/W)^3
//   gm[i] = 1/(s(cent_i) * s(face_{i-1/2})), gp[i] = 1/(s(cent_i) * s(face_{i+1/2}))
static double pml_sigma(double u, double lo, double hi, double W, double smax)
{
    double d = 0;
    if (u < lo)      d = lo - u;
    else if (u > hi) d = u - hi;
    if (d <= 0) return 0;
    if (d > W) d = W;
    return smax * std::pow(d / W, 3.0);
}

static void pml_axis(int N, double d, double W, double smax,
                     std::vector<cplx> &gm, std::vector<cplx> &gp)
{
    const double half = 0.5 * N * d;      // 領域の半幅 (中心原点)
    const double lo = -half + W;
    const double hi = half - W;
    gm.resize(N);
    gp.resize(N);
    for (int i = 0; i < N; i++) {
        const double uc = coord(i, N, d);
        const double um = uc - 0.5 * d;
        const double up = uc + 0.5 * d;
        const cplx sc(1.0, -pml_sigma(uc, lo, hi, W, smax));
        const cplx sm(1.0, -pml_sigma(um, lo, hi, W, smax));
        const cplx sp(1.0, -pml_sigma(up, lo, hi, W, smax));
        gm[i] = 1.0 / (sc * sm);
        gp[i] = 1.0 / (sc * sp);
    }
}

// 傾いたガウシアンを伝搬させ、内部窓 |x|,|y| < win に残る電力を返す。
// 一様媒質 (n = n0) なので導波はなく、ビームは壁へ向かって walk-off する。
static double propagate(int Nx, int Ny, double dx, double dz, int nz,
                        double k0, double n0, double w0, double tiltdeg,
                        double pmlW, double R0, double win, double *pin)
{
    wabpm_params W;
    W.Nx = Nx;
    W.Ny = Ny;
    W.dx = dx;
    W.dy = dx;
    W.dz = dz;
    W.k0 = k0;
    W.n0 = n0;
    W.wideangle = 0;
    W.pol = 0;
    W.symx = 0;
    W.symy = 0;

    std::vector<cplx> gxm, gxp, gym, gyp;
    if (pmlW > 0) {
        // かすめ入射 (kx = k0*n0) で往復振幅反射率が R0 になる sigma_max
        const double smax = -(3.0 + 1.0) * std::log(R0) / (2 * k0 * n0 * pmlW);
        pml_axis(Nx, dx, pmlW, smax, gxm, gxp);
        pml_axis(Ny, dx, pmlW, smax, gym, gyp);
        W.gxm = &gxm[0];
        W.gxp = &gxp[0];
        W.gym = &gym[0];
        W.gyp = &gyp[0];
    }

    const double ktx = k0 * n0 * std::sin(tiltdeg * M_PI / 180.0);
    std::vector<cplx> E((size_t)Nx * Ny), n2((size_t)Nx * Ny, cplx(n0 * n0, 0.0));
    double p0 = 0;
    for (int iy = 0; iy < Ny; iy++) {
        for (int ix = 0; ix < Nx; ix++) {
            const double x = coord(ix, Nx, dx);
            const double y = coord(iy, Ny, dx);
            const double amp = std::exp(-(x * x + y * y) / (w0 * w0));
            E[ix + (size_t)iy * Nx] = amp * cplx(std::cos(ktx * x), -std::sin(ktx * x));
            p0 += amp * amp;
        }
    }
    if (pin != NULL) *pin = p0;

    for (int iz = 0; iz < nz; iz++) {
        wabpm_step(&W, &E[0], &n2[0]);
    }

    double p = 0;
    for (int iy = 0; iy < Ny; iy++) {
        for (int ix = 0; ix < Nx; ix++) {
            const double x = coord(ix, Nx, dx);
            const double y = coord(iy, Ny, dx);
            if ((std::fabs(x) < win) && (std::fabs(y) < win)) {
                p += std::norm(E[ix + (size_t)iy * Nx]);
            }
        }
    }
    return p;
}

int main(void)
{
    const double lambda = 1.55e-6;
    const double k0 = 2 * M_PI / lambda;
    const double n0 = 1.444;
    const double dx = 0.25e-6;
    const double dz = 0.5e-6;

    // --- 1. 後方互換 : 伸長なし (sigma = 0) の PML 配列は「PML 無し」と同一 ---
    {
        const int N = 48;
        wabpm_params Wa, Wb;
        Wa.Nx = Wb.Nx = N;
        Wa.Ny = Wb.Ny = N;
        Wa.dx = Wb.dx = dx;
        Wa.dy = Wb.dy = dx;
        Wa.dz = Wb.dz = dz;
        Wa.k0 = Wb.k0 = k0;
        Wa.n0 = Wb.n0 = n0;
        Wa.wideangle = Wb.wideangle = 0;
        Wa.pol = Wb.pol = 0;
        Wa.symx = Wb.symx = 0;
        Wa.symy = Wb.symy = 0;
        std::vector<cplx> one((size_t)N, cplx(1.0, 0.0));
        Wb.gxm = Wb.gxp = Wb.gym = Wb.gyp = &one[0];

        std::vector<cplx> Ea((size_t)N * N), Eb((size_t)N * N);
        std::vector<cplx> n2((size_t)N * N, cplx(n0 * n0, 0.0));
        for (int iy = 0; iy < N; iy++) {
            for (int ix = 0; ix < N; ix++) {
                const double x = coord(ix, N, dx), y = coord(iy, N, dx);
                const cplx v = std::exp(-(x * x + y * y) / (3e-6 * 3e-6));
                Ea[ix + (size_t)iy * N] = v;
                Eb[ix + (size_t)iy * N] = v;
            }
        }
        for (int iz = 0; iz < 20; iz++) {
            wabpm_step(&Wa, &Ea[0], &n2[0]);
            wabpm_step(&Wb, &Eb[0], &n2[0]);
        }
        double maxdiff = 0;
        for (size_t i = 0; i < Ea.size(); i++) {
            maxdiff = std::max(maxdiff, std::abs(Ea[i] - Eb[i]));
        }
        char det[128];
        std::snprintf(det, sizeof(det), "max|E_pml(sigma=0) - E_nopml| = %.3g", maxdiff);
        check_true("pml identity (sigma = 0)", maxdiff == 0.0, det);
    }

    // --- 2. 境界反射 : 壁の無い広領域を参照にして反射電力を測る ---
    // 小領域 (Nx = 160, ±20um) に 15 度で入射させ、内部窓 |x|,|y| < 16um の
    // 残留電力を、広領域 (Nx = 640, ±80um) の同じ窓と比べる。
    const double w0 = 4e-6, tilt = 15.0, win = 16e-6;
    const int nz = 400;
    double pin = 0, pin_ref = 0;
    const double ref  = propagate(640, 160, dx, dz, nz, k0, n0, w0, tilt, 0.0, 1.0, win, &pin_ref);
    const double pml  = propagate(160, 160, dx, dz, nz, k0, n0, w0, tilt, 4e-6, 1e-20, win, &pin);
    const double hard = propagate(160, 160, dx, dz, nz, k0, n0, w0, tilt, 0.0, 1.0, win, NULL);
    {
        const double r_pml  = std::fabs(pml - ref) / pin;
        const double r_hard = std::fabs(hard - ref) / pin;
        char det[192];
        std::snprintf(det, sizeof(det), "R_pml = %.3g (< 1e-5), R_dirichlet = %.3g (> 0.1)",
                      r_pml, r_hard);
        check_true("pml reflection vs wall-free reference",
                   (r_pml < 1e-5) && (r_hard > 0.1), det);
    }

    // --- 3. 理論との一致 : 弱吸収側で R = R0^(2 sin(theta)) ---
    // R0 を大きく (1e-3) 取ると離散化由来の反射床より連続体の理論反射が支配する。
    {
        const double R0 = 1e-3;
        const double weak = propagate(160, 160, dx, dz, nz, k0, n0, w0, tilt, 4e-6, R0, win, NULL);
        const double r = std::fabs(weak - ref) / pin;
        const double theory = std::pow(R0, 2 * std::sin(tilt * M_PI / 180.0));
        const double ratio = r / theory;
        char det[192];
        std::snprintf(det, sizeof(det), "R = %.3g, theory R0^(2 sin15deg) = %.3g (ratio %.2f)",
                      r, theory, ratio);
        // 連続体理論なので離散化誤差ぶんのずれを許す (係数 3 以内)
        check_true("pml reflection follows R0^(2 sin theta)",
                   (ratio > 1.0 / 3.0) && (ratio < 3.0), det);
    }

    std::printf("\n%s : %d failure(s)\n", g_failures ? "FAILED" : "PASSED", g_failures);
    return g_failures ? 1 : 0;
}
