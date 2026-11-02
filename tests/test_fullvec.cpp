/*
test_fullvec.cpp

完全ベクトル BPM (bpm/wabpm.cpp の wabpm_fv_*) の単体テスト。
外部テストフレームワークに依存しない自己完結ハーネス。

検証内容:
  1. スカラー帰着 : 一様媒質では交差項 Pxy / Pyx が厳密に 0 になり、Ex / Ey は
                    それぞれ独立なスカラー伝搬とビット一致する (近軸・広角の両方)
  2. HE11 の実効屈折率 : 高コントラストのステップインデックスファイバ
                    (n1 = 1.45 / n2 = 1.0 (空気), a = 0.5um, lambda = 1.55um, V = 2.13)
                    で、ベクトル厳密解 HE11 = 1.176447 と比べる。スカラー LP01 = 1.225891
                    との差 0.049 が「偏波補正」で、スカラー・半ベクトルでは再現できない
  3. 偏波補正の大きさ : (neff_FV - neff_scalar) / (HE11 - LP01) が 1 に近いこと
                    (同じ格子の階段近似誤差が差を取ると相殺される)
  4. 交差偏波成分 : HE11 (x 偏波) の Ey は sin(2 phi) 型の四重極で、その電力比は
                    厳密場 (Snyder & Love) から 1.0561e-2。半ベクトルでは 0 になる量で、
                    交差項の大きさを直接検証する。Ey が x, y の両方に奇であることも確認
  5. 縮退 : x 偏波と y 偏波の基本モードの neff が一致する (正方格子の対称性)
  6. 伝搬の整合 (近軸) : 求めた FV モードを n0 = neff で伝搬させると形が保たれ、
                    位相がほぼ進まない (交差項込みの伝搬式がモードを固有状態として保つ)
  7. 広角の伝搬 : 中程度のコントラスト (1.45 / 1.40) のファイバで、参照屈折率を
                    クラッド (n0 = 1.40) に取ってモードを広角 (Pade(1,1)) で伝搬させ、
                    形がほぼ保たれ (overlap > 0.98。既存の半ベクトル広角でも同条件で 0.990)、
                    位相から求めた neff がモードの neff と
                    1e-4 以内で一致すること。近軸では理論どおり +5.3e-5 (= n0 + (neff^2 - n0^2)/(2 n0)
                    - neff) ずれる。広角は Pade の理論誤差 -3e-6 より大きい -5e-5 ずれるが、
                    これは既存の広角 ADI (スカラーも同じ値) の分離誤差
                    (dp^2 - dm^2) Px Py = i dz/(4 k0^3 n0^3) Px Py が毎ステップ O(dz) で
                    残るためで、完全ベクトル固有の誤差ではない (チェックリスト第 24 回参照)

既知の制約 (本テストの範囲外。docs/implementation-checklist.md 参照):
  空気クラッド級の高コントラスト (Δn >~ 0.25) の 2D 構造では、半ベクトル (Stern) と
  同じく横方向演算子が非対称で高周波側に複素固有値を持つため、伝搬で偽の成長モードが
  現れる (広角 Pade では特に顕著)。モード解析 (虚軸伝搬) と近軸のモード伝搬は影響を
  受けない。テスト 6 を近軸、テスト 7 を中程度のコントラストで行うのはこのため。

格子と許容誤差の根拠 (dx = 0.05um = a/10):
  同じ差分式を疎行列の固有値問題として直接解き、サブピクセル平均で階段近似の
  雑音を除いて格子を細かくすると、HE11 との誤差は
    dx = 0.05 / 0.025 / 0.0125 / 0.00833 um で +1.92e-3 / +1.12e-3 / +6.3e-4 / +4.6e-4
  と 1 次収束し (外挿値の誤差 1.1e-4)、離散化が厳密解へ向かうことを確認している。
  斜めの界面 (円周) で法線成分 n^2 E_r の連続を x / y 方向の差分でしか表さないため
  1 次精度になる。本テストの階段近似 (サブピクセル平均なし) では dx = 0.05um で
  +2.7e-3 なので許容を 4e-3 とする (偏波補正 0.049 の 8%、半ベクトルの誤差 5.4e-3 未満)。

失敗が 1 件でもあれば終了コード 1 を返す (CTest 連携用)。
*/

#include "bpm/wabpm.h"
#include "bessel_ref.h"   // std::cyl_bessel_j/k は libc++ 未実装のため自前実装を使用

#include <complex>
#include <vector>
#include <cmath>
#include <cstdio>
#include <string>
#include <algorithm>

typedef std::complex<double> cplx;

static int g_failures = 0;

static void check_true(const std::string &name, bool cond, const std::string &detail)
{
    std::printf("[%s] %-44s %s\n", cond ? " OK " : "FAIL", name.c_str(), detail.c_str());
    if (!cond) g_failures++;
}

// ---------------------------------------------------------------
// ステップインデックスファイバの厳密解 (Snyder & Love, 12 章)
// ---------------------------------------------------------------
struct Fiber {
    double n1, n2, a, k0;
};

static void uw_of(const Fiber &f, double ne, double *U, double *W)
{
    *U = f.a * f.k0 * std::sqrt(f.n1 * f.n1 - ne * ne);
    *W = f.a * f.k0 * std::sqrt(ne * ne - f.n2 * f.n2);
}

// HE / EH (nu = 1) の固有値方程式
//   [J'/(UJ) + K'/(WK)] [J'/(UJ) + (n2/n1)^2 K'/(WK)] = (1/U^2 + 1/W^2)(1/U^2 + (n2/n1)^2/W^2)
static double he1_eq(const Fiber &f, double ne)
{
    double U, W;
    uw_of(f, ne, &U, &W);
    const double J1 = bessel_ref_j(1, U);
    const double K1 = bessel_ref_k(1, W);
    const double Jp = bessel_ref_j(0, U) - J1 / U;            // J1'
    const double Kp = -bessel_ref_k(0, W) - K1 / W;           // K1'
    const double jt = Jp / (U * J1);
    const double kt = Kp / (W * K1);
    const double r = (f.n2 / f.n1) * (f.n2 / f.n1);
    return (jt + kt) * (jt + r * kt)
         - (1 / (U * U) + 1 / (W * W)) * (1 / (U * U) + r / (W * W));
}

// LP01 (弱導波近似) の固有値方程式 U J1/J0 = W K1/K0
static double lp01_eq(const Fiber &f, double ne)
{
    double U, W;
    uw_of(f, ne, &U, &W);
    return U * bessel_ref_j(1, U) / bessel_ref_j(0, U)
         - W * bessel_ref_k(1, W) / bessel_ref_k(0, W);
}

// neff を n1 側から下へ走査し、最初の (= 最大の) 根を二分法で求める。
// J の零点では方程式が発散して符号が変わるので、両端の値が有限で小さい区間だけ採る。
static double largest_root(const Fiber &f, double (*eq)(const Fiber &, double))
{
    const int ns = 20000;
    double hi = f.n1 - 1e-9;
    double fhi = eq(f, hi);
    for (int k = 1; k <= ns; k++) {
        const double lo = f.n1 - 1e-9 - (f.n1 - f.n2 - 2e-9) * k / ns;
        const double flo = eq(f, lo);
        if (((flo > 0) != (fhi > 0)) && (std::fabs(flo) < 1e3) && (std::fabs(fhi) < 1e3)) {
            double a = lo, b = hi, fa = flo;
            for (int it = 0; it < 200; it++) {
                const double m = 0.5 * (a + b);
                const double fm = eq(f, m);
                if ((fm > 0) == (fa > 0)) { a = m; fa = fm; } else { b = m; }
            }
            return 0.5 * (a + b);
        }
        hi = lo;
        fhi = flo;
    }
    return -1;
}

// HE11 (x 偏波) の横電界 (Snyder & Love 表 12-3, nu = 1)。R = r/a
static void he11_field(const Fiber &f, double ne, double R, double phi, double *ex, double *ey)
{
    double U, W;
    uw_of(f, ne, &U, &W);
    const double V2 = U * U + W * W;
    const double b1 = (bessel_ref_j(0, U) - bessel_ref_j(2, U)) / (2 * U * bessel_ref_j(1, U));
    const double b2 = -(bessel_ref_k(0, W) + bessel_ref_k(2, W)) / (2 * W * bessel_ref_k(1, W));
    const double F2 = (V2 / (U * U * W * W)) / (b1 + b2);
    const double a1 = (F2 - 1) / 2, a2 = (F2 + 1) / 2;
    if (R < 1) {
        const double J1 = bessel_ref_j(1, U);
        *ex = -(a1 * bessel_ref_j(0, U * R) + a2 * bessel_ref_j(2, U * R) * std::cos(2 * phi)) / J1;
        *ey = -(a2 * bessel_ref_j(2, U * R) * std::sin(2 * phi)) / J1;
    }
    else {
        const double K1 = bessel_ref_k(1, W);
        *ex = -(U / W) * (a1 * bessel_ref_k(0, W * R) - a2 * bessel_ref_k(2, W * R) * std::cos(2 * phi)) / K1;
        *ey = (U / W) * (a2 * bessel_ref_k(2, W * R) * std::sin(2 * phi)) / K1;
    }
}

// 厳密場の Ey 電力比 ∫|ey|^2 / ∫(|ex|^2 + |ey|^2) (極座標の中点則)
static double he11_ey_fraction(const Fiber &f, double ne)
{
    const int nr = 4000, np = 360;
    const double Rmax = 10.0;   // W = 1.7 なので exp(-2*1.7*10) で無視できる
    double px = 0, py = 0;
    for (int ir = 0; ir < nr; ir++) {
        const double R = (ir + 0.5) * Rmax / nr;
        for (int ip = 0; ip < np; ip++) {
            const double phi = (ip + 0.5) * 2 * M_PI / np;
            double ex, ey;
            he11_field(f, ne, R, phi, &ex, &ey);
            px += ex * ex * R;
            py += ey * ey * R;
        }
    }
    return py / (px + py);
}

// ---------------------------------------------------------------

static void setup(wabpm_params *W, int N, double dx, double k0, double n0)
{
    W->Nx = N;
    W->Ny = N;
    W->dx = dx;
    W->dy = dx;
    W->dz = 0;
    W->k0 = k0;
    W->n0 = n0;
    W->wideangle = 0;
    W->pol = WABPM_POL_FULL;
    W->symx = 0;
    W->symy = 0;
}

static inline double coord(int i, int N, double d) { return (i - (N - 1) / 2.0) * d; }

int main(void)
{
    const double lambda = 1.55e-6;
    const double k0 = 2 * M_PI / lambda;

    // --- 1. スカラー帰着 (一様媒質) ---
    for (int wa = 0; wa < 2; wa++) {
        const int N = 48;
        const double dx = 0.25e-6, n0 = 1.444;
        wabpm_params Wf, Ws;
        setup(&Wf, N, dx, k0, n0);
        Wf.dz = 0.5e-6;
        Wf.wideangle = wa;
        Ws = Wf;
        Ws.pol = 0;
        std::vector<cplx> n2((size_t)N * N, cplx(1.46 * 1.46, -1e-4));
        std::vector<cplx> Ex((size_t)N * N), Ey((size_t)N * N);
        for (int iy = 0; iy < N; iy++) {
            for (int ix = 0; ix < N; ix++) {
                const double x = coord(ix, N, dx), y = coord(iy, N, dx);
                Ex[ix + (size_t)iy * N] = std::exp(-(x * x + y * y) / 9e-12);
                Ey[ix + (size_t)iy * N] = cplx(0.0, 0.3) * std::exp(-((x - 1e-6) * (x - 1e-6) + y * y) / 4e-12);
            }
        }
        std::vector<cplx> Sx = Ex, Sy = Ey;
        for (int iz = 0; iz < 20; iz++) {
            wabpm_fv_step(&Wf, &Ex[0], &Ey[0], &n2[0]);
            wabpm_step(&Ws, &Sx[0], &n2[0]);
            wabpm_step(&Ws, &Sy[0], &n2[0]);
        }
        double md = 0;
        for (size_t i = 0; i < Ex.size(); i++) {
            md = std::max(md, std::max(std::abs(Ex[i] - Sx[i]), std::abs(Ey[i] - Sy[i])));
        }
        char det[128];
        std::snprintf(det, sizeof(det), "max|E_fv - E_scalar| = %.3g", md);
        check_true(wa ? "uniform medium -> scalar (wide-angle)" : "uniform medium -> scalar (paraxial)",
                   md == 0.0, det);
    }

    // --- 2-5. 高コントラストファイバの HE11 ---
    const Fiber fib = {1.45, 1.0, 0.5e-6, k0};
    const double he11 = largest_root(fib, he1_eq);
    const double lp01 = largest_root(fib, lp01_eq);
    const int N = 80;
    const double dx = 0.05e-6;
    std::vector<cplx> n2((size_t)N * N);
    for (int iy = 0; iy < N; iy++) {
        for (int ix = 0; ix < N; ix++) {
            const double x = coord(ix, N, dx), y = coord(iy, N, dx);
            n2[ix + (size_t)iy * N] = ((x * x + y * y) < fib.a * fib.a) ? fib.n1 * fib.n1 : fib.n2 * fib.n2;
        }
    }
    wabpm_params W;
    setup(&W, N, dx, k0, 1.2);
    // 段階 1 (後退オイラー) の虚軸ステップ : 格子スケールの 50 倍
    const double a_grid = 1.0 / (8.0 / (dx * dx));
    W.dz = 50 * a_grid * 4 * k0 * W.n0;

    std::vector<cplx> Ex((size_t)N * N), Ey((size_t)N * N);
    double nfv = 0;
    const int ok = wabpm_fv_find_mode(&W, &n2[0], 1, 20000, 1e-11, &Ex[0], &Ey[0], &nfv);

    // 比較用 : スカラー (LP01 相当) と半ベクトル (x 偏波, Stern)
    double nsc = 0, nsv = 0;
    {
        std::vector<cplx> m((size_t)N * N);
        wabpm_params S = W;
        S.pol = 0;
        wabpm_find_modes(&S, &n2[0], 1, 20000, 1e-11, &m[0], &nsc);
        S.pol = 1;
        wabpm_find_modes(&S, &n2[0], 1, 20000, 1e-11, &m[0], &nsv);
    }
    {
        char det[256];
        std::snprintf(det, sizeof(det),
                      "FV %.6f vs HE11 %.6f (err %+.2e, tol 4e-3); semivector err %+.2e, scalar %.6f (LP01 %.6f)",
                      nfv, he11, nfv - he11, nsv - he11, nsc, lp01);
        check_true("HE11 neff (full-vector vs exact)",
                   ok && (std::fabs(nfv - he11) < 4e-3) && (std::fabs(nfv - he11) < std::fabs(nsv - he11)),
                   det);
        const double ratio = (nfv - nsc) / (he11 - lp01);
        std::snprintf(det, sizeof(det), "(FV - scalar)/(HE11 - LP01) = %.4f (= %.5f / %.5f)",
                      ratio, nfv - nsc, he11 - lp01);
        check_true("polarization correction", (ratio > 0.93) && (ratio < 1.07), det);
    }

    // --- 4. 交差偏波成分 Ey ---
    {
        double px = 0, py = 0, eymax = 0, asym = 0;
        for (int iy = 0; iy < N; iy++) {
            for (int ix = 0; ix < N; ix++) {
                const size_t i = ix + (size_t)iy * N;
                px += std::norm(Ex[i]);
                py += std::norm(Ey[i]);
                eymax = std::max(eymax, std::abs(Ey[i]));
                // Ey(-x, y) = -Ey(x, y), Ey(x, -y) = -Ey(x, y)
                asym = std::max(asym, std::abs(Ey[i] + Ey[(N - 1 - ix) + (size_t)iy * N]));
                asym = std::max(asym, std::abs(Ey[i] + Ey[ix + (size_t)(N - 1 - iy) * N]));
            }
        }
        const double frac = py / (px + py);
        const double exact = he11_ey_fraction(fib, he11);
        char det[192];
        std::snprintf(det, sizeof(det), "Ey power fraction %.4e vs exact %.4e (%+.1f%%, tol 5%%)",
                      frac, exact, 100 * (frac / exact - 1));
        check_true("cross-polarized Ey of HE11", std::fabs(frac / exact - 1) < 0.05, det);
        std::snprintf(det, sizeof(det), "max|Ey + mirror(Ey)| / max|Ey| = %.2e", asym / eymax);
        check_true("Ey is odd in x and y (sin 2phi)", asym / eymax < 1e-6, det);
    }

    // --- 5. 縮退 (y 偏波) ---
    {
        std::vector<cplx> Fx((size_t)N * N), Fy((size_t)N * N);
        double ny = 0;
        const int oky = wabpm_fv_find_mode(&W, &n2[0], 2, 20000, 1e-11, &Fx[0], &Fy[0], &ny);
        char det[128];
        std::snprintf(det, sizeof(det), "neff_x %.9f, neff_y %.9f", nfv, ny);
        // 虚軸伝搬の ADI 分離誤差 (~2e-5, テスト 6 の注参照) は Stern の方向と ADI の
        // 分離順の組み合わせが Ex と Ey で異なるため偏波間で一致しない。直接固有値解では
        // 縮退は厳密 (差 0)。許容は分離誤差の数倍とする (偏波補正 0.049 の 0.2%)
        check_true("x / y polarization degeneracy", oky && (std::fabs(ny - nfv) < 1e-4), det);
    }

    // --- 6. 伝搬の整合 (近軸, モードを伝搬) ---
    // 位相から求めた neff は、同じ差分式を疎行列の固有値問題として直接解いた値
    // 1.1791435 と一致する。モードソルバの値はそれより +2e-5 高い (虚軸伝搬の
    // ADI 分離誤差 O(a^2)) ので、許容はその分を見込んで 5e-5 とする
    {
        wabpm_params P = W;
        P.n0 = nfv;          // 位相基準を neff に取ると、理想的には位相が進まない
        P.wideangle = 0;
        P.dz = 0.1e-6;
        const int nz = 100;  // z = 10um
        std::vector<cplx> Px = Ex, Py = Ey;
        int maxit = 0;
        for (int iz = 0; iz < nz; iz++) {
            const int r = wabpm_fv_step(&P, &Px[0], &Py[0], &n2[0]);
            if ((r < 0) || (r > maxit)) maxit = (r < 0) ? 1000 : r;
        }
        cplx ov = 0;
        double p0 = 0, p1 = 0;
        for (size_t i = 0; i < Px.size(); i++) {
            ov += std::conj(Ex[i]) * Px[i] + std::conj(Ey[i]) * Py[i];
            p0 += std::norm(Ex[i]) + std::norm(Ey[i]);
            p1 += std::norm(Px[i]) + std::norm(Py[i]);
        }
        const double eta = std::norm(ov) / (p0 * p1);
        // 位相 arg<E0, E(z)> = -(beta - k0 n0) z  ->  neff_prop = n0 - arg/(k0 z)
        const double nprop = P.n0 - std::arg(ov) / (k0 * nz * P.dz);
        char det[192];
        std::snprintf(det, sizeof(det),
                      "overlap %.6f (> 0.999), neff from phase %.7f vs mode %.7f (diff %.1e), iters <= %d",
                      eta, nprop, nfv, nprop - nfv, maxit);
        check_true("mode propagation (paraxial)",
                   (eta > 0.999) && (std::fabs(nprop - nfv) < 5e-5) && (maxit < 50), det);
    }

    // --- 7. 広角の伝搬 (中程度のコントラスト) ---
    {
        const Fiber f2 = {1.45, 1.40, 1.0e-6, k0};
        const int M = 80;
        const double d2 = 0.125e-6;
        std::vector<cplx> m2((size_t)M * M);
        for (int iy = 0; iy < M; iy++) {
            for (int ix = 0; ix < M; ix++) {
                const double x = coord(ix, M, d2), y = coord(iy, M, d2);
                m2[ix + (size_t)iy * M] = ((x * x + y * y) < f2.a * f2.a) ? f2.n1 * f2.n1 : f2.n2 * f2.n2;
            }
        }
        wabpm_params Q;
        setup(&Q, M, d2, k0, 1.42);
        Q.dz = 50 * (1.0 / (8.0 / (d2 * d2))) * 4 * k0 * Q.n0;
        std::vector<cplx> Gx((size_t)M * M), Gy((size_t)M * M);
        double nm = 0;
        const int okm = wabpm_fv_find_mode(&Q, &m2[0], 1, 20000, 1e-12, &Gx[0], &Gy[0], &nm);
        const double he = largest_root(f2, he1_eq);
        // 近軸の位相速度の理論値 : beta - k0 n0 = k0 (neff^2 - n0^2)/(2 n0)
        const double npar = f2.n2 + (nm * nm - f2.n2 * f2.n2) / (2 * f2.n2);
        double nph[2] = {0, 0}, etam = 0;
        for (int wa = 0; wa < 2; wa++) {
            wabpm_params P = Q;
            P.n0 = f2.n2;        // 参照屈折率をクラッドに取る (X = (neff^2 - n0^2)/n0^2 ~ 0.029)
            P.wideangle = wa;
            P.dz = 0.25e-6;
            const int nz = 80;   // z = 20um
            std::vector<cplx> Px = Gx, Py = Gy;
            for (int iz = 0; iz < nz; iz++) wabpm_fv_step(&P, &Px[0], &Py[0], &m2[0]);
            cplx ov = 0;
            double p0 = 0, p1 = 0;
            for (size_t i = 0; i < Px.size(); i++) {
                ov += std::conj(Gx[i]) * Px[i] + std::conj(Gy[i]) * Py[i];
                p0 += std::norm(Gx[i]) + std::norm(Gy[i]);
                p1 += std::norm(Px[i]) + std::norm(Py[i]);
            }
            // 位相の巻き戻し : 期待される進み (neff - n0) k0 z は 2 pi を超えるので、
            // 期待値の周りで arg を合わせる
            const double expect = -(nm - P.n0) * k0 * nz * P.dz;
            double ph = std::arg(ov);
            ph += 2 * M_PI * std::round((expect - ph) / (2 * M_PI));
            nph[wa] = P.n0 - ph / (k0 * nz * P.dz);
            if (wa) etam = std::norm(ov) / (p0 * p1);
        }
        char det[256];
        std::snprintf(det, sizeof(det),
                      "mode %.7f (HE11 %.7f): wide-angle %.7f (diff %+.1e, tol 1e-4), paraxial %.7f (diff %+.1e, theory %+.1e), overlap %.6f",
                      nm, he, nph[1], nph[1] - nm, nph[0], nph[0] - nm, npar - nm, etam);
        check_true("wide-angle propagation (moderate contrast)",
                   okm && (std::fabs(nph[1] - nm) < 1e-4) && (etam > 0.98)
                   && (std::fabs(nph[0] - npar) < 5e-6), det);
    }

    std::printf("\n%s : %d failure(s)\n", g_failures ? "FAILED" : "PASSED", g_failures);
    return g_failures ? 1 : 0;
}
