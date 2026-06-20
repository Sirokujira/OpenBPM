/*
 * test_bpmmodel.cpp
 *
 * BPMModel クラス (bpm/model.cpp) の単体テスト
 *
 * テスト一覧:
 *   1. test_dx_dy          – dx = Lx_main/Nx_main, dy = Ly_main/Ny_main
 *   2. test_nx_ny          – Nx = round(padfactor * Lx_main / dx) = 120
 *   3. test_nz             – Nz = max(updates, round(Lz/dz_target)) = 50
 *   4. test_lx_ly          – Lx = dx * Nx = 12.0
 *   5. test_dz             – dz = Lz / Nz = 0.1
 *   6. test_grid_x_length  – getX() の長さが Nx と一致し、間隔が dx
 *   7. test_grid_y_length  – getY() の長さが Ny と一致し、間隔が dy
 *   8. test_name_not_empty – コンストラクタでタイムスタンプ入り name がセットされる
 *
 * デフォルトパラメータ (model.hpp private 初期値):
 *   Lx_main=10, Ly_main=10, Lz=5, dz_target=0.1
 *   padfactor=1.2, Nx_main=100, Ny_main=100, updates=10
 *   xSymmetry=NoSymmetry, ySymmetry=NoSymmetry
 */

#include <cmath>
#include <cstdio>

#include "bpm/model.hpp"

// -----------------------------------------------------------------------
// 軽量テストマクロ
// -----------------------------------------------------------------------

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

#define CHECK_EQ(a, b, msg)                                                    \
    do {                                                                       \
        long _a = (long)(a), _b = (long)(b);                                   \
        if (_a != _b) {                                                        \
            fprintf(stderr,                                                    \
                    "FAIL [%s:%d] %s: %ld != %ld\n",                           \
                    __FILE__, __LINE__, msg, _a, _b);                          \
            return 1;                                                          \
        }                                                                      \
    } while (0)

#define CHECK(cond, msg)                                                       \
    do {                                                                       \
        if (!(cond)) {                                                         \
            fprintf(stderr, "FAIL [%s:%d] %s\n", __FILE__, __LINE__, msg);    \
            return 1;                                                          \
        }                                                                      \
    } while (0)

// -----------------------------------------------------------------------
// Test 1: dx / dy
//
// dx = Lx_main / Nx_main = 10.0 / 100 = 0.1
// -----------------------------------------------------------------------
static int test_dx_dy()
{
    BPMModel m;
    CHECK_NEAR(m.getDx(), 0.1, 1e-12, "dx = Lx_main/Nx_main = 0.1");
    CHECK_NEAR(m.getDy(), 0.1, 1e-12, "dy = Ly_main/Ny_main = 0.1");
    printf("PASS test_dx_dy     dx=%.3f  dy=%.3f\n", m.getDx(), m.getDy());
    return 0;
}

// -----------------------------------------------------------------------
// Test 2: Nx / Ny
//
// Nx = round(padfactor * Lx_main / dx)
//     = round(1.2 * 10.0 / 0.1) = round(120) = 120
// -----------------------------------------------------------------------
static int test_nx_ny()
{
    BPMModel m;
    CHECK_EQ(m.getNx(), 120, "Nx = round(padfactor * Lx_main / dx) = 120");
    CHECK_EQ(m.getNy(), 120, "Ny = round(padfactor * Ly_main / dy) = 120");
    printf("PASS test_nx_ny     Nx=%d  Ny=%d\n", m.getNx(), m.getNy());
    return 0;
}

// -----------------------------------------------------------------------
// Test 3: Nz
//
// Nz = max(updates, round(Lz / dz_target))
//     = max(10, round(5.0 / 0.1)) = max(10, 50) = 50
// -----------------------------------------------------------------------
static int test_nz()
{
    BPMModel m;
    CHECK_EQ(m.getNz(), 50, "Nz = max(updates, round(Lz/dz_target)) = 50");
    printf("PASS test_nz        Nz=%d\n", m.getNz());
    return 0;
}

// -----------------------------------------------------------------------
// Test 4: Lx / Ly (実際の計算領域サイズ)
//
// Lx = getDx() * getNx() = 0.1 * 120 = 12.0
// これは padfactor=1.2 によって Lx_main=10 の 1.2 倍になっている。
// -----------------------------------------------------------------------
static int test_lx_ly()
{
    BPMModel m;
    CHECK_NEAR(m.getLx(), 12.0, 1e-10, "Lx = dx * Nx = 12.0");
    CHECK_NEAR(m.getLy(), 12.0, 1e-10, "Ly = dy * Ny = 12.0");
    printf("PASS test_lx_ly     Lx=%.4f  Ly=%.4f\n", m.getLx(), m.getLy());
    return 0;
}

// -----------------------------------------------------------------------
// Test 5: dz
//
// dz = Lz / Nz = 5.0 / 50 = 0.1
// -----------------------------------------------------------------------
static int test_dz()
{
    BPMModel m;
    CHECK_NEAR(m.getDz(), 0.1, 1e-12, "dz = Lz/Nz = 0.1");
    printf("PASS test_dz        dz=%.6f\n", m.getDz());
    return 0;
}

// -----------------------------------------------------------------------
// Test 6: getX() の長さと間隔
//
// getX() は長さ Nx の配列を返し、
// 各要素は x[i] = i * dx で等間隔。
// -----------------------------------------------------------------------
static int test_grid_x_length()
{
    BPMModel m;
    auto x = m.getX();
    CHECK_EQ((long)x.size(), (long)m.getNx(), "getX().size() == getNx()");
    CHECK_NEAR(x[0], 0.0,       1e-12, "x[0] = 0");
    CHECK_NEAR(x[1], m.getDx(), 1e-12, "x[1] = dx");
    // 末尾要素が (Nx-1)*dx であることも確認
    CHECK_NEAR(x.back(), (m.getNx() - 1) * m.getDx(), 1e-10,
               "x[Nx-1] = (Nx-1)*dx");
    printf("PASS test_grid_x_length  len=%zu  x[0]=%.4f  x[Nx-1]=%.4f\n",
           x.size(), x[0], x.back());
    return 0;
}

// -----------------------------------------------------------------------
// Test 7: getY() の長さと間隔
//
// getY() は長さ Ny の配列を返し、y[i] = i * dy。
// -----------------------------------------------------------------------
static int test_grid_y_length()
{
    BPMModel m;
    auto y = m.getY();
    CHECK_EQ((long)y.size(), (long)m.getNy(), "getY().size() == getNy()");
    CHECK_NEAR(y[0], 0.0,       1e-12, "y[0] = 0");
    CHECK_NEAR(y[1], m.getDy(), 1e-12, "y[1] = dy");
    CHECK_NEAR(y.back(), (m.getNy() - 1) * m.getDy(), 1e-10,
               "y[Ny-1] = (Ny-1)*dy");
    printf("PASS test_grid_y_length  len=%zu  y[0]=%.4f  y[Ny-1]=%.4f\n",
           y.size(), y[0], y.back());
    return 0;
}

// -----------------------------------------------------------------------
// Test 8: コンストラクタで name がセットされる
//
// BPMModel のコンストラクタはタイムスタンプ入りの name をセットする。
// 空でないことを確認する。
// -----------------------------------------------------------------------
static int test_name_not_empty()
{
    BPMModel m;
    CHECK(!m.name.empty(), "BPMModel.name は空であってはならない");
    printf("PASS test_name_not_empty  name=\"%s\"\n", m.name.c_str());
    return 0;
}

// -----------------------------------------------------------------------
// エントリポイント
// -----------------------------------------------------------------------
int main()
{
    int fail = 0;
    fail += test_dx_dy();
    fail += test_nx_ny();
    fail += test_nz();
    fail += test_lx_ly();
    fail += test_dz();
    fail += test_grid_x_length();
    fail += test_grid_y_length();
    fail += test_name_not_empty();

    if (fail == 0)
        printf("\nAll 8 tests passed.\n");
    else
        fprintf(stderr, "\n%d test(s) FAILED.\n", fail);
    return fail;
}
