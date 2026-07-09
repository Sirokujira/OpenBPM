#include <iostream>
#include <functional>
#include <vector>
#include <Eigen/Dense>
#include <cmath>
#include <complex>

#include "bpm/wabpm.h"

using namespace std;
using namespace Eigen;

typedef function<MatrixXf(const MatrixXf&, const MatrixXf&, const vector<float>&)> Func2D_E;

struct P_Struct {
    vector<float> x;
    vector<float> y;
    float n_background;
    float Lz;
    float dx;
    float dy;
    string xSymmetry;
    string ySymmetry;
    float n_0;
    float lambda;
    struct N_Struct {
        MatrixXf n;
        float Lx;
        float Ly;
    } n;
    struct E_Struct {
        MatrixXcf field;   // 複素電界 (位相チルト等を保持できるよう複素型)
        float Lx;
        float Ly;
        string xSymmetry;
        string ySymmetry;
    } E;
    struct Mode_Struct {
        MatrixXcf field;
        float Lx;
        float Ly;
        string xSymmetry;
        string ySymmetry;
        complex<float> neff;
        string label;
    };
    vector<Mode_Struct> modes;
};

P_Struct initializeEfromFunction(P_Struct P, Func2D_E hFunc, vector<float> Eparameters = {}) {
    if (Eparameters.empty()) {
        Eparameters = {};
    }

    MatrixXf X(P.x.size(), P.y.size());
    MatrixXf Y(P.x.size(), P.y.size());

    for (size_t i = 0; i < P.x.size(); ++i) {
        for (size_t j = 0; j < P.y.size(); ++j) {
            X(i, j) = P.x[i];
            Y(i, j) = P.y[j];
        }
    }

    MatrixXf E = hFunc(X, Y, Eparameters); // Call function to initialize E field

    float powerFraction = 1.0f / (1 + (P.xSymmetry != "none")) / (1 + (P.ySymmetry != "none")); // How large a fraction of the total power we are simulating

    P.E.field = (E / sqrt((E.array().abs2().sum()) / powerFraction)).cast<complex<float>>();

    P.E.Lx = P.n.Lx;
    P.E.Ly = P.n.Ly;
    P.E.xSymmetry = P.xSymmetry;
    P.E.ySymmetry = P.ySymmetry;

    return P;
}

P_Struct tiltField(P_Struct P, float direction, float angle) {
    if (typeid(P.E.field) == typeid(function<void()>)) {
        throw invalid_argument("Error: BPM-Matlab currently doesn't support field tilts in segments where the field is specified as a function handle. You can avoid this by adding a short initial propagation segment and applying the field tilt between the first and second segments.");
    }

    MatrixXf X(P.x.size(), P.y.size());
    MatrixXf Y(P.x.size(), P.y.size());

    for (size_t i = 0; i < P.x.size(); ++i) {
        for (size_t j = 0; j < P.y.size(); ++j) {
            X(i, j) = P.x[i];
            Y(i, j) = P.y[j];
        }
    }

    float k = P.n_0 * 2 * M_PI / P.lambda * angle / 180 * M_PI; // Refractive index in wedge in between fiber segments assumed to be uniformly P.n_0
    ArrayXXf phi = cos(direction * M_PI / 180) * X.array() + sin(direction * M_PI / 180) * Y.array();
    P.E.field = (P.E.field.array() * (complex<float>(0.0f, -k) * phi.cast<complex<float>>()).exp()).matrix();

    return P;
}

P_Struct offsetField(P_Struct P, float direction, float distance) {
    if (typeid(P.E.field) == typeid(function<void()>)) {
        throw invalid_argument("Error: BPM-Matlab currently doesn't support field offsets in segments where the field is specified as a function handle. You can avoid this by adding a short initial propagation segment and applying the field offset between the first and second segments.");
    }

    MatrixXf X(P.x.size(), P.y.size());
    MatrixXf Y(P.x.size(), P.y.size());

    for (size_t i = 0; i < P.x.size(); ++i) {
        for (size_t j = 0; j < P.y.size(); ++j) {
            X(i, j) = P.x[i];
            Y(i, j) = P.y[j];
        }
    }

    // 界を +distance 方向へ平行移動する : 移動先の各格子点で
    // 元の界を (x - ox, y - oy) の位置からインデックス空間の双一次補間で取得する
    // (格子外は 0 とみなす)
    const float ox = distance * cos(direction * M_PI / 180);
    const float oy = distance * sin(direction * M_PI / 180);
    const int Nx = (int)P.x.size();
    const int Ny = (int)P.y.size();
    const float dxs = (Nx > 1) ? (P.x[1] - P.x[0]) : P.dx;
    const float dys = (Ny > 1) ? (P.y[1] - P.y[0]) : P.dy;

    const MatrixXcf src = P.E.field;
    P.E.field = MatrixXcf::Zero(Nx, Ny);
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            const float fi = (P.x[i] - ox - P.x[0]) / dxs;   // 元格子上の実数インデックス
            const float fj = (P.y[j] - oy - P.y[0]) / dys;
            const int i0 = (int)floor(fi);
            const int j0 = (int)floor(fj);
            if (i0 >= 0 && i0 + 1 < Nx && j0 >= 0 && j0 + 1 < Ny) {
                const float tx = fi - i0;
                const float ty = fj - j0;
                P.E.field(i, j) = (1 - tx) * (1 - ty) * src(i0,     j0)
                                + (    tx) * (1 - ty) * src(i0 + 1, j0)
                                + (1 - tx) * (    ty) * src(i0,     j0 + 1)
                                + (    tx) * (    ty) * src(i0 + 1, j0 + 1);
            }
        }
    }

    return P;
}

P_Struct modeSuperposition(P_Struct P, vector<int> modeIdxs, vector<float> coeffs = {}) {
    if (coeffs.empty()) {
        coeffs = vector<float>(modeIdxs.size(), 1.0f);
    }

    P_Struct::E_Struct result;
    result.Lx = P.modes[modeIdxs[0]].Lx;
    result.Ly = P.modes[modeIdxs[0]].Ly;
    result.field = MatrixXcf::Zero(P.modes[0].field.rows(), P.modes[0].field.cols());
    result.xSymmetry = P.modes[modeIdxs[0]].xSymmetry;
    result.ySymmetry = P.modes[modeIdxs[0]].ySymmetry;

    for (size_t i = 0; i < modeIdxs.size(); ++i) {
        result.field += complex<float>(coeffs[i], 0.0f) * P.modes[modeIdxs[i]].field;
    }

    P.E = result;
    return P;
}

// 導波モードを実効屈折率の降順で求める (虚軸伝搬法, bpm/modes.cpp)
// P.n.n は field(ix, iy) 配置 (行 = x, 列 = y) の実屈折率分布。
// 使用時は bpm/modes.cpp と bpm/wabpm.cpp をリンクすること。
// singleCoreModes / sortByLoss / plotModes は未サポート (無視される)。
P_Struct findModes(P_Struct P, int nModes, bool singleCoreModes = false, bool sortByLoss = false, bool plotModes = true) {
    (void)singleCoreModes; (void)sortByLoss; (void)plotModes;
    cout << "Finding modes..." << endl;

    const int Nx = (int)P.n.n.rows();
    const int Ny = (int)P.n.n.cols();
    const long N = (long)Nx * Ny;

    wabpm_params W;
    W.Nx = Nx;
    W.Ny = Ny;
    W.dx = P.dx;
    W.dy = P.dy;
    W.k0 = 2.0 * M_PI / P.lambda;
    W.n0 = P.n_0;
    W.wideangle = 0;
    W.pol = 0;

    // 屈折率分布 (実数) -> 複素比誘電率
    vector<complex<double>> n2((size_t)N);
    double nmax = P.n_0;
    for (int iy = 0; iy < Ny; ++iy) {
        for (int ix = 0; ix < Nx; ++ix) {
            const double n = P.n.n(ix, iy);
            n2[ix + (long)iy * Nx] = complex<double>(n * n, 0.0);
            if (n > nmax) nmax = n;
        }
    }

    // 虚軸ステップ幅 : 最大固有値 mu_max = k0^2(nmax^2 - n0^2) に対し
    // a*mu_max ~ 0.5 となるよう選ぶ (a = dz/(4 k0 n0))
    const double mu_max = W.k0 * W.k0 * (nmax * nmax - W.n0 * W.n0);
    W.dz = (mu_max > 0) ? (2.0 * W.k0 * W.n0 / mu_max) : (10 * P.dx);

    vector<complex<double>> modes((size_t)nModes * N);
    vector<double> neff((size_t)nModes);
    const int nFound = wabpm_find_modes(&W, n2.data(), nModes, 20000, 1e-12,
                                        modes.data(), neff.data());

    for (int m = 0; m < nFound; ++m) {
        P_Struct::Mode_Struct mode;
        mode.Lx = P.n.Lx;
        mode.Ly = P.n.Ly;
        mode.xSymmetry = P.xSymmetry;
        mode.ySymmetry = P.ySymmetry;
        mode.field = MatrixXcf(Nx, Ny);
        for (int iy = 0; iy < Ny; ++iy) {
            for (int ix = 0; ix < Nx; ++ix) {
                const complex<double> v = modes[(size_t)m * N + ix + (long)iy * Nx];
                mode.field(ix, iy) = complex<float>((float)v.real(), (float)v.imag());
            }
        }
        mode.neff = complex<float>((float)neff[m], 0.0f);
        mode.label = "Mode " + to_string(m + 1);
        P.modes.push_back(mode);
    }
    cout << "Found " << nFound << " modes." << endl;
    return P;
}

int getLabeledModeIndex(const P_Struct& P, const string& label) {
    for (size_t i = 0; i < P.modes.size(); ++i) {
        if (P.modes[i].label.find(label) != string::npos) {
            return static_cast<int>(i);
        }
    }
    throw invalid_argument("Mode " + label + " not found");
}

//int main() {
//    // Example usage of the initializeEfromFunction function
//    P_Struct P;
//    P.x = {0.0f, 1.0f, 2.0f};
//    P.y = {0.0f, 1.0f, 2.0f};
//    P.n_background = 1.0f;
//    P.Lz = 10.0f;
//    P.dx = 0.1f;
//    P.dy = 0.1f;
//    P.xSymmetry = "none";
//    P.ySymmetry = "none";
//    P.n_0 = 1.5f;
//    P.lambda = 0.633f;
//
//    vector<float> Eparameters = {1.5f};
//
//    try {
//        P = initializeEfromFunction(P, [](const MatrixXf& X, const MatrixXf& Y, const vector<float>& Eparameters) -> MatrixXf {
//            MatrixXf E = X + Y;
//            return E;
//        }, Eparameters);
//
//        cout << "P.E.field: \n" << P.E.field << endl;
//
//        // Apply tilt to the field
//        P = tiltField(P, 45.0f, 30.0f);
//
//        cout << "P.E.field after tilt: \n" << P.E.field << endl;
//
//        // Example usage of modeSuperposition
//        P_Struct::Mode_Struct mode1 = {MatrixXf::Ones(3, 3), 1.0f, 1.0f, "none", "none"};
//        P_Struct::Mode_Struct mode2 = {MatrixXf::Constant(3, 3, 2.0f), 1.0f, 1.0f, "none", "none"};
//        P.modes.push_back(mode1);
//        P.modes.push_back(mode2);
//
//        vector<int> modeIdxs = {0, 1};
//        vector<float> coeffs = {0.5f, 1.0f};
//
//        P = modeSuperposition(P, modeIdxs, coeffs);
//
//        cout << "P.E.field after mode superposition: \n" << P.E.field << endl;
//
//        // Example usage of findModes
//        P = findModes(P, 3);
//        for (const auto& mode : P.modes) {
//            cout << mode.label << " with neff = " << mode.neff << endl;
//        }
//
//        // Example usage of getLabeledModeIndex
//        int idx = getLabeledModeIndex(P, "Mode 2");
//        cout << "Index of 'Mode 2': " << idx << endl;
//
//        // Example usage of offsetField
//        P = offsetField(P, 90.0f, 1.0f);
//        cout << "P.E.field after offset: \n" << P.E.field << endl;
//    } catch (const invalid_argument& e) {
//        cerr << "Error: " << e.what() << endl;
//    }
//
//    return 0;
//}



