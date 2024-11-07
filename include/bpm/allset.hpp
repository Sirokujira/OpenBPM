#include <iostream>
#include <functional>
#include <vector>
#include <Eigen/Dense>
#include <cmath>
#include <complex>

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
        MatrixXf field;
        float Lx;
        float Ly;
        string xSymmetry;
        string ySymmetry;
    } E;
    struct Mode_Struct {
        MatrixXf field;
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

    P.E.field = E / sqrt((E.array().abs2().sum()) / powerFraction);

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
    P.E.field = P.E.field.cwiseProduct((-1i * k * (cos(direction * M_PI / 180) * X.array() + sin(direction * M_PI / 180) * Y.array())).exp());

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

    MatrixXf X_offset = X.array() + distance * cos(direction * M_PI / 180);
    MatrixXf Y_offset = Y.array() + distance * sin(direction * M_PI / 180);

    // Interpolate the field with offset
    P.E.field = MatrixXf::Zero(X.rows(), X.cols());
    for (int i = 0; i < X.rows(); ++i) {
        for (int j = 0; j < X.cols(); ++j) {
            float x = X_offset(i, j);
            float y = Y_offset(i, j);

            // Bilinear interpolation (assuming field values outside the grid are zero)
            int x1 = floor(x);
            int x2 = ceil(x);
            int y1 = floor(y);
            int y2 = ceil(y);

            if (x1 >= 0 && x2 < X.cols() && y1 >= 0 && y2 < Y.rows()) {
                float Q11 = P.E.field(y1, x1);
                float Q12 = P.E.field(y2, x1);
                float Q21 = P.E.field(y1, x2);
                float Q22 = P.E.field(y2, x2);

                float R1 = ((x2 - x) / (x2 - x1)) * Q11 + ((x - x1) / (x2 - x1)) * Q21;
                float R2 = ((x2 - x) / (x2 - x1)) * Q12 + ((x - x1) / (x2 - x1)) * Q22;

                P.E.field(i, j) = ((y2 - y) / (y2 - y1)) * R1 + ((y - y1) / (y2 - y1)) * R2;
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
    result.field = MatrixXf::Zero(P.modes[0].field.rows(), P.modes[0].field.cols());
    result.xSymmetry = P.modes[modeIdxs[0]].xSymmetry;
    result.ySymmetry = P.modes[modeIdxs[0]].ySymmetry;

    for (size_t i = 0; i < modeIdxs.size(); ++i) {
        result.field += coeffs[i] * P.modes[modeIdxs[i]].field;
    }

    P.E = result;
    return P;
}

P_Struct findModes(P_Struct P, int nModes, bool singleCoreModes = false, bool sortByLoss = false, bool plotModes = true) {
    // Placeholder implementation to find modes
    cout << "Finding modes..." << endl;
    // For now, we'll simulate finding some modes
    for (int i = 0; i < nModes; ++i) {
        P_Struct::Mode_Struct mode;
        mode.Lx = P.n.Lx;
        mode.Ly = P.n.Ly;
        mode.xSymmetry = P.xSymmetry;
        mode.ySymmetry = P.ySymmetry;
        mode.field = MatrixXf::Random(P.n.n.rows(), P.n.n.cols());
        mode.neff = complex<float>(1.5f, 0.01f * i);
        mode.label = "Mode " + to_string(i + 1);
        P.modes.push_back(mode);
    }
    cout << "Found " << nModes << " modes." << endl;
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



