#include <iostream>
#include <vector>
#include <complex>
#include <stdexcept>
#include <string>
#include <functional>
#include <algorithm>
#include <Eigen/Dense>

#include "define.h"

using namespace Eigen;

class RefractiveIndexProfile {
public:
    // Constructor
    RefractiveIndexProfile(double Lx = 1.0, double Ly = 1.0)
        : Lx(Lx), Ly(Ly), xSymmetry(Symmetry::NoSymmetry), ySymmetry(Symmetry::NoSymmetry), label(""), neff(NAN) {
        if (Lx <= 0 || Ly <= 0) {
            throw std::invalid_argument("Lx and Ly must be positive.");
        }
    }

    // Enums to represent symmetry types
    enum class Symmetry : uint8_t {
        NoSymmetry = 0,
        Symmetric = 1,
        AntiSymmetric = 2
    };

    // Enums to represent colormap types
    enum class Colormap {
        GPBGYR,
        HSV,
        Parula,
        Gray,
        Cividis
    };

    // Function to initialize refractive index from a function
    void initializeRIfromFunction(const std::function<std::complex<float>(float, float, float, const std::vector<std::complex<float>>&)>& hFunc, const std::vector<std::complex<float>>& nParameters = {}, int Nz = 1) {
        if (Nz <= 0) {
            throw std::invalid_argument("Nz must be a positive integer.");
        }

        // Generate the grid
        if (Nz > 1) {
            float dz_n = Lz / (Nz - 1);
            for (int k = 0; k < Nz; ++k) {
                float z = k * dz_n;
                for (size_t i = 0; i < x.size(); ++i) {
                    for (size_t j = 0; j < y.size(); ++j) {
                        field[i][j] = hFunc(x[i], y[j], z, nParameters);
                    }
                }
            }
        } else { // 2D case
            for (size_t i = 0; i < x.size(); ++i) {
                for (size_t j = 0; j < y.size(); ++j) {
                    field[i][j] = hFunc(x[i], y[j], 0.0f, nParameters);
                }
            }
        }
    }

	// 任意の型のベクターを処理するテンプレート関数
	template <typename T>
	std::vector<T> reverseAndNegate(const std::vector<T>& vec) {
	    std::vector<T> result;
	    if (vec.empty()) {
	        return result;
	    }

	    result.reserve(vec.size() - 1);

	    // vecの要素を逆順に走査し、反転してresultに追加
	    std::transform(vec.rbegin() + 1, vec.rend(), std::back_inserter(result), [](const T& val) {
	        return -val;
	    });

	    return result;
	}

	// 使い方
	//x.insert(x.begin(), reverseAndNegate(x).begin(), reverseAndNegate(x).end());
	
    // Function to calculate the full refractive index based on symmetry
    void calcFullRI() {
        if (x.empty() || y.empty() || field.empty()) {
            throw std::runtime_error("x, y, or field data is missing.");
        }

        // Handle xSymmetry
        switch (xSymmetry) {
            case Symmetry::NoSymmetry:
                // Do nothing
                break;
            case Symmetry::AntiSymmetric:
                x.insert(x.begin(), reverseAndNegate(x).begin(), reverseAndNegate(x).end());
                for (auto& row : field) {
                    row.insert(row.begin(), reverseAndNegate(y).begin(), reverseAndNegate(y).end());
                }
                break;
            case Symmetry::Symmetric:
                x.insert(x.begin(), std::vector<double>(x.rbegin(), x.rend()).begin(), std::vector<double>(x.rbegin(), x.rend()).end());
                for (auto& row : field) {
                    //row.insert(row.begin(), std::vector<std::complex<float>>(row.rbegin(), row.rend()));
                }
                break;
        }

        // Handle ySymmetry
        switch (ySymmetry) {
            case Symmetry::NoSymmetry:
                // Do nothing
                break;
            case Symmetry::AntiSymmetric:
                y.insert(y.begin(), reverseAndNegate(y).begin(), reverseAndNegate(y).end());
                field.insert(field.begin(), std::vector<std::complex<float>>(field[0].size(), {0, 0}));
                for (size_t i = 1; i < field.size(); ++i) {
                    for (auto& element : field[i]) {
                        element = -element;
                    }
                }
                break;
            case Symmetry::Symmetric:
                y.insert(y.begin(), std::vector<double>(y.rbegin(), y.rend()).begin(), std::vector<double>(y.rbegin(), y.rend()).end());
                field.insert(field.begin(), std::vector<std::complex<float>>(field[0].size(), {0, 0}));
                break;
        }
    }

    // Function to calculate the full electric field based on symmetry
    void calcFullField() {
        if (x.empty() || y.empty() || field.empty()) {
            throw std::runtime_error("x, y, or field data is missing.");
        }

        // Handle ySymmetry
        switch (ySymmetry) {
            case Symmetry::NoSymmetry:
                // Do nothing
                break;
            case Symmetry::AntiSymmetric:
                x.insert(x.begin(), reverseAndNegate(x).begin(), reverseAndNegate(x).end());
                for (auto& row : field) {
                    row.insert(row.begin(), reverseAndNegate(row).begin(), reverseAndNegate(row).end());
                }
                break;
            case Symmetry::Symmetric:
                x.insert(x.begin(), std::vector<double>(x.rbegin(), x.rend()).begin(), std::vector<double>(x.rbegin(), x.rend()).end());
                for (auto& row : field) {
                    //row.insert(row.begin(), std::vector<std::complex<float>>(row.rbegin(), row.rend()));
                }
                break;
        }

        // Handle xSymmetry
        switch (xSymmetry) {
            case Symmetry::NoSymmetry:
                // Do nothing
                break;
            case Symmetry::AntiSymmetric:
                y.insert(y.begin(), reverseAndNegate(y).begin(), reverseAndNegate(y).end());
                field.insert(field.begin(), std::vector<std::complex<float>>(field[0].size(), {0, 0}));
                for (size_t i = 1; i < field.size(); ++i) {
                    for (auto& element : field[i]) {
                        element = -element;
                    }
                }
                break;
            case Symmetry::Symmetric:
                y.insert(y.begin(), std::vector<double>(y.rbegin(), y.rend()).begin(), std::vector<double>(y.rbegin(), y.rend()).end());
                field.insert(field.begin(), std::vector<std::complex<float>>(field[0].size(), {0, 0}));
                break;
        }
    }

    // Getters and setters
    double getLx() const { return Lx; }
    void setLx(double value) {
        if (value <= 0) {
            throw std::invalid_argument("Lx must be positive.");
        }
        Lx = value;
    }

    double getLy() const { return Ly; }
    void setLy(double value) {
        if (value <= 0) {
            throw std::invalid_argument("Ly must be positive.");
        }
        Ly = value;
    }

    const std::vector<std::vector<std::complex<float>>>& getField() const { return field; }
    void setField(const std::vector<std::vector<std::complex<float>>>& value) {
        field = value;
    }

    Symmetry getXSymmetry() const { return xSymmetry; }
    void setXSymmetry(Symmetry value) { xSymmetry = value; }

    Symmetry getYSymmetry() const { return ySymmetry; }
    void setYSymmetry(Symmetry value) { ySymmetry = value; }

    std::string getLabel() const { return label; }
    void setLabel(const std::string& value) { label = value; }

    double getNeff() const { return neff; }
    void setNeff(double value) { neff = value; }

private:
    double Lx;
    double Ly;
    double Lz;
    std::vector<double> x;
    std::vector<double> y;
    std::vector<std::vector<std::complex<float>>> field; // Usually complex
    Symmetry xSymmetry;
    Symmetry ySymmetry;
    std::string label;
    double neff;
};

//int main() {
//    try {
//        RefractiveIndexProfile profile;
//        profile.setLx(2.0);
//        profile.setLy(3.0);
//        profile.setLabel("Test Label");
//    } catch (const std::exception& e) {
//        std::cerr << e.what() << std::endl;
//    }
//
//    return 0;
//}


