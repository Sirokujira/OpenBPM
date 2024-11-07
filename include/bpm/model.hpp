#ifndef BPM_MODEL_HPP
#define BPM_MODEL_HPP

#include <string>
#include <vector>
#include <ctime>
#include <iostream>
#include <stdexcept>
#include <cmath>
#include <algorithm>

class BPMModel {
public:
    // Constructor
    BPMModel();

    // Getters for dimensions
    int getNx() const;
    int getNy() const;
    int getNz() const;
    double getLx() const;
    double getLy() const;
    double getDx() const;
    double getDy() const;
    double getDz() const;

    // Getters for grid points
    std::vector<double> getX() const;
    std::vector<double> getY() const;

    // Other member functions
    void setNCladding();
    void setShapes();
    void setDisplayScaling();
    void finalizeVideo();

    // Public member variables
    std::string name;

private:
    // Private member functions
    std::vector<double> getGridArray(int size, double spacing, const std::string& symmetry) const;

    // Private member variables
    double Lx_main = 10.0;       // Example value, can be set as needed
    double Ly_main = 10.0;       // Example value, can be set as needed
    double Lz = 5.0;             // Example value, can be set as needed
    double dz_target = 0.1;      // Target spacing in z-direction
    double padfactor = 1.2;      // Padding factor
    int Nx_main = 100;           // Example value, can be set as needed
    int Ny_main = 100;           // Example value, can be set as needed
    int updates = 10;            // Example value, can be set as needed
    std::string xSymmetry = "NoSymmetry";  // Symmetry setting for x-direction
    std::string ySymmetry = "NoSymmetry";  // Symmetry setting for y-direction
};

#endif // BPM_MODEL_HPP



