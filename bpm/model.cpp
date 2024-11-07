#include "bpm/model.hpp"

BPMModel::BPMModel() {
    time_t now = time(0);
    char timestamp[20];
    strftime(timestamp, sizeof(timestamp), "%d-%b-%Y %H.%M.%S", localtime(&now));
    name = "BPM-Matlab model " + std::string(timestamp);
}

int BPMModel::getNx() const {
    double targetLx = padfactor * Lx_main;
    int Nx = static_cast<int>(std::round(targetLx / getDx()));
    if (ySymmetry != "NoSymmetry") {
        Nx += (Nx % 2 != Nx_main % 2);  // Ensure symmetry matches Nx_main if it was set odd
    }
    return Nx;
}

int BPMModel::getNy() const {
    double targetLy = padfactor * Ly_main;
    int Ny = static_cast<int>(std::round(targetLy / getDy()));
    if (xSymmetry != "NoSymmetry") {
        Ny += (Ny % 2 != Ny_main % 2);  // Ensure symmetry matches Ny_main if it was set odd
    }
    return Ny;
}

int BPMModel::getNz() const {
    return std::max(updates, static_cast<int>(std::round(Lz / dz_target)));
}

double BPMModel::getLx() const {
    return getDx() * getNx();
}

double BPMModel::getLy() const {
    return getDy() * getNy();
}

double BPMModel::getDx() const {
    return Lx_main / Nx_main;
}

double BPMModel::getDy() const {
    return Ly_main / Ny_main;
}

double BPMModel::getDz() const {
    return Lz / getNz();
}

std::vector<double> BPMModel::getX() const {
    return getGridArray(getNx(), getDx(), ySymmetry);
}

std::vector<double> BPMModel::getY() const {
    return getGridArray(getNy(), getDy(), xSymmetry);
}

std::vector<double> BPMModel::getGridArray(int size, double spacing, const std::string& symmetry) const {
    std::vector<double> grid;
    grid.reserve(size);
    for (int i = 0; i < size; ++i) {
        grid.push_back(i * spacing);
    }
    return grid;
}

void BPMModel::setNCladding() {
    throw std::runtime_error("Error: n_cladding has been renamed n_background");
}

void BPMModel::setShapes() {
    throw std::runtime_error("Error: The shapes field has been deprecated. Use the n property to define the refractive index instead.");
}

void BPMModel::setDisplayScaling() {
    throw std::runtime_error("Error: displayScaling has been renamed plotZoom.");
}

void BPMModel::finalizeVideo() {
    // Placeholder for closing video handle
    std::cout << "Finalizing video..." << std::endl;
}

//int main() {
//    BPMModel model;
//    std::cout << "Model name: " << model.name << std::endl;
//    return 0;
//}



