#pragma once

#include <cmath>
#include <string>
#include <vector>

struct WeightSystematic {
    std::string suffix;
    double weight;
};

struct MetSystematic {
    std::string suffix;
    double px;
    double py;
};

inline void cartesianToPolar(double px, double py, double& met, double& metphi) {
    met = std::hypot(px, py);
    metphi = std::atan2(py, px);
}
