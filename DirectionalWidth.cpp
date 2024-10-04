#define _USE_MATH_DEFINES
#include "DirectionalWidth.h"
#include <cmath>
#include <limits>
using namespace std;
using namespace Eigen;

Eigen::MatrixXf fibonacci_lattice(int N) {
    Eigen::MatrixXf directions(N, 3);  // N rows, 3 columns (for x, y, z)
    const float golden_ratio = (1 + std::sqrt(5)) / 2;
    
    for (int i = 0; i < N; ++i) {
        float theta = 2.0 * M_PI * i / golden_ratio;  // Azimuth angle
        float phi = std::acos(1 - 2.0 * (i + 0.5) / N);  // Polar angle
        
        // Convert spherical coordinates to Cartesian coordinates (x, y, z)
        float x = std::sin(phi) * std::cos(theta);
        float y = std::sin(phi) * std::sin(theta);
        float z = std::cos(phi);
        
        directions(i, 0) = x;
        directions(i, 1) = y;
        directions(i, 2) = z;
    }

    return directions;
}
DirectionalWidth::DirectionalWidth(int n_dirs): n_dirs(n_dirs)
{
    dirs.resize(n_dirs + 3, 3);
    dirs.block<3, 3>(0, 0) = Eigen::Matrix3f::Identity();
    dirs.block(3, 0, n_dirs, 3) = fibonacci_lattice(n_dirs);
}

float DirectionalWidth::W(Eigen::MatrixXf &bounds) const {
    auto b = bounds.col(1) - bounds.col(0);
    return b.minCoeff();
} 

MatrixXf DirectionalWidth::eval(const MatrixXf &P) const {
    MatrixXf ret(n_dirs + 3, 2);
    for (int j = 0; j < n_dirs + 3; j++) {
        float min = numeric_limits<float>::max(), max = numeric_limits<float>::lowest();
        for (int i = 0; i < P.rows(); i ++) {
            Vector3f p = P.row(i);
            float dot = dirs.row(j).dot(p);
            if (dot < min) min = dot;
            if (dot > max) max = dot;
        }
        ret(j, 0) = min;
        ret(j, 1) = max;
    }
    return ret;
}