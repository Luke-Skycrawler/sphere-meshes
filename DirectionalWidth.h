#pragma once
#include <Eigen/Dense>

struct DirectionalWidth {
    Eigen::MatrixXf dirs;
    DirectionalWidth(int n_dirs = 30);
    float W(Eigen::MatrixXf &P) const;
    int n_dirs;
};