#pragma once
#include <Eigen/Dense>

struct DirectionalWidth {
    Eigen::MatrixXf dirs;
    DirectionalWidth(int n_dirs = 30);
    float W(Eigen::MatrixXf &bounds) const;
    int n_dirs;
    Eigen::MatrixXf eval(const Eigen::MatrixXf &polygons) const;
};