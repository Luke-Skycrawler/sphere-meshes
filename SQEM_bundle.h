#pragma once
#include <vector>
#include <Eigen/Dense>
#include "SQEM.h"


struct SQEMBundle{
    std::vector<SQEM> sqem;
    double minimize(Eigen::Vector3f &center, float &r, const Eigen::Vector3f &pa, const Eigen::Vector3f &pb) const;
    double evaluate(Eigen::Vector3f &center, float r);
    
};