#include "SQEM_bundle.h"
using namespace std;
using namespace Eigen;
double SQEMBundle::minimize(Eigen::Vector3f &center, float &r, const Eigen::Vector3f &pa, const Eigen::Vector3f &pb) const
{
    double min = std::numeric_limits<double>::max();
    for (const SQEM &s : sqem)
    {

    }
    return min;
}
double SQEMBundle::evaluate(Eigen::Vector3f &center, float r)
{
    double sum = 0.0;
    for (const SQEM &s : sqem)
    {
        Vector3f c1 {deform(c)}; 
        sum += s.evaluate(center, r);
    }
    return sum;
}