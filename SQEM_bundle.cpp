#include "SQEM_bundle.h"
using namespace std;
using namespace Eigen;
double SQEMBundle::minimize(Eigen::Vector3f &center, float &radius, const Eigen::Vector3f &pa, const Eigen::Vector3f &pb, float radius_bound)
{
    // constructs a lower bound for sum of SQEM under all possible deformations 
    double sum = 0.0;
    // naive imple
    for (int i = 0; i < sqem.size(); i++){
        auto &s {sqem[i]};

        s.minimize(this->center[i], this->radius[i], pa, pb, radius_bound);
        eval[i] = s.evaluate(this->center[i], this->radius[i]);
        sum += eval[i];
    }
    return sum;
}

double SQEMBundle::evaluate(Eigen::Vector3f &center, float r)
{
    double sum = 0;
    for (auto f: eval) {
        sum += f;
    }
    return sum;
}


SQEMBundle BundledSphereMesh::Qbundle(int u) const {
    SQEMBundle q(n_deforms);
    for (int i = 0; i < n_deforms; i++){
        Vector3f Vu{V_deform[i].row(u)};
        assert(valid(u));
        for (int j = NI(u); j < NI(u + 1); j ++) if (Fvalid(VF(j))){
            int fu = VF(j);
            Vector3f nu {N_deform[i].row(fu)};
            Vector3f _Vu = Vu + nu * rest_shape.R(u);
            q.sqem[i] += SQEM(_Vu, nu) * (this->area(F.row(fu), i) / 3.0);
        }
    }
    return q;
}

float BundledSphereMesh::area(const Vector3i &f, int i) const {
    return SphereMeshBase::area(V_deform[i], f);
}