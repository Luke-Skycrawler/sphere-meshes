#include "../SphereMeshes.h"
#include <iostream>
using namespace std;
using namespace Eigen;
static const float tol = 1e-10f;
Vector3f compute_normal(const Vector3f &v0, const Vector3f &v1, const Vector3f &v2, float r0, float r1, float r2, const Vector3f &n0) {
   
    auto r01 = v1 - v0;
    auto r12 = v2 - v1;

    Vector3f b{r0 - r1, r1 - r2, 1.0f};
    Matrix3f A;
    A.row(0) = r01;
    A.row(1) = r12;
    const auto residue = [&](const Vector3f &n) -> Vector3f {
        A.row(2) = n;
        return A * n - b;
    };
    // f(n) = (r01 ^T, r12^T, n^T) n - (r0 - r1, r1 - r2, 1)^T
    // solve for f(n) = 0
    Vector3f n = n0;
    auto res = residue(n);
    while (res.squaredNorm() > tol) {
        A.row(2) = 2 * n;
        auto dn = A.inverse() * res;
        n -= dn;
        res = residue(n);
    }
    assert(n.dot(n0) > 0.0f);
    return n.normalized();
}
int main() {
    // SphereMesh mesh("bar.obj");
    // mesh.compute_normal(n0);
    Vector3f vt0 = Vector3f::Random();
    Vector3f vt1 = Vector3f::Random();
    Vector3f vt2 = Vector3f::Random();

    vt0[2] = 0.0f;
    vt1[2] = 0.0f;
    vt2[2] = 0.0f;

    Vector3f n = Vector3f::Random() * 0.25;
    n[2] = -1.0f;
    n = n.normalized();
    Vector3f r = (Vector3f::Random() + Vector3f{1.0f, 1.0f, 1.0f}) * 0.1;
    Vector3f nref{0.0f, 0.0f, 1.0f};

    vt0 += r[0] * nref;
    vt1 += r[1] * nref;
    vt2 += r[2] * nref;

    Vector3f n1 = compute_normal(vt0, vt1, vt2, r[0], r[1], r[2], n);
    cout << n1;

    return 0;
}