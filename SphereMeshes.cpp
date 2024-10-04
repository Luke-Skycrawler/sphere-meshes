#include "SQEM.h"
#include <igl/adjacency_list.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/readOBJ.h>
#include <igl/edges.h>
#include <igl/per_face_normals.h>
#include <Eigen/Dense>
#include <iostream>
#include <algorithm>
#include <assert.h>
#include "SphereMeshes.h"

using namespace std;
using namespace Eigen;


void remove_duplicates(vector<int> &a) {
    sort(a.begin(), a.end());
    auto last = std::unique(a.begin(), a.end());
    a.erase(last, a.end());
}

MatrixXf SphereMesh::polygon_fan(int u) const{
    Vector3f pu {V.row(u)};
    int len = (NI(u + 1) - NI(u)) * 3 + 1;
    MatrixXf ret(len, 3);   
    ret.row(0) = pu;
    for (int i = NI(u); i < NI(u + 1); i ++) {
        Vector3i face = F.row(VF(i, 0));

        int j = 0;
        for (int k = 0; k < 3; k ++) if (face[k] == u) j = k;
        int offset = 1 + (i - NI(u)) * 3;
        Vector3f p0 = V.row(face[0]); 
        Vector3f p1 = V.row(face[1]);
        Vector3f p2 = V.row(face[2]);
        Vector3f barycenter = (p0 + p1 + p2) / 3.0f;
        ret.row(offset ++) = barycenter;
        for (int ii = 0; ii < 3; ii++)
            if (ii != j) {
                ret.row(offset ++) = (V.row(face[ii]) + V.row(face[j])) * 0.5f;
            }
    }
    return ret;
}
SphereMesh::SphereMesh(const std::string &filename){
    igl::readOBJ(filename, V, F);
    igl::edges(F, E);
    nv = V.rows();
    ne = E.rows(); 
    nf = F.rows();

    Fvalid.resize(nf);
    Fvalid.fill(true);
    igl::per_face_normals_stable(V, F, N);
    assert(N.rows() == nf);

    valid.resize(nv);
    valid.fill(true);
    R.resize(nv);
    R.fill(0.0f);
    igl::vertex_triangle_adjacency(F, nv, VF, NI);
    igl::adjacency_list(F, adj);
    assert(adj.size() == nv);
    // does not include isolate vertices
    assert(NI.rows() == nv + 1);
    assert(VF.rows() == nv);

}

Vector3f SphereMesh::compute_normal(const Vector3i &f, const Vector3f &n0) const{
    const float tol = 1e-6f;
    Vector3f v0 = V.row(f[0]);
    Vector3f v1 = V.row(f[1]);
    Vector3f v2 = V.row(f[2]);

    float r0 = R(f[0]);
    float r1 = R(f[1]);
    float r2 = R(f[2]);
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

SQEM SphereMesh::Q(int u) const {
    // computes spherical quadric error metric for vertex u
    SQEM q;
    Vector3f Vu{V.row(u)};
    for (int i = NI(u); i < NI(u + 1); i ++) {
        int fu = VF(i);
        Vector3f nu {N.row(fu)};
        assert(Fvalid(fu));
        Vector3f _Vu = Vu + nu * R(u);
        q += SQEM(_Vu, nu);
    }
    return q;
}

ColapsedEdge SphereMesh::argmin_sqe(int u, int v) const {
    SQEM sqem {Q(u) + Q(v)};
    Vector3f center, Vu(V.row(u)), Vv(V.row(v));
    float r;
    auto fanu = polygon_fan(u), fanv = polygon_fan(v);
    MatrixX3f fanuv(fanu.rows() + fanv.rows(), 3);
    fanuv.block(0, 0, fanu.rows(), 3) = fanu;
    fanuv.block(fanu.rows(), 0, fanv.rows(), 3) = fanv;

    float radius_bound = dw.W(fanu);
    sqem.minimize(center, r, Vu, Vv, radius_bound);

    return {static_cast<float>(sqem.evaluate(center, r)), u, v, {center, r}};
}

void SphereMesh::reconnect_triangles(int u, int v, int w){
    /****************************************************
    collapse edge uv to w
    reconnect all vertices adjacent to u to w
    NOTE: often need to call reconnect(v, u, w) as well
    *****************************************************/
    valid(u) = false;
    for (auto fu: adj[u]) {
        assert(valid(fu));
        Vector3i Fu = F.row(fu);
        if (Fu[0] != v && Fu[1] != v && Fu[2] != v) {
            // one of the vertex in {u, v}
            // face still valid but replace the vertex with w

            
            // int iw = 0;
            // for (int i = 0; i < 3; i++) if (Fu(i) == u) {
            //     Fu(i) = w;
            //     iw = i;
            // }

            F.row(fu) = Fu; 
            N.row(fu) = compute_normal(Fu, N.row(fu));

            // construct face adj list for new node w
            VF.conservativeResize(VF.rows() + 1);
            VF(VF.rows() - 1) = fu;
            NI(w + 1) += 1;

            for (int i = 0; i < 3; i ++) {  // assuming triangle mesh
                int j = Fu(i);
                assert(j != v && valid(j));
                if (j != w) {
                    // delete u from adjacent list and replace them with w
                    auto it = find(adj[j].begin(), adj[j].end(), u);
                    assert(it != adj[j].end());
                    *it = w;

                    // construct vertex adj list for new node w
                    adj[w].push_back(j);

                }
            }


        } else {
            // face with collapsed edge uv, invalidate it
            Fvalid(fu) = false;
        }
    }
}
void SphereMesh::simplify(int nv_target) {
    for (int i = 0; i < ne; i ++) {
        int u = E(i, 0), v = E(i, 1);
        auto cuv = argmin_sqe(u, v);
        q.push(cuv);
    }
    while (!q.empty() && nv_valid > nv_target) {
        auto cuv = q.top();
        q.pop();
        
        if (valid(cuv.u) && valid(cuv.v)) {
            int u = cuv.u, v = cuv.v;           
            int w = add_sphere(cuv.s);

            reconnect_triangles(u, v, w);
            reconnect_triangles(v, u, w);
            remove_duplicates(adj[w]);
            for (auto i: adj[w]) {
                assert(valid(i));
                auto ciw = argmin_sqe(i, w);
                q.push(ciw);
            }

            nv_valid --;
            // two vertices are collapsed into one
        }
        
    }

}
