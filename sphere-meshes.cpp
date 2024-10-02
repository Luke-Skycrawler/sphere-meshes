#include "SQEM.h"
#include <queue> 
#include <vector>
#include <igl/adjacency_list.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/readOBJ.h>
#include <igl/edges.h>
#include <igl/per_face_normals.h>
#include <Eigen/Dense>
#include <string>
#include <iostream>
#include <algorithm>
using namespace std;
using namespace Eigen;

struct Sphere {
    Vector3f q; 
    float r;
};
struct ColapsedEdge {
    float c;
    int u, v;
    Sphere s;

    bool operator < (const ColapsedEdge &e) const {
        return c < e.c;
    }
    

};

void remove_duplicates(vector<int> &a) {
    sort(a.begin(), a.end());
    auto last = std::unique(a.begin(), a.end());
    a.erase(last, a.end());
}
ColapsedEdge argmin_sqe(int u, int v, const MatrixXf &V, const MatrixXi &VF, const VectorXi &NI, const MatrixXf &N) {
    SQEM sqem;
    Vector3f Vu(V(u)), Vv(V(v));
    for (int i = NI(u); i < NI(u + 1); i ++) {
        int fu = VF(i);
        Vector3f nu {N(fu)};
        sqem += SQEM(Vu, nu);
    }
    for (int i = NI(v); i < NI(v + 1); i ++) {
        int fv = VF(i);
        Vector3f nv {N(fv)};
        sqem += SQEM(Vv, nv);
    }
    Vector3f center;
    float r;
    sqem.minimize(center, r, Vu, Vv);
    return {sqem.evaluate(center, r), u, v, {center, r}};
}
float compute_sqe(Vector3f v, MatrixXi F) {

}

struct SphereMesh {
    Eigen::Vector<bool, -1> valid; 
    Eigen::MatrixXf V, N;
    Eigen::MatrixXi F, E;
    int nv, ne, nf, cnt;
    vector<vector<int>> adj;
    Eigen::MatrixXi VF;
    Eigen::VectorXi NI;
    SphereMesh(const string &filename = "bar.obj") {
        cnt = 0;
        igl::readOBJ(filename, V, F);
        igl::edges(F, E);
        nv = V.rows();
        ne = E.rows(); 
        nf = F.rows();
        igl::per_face_normals_stable(V, F, N);

        igl::vertex_triangle_adjacency(F, nv, VF, NI);
        igl::adjacency_list(F, adj);
    }
    void simplify(int nv_target);
private: 
    priority_queue<ColapsedEdge> q;
    void connect_triangles(int u, int v, int w);
};

void SphereMesh::connect_triangles(int u, int v, int w){

    for (auto fu: adj[u]) {
        Vector3i Fu = F.row(fu);
        if (Fu[0] != v && Fu[1] != v && Fu[2] != v) {
            Vector3i Fuw {Fu};
            for (int i = 0; i < 3; i ++) {
                if (Fuw[i] == u) Fuw[i] = w;
                else {
                    adj[w].push_back(Fuw[i]);
                    // might include invalid vertices

                    // TODO:
                    // disconnect from invalid u and v
                }
            }
            
            F.conservativeResize(V.rows() + 1, NoChange);
            F.row(nf ++) = Fuw;
            // keeping nf = F.rows() 
        }
    }
}
void SphereMesh::simplify(int nv_target) {
    for (int i = 0; i < ne; i ++) {
        int u = E(i, 0), v = E(i, 1);
        auto cuv = argmin_sqe(u, v, V, VF, NI, N);
        q.push(cuv);
    }
    while (!q.empty() && nv > nv_target) {
        auto cuv = q.top();
        q.pop();
        
        if (valid[cuv.u] && valid[cuv.v]) {
            int u = cuv.u, v = cuv.v;
            valid[u] = false;
            valid[v] = false;

            int w = cnt ++;
            V.conservativeResize(V.rows() + 1, NoChange);
            V.row(w) = cuv.s.q;
            adj.push_back({});
            
            connect_triangles(u, v, w);
            connect_triangles(v, u, w);
            remove_duplicates(adj[w]);
            for (auto i: adj[w]) {
                if (valid[i]) {
                    auto ciw = argmin_sqe(i, w, )
                }
            }

            nv --;
            // two vertices become one
        }
    }

}
int main() {


    SphereMesh mesh("bar.obj");
    mesh.simplify(10);
    
    return 0;
}