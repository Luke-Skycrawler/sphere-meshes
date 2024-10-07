#include "SQEM.h"
#include <igl/adjacency_list.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/readOBJ.h>
#include <igl/edges.h>
#include <igl/per_face_normals.h>
#include <igl/readOFF.h>
#include <Eigen/Dense>
#include <iostream>
#include <algorithm>
#include <assert.h>
#include <iomanip>
#include "SphereMeshes.h"

using namespace std;
using namespace Eigen;
static const int INVALID = -1;


MatrixXf SphereMeshBase::polygon_fan(int u, const MatrixXf &V) const{
    
    Vector3f pu {V.row(u)};
    int len = (NI(u + 1) - NI(u)) * 3 + 1;

    MatrixXf ret(len, 3);   
    int offset = 0;
    ret.row(offset++) = pu;
    for (int i = NI(u); i < NI(u + 1); i ++) if (Fvalid(VF(i))) {
        Vector3i face = F.row(VF(i));

        int j = -1;
        for (int k = 0; k < 3; k ++) if (face[k] == u) j = k;
        assert(j != -1);
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
    ret.conservativeResize(offset, NoChange);
    return ret;
}

Positions::Positions(const string &filename, MatrixXi &F) {
    igl::readOFF(filename, V, F);
    int nv = V.rows(), nf = F.rows();
    igl::per_face_normals_stable(V, F, N);
    assert(N.rows() == nf);
    R.resize(nv);
    R.fill(0.0f);
    node_q.resize(nv);
    dir_widths.resize(nv);

}
void SphereMeshBase::init_connectivity(){
    igl::edges(F, E);
    ne = E.rows(); 
    nf = F.rows();
    nv_valid = nv;
    nf_valid = nf;

    Fvalid.resize(nf);
    Fvalid.fill(true);

    valid.resize(nv);
    valid.fill(true);
    igl::vertex_triangle_adjacency(F, nv, VF, NI);
    igl::adjacency_list(F, adj);
    // TODO: assert adj is sorted
    assert(adj.size() == nv);
    // does not include isolate vertices
    assert(NI.rows() == nv + 1);
    assert(VF.rows() == NI(nv));
}

void SphereMesh::init_nodes() {
    for (int i = 0; i < nv; i++) {
        for (auto &shape: deformed_shapes) {
            auto pf{ polygon_fan(i, shape.V)};
            shape.node_q[i] = Q(i);
            shape.dir_widths[i] = dw.eval(pf);
        }
    }
}
SQEM SphereMeshBase::Q(int u) const {
    // computes spherical quadric error metric for vertex u
    SQEM q;
    q.setZero();
    Vector3f Vu{rest_shape.V.row(u)};
    assert(valid(u));
    for (int i = NI(u); i < NI(u + 1); i ++) if (Fvalid(VF(i))){
        int fu = VF(i);
        Vector3f nu {rest_shape.N.row(fu)};
        Vector3f _Vu = Vu + nu * rest_shape.R(u);
        q += SQEM(_Vu, nu) * (area(F.row(fu)) / 3.0);
    }
    return q;
}

float SphereMeshBase::area(const Vector3i &f) const {
    return area(rest_shape.V, f);
}

float SphereMeshBase::area(const MatrixXf &V, const Vector3i &f) const {
    
    if (f[0] == INVALID || f[1] == INVALID || f[2] == INVALID) {
        return 0.0f;
    }
    Vector3f v0 = V.row(f[0]);
    Vector3f v1 = V.row(f[1]);
    Vector3f v2 = V.row(f[2]);
    Vector3f e0 = v1 - v0;
    Vector3f e1 = v2 - v0;
    return 0.5f * e0.cross(e1).norm();
}

ColapsedEdge SphereMeshBase::argmin_sqe(int u, int v, const Positions &p) const {
    // SQEM sqem {Q(u) + Q(v)};
    auto &node_q {p.node_q};
    auto &dir_widths {p.dir_widths};
    auto &V {p.V};

    SQEM sqem{node_q[u] + node_q[v]};
    Vector3f center, Vu(V.row(u)), Vv(V.row(v));
    float r;
    auto boundu = dir_widths[u], boundv = dir_widths[v];
    int n_dirs = dw.n_dirs + 3;
    VectorXf lu = boundu.col(0), lv = boundv.col(0);
    VectorXf uu = boundu.col(1), uv = boundv.col(1);
    auto boundw = MatrixXf(n_dirs, 2);
    boundw.col(0) = lu.cwiseMin(lv);
    boundw.col(1) = uu.cwiseMax(uv);
    float radius_bound = dw.W(boundw);

    double c = sqem.minimize(center, r, Vu, Vv, radius_bound);
    // double c = sqem.evaluate(center, r);
    return {static_cast<float>(c), u, v, {center, r}, sqem, boundw};
}

void SphereMeshBase::reconnect_triangles(int u, int v, int w){
    /****************************************************
    collapse edge uv to w
    reconnect all vertices adjacent to u to w
    NOTE: often need to call reconnect(v, u, w) as well
    *****************************************************/
    valid(u) = false;
    for (int jj = NI(u); jj < NI(u + 1); jj ++) {
        int fu = VF(jj);
        if (Fvalid(fu)) {

            Vector3i Fu = F.row(fu);
            if (Fu[0] != v && Fu[1] != v && Fu[2] != v) {
                // one of the vertex in {u, v}
                // still a triangle, replace the vertex with w
                int iw = 0;
                for (int i = 0; i < 3; i++) if (Fu(i) == u) {
                    Fu(i) = w;
                    iw = i;
                }

                F.row(fu) = Fu; 
                // N.row(fu) = compute_normal(Fu, N.row(fu));

                // construct face adj list for new node w
                VF.conservativeResize(VF.rows() + 1);
                VF(VF.rows() - 1) = fu;
                NI(w + 1) += 1;
            } else {
                // face with collapsed edge uv
                for (int i = 0; i < 3; i ++) {
                    if (Fu(i) != u && Fu(i) != v && Fu(i) != INVALID) {
                        // triangle collapsed to an edge 
                        int j = Fu(i);
                        sort(adj[u].begin(), adj[u].end());
                        sort(adj[j].begin(), adj[j].end());
                        sort(adj[v].begin(), adj[v].end());
                        vector<int> uj_inter, vj_inter;

                        set_intersection(adj[u].begin(), adj[u].end(), adj[j].begin(), adj[j].end(), back_inserter(uj_inter));
                        set_intersection(adj[v].begin(), adj[v].end(), adj[j].begin(), adj[j].end(), back_inserter(vj_inter));

                        const auto filter = [&](int x) {
                            return valid(x) && x != v;
                        };

                        int countu = count_if(adj[u].begin(), adj[u].end(), filter);
                        int countv = count_if(adj[v].begin(), adj[v].end(), filter);
                        

                        if (countu || countv) {
                            delete_face(fu);
                        }
                        else {
                            F.row(fu) = Vector3i{w, j, INVALID};
                        }
                        break;       
                    }
                    else if (i == 2) {
                        // edge collapsed to a point
                        delete_face(fu);
                    }
                }
            }
        
        }
    }

    for (auto i: adj[u]) {
        if (valid(i) && i != v){
            auto it = find(adj[i].begin(), adj[i].end(), u);
            if (it != adj[i].end()) {
                *it = w;
            }
            adj[w].push_back(i);
        }
    }


}

void SphereMeshBase::delete_face(int f) {
    // lazy delete
    nf_valid --;
    Fvalid(f) = false;
}
void SphereMeshBase::simplify(int nv_target) {
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
            int w = add_sphere(cuv);

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
    lazy_delete();
}

void SphereMeshBase::lazy_delete_vertex(Positions &p) {
    MatrixXf Vnew(nv_valid, 3);
    VectorXf Rnew(nv_valid);
    VectorXi Vmap(nv);
    int v = 0;
    for (int i = 0; i < nv;  i ++ ) {
        if (valid(i)) {
            Rnew(v) = p.R(i);
            Vnew.row(v++) = p.V.row(i);
        }
    }
    p.V = Vnew;
    p.R = Rnew;
}

void SphereMeshBase::lazy_delete_face(const VectorXi &Vmap) {
    MatrixXi Fnew(nf, 3);
    const auto map = [&](Vector3i f) -> Vector3i {
        return {Vmap(f(0)), Vmap(f(1)), Vmap(f(2))};
    };
    int f = 0;
    for (int i = 0; i < nf; i ++) {
        if (Fvalid(i)) {
            Vector3i fi{F.row(i)};
            Fnew.row(f++) = map(fi);
        }
    }
    F = Fnew;
    F.conservativeResize(f, NoChange);
}
void SphereMeshBase::export_ply(const string &fname) const {

    std::string plyname = fname;
    std::ofstream fout(plyname);

    //	GraphVertexIterator gvi,gvi_end;
    fout << "ply" << std::endl;
    fout << "format ascii 1.0" << std::endl;
    fout << "element vertex " << nv_valid << std::endl;
    fout << "property float x" << endl;
    fout << "property float y" << endl;
    fout << "property float z" << endl;
    fout << "property float r" << endl;
    fout << "element face " << nf_valid << endl;
    fout << "property list uchar uint vertex_indices" << endl;
    fout << "end_header" << endl;
    
    auto &V {rest_shape.V};
    for(int i = 0; i < V.rows(); i ++)
        fout << setiosflags(ios::fixed) << setprecision(15) << V(i, 0) << " " << V(i, 1) << " " << V(i, 2) << " " << rest_shape.R[i] << std::endl;


    for (int i = 0; i < F.rows(); i++) {
        if (F(i, 2) == INVALID) {
            fout << "2 " << F(i, 0) << " " << F(i, 1) << std::endl;
        }
        else {
            fout << "3 " << F(i, 0) << " " << F(i, 1) << " " << F(i, 2) << std::endl;
        }
    }
    
    fout.close();

}