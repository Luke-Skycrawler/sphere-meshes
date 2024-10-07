#pragma once
#include <Eigen/Dense>
#include <queue>
#include <vector>
#include "SQEM.h"
#include "DirectionalWidth.h"

inline void remove_duplicates(std::vector<int> &a) {
    std::sort(a.begin(), a.end());
    auto last = std::unique(a.begin(), a.end());
    a.erase(last, a.end());
}

struct Sphere {
    Eigen::Vector3f q; 
    float r;
};
struct ColapsedEdge {
    float c;
    int u, v;
    Sphere s;
    SQEM q;
    Eigen::MatrixXf dir_width;
    bool operator < (const ColapsedEdge &e) const {
        return c > e.c;
    }
};

struct Positions {
    Eigen::MatrixXf V;
    Eigen::VectorXf R;
    std::vector<SQEM> node_q; 
    std::vector<Eigen::MatrixXf> dir_widths;    

    Eigen::MatrixXf N;
    inline void new_sphere_node(const Sphere &s, const SQEM &sqe, const Eigen::MatrixXf &dir_width){
        int nv = V.rows(); 

        int w = nv ++;
        V.conservativeResize(nv, Eigen::NoChange);
        R.conservativeResize(nv);

        V.row(w) = s.q;
        R(w) = s.r;
        node_q.push_back(sqe);
        dir_widths.push_back(dir_width);
    }
    Positions(const std::string& file, Eigen::MatrixXi &F);
    Positions() {}
};

struct SphereMeshBase {
    // per-vertex attributes
    Eigen::Vector<bool, -1> valid;
    std::vector<std::vector<int>> adj;
    Eigen::VectorXi NI;
    Eigen::VectorXi VF;

    void SphereMeshBase::init_connectivity();


    // per-face attributes
    Eigen::MatrixXi F;

    // vertex attributes vary by deformation
    // must init after F
    std::vector<Positions> deformed_shapes;
    Positions& rest_shape;

    Eigen::Vector<bool, -1> Fvalid; 


    Eigen::MatrixXi E;
    int nv, ne, nf, nv_valid, nf_valid;
    virtual void simplify(int nv_target);
    // Eigen::Vector3f compute_normal(const Eigen::Vector3i &f, const Eigen::Vector3f &n0) const;
    void export_ply(const std::string &filename) const;

    inline Eigen::VectorXi construct_vmap(){
        Eigen::VectorXi Vmap(nv);
        int v = 0;
        for (int i = 0; i < nv; i++) if (valid(i)) {
            Vmap(i) = v ++;
        }
        return Vmap;
    }
    inline void lazy_delete() {
        auto vmap = construct_vmap();
        for (auto &p: deformed_shapes) {
            lazy_delete_vertex(p);
        }
        lazy_delete_face(vmap);
    }
    void lazy_delete_vertex(Positions &p);
    void lazy_delete_face(const Eigen::VectorXi &vmap);
// private: 
    DirectionalWidth dw;
    Eigen::MatrixXf polygon_fan(int u, const Eigen::MatrixXf &V) const;
    ColapsedEdge argmin_sqe(int u, int v, const Positions &p) const;
    ColapsedEdge argmin_sqe(int u, int v) const {
        return argmin_sqe(u, v, rest_shape);
    }
    std::priority_queue<ColapsedEdge> q;
    void reconnect_triangles(int u, int v, int w);
    SQEM Q(int u) const;
    SphereMeshBase(const std::string &filename = "bar.obj"): deformed_shapes({Positions(filename, F)}), rest_shape(deformed_shapes[0]) {
    // SphereMeshBase(const std::string &filename = "bar.obj"): rest_shape(filename, F) {
        nv = rest_shape.V.rows();
        init_connectivity();
    }

    float area(const Eigen::Vector3i &f) const;
    float area(const Eigen::MatrixXf &V, const Eigen::Vector3i &f) const;
    void delete_face(int f);
    inline int new_sphere_connectivity() {
        int w = nv ++;
        valid.conservativeResize(nv);
        valid(w) = true;
        adj.push_back({});
        NI.conservativeResize(nv + 1);
        NI(w + 1) = NI(w);
        // no need to resize VF for now
        return w;
    }

    virtual inline int add_sphere(ColapsedEdge &cuv) {
        rest_shape.new_sphere_node(cuv.s, cuv.q, cuv.dir_width);
        return new_sphere_connectivity();
    }
    virtual void init_nodes() = 0;
};

struct SphereMesh: public SphereMeshBase {
    SphereMesh(const std::string &filename = "bar.obj"): SphereMeshBase(filename) {
        init_nodes();
    }
    void init_nodes();
    
};

