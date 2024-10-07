#pragma once
#include <Eigen/Dense>
#include <queue>
#include <vector>
#include "SQEM.h"
#include "DirectionalWidth.h"
struct Sphere {
    Eigen::Vector3f q; 
    float r;
};
struct ColapsedEdge {
    float c;
    int u, v;
    Sphere s;
    SQEM q;
    Eigen::MatrixXf fan;
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
    inline void new_sphere_node(const ColapsedEdge &cuv){
        int nv = V.rows(); 

        int w = nv ++;
        V.conservativeResize(nv, Eigen::NoChange);
        R.conservativeResize(nv);

        auto s {cuv.s};
        auto sqe {cuv.q};
        V.row(w) = s.q;
        R(w) = s.r;
        node_q.push_back(cuv.q);
        dir_widths.push_back(cuv.fan);
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
    // std::vector<Positions> deformed_shapes;
    Positions rest_shape;

    Eigen::Vector<bool, -1> Fvalid; 


    Eigen::MatrixXi E;
    int nv, ne, nf, nv_valid, nf_valid;
    void simplify(int nv_target);
    // Eigen::Vector3f compute_normal(const Eigen::Vector3i &f, const Eigen::Vector3f &n0) const;
    void export_ply(const std::string &filename) const;

// private: 
    DirectionalWidth dw;
    Eigen::MatrixXf polygon_fan(int u) const;
    ColapsedEdge argmin_sqe(int u, int v) const;
    std::priority_queue<ColapsedEdge> q;
    void reconnect_triangles(int u, int v, int w);
    SQEM Q(int u) const;
    // SphereMeshBase(const std::string &filename = "bar.obj"): deformed_shapes({Positions(filename, F)}), rest_shape(deformed_shapes[0]) {
    SphereMeshBase(const std::string &filename = "bar.obj"): rest_shape(filename, F) {
        nv = rest_shape.V.rows();
        init_connectivity();
        init_nodes();
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
        rest_shape.new_sphere_node(cuv);
        return new_sphere_connectivity();
    }
    virtual void init_nodes();
};


struct SphereMesh: public SphereMeshBase {
    SphereMesh(const std::string &filename = "bar.obj"): SphereMeshBase(filename) {}
};