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

    bool operator < (const ColapsedEdge &e) const {
        return c < e.c;
    }
    

};


struct SphereMesh {
    // per-vertex attributes
    Eigen::MatrixXf V;
    Eigen::VectorXf R;
    Eigen::Vector<bool, -1> valid;
    std::vector<std::vector<int>> adj;
    Eigen::VectorXi NI;
    Eigen::VectorXi VF;

    // per-face attributes
    Eigen::MatrixXi F;
    Eigen::MatrixXf N;
    Eigen::Vector<bool, -1> Fvalid; 


    Eigen::MatrixXi E;
    int nv, ne, nf, nv_valid, nf_valid;
    void simplify(int nv_target);
    Eigen::Vector3f compute_normal(const Eigen::Vector3i &f, const Eigen::Vector3f &n0) const;

// private: 
    DirectionalWidth dw;
    Eigen::MatrixXf polygon_fan(int u) const;
    ColapsedEdge argmin_sqe(int u, int v) const;
    std::priority_queue<ColapsedEdge> q;
    void reconnect_triangles(int u, int v, int w);
    SQEM Q(int u) const;
    SphereMesh::SphereMesh(const std::string &filename = "bar.obj");

    inline int add_sphere(const Sphere& s) {
        int w = nv ++;
        V.conservativeResize(nv, Eigen::NoChange);
        V.row(w) = s.q;
        R.conservativeResize(nv);
        R(w) = s.r;
        valid.conservativeResize(nv);
        valid(w) = true;
        adj.push_back({});
        NI.conservativeResize(nv + 1);
        NI(w) = NI(w - 1);
        // no need to resize VF for now
        return w;
    }
};
