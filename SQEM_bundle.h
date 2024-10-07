#pragma once
#include <vector>
#include <Eigen/Dense>
#include "SQEM.h"
#include "SphereMeshes.h"

struct SQEMBundle{
    std::vector<SQEM> sqem;
    std::vector<Eigen::Vector3f> center; 
    std::vector<float> radius;
    std::vector<float> eval;
    int n;

    double minimize(Eigen::Vector3f &center, float &r, const Eigen::Vector3f &pa, const Eigen::Vector3f &pb, float radius_bound);
    double evaluate(Eigen::Vector3f &center, float r);

    SQEMBundle(int n = 1): n(n), sqem(n), center(n), radius(n), eval(n) {
        for (int i = 0; i < n; i++) {
            sqem[i].setZero();
        }
    }
};
struct CollapsedEdgeDeformed {
    float c;
    int u, v;
    SQEMBundle q;
    std::vector<Eigen::MatrixXf> fan;
    bool operator < (const CollapsedEdgeDeformed &e) const {
        return c > e.c;
    }
};

struct BundledSphereMesh: public SphereMesh {
    SQEMBundle Qbundle(int u) const;
    int n_deforms;
    std::vector<Eigen::MatrixXf> V_deform;
    std::vector<SQEMBundle> node_q;
    std::vector<std::vector<Eigen::MatrixXf>> node_fan;
    std::vector<Eigen::MatrixXf> N_deform;
    std::vector<Eigen::VectorXf> R_deform;
    
    float area(const Eigen::Vector3i &f, int i) const;
    inline int add_sphere(const CollapsedEdgeDeformed &cuv) {
        int w = nv ++;
        auto sqe {cuv.q};

        
        for (int i = 0; i < n_deforms; i++) {
            // per deform config 
            V_deform[i].conservativeResize(nv, Eigen::NoChange);
            Sphere s{sqe.center[i], sqe.radius[i]};
            V_deform[i].row(w) = s.q;
            R_deform[i].conservativeResize(nv);
            R_deform[i](w) = s.r;
        }

        valid.conservativeResize(nv);
        valid(w) = true;
        node_q.push_back(sqe);
        node_fan.push_back(cuv.fan);

        adj.push_back({});
        NI.conservativeResize(nv + 1);
        NI(w + 1) = NI(w);
        return w;
    }

    BundledSphereMesh(
        const std::string &filename = "bar.obj",
        int n_deforms = 1
        ): SphereMesh(filename), 
        n_deforms(n_deforms),
        V_deform(n_deforms), N_deform(n_deforms), R_deform(n_deforms) {
            for (auto &v: node_fan) {
                v.resize(n_deforms);
            }
        }
    
};