#pragma once
#include <vector>
#include <Eigen/Dense>
#include "SQEM.h"
#include "SphereMeshes.h"

struct SQEMBundle{
    std::vector<SQEM> sqem;
    std::vector<Eigen::Vector3f> center; 
    std::vector<float> radius;
    int n;


    SQEMBundle(int n = 1): n(n), sqem(n), center(n), radius(n) {
        for (int i = 0; i < n; i++) {
            sqem[i].setZero();
        }
    }

    
};
struct CollapsedEdgeDeformed {
    float c;
    int u, v;
    SQEMBundle q;
    std::vector<Eigen::MatrixXf> dir_width;
    bool operator < (const CollapsedEdgeDeformed &e) const {
        return c > e.c;
    }
    inline Sphere get_sphere(int i) const {
        return Sphere{q.center[i], q.radius[i]};
    }
    CollapsedEdgeDeformed(int n): q(n), dir_width(n), c(0.0f) {}
};

struct DeformedMesh: public SphereMeshBase {
    int n_deforms;
    inline int add_sphere(CollapsedEdgeDeformed &cuv) {

        for (int i = 0; i < n_deforms; i++) {
            Sphere s {cuv.get_sphere(i)};
            deformed_shapes[i].new_sphere_node(s, cuv.q.sqem[i], cuv.dir_width[i]);
        }
        return new_sphere_connectivity();
    }
    void simplify(int nv_target);
    CollapsedEdgeDeformed argmin_sqe_bundle(int u, int v);
    std::priority_queue<CollapsedEdgeDeformed> q;

    DeformedMesh(int n_deforms = 1, const std::string &file = "bar.obj"): SphereMeshBase(file), n_deforms(n_deforms) {
        init_nodes();
    }

    void init_nodes();
};