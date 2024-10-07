#include "DeformedMesh.h"
#include <iostream>
using namespace std;
using namespace Eigen;

void DeformedMesh::init_nodes() {
    for (int i = 0; i < nv; i++) {
        for (auto &shape: deformed_shapes) {
            auto pf{ polygon_fan(i, shape.V)};
            shape.node_q[i] = Q(i);
            shape.dir_widths[i] = dw.eval(pf);
        }
    }
}

void DeformedMesh::simplify(int nv_target) {
    for (int i = 0; i < ne; i ++) {
        int u{E(i, 0)}, v{E(i, 1)};
        auto cuv = argmin_sqe_bundle(u, v);
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
            for (auto i: adj[w]){
                assert(valid(i));
                auto ciw = argmin_sqe_bundle(i, w);
                q.push(ciw);
            }

            nv_valid --;
        }
    }
    lazy_delete();
}

CollapsedEdgeDeformed DeformedMesh::argmin_sqe_bundle(int u, int v) {
    CollapsedEdgeDeformed ret(n_deforms);
    ret.u = u; 
    ret.v = v;
    for (int i = 0; i < n_deforms; i++) {
        auto cuv_i = argmin_sqe(u, v, deformed_shapes[i]);
        ret.q.radius[i] = cuv_i.s.r;
        ret.q.center[i] = cuv_i.s.q;
        ret.q.sqem[i] = cuv_i.q;
        ret.dir_width[i] = cuv_i.dir_width;
        ret.c += cuv_i.c;
    }
    return ret;
}
