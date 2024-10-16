#include "SphereMeshes.h"
#include <iostream>
#include <iomanip>
#include <igl/readOBJ.h>
#include <igl/readOFF.h>
#include <igl/per_face_normals.h>
#include <algorithm>
using namespace std;
using namespace Eigen;
static const int INVALID = -1;

float SphereCloud::area(const Vector3i &f) const {
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

Vector3f SphereCloud::com(const Vector3i &f) const {
    Vector3f v0 = V.row(f[0]);
    Vector3f v1 = V.row(f[1]);
    Vector3f v2 = V.row(f[2]);
    
    return (1.0f / 3.0f) * (v0 + v1 + v2);
}



SQEM SphereCloud::Q(int f) const {
    assert(Fvalid(f));
    Vector3f Vu(V.row(F(f, 0)));
    Vector3f n = N.row(f);
    SQEM q(Vu, n);
    q *= (area(F.row(f)) / 3.0f);
    return q;
}

ColapsedEdge SphereCloud:: argmin_sqe(int u, int v) const {
    SQEM sqem{tile_q[u] + tile_q[v]};
    Vector3f center, pa{-10.0f, -10.0f, -10.0f}, pb {10.0f, 10.0f, 10.0f}; 
    float r;
    double c = sqem.minimize(center, r, pa, pb);

    pa = tile_com.row(u);
    pb = tile_com.row(v);

    return {static_cast<float>(c), u, v, {center, r}, sqem, {}};
}

void SphereCloud::init_queue() {
    for (int u = 0; u < nf_valid; u++) {
        for (int v = u + 1; v < nf_valid; v++) {
            q.push(argmin_sqe(u, v));
        }
    }
}


int SphereCloud::new_patch(ColapsedEdge &cuv) {
    int w = nf ++;
    int u = cuv.u, v = cuv.v;
    auto s {cuv.s};
    auto sqe {cuv.q};
    Fvalid.conservativeResize(nf);
    Fvalid[w] = true;
    

    tile_q.push_back(sqe);
    tile_R.conservativeResize(nf); 
    tile_R(w) = s.r;
    tile_C.conservativeResize(nf, NoChange);
    tile_C.row(w) = s.q;
    
    merged.push_back({});
    
    std::merge(merged[u].begin(), merged[u].end(), merged[v].begin(), merged[v].end(), std::back_inserter(merged[w]));

    tile_area.conservativeResize(nf);
    tile_area(w) = tile_area(u) + tile_area(v);
    tile_com.conservativeResize(nf, NoChange);
    tile_com.row(w) = (tile_com.row(u) * tile_area(u) + tile_com.row(v) * tile_area(v)) / tile_area(w);
    
    return w;


}
void SphereCloud::merge(int nv_target) {
    while (!q.empty() && nf_valid > nv_target) {
        auto cuv = q.top();
        q.pop();
        if (Fvalid(cuv.u) && Fvalid(cuv.v)) {
            int u = cuv.u, v = cuv.v;
            int w = new_patch(cuv);

            Fvalid(u) = false;
            Fvalid(v) = false;

            valid_set.erase(u);
            valid_set.erase(v);


            for (int i : valid_set) {
                assert(Fvalid(i));
                q.push(argmin_sqe(i, w));
            }
            valid_set.insert(w);
            nf_valid --;
        }
    }
}


SphereCloud::SphereCloud(const std::string &filename){
    string format = filename.substr(filename.size() - 3);
    if (format == "obj") {
        igl::readOBJ(filename, V, F);
    } else if (format == "off") {
        igl::readOFF(filename, V, F);
    } else {
        cout << "unsupported format" << endl;
        exit(1);
    }
    nf = F.rows();
    nf_valid = nf;

    Fvalid.resize(nf);
    Fvalid.fill(true);
    igl::per_face_normals_stable(V, F, N);
    assert(N.rows() == nf);

    merged.resize(nf);

    tile_q.resize(nf);
    tile_C.resize(nf, 3);
    tile_R.resize(nf);
    tile_com.resize(nf, 3);  
    tile_area.resize(nf);
    // valid_set.resize(nf);
    for (int i = 0; i < nf; i++) {
        tile_q[i] = Q(i);
        valid_set.insert(i);
        merged[i] = {i};
        tile_com.row(i) = com(F.row(i));
        tile_area(i) = area(F.row(i));
    }
    
}


void SphereCloud::export_ply(const string &fname) const {
    std::string plyname = fname;
    std::ofstream fout(plyname);

    // int nv_valid = tile_C.rows();
    int nv_valid = valid_set.size();

    fout << "ply" << std::endl;
    fout << "format ascii 1.0" << std::endl;
    fout << "element vertex " << nv_valid << std::endl;
    fout << "property float x" << endl;
    fout << "property float y" << endl;
    fout << "property float z" << endl;
    fout << "property float r" << endl;
    // fout << "element face " << nf_valid << endl;
    // fout << "property list uchar uint vertex_indices" << endl;
    fout << "end_header" << endl;
    
    
    for (auto i: valid_set) {
        fout << setiosflags(ios::fixed) << setprecision(15) << tile_C(i, 0) << " " << tile_C(i, 1) << " " << tile_C(i, 2) << " " << tile_R[i] << std::endl;
    }
        
}