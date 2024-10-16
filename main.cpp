#include "SphereMeshes.h"
#include <iostream>
#include "cnpy.h"
#include <igl/writeOBJ.h>
#include <nlohmann/json.hpp>


using namespace std;
using namespace Eigen;
using json = nlohmann::json;


MatrixXd load_Phi() {
    auto Phi = cnpy::npy_load("../data/Q.npy");
    assert(Phi.word_size == sizeof(double));
    int rows = Phi.shape[0];
    int cols = Phi.shape[1];
    // cout<<rows<<" "<<cols<<endl;
    Map<Matrix<double, Dynamic, Dynamic, RowMajor>> Q(Phi.data<double>(), rows, cols);
    // cout << Q.col(1244) << endl;
    return Q;
}
void displace(int col, double magnitude, const MatrixXd &Q, MatrixXf &V) {

    VectorXd q = Q.row(col) * magnitude;
    // assert(V.rows() * 3 == q.rows());
    for (int i = 0; i < V.rows(); i++) {
        V.row(i) += q.segment<3>(3 * i).cast<float>();
    }
}
int main() {

    std::ifstream config("../config.json");
    json data = json::parse(config);
    string input_file = data["input_file"];
    string output_file = data["output_file"];
    string basename = output_file.substr(0, output_file.find_last_of("."));
    int nv_end = data["nv_end"];
    int nv_start = data["nv_start"];
    int nv_step = data["nv_step"];
    
    // SphereMesh mesh("../assets/hand.off");
    // DeformedMesh mesh(1, "../assets/hand.off");
    // DeformedMesh mesh(1, input_file);
    // DeformedMesh mesh(1, "../output/bar2_deformed.obj");
    SphereMesh mesh(input_file);

    // auto &V {mesh.V};
    // auto &R {mesh.R};
    
    // MatrixXd Q = load_Phi();
    // displace(1238, 1.0, Q, V);
    // igl::writeOBJ("../output/bar2_deformed.obj", V, mesh.F);

    int frame = 0;
    mesh.init_queue();
    for (int nv_target = nv_start; nv_target >= nv_end; nv_target --) {
        mesh.simplify(nv_target); 
        if (nv_target % nv_step == 0) {
            auto [V, R, F] = mesh.lazy_delete();
            string seq = basename + to_string(frame ++) + ".ply";
            mesh.export_ply(seq, V, R, F);
        }
    }
    

    auto [V, R, F] = mesh.lazy_delete();

    // mesh.export_ply("../output/hand.ply");
    string seq = basename + to_string(frame ++) + ".ply";
    cout << "final mesh: " << seq << endl;
    mesh.export_ply(seq, V, R, F);

    cout << "#V = " << V.rows() << endl << V.transpose() << endl;
    cout << "#R = " << R.rows() << endl << R.transpose() << endl;
    cout << "#F = " << F.rows() << endl << F << endl;

    return 0;
}