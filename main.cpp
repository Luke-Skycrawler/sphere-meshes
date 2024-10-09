// #include "SphereMeshes.h"
#include "DeformedMesh.h"
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
    assert(V.rows() * 3 == q.rows());
    for (int i = 0; i < V.rows(); i++) {
        V.row(i) += q.segment<3>(3 * i).cast<float>();
    }
}
int main() {

    std::ifstream config("../config.json");
    json data = json::parse(config);
    string input_file = data["input_file"];
    string output_file = data["output_file"];
    int nv_target = data["nv_target"];

    // SphereMesh mesh("../assets/hand.off");
    // DeformedMesh mesh(1, "../assets/hand.off");
    // DeformedMesh mesh(1, input_file);
    // DeformedMesh mesh(1, "../output/bar2_deformed.obj");
    SphereMesh mesh(input_file);

    auto &V {mesh.rest_shape.V};
    auto &R {mesh.rest_shape.R};
    
    // MatrixXd Q = load_Phi();
    // displace(1238, 1.0, Q, V);
    // igl::writeOBJ("../output/bar2_deformed.obj", V, mesh.F);

    mesh.simplify(nv_target);
    cout << "#V = " << V.rows() << endl << V.transpose() << endl;
    cout << "#R = " << R.rows() << endl << R.transpose() << endl;
    cout << "#F = " << mesh.F.rows() << endl << mesh.F << endl;

    // mesh.export_ply("../output/hand.ply");
    mesh.export_ply(output_file);

    return 0;
}
