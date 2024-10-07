#include "SphereMeshes.h"
#include <iostream>
using namespace std;
using namespace Eigen;
int main() {
    SphereMeshBase mesh("../assets/hand.off");
    //SphereMesh mesh("teapot.obj");
    
    auto &V {mesh.rest_shape.V};
    auto &R {mesh.rest_shape.R};
    mesh.simplify(35);
    cout << "#V = " << V.rows() << endl << V.transpose() << endl;
    cout << "#R = " << R.rows() << endl << R.transpose() << endl;
    cout << "#F = " << mesh.F.rows() << endl << mesh.F << endl;

    mesh.export_ply("../output/hand.ply");
    return 0;
}