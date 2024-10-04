#include "SphereMeshes.h"
#include <iostream>
using namespace std;
using namespace Eigen;
int main() {
    SphereMesh mesh("../assets/hand.off");
    //SphereMesh mesh("teapot.obj");

    mesh.simplify(35);
    cout << "#V = " << mesh.V.rows() << endl << mesh.V.transpose() << endl;
    cout << "#R = " << mesh.R.rows() << endl << mesh.R.transpose() << endl;
    cout << "#F = " << mesh.F.rows() << endl << mesh.F << endl;

    mesh.export_ply("../output/hand.ply");
    return 0;
}