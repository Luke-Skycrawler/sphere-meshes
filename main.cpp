#include "SphereMeshes.h"
#include <iostream>
using namespace std;
using namespace Eigen;
int main() {
    SphereMesh mesh("bar.obj");
    mesh.simplify(2);
    cout << "#V = " << mesh.V.rows() << endl << mesh.V.transpose() << endl;
    cout << "#R = " << mesh.R.rows() << endl << mesh.R.transpose() << endl;
    cout << "#F = " << mesh.F.rows() << endl << mesh.F << endl;
    return 0;
}