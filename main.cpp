#include "SphereMeshes.h"
int main() {
    SphereMesh mesh("bar.obj");
    mesh.simplify(10);
    return 0;
}