#include "SphereMeshes.h"
#include "cnpy.h"
#include <iostream>
#include <nlohmann/json.hpp>
#include <fstream>

using namespace std;
using namespace Eigen;
using json = nlohmann::json;
int main(){
    std::ifstream config("../config_tile.json");
    json data = json::parse(config);
    string input_file = data["input_file"];
    string output_file = data["output_file"];
    string basename = output_file.substr(0, output_file.find_last_of("."));
    int nv_end = data["nv_end"];
    int nv_start = data["nv_start"];
    int nv_step = data["nv_step"];


    SphereCloud tiles_soup(input_file);
    tiles_soup.init_queue();

    int frame = 0;
    for (int nv_target = nv_start; nv_target >= nv_end; nv_target --) {
        tiles_soup.merge(nv_target);
        assert(tiles_soup.nf_valid == nv_target);
        if (nv_target % nv_step == 0) {
            tiles_soup.export_ply(basename + "_" + to_string(frame) + ".ply");
            frame ++;
        }
    }

    return 0;
}