import igl
import numpy as np

def tet2surface():
    file_surface = "assets/bar.obj"
    file_deformed = "assets/bar_deformed.obj"
    file_output = "assets/bar_deformed_surface.obj"

    surface_V, _, _, surface_F, _, _ = igl.read_obj(file_surface)
    deformed_V, _, _, deformed_F, _, _ = igl.read_obj(file_deformed)
    output_V = deformed_V[:surface_V.shape[0], :]
    igl.write_obj(file_output, output_V, surface_F)

if __name__ == "__main__":
    tet2surface()
