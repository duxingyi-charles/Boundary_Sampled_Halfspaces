#include <igl/copyleft/marching_cubes.h>
#include <igl/write_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>

#include <CLI/CLI.hpp>
#include <Eigen/Core>
#include <array>
#include <string>
#include <map>

#include "PyMesh/CellPartition.h"
#include "ScopedTimer.h"
#include "config.h"
#include "mesh_arrangement_utils.h"

int main(int argc, char** argv) {
    struct {
        std::string config_file;
        std::string grid_file;
    } args;

    CLI::App app{"Piecewise implicit surface demo"};
    app.add_option("-G,--grid", args.grid_file, "Grid spec file")->required();
    app.add_option("config_file", args.config_file, "Configuration file")
        ->required();
    CLI11_PARSE(app, argc, argv);

    auto implicit_functions =
        initialize_sampled_implicit_functions(args.config_file);
    auto grid_spec = parse_grid_spec(args.grid_file);

    std::vector<IGL_Mesh> meshes;
    meshes.reserve(implicit_functions.size() + 1);

    meshes.push_back(generate_cube(grid_spec));
    for (const auto& fn : implicit_functions) {
        meshes.push_back(marching_cubes(*fn, grid_spec));
    }

    auto cells = compute_arrangement(meshes);
    auto num_cells = cells.size();

    igl::opengl::glfw::Viewer viewer;

    std::map<int, Eigen::RowVector3d> colors;
    for (size_t i=1; i<num_cells; i++) {
        viewer.append_mesh();
        viewer.data(i).set_mesh(cells[i].vertices, cells[i].faces);
        colors.emplace(viewer.data().id, 0.5*Eigen::RowVector3d::Random().array() + 0.5);
    }
    viewer.launch();

    return 0;
}
