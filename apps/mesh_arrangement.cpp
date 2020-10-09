#include <igl/copyleft/marching_cubes.h>
#include <igl/write_triangle_mesh.h>

#include <CLI/CLI.hpp>
#include <Eigen/Core>
#include <array>
#include <string>

#include "PyMesh/CellPartition.h"
#include "ScopedTimer.h"
#include "config.h"
#include "mesh_arrangement_utils.h"

int main(int argc, char** argv) {
    struct {
        std::string config_file;
        std::string output_mesh;
        std::string grid_file;
    } args;

    CLI::App app{"Piecewise implicit surface demo"};
    app.add_option("-G,--grid", args.grid_file, "Grid spec file")->required();
    app.add_option("config_file", args.config_file, "Configuration file")
        ->required();
    app.add_option("output_mesh", args.output_mesh, "Output mesh file")
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

    auto dot_pos = args.output_mesh.find_last_of(".");
    auto basename = args.output_mesh.substr(0, dot_pos);
    auto ext = args.output_mesh.substr(dot_pos);

    size_t num_meshes = meshes.size();
    for (size_t i = 0; i < num_meshes; i++) {
        std::string out_name = basename + ".input_" + std::to_string(i) + ext;
        igl::write_triangle_mesh(out_name, meshes[i].vertices, meshes[i].faces);
    }

    auto cells = compute_arrangement(meshes);
    auto num_cells = cells.size();
    for (size_t i = 0; i < num_cells; i++) {
        std::string out_name = basename + ".comp_" + std::to_string(i) + ext;
        igl::write_triangle_mesh(out_name, cells[i].vertices, cells[i].faces);
    }

    return 0;
}
