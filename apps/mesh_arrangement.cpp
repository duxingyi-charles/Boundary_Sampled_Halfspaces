#include <igl/copyleft/marching_cubes.h>
#include <igl/write_triangle_mesh.h>

#include <CLI/CLI.hpp>
#include <Eigen/Core>
#include <array>
#include <string>

#include "PyMesh/CellPartition.h"
#include "config.h"
#include "ScopedTimer.h"

struct IGL_Mesh {
    Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> vertices;
    Eigen::Matrix<int, Eigen::Dynamic, 3, Eigen::RowMajor> faces;
};

IGL_Mesh generate_cube(const GridSpec& grid_spec) {
    ScopedTimer<> timer("Generate cube");
    constexpr double eps = 0;
    auto bbox_min = grid_spec.bbox_min.array() + eps;
    auto bbox_max = grid_spec.bbox_max.array() - eps;
    IGL_Mesh cube;
    cube.vertices.resize(8, 3);
    cube.vertices <<
        bbox_min[0], bbox_min[1], bbox_max[2],
        bbox_min[0], bbox_min[1], bbox_min[2],
        bbox_max[0], bbox_min[1], bbox_min[2],
        bbox_max[0], bbox_min[1], bbox_max[2],
        bbox_min[0], bbox_max[1], bbox_max[2],
        bbox_max[0], bbox_max[1], bbox_max[2],
        bbox_max[0], bbox_max[1], bbox_min[2],
        bbox_min[0], bbox_max[1], bbox_min[2];

    cube.faces.resize(12, 3);
    cube.faces << 0, 1, 2, 0, 2, 3, 4, 5, 6, 4, 6, 7, 0, 3, 5, 0, 5, 4, 1, 7, 6,
        1, 6, 2, 2, 6, 5, 2, 5, 3, 0, 4, 7, 0, 7, 1;

    return cube;
}

IGL_Mesh merge_meshes(const std::vector<IGL_Mesh>& meshes) {
    ScopedTimer<> timer("Merge meshes");
    size_t num_meshes = meshes.size();
    std::vector<size_t> num_vertices(num_meshes + 1, 0);
    std::vector<size_t> num_faces(num_meshes + 1, 0);

    for (size_t i = 0; i < num_meshes; i++) {
        num_vertices[i + 1] = meshes[i].vertices.rows();
        num_faces[i + 1] = meshes[i].faces.rows();
    }

    std::partial_sum(num_vertices.begin(), num_vertices.end(),
                     num_vertices.begin());
    std::partial_sum(num_faces.begin(), num_faces.end(), num_faces.begin());

    IGL_Mesh merged_mesh;

    merged_mesh.vertices.resize(num_vertices.back(), 3);
    merged_mesh.faces.resize(num_faces.back(), 3);

    for (size_t i = 0; i < num_meshes; i++) {
        merged_mesh.vertices.block(num_vertices[i], 0,
                                   num_vertices[i + 1] - num_vertices[i], 3) =
            meshes[i].vertices;
        merged_mesh.faces.block(num_faces[i], 0,
                                num_faces[i + 1] - num_faces[i], 3) =
            meshes[i].faces.array() + num_vertices[i];
    }

    return merged_mesh;
}

IGL_Mesh marching_cubes(Sampled_Implicit& fn, const GridSpec& grid_spec) {
    ScopedTimer<> timer("Marching cube");
    size_t num_grid_pts = (grid_spec.resolution[0] + 1) *
                          (grid_spec.resolution[1] + 1) *
                          (grid_spec.resolution[2] + 1);

    Eigen::Matrix<double, Eigen::Dynamic, 1> values(num_grid_pts);
    Eigen::Matrix<double, Eigen::Dynamic, 3> grid(num_grid_pts, 3);
    for (int i = 0; i <= grid_spec.resolution[2]; i++) {
        for (int j = 0; j <= grid_spec.resolution[1]; j++) {
            for (int k = 0; k <= grid_spec.resolution[0]; k++) {
                int idx = i * (grid_spec.resolution[0] + 1) *
                              (grid_spec.resolution[1] + 1) +
                          j * (grid_spec.resolution[0] + 1) + k;
                double x = double(k) / double(grid_spec.resolution[0]) *
                               (grid_spec.bbox_max[0] - grid_spec.bbox_min[0]) +
                           grid_spec.bbox_min[0];
                double y = double(j) / double(grid_spec.resolution[1]) *
                               (grid_spec.bbox_max[1] - grid_spec.bbox_min[1]) +
                           grid_spec.bbox_min[1];
                double z = double(i) / double(grid_spec.resolution[2]) *
                               (grid_spec.bbox_max[2] - grid_spec.bbox_min[2]) +
                           grid_spec.bbox_min[2];
                grid.row(idx) << x, y, z;
                values[idx] = fn.function_at({x, y, z});
            }
        }
    }

    IGL_Mesh mesh;
    igl::copyleft::marching_cubes(
        values, grid, grid_spec.resolution[0]+1, grid_spec.resolution[1]+1,
        grid_spec.resolution[2]+1, mesh.vertices, mesh.faces);
    return mesh;
}

std::vector<IGL_Mesh> compute_arrangement(const std::vector<IGL_Mesh>& meshes) {
    ScopedTimer<> timer("mesh arrangement");
    auto merged_mesh = merge_meshes(meshes);
    auto engine = PyMesh::CellPartition::create_raw(merged_mesh.vertices, merged_mesh.faces);
    engine->run();

    size_t num_cells = engine->get_num_cells();
    std::vector<IGL_Mesh> cells(num_cells);
    for (size_t i=0; i<num_cells; i++) {
        cells[i].vertices = engine->get_vertices();
        cells[i].faces = engine->get_cell_faces(i);
    }

    return cells;
}

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
    meshes.reserve(implicit_functions.size()+1);

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
    for (size_t i=0; i<num_cells; i++) {
        std::string out_name = basename + ".comp_" + std::to_string(i) + ext;
        igl::write_triangle_mesh(out_name, cells[i].vertices, cells[i].faces);
    }

    return 0;
}
