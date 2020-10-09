#include <Eigen/Core>
#include "config.h"

struct IGL_Mesh {
    Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> vertices;
    Eigen::Matrix<int, Eigen::Dynamic, 3, Eigen::RowMajor> faces;
};

IGL_Mesh generate_cube(const GridSpec& grid_spec);

IGL_Mesh merge_meshes(const std::vector<IGL_Mesh>& meshes);

IGL_Mesh marching_cubes(Sampled_Implicit& fn, const GridSpec& grid_spec);

std::vector<IGL_Mesh> compute_arrangement(const std::vector<IGL_Mesh>& meshes);

