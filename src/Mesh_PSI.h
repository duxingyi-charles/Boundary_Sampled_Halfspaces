//
// Created by Charles Du on 11/5/20.
//

#ifndef PSI_MESH_PSI_H
#define PSI_MESH_PSI_H

#include "PSI.h"

struct IGL_Mesh {
    Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> vertices;
    Eigen::Matrix<int, Eigen::Dynamic, 3, Eigen::RowMajor> faces;
};

class Mesh_PSI : public PSI {

public:
    Mesh_PSI() : PSI() {};

    ~Mesh_PSI() override = default;

private:
    void compute_arrangement_for_graph_cut(
            const GridSpec &grid,
            const std::vector<std::unique_ptr<Sampled_Implicit>> &implicits) override;

    void compute_arrangement(
            const GridSpec &grid,
            const std::vector<std::unique_ptr<Sampled_Implicit>> &implicits);

    static IGL_Mesh generate_cube(const GridSpec& grid_spec);

    static IGL_Mesh marching_cubes(Sampled_Implicit& fn, const GridSpec& grid_spec);

    static void merge_meshes(const std::vector<IGL_Mesh>& meshes,
            // output
            IGL_Mesh &merged_mesh,
            Eigen::VectorXi &face_to_mesh);

    IGL_Mesh generate_plane(const GridSpec &grid,
                            const Point &p, const Eigen::Vector3d &normal);

    IGL_Mesh generate_cylinder(const GridSpec &grid, int n,
                               const Point& axis_point, const Eigen::Vector3d &axis_unit_vector,
                               double radius, bool is_flipped);

    IGL_Mesh generate_cone(const GridSpec &grid, int n,
                           const Point& apex, const Eigen::Vector3d& axis_unit_vector,
                           double apex_angle, bool is_flipped);
//
//    IGL_Mesh generate_sphere(const GridSpec &grid, int num_subdivide,
//                             const Point& center, double radius, bool is_flipped);
//
    IGL_Mesh generate_torus(const GridSpec &grid, int n_major, int n_minor,
                            const Point& center, const Eigen::Vector3d& axis_unit_vector,
                            double major_radius, double minor_radius,
                            bool is_flipped);

    IGL_Mesh generate_random_plane(const GridSpec &grid);

};


#endif //PSI_MESH_PSI_H
