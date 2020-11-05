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

    static IGL_Mesh generate_cube(const GridSpec& grid_spec);

    static IGL_Mesh marching_cubes(Sampled_Implicit& fn, const GridSpec& grid_spec);

    static void merge_meshes(const std::vector<IGL_Mesh>& meshes,
            // output
            IGL_Mesh &merged_mesh,
            Eigen::VectorXi &face_to_mesh);

};


#endif //PSI_MESH_PSI_H
