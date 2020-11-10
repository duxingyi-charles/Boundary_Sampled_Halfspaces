//
// Created by Charles Du on 11/5/20.
//

#include "Mesh_PSI.h"

#include <igl/copyleft/marching_cubes.h>

#include <limits>

#include <numeric> //std::partial_sum

#include "PyMesh/CellPartition.h"
#include "ScopedTimer.h"


IGL_Mesh Mesh_PSI::generate_cube(const GridSpec& grid_spec) {
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

IGL_Mesh Mesh_PSI::marching_cubes(Sampled_Implicit& fn, const GridSpec& grid_spec) {
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

void Mesh_PSI::merge_meshes(const std::vector<IGL_Mesh>& meshes,
                      // output
                      IGL_Mesh &merged_mesh,
                      Eigen::VectorXi &face_to_mesh) {
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

    merged_mesh.vertices.resize(num_vertices.back(), 3);
    merged_mesh.faces.resize(num_faces.back(), 3);

    // face id to mesh id
    face_to_mesh = Eigen::VectorXi::Ones(num_faces.back());


    for (size_t i = 0; i < num_meshes; i++) {
        merged_mesh.vertices.block(num_vertices[i], 0,
                                   num_vertices[i + 1] - num_vertices[i], 3) =
                meshes[i].vertices;
        merged_mesh.faces.block(num_faces[i], 0,
                                num_faces[i + 1] - num_faces[i], 3) =
                meshes[i].faces.array() + num_vertices[i];
        //
        face_to_mesh.segment(num_faces[i],num_faces[i+1]-num_faces[i]) *= i;
    }

}



void Mesh_PSI::compute_arrangement_for_graph_cut(
        const GridSpec &grid,
        const std::vector<std::unique_ptr<Sampled_Implicit>> &implicits)
{
    // generate meshes
    std::vector<IGL_Mesh> meshes;
    meshes.reserve(implicits.size()+1);

    meshes.push_back(generate_cube(grid));
    for (const auto& fn : implicits) {
        meshes.push_back(marching_cubes(*fn, grid));
    }


    // merge meshes
    IGL_Mesh merged_mesh;
    Eigen::VectorXi face_to_mesh;
    merge_meshes(meshes,merged_mesh,face_to_mesh);


    // compute arrangement
    ScopedTimer<> timer("mesh arrangement for graph-cut");
    auto engine = PyMesh::CellPartition::create_raw(merged_mesh.vertices, merged_mesh.faces);
    engine->run();

    auto vertices = engine->get_vertices();
    auto faces = engine->get_faces();

    //copy vertices to V
    V.clear();
    for (int i = 0; i < vertices.rows(); ++i) {
        V.emplace_back(vertices(i,0),vertices(i,1),vertices(i,2));
    }
    //copy faces to F
    F.clear();
    F.resize(faces.rows());
    for (int i = 0; i < faces.rows(); ++i) {
        for (int j = 0; j < faces.cols(); ++j) {
            F[i].push_back(faces(i,j));
        }
    }

    // compute map: result face id -> input implicit id
    auto source_faces = engine->get_source_faces();
    Eigen::VectorXi result_face_to_implicit(source_faces.size());
    for (int i=0; i<result_face_to_implicit.size(); ++i) {
        // -1: grid bounding cube
        result_face_to_implicit(i) = face_to_mesh(source_faces(i)) - 1;
    }


    // compute distance weighted area, as well as sample points of patches
    std::vector<std::vector<double>> sample_min_dist;
    std::vector<std::vector<int>> sample_nearest_patch;
    double infinity = std::numeric_limits<double>::infinity();
    for (const auto& fn : implicits) {
        int num_samples = fn->get_sample_points().size();
        sample_min_dist.emplace_back(num_samples,infinity);
        sample_nearest_patch.emplace_back(num_samples,-1);
    }

    // compute distance weighted area of patches
    auto cells = engine->get_cells();
    int num_patch = cells.rows();

    //distance weighted area of each patch
    std::vector<double> patch_dist(num_patch, 0);
    //index of implicit for each patch
    std::vector<int> patch_impl(num_patch, -1);

    auto face_to_patch = engine->get_patches();


    for (int i=0; i<F.size(); i++) {
        int implicit_id = result_face_to_implicit(i);
        if (implicit_id == -1) continue; // the current face belongs to grid cube
        int patch_id = face_to_patch(i);
        patch_impl[patch_id] = implicit_id;

        auto samples = implicits[implicit_id]->get_sample_points();

        Point p1 = V[F[i][0]];
        Point p2 = V[F[i][1]];
        Point p3 = V[F[i][2]];

        auto face_center = (p1 + p2 + p3)/3;
        // distance weighted area
        double tri_area = (p2-p1).cross(p3-p1).norm()/2;
        double min_distance = infinity;
        for (int j = 0; j < samples.size(); ++j) {
            double distance = (face_center - samples[j]).norm();
            if (distance < sample_min_dist[implicit_id][j]) {
                sample_min_dist[implicit_id][j] = distance;
                sample_nearest_patch[implicit_id][j] = patch_id;
            }
            if (distance < min_distance) {
                min_distance = distance;
            }
        }
        patch_dist[patch_id] += min_distance * tri_area;
    }

    // extract sample points on patches
    std::vector<std::vector<int>> patch_samples(num_patch);

    for (auto& nearest_patches : sample_nearest_patch) {
        for (int i = 0; i < nearest_patches.size(); ++i) {
            int nearest_patch = nearest_patches[i];
            if (nearest_patch != -1) {
                patch_samples[nearest_patch].push_back(i);
            }
        }
    }

    //blocks incident to each patch
    std::vector<std::vector<int>> patch_block;
    //sign of implicit on blocks incident to each patch
    std::vector<std::vector<int>> patch_sign;

    for (int i=0; i < num_patch; ++i) {
        patch_block.emplace_back(std::vector<int> {cells(i, 0), cells(i, 1)});
        patch_sign.emplace_back(std::vector<int> {1, -1});
    }

    // remove patches from the grid cube
    std::vector<int> P_old_to_new_index(patch_impl.size());

    P_dist.clear();
    P_samples.clear();
    P_block.clear();
    P_sign.clear();
    P_Impl.clear();

    int patch_count = 0;
    for (int i=0; i < patch_impl.size(); ++i) {
        if (patch_impl[i] != -1) {
            P_old_to_new_index[i] = patch_count;
            P_dist.push_back(patch_dist[i]);
            P_samples.emplace_back(patch_samples[i]);
            P_block.emplace_back(patch_block[i]);
            P_sign.emplace_back(patch_sign[i]);
            P_Impl.push_back(patch_impl[i]);
            patch_count++;
        }
        else {
            P_old_to_new_index[i] = -1; //patch deleted
        }
    }
    num_patch = P_Impl.size();

    // P_block: make sure the minimal index is 0
    int min_block_index = std::numeric_limits<int>::max();
    int max_block_index = -1;
    for (auto &blocks : P_block) {
        for (auto b : blocks) {
            min_block_index = (b < min_block_index) ? b : min_block_index;
            max_block_index = (b > max_block_index) ? b : max_block_index;
        }
    }
    if (min_block_index > 0) {
        for (auto &blocks : P_block) {
            for (auto &b : blocks) {
                b -= min_block_index;
            }
        }
    }
    int num_block = max_block_index - min_block_index + 1;


    // collect faces belonging to each patches
    P.clear();
    P.resize(num_patch);
    F_Impl.clear();
    F_Impl.resize(F.size(), -1);  // -1 for faces on grid cube
    for (int i=0; i<F.size(); ++i) {
        int patch_id = P_old_to_new_index[face_to_patch(i)];
        if (patch_id != -1) {
            P[patch_id].push_back(i);
            F_Impl[i] = P_Impl[patch_id];
        }
    }

    // collect patches belonging to each block
    B_patch.clear();
    B_patch.resize(num_block);
    for (int i=0; i < P_block.size(); ++i) {
        for (int j=0; j < P_block[i].size(); ++j) {
            B_patch[P_block[i][j]].push_back(i);
        }
    }

}