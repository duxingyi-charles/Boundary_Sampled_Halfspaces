#include <igl/copyleft/marching_cubes.h>
#include <igl/write_triangle_mesh.h>

#include <CLI/CLI.hpp>
#include <Eigen/Core>
#include <array>
#include <string>

#include <limits>
#include <Eigen/Geometry>

#include "PyMesh/Arrangement.h"
#include "config.h"
#include "ScopedTimer.h"

#include "GridSpec.h"

// graph cut
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/boykov_kolmogorov_max_flow.hpp>

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
    PyMesh::VectorI face_labels(merged_mesh.faces.rows());
    face_labels.setZero(); // Use dummy labels.
    auto engine = PyMesh::Arrangement::create_mesh_arrangement(merged_mesh.vertices, merged_mesh.faces, face_labels);
    engine->run();

    size_t num_cells = engine->get_num_cells();
    std::vector<IGL_Mesh> cells(num_cells);
    for (size_t i=0; i<num_cells; i++) {
        cells[i].vertices = engine->get_vertices();
        cells[i].faces = engine->get_cell_faces(i);
    }

    return cells;
}

//
void prepare_graph_cut(const std::vector<std::unique_ptr<Sampled_Implicit>> &implicit_functions,
                       const std::vector<IGL_Mesh>& meshes,
                       // output
                       std::vector<Point> &V,
                       std::vector<std::vector<int>> &F,
                       std::vector<int> &F_Impl,
                       std::vector<std::vector<int>> &P,
                       std::vector<std::vector<int>> &B_patch,
                       std::vector<int> &P_Impl,
                       // util: for latter use in graph cut
                       std::vector<double> &P_dist,
                       std::vector<std::vector<int>> &P_samples,
                       std::vector<std::vector<int>> &P_block,
                       std::vector<std::vector<int>> &P_sign
                       ) {
    // merge meshes
//    ScopedTimer<> timer("Merge meshes");
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

    // face id to mesh id
    Eigen::VectorXi face_to_mesh = Eigen::VectorXi::Ones(num_faces.back());


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

    // compute arrangement
//    ScopedTimer<> timer("mesh arrangement");
    auto engine = PyMesh::Arrangement::create_mesh_arrangement(merged_mesh.vertices, merged_mesh.faces, face_to_mesh);
    engine->run();

    auto vertices = engine->get_vertices();
    auto faces = engine->get_faces();

    //vertices
    V.clear();
    for (int i = 0; i < vertices.rows(); ++i) {
        V.emplace_back(vertices(i,0),vertices(i,1),vertices(i,2));
    }
    //faces
    F.clear();
    F.resize(faces.rows());
    for (int i = 0; i < faces.rows(); ++i) {
        for (int j = 0; j < faces.cols(); ++j) {
            F[i].push_back(faces(i,j));
        }
    }

    // compute map: result face id -> input implicit id
    Eigen::VectorXi result_face_to_implicit = engine->get_out_face_labels().array()-1;


    // compute distance weighted area, as well as sample points of patches
    std::vector<std::vector<double>> sample_min_dist;
    std::vector<std::vector<int>> sample_nearest_patch;
    double infinity = std::numeric_limits<double>::infinity();
    for (const auto& fn : implicit_functions) {
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


    for (int i=0; i<faces.rows(); i++) {
        int implicit_id = result_face_to_implicit(i);
        if (implicit_id == -1) continue; // the current face belongs to grid cube
        int patch_id = face_to_patch(i);
        patch_impl[patch_id] = implicit_id;

        auto samples = implicit_functions[implicit_id]->get_sample_points();


        Point p1(vertices(faces(i,0),0),vertices(faces(i,0),1),vertices(faces(i,0),2));
        Point p2(vertices(faces(i,1),0),vertices(faces(i,1),1),vertices(faces(i,1),2));
        Point p3(vertices(faces(i,2),0),vertices(faces(i,2),1),vertices(faces(i,2),2));

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

//    std::vector<double> P_dist;
//    std::vector<std::vector<int>> P_samples;
//    std::vector<std::vector<int>> P_block;
//    std::vector<std::vector<int>> P_sign;
//    std::vector<int> P_Impl;
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
    for (int i=0; i<P_block.size(); ++i) {
        for (int j=0; j<P_block[i].size(); ++j) {
            min_block_index = (P_block[i][j] < min_block_index) ? P_block[i][j] : min_block_index;
            max_block_index = (P_block[i][j] > max_block_index) ? P_block[i][j] : max_block_index;
        }
    }
    if (min_block_index > 0) {
        for (int i=0; i<P_block.size(); ++i) {
            for (int j=0; j<P_block[i].size(); ++j) {
                P_block[i][j] -= min_block_index;
            }
        }
    }
    int num_block = max_block_index - min_block_index + 1;


    // collect faces belonging to each patches
//    std::vector<std::vector<int>> P(num_patch);
    P.clear();
    P.resize(num_patch);
//    std::vector<int> F_Impl(faces.rows(),-1); // -1 for faces on grid cube
    F_Impl.clear();
    F_Impl.resize(faces.rows(), -1);  // -1 for faces on grid cube
    for (int i=0; i<faces.rows(); ++i) {
        int patch_id = P_old_to_new_index[face_to_patch(i)];
        if (patch_id != -1) {
            P[patch_id].push_back(i);
//            F_Impl[i] = patch_impl[patch_id];
            F_Impl[i] = P_Impl[patch_id];
        }
    }

    // collect patches belonging to each block
    //(unordered) list of boundary patches of each block
//    int num_block = engine->get_num_cells();
//    num_block += 1; //some bug here
//    std::vector<std::vector<int>> B_patch(num_block);
    B_patch.clear();
    B_patch.resize(num_block);
    for (int i=0; i < P_block.size(); ++i) {
        for (int j=0; j < P_block[i].size(); ++j) {
            B_patch[P_block[i][j]].push_back(i);
        }
    }

}

void graph_cut(const std::vector<double> &P_dist,
               const std::vector<std::vector<int>> &P_samples,
               const std::vector<std::vector<int>> &P_block,
               const std::vector<std::vector<int>> &P_sign,
               const std::vector<std::vector<int>> &B_patch,
               // output
               // block labels: object -> true, background -> false
               std::vector<bool> &B_label,
               // patch labels: surface -> true, not surface -> false
               std::vector<bool> &P_label
               ) {
    typedef std::pair<int,int> Edge;

    // constants
    double inf = std::numeric_limits<double>::infinity();
    double Delta = 0;
    for (auto d : P_dist) Delta += d;
    Delta *= 2;

    // define per-cell costs: hPos, hNeg
    int nBlock = B_patch.size();
    std::vector<double> hPos(nBlock);
    std::vector<double> hNeg(nBlock);

    int nPatch = P_samples.size();
    for (int p = 0; p < nPatch; ++p) {
        double cost = Delta * P_samples[p].size();
        //
        auto b = P_block[p][0];
        if (P_sign[p][0] > 0) hNeg[b] += cost;
        else hPos[b] += cost;
        //
        b = P_block[p][1];
        if (P_sign[p][1] > 0) hNeg[b] += cost;
        else hPos[b] += cost;
    }

    // define pair-cell costs: hPair
    std::map<Edge, double> hPair;
    // init hPair
    for (int p = 0; p < nPatch; ++p) {
        int b1 = P_block[p][0];
        int b2 = P_block[p][1];
        hPair[Edge(b1,b2)] = 0;
        hPair[Edge(b2,b1)] = 0;
    }
    // compute hPair
    for (int p = 0; p < nPatch; ++p) {
        int b1 = P_block[p][0];
        int b2 = P_block[p][1];
        //
        if (P_sign[p][0] == 1) {
            if (hPair[Edge(b1,b2)] != inf) {
                hPair[Edge(b1,b2)] += P_dist[p];
            }
        }
        else {
            hPair[Edge(b1,b2)] = inf;
        }
        //
        if (P_sign[p][1] == 1) {
            if (hPair[Edge(b2,b1)] != inf) {
                hPair[Edge(b2,b1)] += P_dist[p];
            }
        }
        else {
            hPair[Edge(b2,b1)] = inf;
        }
    }

    // create the graph
    typedef boost::adjacency_list_traits<boost::vecS, boost::vecS, boost::directedS> Traits;
    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS,
    boost::property<boost::vertex_color_t, boost::default_color_type,
            boost::property<boost::vertex_predecessor_t,Traits::edge_descriptor,
            boost::property<boost::vertex_distance_t, double,
            boost::property<boost::vertex_index_t, long> > >
    >,

    boost::property<boost::edge_capacity_t, double,
            boost::property<boost::edge_residual_capacity_t, double,
            boost::property<boost::edge_reverse_t, Traits::edge_descriptor > > > >
                                                   Graph;

    Graph g(nBlock);

    boost::property_map<Graph,boost::edge_capacity_t>::type
            e_weights = get(boost::edge_capacity,g);
    boost::property_map<Graph,boost::edge_reverse_t>::type
            e_reverse = get(boost::edge_reverse,g);

    // add edges between blocks
    for (int p = 0; p < nPatch; ++p) {
        int b1 = P_block[p][0];
        int b2 = P_block[p][1];
        auto e = add_edge(b1,b2,g).first;
        e_weights[e] = hPair[Edge(b1,b2)];
        auto re = add_edge(b2,b1,g).first;
        e_weights[re] = hPair[Edge(b2,b1)];
        e_reverse[e] = re;
        e_reverse[re] = e;
    }

    // add edges between blocks and terminals (s,t)
    auto sid = add_vertex(g);
    auto tid = add_vertex(g);
    for (int b = 0; b < nBlock; ++b) {
        auto e = add_edge(sid,b,g).first;
        e_weights[e] = hNeg[b];
        auto re = add_edge(b,sid,g).first;
        e_weights[re] = 0;
        e_reverse[e] = re;
        e_reverse[re] = e;
        //
        e = add_edge(b,tid,g).first;
        e_weights[e] = hPos[b];
        re = add_edge(tid,b,g).first;
        e_weights[re] = 0;
        e_reverse[e] = re;
        e_reverse[re] = e;
    }

    // max-flow-min-cut
    /*double flow = */boykov_kolmogorov_max_flow(g ,sid, tid);

    // print max-flow result
//    std::cout << "c  The total flow:" << std::endl;
//    std::cout << "s " << flow << std::endl << std::endl;
//
//    std::cout << "c flow values:" << std::endl;
//    boost::graph_traits<Graph>::vertex_iterator u_iter, u_end;
//    boost::graph_traits <Graph>::out_edge_iterator ei, e_end;
//    for (boost::tie(u_iter, u_end) = vertices(g); u_iter != u_end; ++u_iter)
//        for (boost::tie(ei, e_end) = out_edges(*u_iter, g); ei != e_end; ++ei)
//            if (e_weights[*ei] > 0)
//                std::cout << "f " << *u_iter << " " << target(*ei, g) << " "
//                << (e_weights[*ei]) << " "
//                << (e_weights[*ei] - e_residual[*ei]) << std::endl;

    // print block labels
//    std::cout << "vertex color:" << std::endl;
//    boost::property_map<Graph ,boost::vertex_color_t>::type v_color = get(boost::vertex_color,g);
//    for (boost::tie(u_iter, u_end) = vertices(g); u_iter != u_end; ++u_iter) {
//        std::cout << *u_iter << " " << v_color[*u_iter] << std::endl;
//    }
    //

    // get block and patch labels
    boost::property_map<Graph,boost::vertex_color_t>::type
            block_labels = get(boost::vertex_color,g);
    auto s_label = block_labels[sid];

    B_label.clear();
    P_label.clear();

    B_label.resize(nBlock);
    for (int b = 0; b < nBlock; ++b) {
        B_label[b] = (block_labels[b] == s_label);
    }

    P_label.resize(nPatch);
    for (int p = 0; p < nPatch; ++p) {
        P_label[p] = (B_label[P_block[p][0]] != B_label[P_block[p][1]]);
    }

}

bool export_grid(const std::string &filename,
                 const std::vector<std::unique_ptr<Sampled_Implicit>> &implicit_functions,
                 const std::vector<Point> &V,
                 const std::vector<std::vector<int>> &F,
                 const std::vector<int> &F_Impl,
                 const std::vector<std::vector<int>> &P,
                 const std::vector<std::vector<int>> &B_patch,
                 const std::vector<int> &P_Impl,
                 const std::vector<bool> &B_label,
                 const std::vector<bool> &P_label
)
{
    std::ofstream fout(filename, std::ofstream::out);
    if (!fout.good()) {
        std::cout << "Can not create output file " << filename << std::endl;
        return false;
    }

    // V
    fout << "vert ";
    fout << V.size() << " " << V[0].size() << std::endl;
    for (auto &v : V) {
        fout << v.transpose() << std::endl;
    }

    // E
//    fout << "edge ";
//    fout << E.size() << " " << "2" << std::endl;
//    for (auto &e : E) {
//        fout << e.first << " " << e.second << std::endl;
//    }

    // F
    fout << "face ";
    fout << F.size() << std::endl;
    for (auto &f : F) {
        for (auto  &e : f) {
            fout << e << " ";
        }
        fout << std::endl;
    }

    // C
//    fout << "cell ";
//    fout << C.size() << std::endl;
//    for (auto &c : C) {
//        for (auto &f : c) {
//            fout << f << " ";
//        }
//        fout << std::endl;
//    }

    // V_Impl
//    fout << "vert_implicit ";
//    fout << V_Impl.size() << std::endl;
//    for (auto &v_impl : V_Impl) {
//        for (auto &impl : v_impl) {
//            fout << impl << " ";
//        }
//        fout << std::endl;
//    }

    // E_Impl
//    fout << "edge_implicit ";
//    fout << E_Impl.size() << std::endl;
//    for (auto &e_impl : E_Impl) {
//        for (auto &impl : e_impl) {
//            fout << impl << " ";
//        }
//        fout << std::endl;
//    }

    // F_Impl
    fout << "face_implicit ";
    fout << 1 << std::endl; // row vector
    for (auto &impl : F_Impl) {
        fout << impl << " ";
    }
    fout << std::endl;

    // P
    fout << "patch_faces ";
    fout << P.size() << std::endl;
    for (auto &patch : P) {
        for (auto & face : patch) {
            fout << face << " ";
        }
        fout << std::endl;
    }

    // P_Impl
    fout << "patch_implicit ";
    fout << 1 << std::endl; // row vector
    for (auto & impl: P_Impl) {
        fout << impl << " ";
    }
    fout << std::endl;

    // samples
    fout << "samples ";
    fout << implicit_functions.size() << std::endl;
    for (const auto &impl : implicit_functions) {
        for (auto &point : impl->get_sample_points()) {
            fout << point.transpose() << " ";
        }
        fout << std::endl;
    }

    // P_samples
//    fout << "patch_samples ";
//    fout << P_samples.size() << std::endl;
//    for (auto &samples : P_samples) {
//        for (auto & sample : samples) {
//            fout << sample << " ";
//        }
//        fout << std::endl;
//    }

    // P_dist
//    fout << "patch_distance_area ";
//    fout << 1 << std::endl; // row vector
//    for (auto & dist: P_dist) {
//        fout << dist << " ";
//    }
//    fout << std::endl;

    // B_patch
    fout << "block_patches ";
    fout << B_patch.size() << std::endl;
    for (auto &patches : B_patch) {
        for (auto & patch : patches) {
            fout << patch << " ";
        }
        fout << std::endl;
    }

    // B_cell
//    fout << "block_cells ";
//    fout << B_cell.size() << std::endl;
//    for (auto &cells : B_cell) {
//        for (auto & cell : cells) {
//            fout << cell << " ";
//        }
//        fout << std::endl;
//    }

    // P_block
//    fout << "patch_blocks ";
//    fout << P_block.size() << std::endl;
//    for (auto &blocks : P_block) {
//        for (auto & block : blocks) {
//            fout << block << " ";
//        }
//        fout << std::endl;
//    }

    // P_sign
//    fout << "patch_signs ";
//    fout << P_sign.size() << std::endl;
//    for (auto &signs : P_sign) {
//        for (auto & sign : signs) {
//            fout << sign << " ";
//        }
//        fout << std::endl;
//    }

//  P_label
    fout << "patch_labels ";
    fout << 1 << std::endl; // row vector
    for (auto label : P_label) {
        fout << label << " ";
    }
    fout << std::endl;


    // B_label
    fout << "block_labels ";
    fout << 1 << std::endl; // row vector
    for (auto label : B_label) {
        fout << label << " ";
    }
    fout << std::endl;


    fout.close();
    std::cout << "export_grid finish: " << filename << std::endl;
    return true;
}



int main(int argc, char** argv) {
    struct {
        std::string config_file;
//        std::string output_mesh;
        std::string grid_file;
        std::string output_grid_file;
    } args;

    CLI::App app{"Piecewise implicit surface demo"};
    app.add_option("-G,--grid", args.grid_file, "Grid spec file")->required();
    app.add_option("config_file", args.config_file, "Configuration file")
        ->required();
//    app.add_option("output_mesh", args.output_mesh, "Output mesh file")
//        ->required();
    app.add_option("output_grid_file", args.output_grid_file,
                   "Output grid file")
            ->required();
    CLI11_PARSE(app, argc, argv);

    auto implicit_functions =
        initialize_sampled_implicit_functions(args.config_file);
    auto grid_spec = GridSpec::parse_grid_spec(args.grid_file);



    std::vector<IGL_Mesh> meshes;
    meshes.reserve(implicit_functions.size()+1);

    meshes.push_back(generate_cube(grid_spec));
    for (const auto& fn : implicit_functions) {
        meshes.push_back(marching_cubes(*fn, grid_spec));
    }

//    auto dot_pos = args.output_mesh.find_last_of(".");
//    auto basename = args.output_mesh.substr(0, dot_pos);
//    auto ext = args.output_mesh.substr(dot_pos);

//    size_t num_meshes = meshes.size();
//    for (size_t i = 0; i < num_meshes; i++) {
//        std::string out_name = basename + ".input_" + std::to_string(i) + ext;
//        igl::write_triangle_mesh(out_name, meshes[i].vertices, meshes[i].faces);
//    }


//    auto cells = compute_arrangement(meshes);
//    auto num_cells = cells.size();
//    for (size_t i=0; i<num_cells; i++) {
//        std::string out_name = basename + ".comp_" + std::to_string(i) + ext;
//        igl::write_triangle_mesh(out_name, cells[i].vertices, cells[i].faces);
//    }

    // prepare for graph cut
    std::vector<Point> V;
    std::vector<std::vector<int>> F;
    std::vector<int> F_Impl;
    std::vector<std::vector<int>> P;
    std::vector<std::vector<int>> B_patch;
    std::vector<int> P_Impl;
    // util: for latter use in graph cut
    std::vector<double> P_dist;
    std::vector<std::vector<int>> P_samples;
    std::vector<std::vector<int>> P_block;
    std::vector<std::vector<int>> P_sign;

    prepare_graph_cut(implicit_functions, meshes,
    // output
    V, F, F_Impl,P, B_patch,P_Impl,
    // util: for latter use in graph cut
    P_dist,P_samples,P_block,P_sign
    );

    // graph cut
    // block labels: object -> true, background -> false
    std::vector<bool> B_label;
    // patch labels: surface -> true, not surface -> false
    std::vector<bool> P_label;
    graph_cut(P_dist,P_samples,P_block,P_sign,B_patch,
    // output
   B_label,P_label);

    // export
    export_grid(args.output_grid_file,
                     implicit_functions,
                     V,
                     F,
                     F_Impl,
                     P,
                     B_patch,
                     P_Impl,
                     B_label,
                     P_label
    );



    return 0;
}
