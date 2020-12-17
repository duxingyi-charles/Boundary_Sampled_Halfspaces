//
// Created by Charles Du on 11/5/20.
//

#include "PSI.h"

#include <fstream>
#include <iostream>

#include <limits>
#include <map>

#include <cmath>

#include "ScopedTimer.h"

// cross product
#include <Eigen/Geometry>

// graph cut
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/boykov_kolmogorov_max_flow.hpp>

void PSI::run(const GridSpec &grid, std::vector<std::unique_ptr<Sampled_Implicit>> &implicits){
    // make a pointer to the input implicits
    Impl_ptr = &implicits;

    compute_arrangement_for_graph_cut(grid, implicits);
    arrangement_ready = true;
    ready_for_graph_cut = true;
    graph_cut();
}

void PSI::process_samples() {
    if (arrangement_ready) {
        process_samples(V,F,P,P_Impl,Impl_ptr,use_distance_weighted_area,P_samples,P_dist);
        ready_for_graph_cut = true;
    } else {
        std::cout << "Error: you should compute arrangement before handling samples." << std::endl;
    }
}

// static process_samples function
void PSI::process_samples(
        //input
        const std::vector<Point> &V,
        const std::vector<std::vector<int>> &F,
        const std::vector<std::vector<int>> &P,
        const std::vector<int> &P_Impl,
        const std::vector<std::unique_ptr<Sampled_Implicit>> *Impl_ptr,
        bool use_distance_weighted_area,
        //output
        std::vector<std::vector<int>> &P_samples,
        std::vector<double> &P_dist
){
    // compute distance weighted area, as well as sample points of patches
    std::vector<std::vector<double>> sample_min_dist;
    std::vector<std::vector<int>> sample_nearest_patch;
    double infinity = std::numeric_limits<double>::infinity();
    for (const auto &impl : (*Impl_ptr)) {
        int num_samples = impl->get_sample_points().size();
        sample_min_dist.emplace_back(num_samples,infinity);
        sample_nearest_patch.emplace_back(num_samples,-1);
    }

    // compute distance weighted area of patches
    P_dist.clear();
    P_dist.resize(P.size(), 0);

    for (int i = 0; i < P.size(); ++i) {
        auto &patch = P[i];
        int impl = P_Impl[i];
        auto samples = (*Impl_ptr)[impl]->get_sample_points();
        for (int f : patch) {
            auto &face = F[f];
            // compute center of the face
            Point face_center(0,0,0);
            for (int v : face) {
                face_center += V[v];
            }
            face_center /= face.size();
            // compute distance between face center and nearest sample point
            double min_distance = infinity;
            for (int j = 0; j < samples.size(); ++j) {
                double distance = (face_center - samples[j]).norm();
                if (distance < sample_min_dist[impl][j]) {
                    sample_min_dist[impl][j] = distance;
                    sample_nearest_patch[impl][j] = i;
                }
                if (distance < min_distance) {
                    min_distance = distance;
                }
            }

            // compute distance weighted area of the face
            double weighted_area = 0;
            if (face.size() == 3) {  //triangle
                const Point &p1 = V[face[0]];
                const Point &p2 = V[face[1]];
                const Point &p3 = V[face[2]];
                double tri_area = (p2-p1).cross(p3-p1).norm()/2;
                if (use_distance_weighted_area) {
                    weighted_area = min_distance * tri_area;
                } else {
                    weighted_area = tri_area;
                }

            }
            else {   // general polygon
                for (size_t vi=0; vi < face.size(); ++vi) {
                    size_t vj = (vi + 1) % face.size();
                    const Point &p1 = V[face[vi]];
                    const Point &p2 = V[face[vj]];
                    double area = ((p1 - face_center).cross(p2 - face_center)).norm() / 2;
                    if (use_distance_weighted_area) {
                        weighted_area += min_distance * area;
                    } else {
                        weighted_area += area;
                    }
                }
            }
            //
            P_dist[i] += weighted_area;
        }
    }
    for (auto &d : P_dist) {
        if (!isfinite((d))) d = infinity;
    }

    // extract sample points on patches
    P_samples.clear();
    P_samples.resize(P.size());
    for (auto& nearest_patches : sample_nearest_patch) {
        for (int i = 0; i < nearest_patches.size(); ++i) {
            int nearest_patch = nearest_patches[i];
            if (nearest_patch != -1) {
                P_samples[nearest_patch].push_back(i);
            }
        }
    }
}


void PSI::graph_cut() {
    if (ready_for_graph_cut) {
        ScopedTimer<> timer("graph cut");
        graph_cut(P_dist,P_samples,P_block,P_sign,B_patch,
                  B_label,P_label);
        graph_cut_finished = true;
    } else {
        graph_cut_finished = false;
    }
}

// static graph_cut function
void PSI::graph_cut(
        // input
        const std::vector<double> &P_dist,
        const std::vector<std::vector<int>> &P_samples,
        const std::vector<std::vector<int>> &P_block,
        const std::vector<std::vector<int>> &P_sign,
        const std::vector<std::vector<int>> &B_patch,
        // output
        std::vector<bool> &B_label,
        std::vector<bool> &P_label) {
    // constants
    double inf = std::numeric_limits<double>::infinity();
    double Delta = 0;
//    for (auto d : P_dist) Delta += d;
    for (auto &d: P_dist) {
        if (isfinite(d)) Delta += d;
    }
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
    /*double flow =*/ boykov_kolmogorov_max_flow(g ,sid, tid);

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
//                << (e_weights[*ei]) << " " << std::endl;
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

bool PSI::export_data(const std::string &filename) const
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

    // Face_edge
//    fout << "face_edge ";
//    fout << Face_edges.size() << std::endl;
//    for (auto &f : Face_edges) {
//        for (auto  &e : f) {
//            fout << e << " ";
//        }
//        fout << std::endl;
//    }

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
    fout << Impl_ptr->size() << std::endl;
    for (const auto &impl : (*Impl_ptr)) {
        for (auto &point : impl->get_sample_points()) {
            fout << point.transpose() << " ";
        }
        fout << std::endl;
    }

    // P_samples
    fout << "patch_samples ";
    fout << P_samples.size() << std::endl;
    for (auto &samples : P_samples) {
        for (auto & sample : samples) {
            fout << sample << " ";
        }
        fout << std::endl;
    }

    // P_dist
    fout << "patch_distance_area ";
    fout << 1 << std::endl; // row vector
    for (auto & dist: P_dist) {
        fout << dist << " ";
    }
    fout << std::endl;

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
    std::cout << "export_data finish: " << filename << std::endl;
    return true;
}


// Algorithms for reverse engineering

void PSI::reduce_samples() {
    // should first run PSI on dense samples
    if (!graph_cut_finished) {
        return;
    }
    //
    // Note: current implementation is only suitable for non-distance-weighted patch area
    // If the patch area is not weighted by distance to samples,
    // then there is no need to update P_dist

    // init samples
    int num_patch = P_samples.size();
    std::vector<std::vector<int>> P_samples_dense = P_samples;
    std::vector<std::vector<int>> P_samples_sparse(num_patch);  // no samples at all

    std::vector<bool> B_label_sparse;
    std::vector<bool> P_label_sparse;

    graph_cut(P_dist,P_samples_sparse,P_block,P_sign,B_patch,B_label_sparse,P_label_sparse);
    int num_diff_label = 0;
    for (int p=0; p<num_patch; ++p) {
        if (P_label[p] != P_label_sparse[p]) {
            ++num_diff_label;
        }
    }

    int num_iter = 0;
    while (num_diff_label != 0) {
        ++num_iter;
        std::cout << "iter " << num_iter << ": " << num_diff_label << " different patch labels." << std::endl;

        // find missing patch with largest area and some samples on it
        double max_area = -1;
        double max_area_patch = -1;
        for (int p = 0; p < num_patch; ++p) {
            if (P_label[p] && !P_label_sparse[p] && !P_samples_dense[p].empty()) {
                if (P_dist[p] > max_area) {
                    max_area = P_dist[p];
                    max_area_patch = p;
                }
            }
        }

        if (max_area_patch == -1) { // no patch found but graph-cut result is not the same
            std::cout << "sample reduction failed." << std::endl;
            return;
        }

        // find the most centered sample on the patch
        double min_squared_dist = std::numeric_limits<double>::infinity();
        double min_i = -1;

        int impl = P_Impl[max_area_patch];
        const auto &samples = (*Impl_ptr)[impl]->get_sample_points();

        const auto &sampleIds = P_samples_dense[max_area_patch];
        for (int i = 0; i < sampleIds.size(); ++i) {
            double square_dist_i = 0;
            for (int j = 0; j < sampleIds.size(); ++j) {
                square_dist_i += (samples[sampleIds[i]] - samples[sampleIds[j]]).squaredNorm();
            }
            if (square_dist_i < min_squared_dist) {
                min_squared_dist = square_dist_i;
                min_i = i;
            }
        }

        // move the most centered sample to sparse samples
        P_samples_sparse[max_area_patch].push_back(P_samples_dense[max_area_patch][min_i]);
        P_samples_dense[max_area_patch].erase(P_samples_dense[max_area_patch].begin() + min_i);

        // re-compute graph-cut
        graph_cut(P_dist,P_samples_sparse,P_block,P_sign,B_patch,B_label_sparse,P_label_sparse);
        num_diff_label = 0;
        for (int p=0; p<num_patch; ++p) {
            if (P_label[p] != P_label_sparse[p]) {
                ++num_diff_label;
            }
        }
    }

    if (num_diff_label == 0) {
        P_samples = P_samples_sparse;
    } else {
        std::cout << "sample reduction failed." << std::endl;
    }

}




