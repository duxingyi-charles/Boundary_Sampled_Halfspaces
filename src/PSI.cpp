//
// Created by Charles Du on 11/5/20.
//

#include "PSI.h"

#include <fstream>
#include <iostream>

#include <limits>
#include <map>

// graph cut
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/boykov_kolmogorov_max_flow.hpp>

void PSI::run(const GridSpec &grid, std::vector<std::unique_ptr<Sampled_Implicit>> &implicits){
    // make a pointer to the input implicits
    Impl_ptr = &implicits;

    compute_arrangement_for_graph_cut(grid, implicits);
    ready_for_graph_cut = true;
    graph_cut();
}

void PSI::graph_cut() {
    if (ready_for_graph_cut) {
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
    std::cout << "export_data finish: " << filename << std::endl;
    return true;
}
