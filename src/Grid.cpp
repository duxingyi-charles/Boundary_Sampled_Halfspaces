//
// Created by Charles Du on 8/26/20.
//

#include "Grid.h"

#include <map>
#include <queue>
#include <set>

#include <fstream>
#include <iostream>

#include <limits>
#include <Eigen/Geometry>

// graph cut
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/boykov_kolmogorov_max_flow.hpp>


Grid::Grid(const Point &p_min, const Point &p_max, int n_cell_x, int n_cell_y, int n_cell_z) {
    init_grid(p_min, p_max, n_cell_x, n_cell_y, n_cell_z);
}

Grid::Grid(const Point &p_min, const Point &p_max) {
    init_grid(p_min, p_max, 1, 1, 1);
}




void Grid::compute_arrangement(const Sampled_Implicit &sImplicit) {
    // record the implicit
    int cur_Impl = Impl.size();
    Impl.push_back(&sImplicit);

    // initialize
    std::vector<Point> Vertices = V;
    std::vector<Edge>  Edges = E;
    std::vector<std::vector<int>> Faces = F;
    std::vector<std::vector<int>> Cells = C;

//    std::vector<std::vector<int>> Vert_Implicits = V_Impl;
    std::vector<std::vector<int>> Edge_Implicits = E_Impl;
    std::vector<int> Face_Implicit = F_Impl;

    std::vector<int> E_vert(E.size(), -1); // -1 for null
    std::vector<std::pair<int,int>> E_new(E.size());
    std::vector<bool> is_split_edge(E.size(), false);

    std::vector<std::vector<int>> F_edge(F.size());
    std::vector<std::vector<int>> F_new(F.size());

    std::vector<bool> C_split(C.size(), false);

    // step 1: get implicit function value at each vertex
    double zero_threshold = 1e-8;
    std::vector<double> V_value(V.size());
    for (int i = 0; i < V.size(); ++i) {
        V_value[i] = sImplicit.function_at(V[i]);
        // if function at vertex is close to zero, slightly perturb it away from zero
        if (abs(V_value[i])<zero_threshold) {
            V_value[i] = (V_value[i] >=0 ? zero_threshold : -zero_threshold);
        }
    }

    // step 2: intersect edges with 0 level-set of the implicit function
    for (int e = 0; e < E.size(); ++e) {
        int vi = E[e].first;
        int vj = E[e].second;

        if (V_value[vi] * V_value[vj] < 0) {
            // linear interpolate
            double s = V_value[vi] / (V_value[vi] - V_value[vj]);
            Point v = (1-s) * V[vi] + s * V[vj];
            int v_index = Vertices.size();
            Vertices.push_back(v);
            if (E_Impl[e].size() == 1 && E_Impl[e][0] == -1) { // e is a grid edge
                V_Impl.emplace_back(std::vector<int>{cur_Impl});
            } else {
                V_Impl.emplace_back(E_Impl[e]);
                V_Impl.back().push_back(cur_Impl);
            }
            // split edge
            Edges.emplace_back(vi, v_index);
            Edges.emplace_back(v_index, vj);
            is_split_edge.push_back(false);
            is_split_edge.push_back(false);
            Edge_Implicits.emplace_back(E_Impl[e]);
            Edge_Implicits.emplace_back(E_Impl[e]);
            // record
            E_vert[e] = v_index;
            E_new[e].first  = Edges.size() - 1;
            E_new[e].second = Edges.size() - 2;
        }
    }

    // step 3: split faces with 0 level-set of the implicit function

    auto orient_edge = [&] (int e, int e_neighbor) {
        int e2 = E_new[e].second;
        if (Edges[e2].first  == Edges[e_neighbor].first  ||
            Edges[e2].first  == Edges[e_neighbor].second ||
            Edges[e2].second == Edges[e_neighbor].first  ||
            Edges[e2].second == Edges[e_neighbor].second)
        {   // if e2 intersects e_neighbor, reverse the pair E_new[e]
            E_new[e].second = E_new[e].first;
            E_new[e].first  = e2;
        }
    };

    for (int f = 0; f < F.size(); ++f) {
        int v_start = -1;
        int v_end   = -1;
        std::vector<int> remain_face;    // the remaining face
        std::vector<int> split_face;  // the split-off face
        //  1: side of the remaining face
        // -1: side of the split-off face
        int side = 1;
        bool is_split = false;  // is face f split

        const auto &face = F[f];
        int n_edge = face.size();
        for (int i = 0; i < n_edge; ++i) {
            int e = face[i];
            if (E_vert[e] != -1) { // edge e is split
                is_split = true;
                if (side == 1) {
                    // start a new split face
                    v_start = E_vert[e];
                    int e_next = face[(i+1)%n_edge];
                    orient_edge(e, e_next);
                    remain_face.push_back(E_new[e].second);
                    split_face.push_back(E_new[e].first);
                    side = -1;
                }
                else { // side == -1
                    // end a split face
                    v_end = E_vert[e];
                    int e_prev = face[(i+n_edge-1)%n_edge];
                    orient_edge(e, e_prev);
                    split_face.push_back(E_new[e].first);
                    //
                    int e_split = Edges.size();
                    Edges.emplace_back(v_start, v_end);
                    if (F_Impl[f] == -1) {
                        Edge_Implicits.emplace_back(std::vector<int>{cur_Impl});
                    }
                    else {
                        Edge_Implicits.emplace_back(std::vector<int>{F_Impl[f],cur_Impl});
                    }
                    is_split_edge.push_back(true);
                    F_edge[f].push_back(e_split);
                    split_face.push_back(e_split);
                    //
                    int f_split = Faces.size();
                    Faces.push_back(split_face);
                    Face_Implicit.push_back(F_Impl[f]);
                    F_new[f].push_back(f_split);
                    split_face.clear();
                    //
                    remain_face.push_back(e_split);
                    remain_face.push_back(E_new[e].second);
                    //
                    side = 1;
                }
            }
            else { // edge e is not split
                if (side == 1) {
                    remain_face.push_back(e);
                }
                else {
                    split_face.push_back(e);
                }
            }
        }
        // if face f is split, record the remaining face
        if (is_split) {
            int f_remain = Faces.size();
            Faces.push_back(remain_face);
            Face_Implicit.push_back(F_Impl[f]);
            F_new[f].push_back(f_remain);
        }
    }

    // step 4: split cells

    auto extract_loops = [&] (const std::vector<int> &edges) -> std::vector<std::vector<int>> {
        // create map: vert -> incident edges
        std::map<int,std::vector<int>> vert_to_edge;
        for (auto e : edges) {
            vert_to_edge[Edges[e].first].push_back(e);
            vert_to_edge[Edges[e].second].push_back(e);
        }
        // trace loops
        std::map<int,bool> visited;
        for (auto e : edges) visited[e] = false;

        std::vector<std::vector<int>> loops;
        for (auto e : edges) {
            if (! visited[e]) {
                int e_start = e;
                visited[e] = true;
                loops.emplace_back();
                auto &new_loop = loops.back();
                new_loop.push_back(e_start);
                int v = Edges[e_start].first;
                auto neighbor_edges = vert_to_edge[v];
                int e_next = (neighbor_edges[0] == e_start ? neighbor_edges[1] : neighbor_edges[0]);
                v = (Edges[e_next].first == v ? Edges[e_next].second : Edges[e_next].first);

                while (e_next != e_start) {
                    visited[e_next] = true;
                    new_loop.push_back(e_next);
                    neighbor_edges = vert_to_edge[v];
                    //
                    e_next = (neighbor_edges[0] == e_next ? neighbor_edges[1] : neighbor_edges[0]);
                    v = (Edges[e_next].first == v ? Edges[e_next].second : Edges[e_next].first);
                }
            }
        }

        return loops;
    };

    for (int c = 0; c < C.size(); ++c) {
        auto &cell = C[c];
        // scan faces to collect split edges
        std::vector<int> split_edges;
        for (auto f : cell) {
            if (! F_edge[f].empty()) {
                split_edges.insert(split_edges.end(), F_edge[f].begin(), F_edge[f].end());
            }
        }

        //
        if (! split_edges.empty()) {
            C_split[c] = true;
            // extract loops from split edges
            auto loops = extract_loops(split_edges);
            // create loop faces and the (map: edge -> belonging loop)
            std::map<int,int> split_edge_to_loop;
            for (auto &loop : loops) {
                int f_loop = Faces.size();
                Faces.push_back(loop);
                Face_Implicit.push_back(cur_Impl);
                for (auto &e : loop) {
                    split_edge_to_loop[e] = f_loop;
                }
            }

            // collect all the faces of cell
            std::vector<int> faces;
            for (auto f : cell) {
                if (F_new[f].empty()) { // f is not split
                    faces.push_back(f);
                }
                else { // f is split
                    faces.insert(faces.end(), F_new[f].begin(), F_new[f].end());
                }
            }
            // create map: edge -> face
            std::map<int,std::vector<int>> edge_to_face;
            for (auto f : faces) {
                for (auto e: Faces[f]) {
                    edge_to_face[e].push_back(f);
                }
            }

            // group connected faces
            std::map<int,bool> visited;
            for (auto f : faces) visited[f] = false;

            std::queue<int> Q;
            for (auto f : faces) {
                if (! visited[f]) {
                    Cells.emplace_back();
                    auto &new_cell = Cells.back();
                    std::set<int> boundary_loops;
                    visited[f] = true;
                    Q.push(f);
                    while (! Q.empty()) {
                        int f_front = Q.front();
                        Q.pop();
                        new_cell.push_back(f_front);
                        for (auto e : Faces[f_front]) {
                            if (is_split_edge[e]) {
                                boundary_loops.insert(split_edge_to_loop[e]);
                            }
                            else {
                                const auto &neighbor_faces = edge_to_face[e];
                                int f_opposite = (neighbor_faces[0] == f_front ? neighbor_faces[1] : neighbor_faces[0]);
                                if (! visited[f_opposite]) {
                                    visited[f_opposite] = true;
                                    Q.push(f_opposite);
                                }
                            }
                        }
                    }
                    new_cell.insert(new_cell.end(), boundary_loops.begin(), boundary_loops.end());
                }
            }
        }
    }

    // step 5: remove the original edges, faces and cells that are split
    // update vertices V
    V = Vertices;

    // update edges E
    std::vector<int> E_new_index(Edges.size());
    int n_E = E.size();
    E.clear();
    E_Impl.clear();
    int cur_index = 0;
    for (int e = 0; e < n_E; ++e) {
        if (E_vert[e] == -1) { // edge e is not split
            E.emplace_back(Edges[e]);
            E_Impl.emplace_back(Edge_Implicits[e]);
            E_new_index[e] = cur_index;
            ++cur_index;
        }
    }
    auto eit = Edges.begin() + n_E;
    E.insert(E.end(), eit, Edges.end());
    auto eiIt = Edge_Implicits.begin() + n_E;
    E_Impl.insert(E_Impl.end(), eiIt, Edge_Implicits.end());
    for (int e = n_E; e < Edges.size(); ++e) {
        E_new_index[e] = cur_index;
        ++cur_index;
    }


    // update faces F
    std::vector<int> F_new_indices(Faces.size());
    int n_F = F.size();
    F.clear();
    F_Impl.clear();
    cur_index = 0;
    for (int f = 0; f < n_F; ++f) {
        if (F_edge[f].empty()) { // face f is not split
            F.emplace_back(Faces[f]);
            F_Impl.emplace_back(Face_Implicit[f]);
            F_new_indices[f] = cur_index;
            ++cur_index;
        }
    }
    auto fit = Faces.begin() + n_F;
    F.insert(F.end(), fit, Faces.end());
    auto fiIt = Face_Implicit.begin() + n_F;
    F_Impl.insert(F_Impl.end(), fiIt, Face_Implicit.end());
    for (int f = n_F; f < Faces.size(); ++f) {
        F_new_indices[f] = cur_index;
        ++cur_index;
    }

    // update edge indices stored in faces
    for (auto &face : F) {
        for (auto &e : face) {
            e = E_new_index[e];
        }
    }

    // update cells C
    int n_C = C.size();
    C.clear();
    for (int c = 0; c < n_C; ++c) {
        if (! C_split[c]) {
            C.emplace_back(Cells[c]);
        }
    }

    auto cit = Cells.begin() + n_C;
    C.insert(C.end(), cit, Cells.end());

    // update face indices stored in cells
    for (auto &cell : C) {
        for (auto &f : cell) {
            f = F_new_indices[f];
        }
    }

    ///
}

void Grid::prepare_graph_data() {
    /// step 1: group faces into patches

    // collect faces on implicit surfaces
    std::vector<int> implicit_faces;
    for (int i = 0; i < F.size(); ++i) {
        if (F_Impl[i] != -1) {
            implicit_faces.push_back(i);
        }
    }

    // for implicit faces, create map: edge -> incident faces
    std::map<int, std::vector<int>> edge_to_faces;
    for (int i = 0; i < implicit_faces.size(); ++i) {
        int f = implicit_faces[i];
        for (int e : F[f]) {
            edge_to_faces[e].push_back(i);
        }
    }

    // group implicit faces into patches
    std::vector<int> F_patch(F.size(),-1);
    std::vector<bool> visited(implicit_faces.size(), false);

    // Patches
    P.clear();

    std::queue<int> Q;
    for (int i = 0; i < implicit_faces.size(); ++i) {
        if (!visited[i]) {
            int cur_patch_index = P.size();
            P.emplace_back();
            auto& patch = P.back();

            visited[i] = true;
            Q.push(i);
            while (!Q.empty()) {
                int i_front = Q.front();
                Q.pop();
                int f_front = implicit_faces[i_front];
                patch.push_back(f_front);
                F_patch[f_front] = cur_patch_index;
                for (int e : F[f_front]) {
                    if (E_Impl[e].size() == 1 && edge_to_faces[e].size()==2) { // edge only belongs to one implicit
                        int i_neighbor = (edge_to_faces[e][0] == i_front) ? edge_to_faces[e][1] : edge_to_faces[e][0];
                        if (!visited[i_neighbor]) {
                            visited[i_neighbor] = true;
                            Q.push(i_neighbor);
                        }
                    }
                }
            }
        }
    }

    // the implicit index of each patch
    P_Impl.clear();
    P_Impl.reserve(P.size());
    for (auto& patch : P) {
        P_Impl.push_back(F_Impl[patch.front()]);
    }

    // compute distance weighted area, as well as sample points of patches
    std::vector<std::vector<double>> sample_min_dist;
    std::vector<std::vector<int>> sample_nearest_patch;
    double infinity = std::numeric_limits<double>::infinity();
    for (const auto *impl : Impl) {
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
        auto samples = Impl[impl]->get_sample_points();
        for (int f : patch) {
            // compute center of the face
            auto &face = F[f];
            Point face_center(0,0,0);
            for (int e : face) {
                face_center += V[E[e].first];
                face_center += V[E[e].second];
            }
            face_center /= (2 * face.size());
            // distance weighted area
            double weighted_area = 0;
            for (int e : face) {
                const Point &p1 = V[E[e].first];
                const Point &p2 = V[E[e].second];
                double area = ((p1-face_center).cross(p2-face_center)).norm()/2;
                Point tri_center = (p1+p2+face_center)/3;
                double min_distance = infinity;
                for (int j = 0; j < samples.size(); ++j) {
                    double distance = (tri_center - samples[j]).norm();
                    if (distance < sample_min_dist[impl][j]) {
                        sample_min_dist[impl][j] = distance;
                        sample_nearest_patch[impl][j] = i;
                    }
                    if (distance < min_distance) {
                        min_distance = distance;
                    }
                }
                weighted_area += min_distance * area;
            }
            //
            P_dist[i] += weighted_area;
        }
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

    /// step 2: group cells into blocks

    // create map: face -> cells
    std::map<int, std::vector<int>> face_to_cells;
    for (int i = 0; i < C.size(); ++i) {
        for (int f : C[i]) {
            face_to_cells[f].push_back(i);
        }
    }

    // group cells into blocks
    B_cell.clear();
    B_patch.clear();

    visited.clear();
    visited.resize(C.size(),false);
    assert(Q.empty());
    for (int i = 0; i < C.size(); ++i) {
        if (!visited[i]) {
            std::set<int> boundary_patches;
            B_cell.emplace_back();
            auto &block_cells = B_cell.back();
            visited[i] = true;
            Q.push(i);
            while (!Q.empty()) {
                int front_cell = Q.front();
                Q.pop();
                block_cells.push_back(front_cell);
                auto &cell = C[front_cell];
                for (int f : cell) {
                    if (F_Impl[f] == -1) { // f is not on implicit surface
                        for (int c_f : face_to_cells[f]) {
                            if (c_f != front_cell && !visited[c_f]) {
                                visited[c_f] = true;
                                Q.push(c_f);
                            }
                        }
                    }
                    else { // f is on some implicit surface
                        boundary_patches.insert(F_patch[f]);
                    }
                }
            }
            //
            B_patch.emplace_back(boundary_patches.begin(),boundary_patches.end());
        }
    }

    // P_block: neighboring blocks for patches
    P_block.clear();
    P_block.resize(P.size());
    // P_sign:  signs of implicit function on neighboring blocks
    P_sign.clear();
    P_sign.resize(P.size());

    // find an interior point for each block
    std::vector<Point> B_interior;
    for (auto &cells : B_cell) {
        bool interior_found = false;
        for (int c : cells) {
            for (int f : C[c]) {
                for (int e : F[f]) {
                    if (E_Impl[e].size()==1 && E_Impl[e][0]==-1) {
                        B_interior.emplace_back(0.5 * (V[E[e].first]+V[E[e].second]));
                        interior_found = true;
                        break;
                    }
                }
                if (interior_found) break;
            }
            if (interior_found) break;
        }
        assert(interior_found);
    }

    for (int i = 0; i < B_patch.size(); ++i) {
        for (int p : B_patch[i]) {
            P_block[p].push_back(i);
            // find sign of p's implicit in block i
            int impl_index = P_Impl[p];
            const auto *impl = Impl[impl_index];
            if (impl->function_at(B_interior[i]) > 0) {
                P_sign[p].push_back(1);
            }
            else {
                P_sign[p].push_back(-1);
            }
        }
    }


}

void Grid::graph_cut() {
    // constants
    double inf = std::numeric_limits<double>::infinity();
    double Delta = 0;
    for (auto d : P_dist) Delta += d;
    Delta *= 2;

    // define per-cell costs: hPos, hNeg
    int nBlock = B_patch.size();
    std::vector<double> hPos(nBlock);
    std::vector<double> hNeg(nBlock);

    int nPatch = P.size();
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
    boost::property_map<Graph,boost::edge_residual_capacity_t>::type
            e_residual = get(boost::edge_residual_capacity,g);
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
    double flow = boykov_kolmogorov_max_flow(g ,sid, tid);

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

}

bool Grid::export_grid(const std::string &filename) const {
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
    fout << "edge ";
    fout << E.size() << " " << "2" << std::endl;
    for (auto &e : E) {
        fout << e.first << " " << e.second << std::endl;
    }

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
    fout << "cell ";
    fout << C.size() << std::endl;
    for (auto &c : C) {
        for (auto &f : c) {
            fout << f << " ";
        }
        fout << std::endl;
    }

    // V_Impl
    fout << "vert_implicit ";
    fout << V_Impl.size() << std::endl;
    for (auto &v_impl : V_Impl) {
        for (auto &impl : v_impl) {
            fout << impl << " ";
        }
        fout << std::endl;
    }

    // E_Impl
    fout << "edge_implicit ";
    fout << E_Impl.size() << std::endl;
    for (auto &e_impl : E_Impl) {
        for (auto &impl : e_impl) {
            fout << impl << " ";
        }
        fout << std::endl;
    }

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
    fout << Impl.size() << std::endl;
    for (const auto *impl : Impl) {
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
    fout << "block_cells ";
    fout << B_cell.size() << std::endl;
    for (auto &cells : B_cell) {
        for (auto & cell : cells) {
            fout << cell << " ";
        }
        fout << std::endl;
    }

    // P_block
    fout << "patch_blocks ";
    fout << P_block.size() << std::endl;
    for (auto &blocks : P_block) {
        for (auto & block : blocks) {
            fout << block << " ";
        }
        fout << std::endl;
    }

    // P_sign
    fout << "patch_signs ";
    fout << P_sign.size() << std::endl;
    for (auto &signs : P_sign) {
        for (auto & sign : signs) {
            fout << sign << " ";
        }
        fout << std::endl;
    }


    fout.close();
    std::cout << "export_grid finish: " << filename << std::endl;
    return true;
}




void Grid::init_grid(const Point &p_min, const Point &p_max,
                     int n_cell_x, int n_cell_y, int n_cell_z) {
    // make sure grid resolution is non-negative
    n_cell_x = abs(n_cell_x);
    n_cell_y = abs(n_cell_y);
    n_cell_z = abs(n_cell_z);

    // clear grid data
    V.clear(); E.clear(); F.clear(); C.clear();

    /// step 1:  create vertices
    double dx = (p_max.x() - p_min.x()) / n_cell_x;
    double dy = (p_max.y() - p_min.y()) / n_cell_y;
    double dz = (p_max.z() - p_min.z()) / n_cell_z;

    double cur_x = p_min.x(), cur_y = p_min.y(), cur_z = p_min.z();
    for (int i = 0; i < n_cell_z+1; ++i) {
        cur_x = p_min.x();
        cur_y = p_min.y();
        for (int j = 0; j < n_cell_y+1; ++j) {
            cur_x = p_min.x();
            for (int k = 0; k < n_cell_x+1; ++k) {
                V.emplace_back(cur_x, cur_y, cur_z);
                cur_x += dx;
            }
            cur_y += dy;
        }
        cur_z += dz;
    }

    /// step 2: create edges
    std::vector<int> V_E_x(V.size(), -1);
    std::vector<int> V_E_y(V.size(), -1);
    std::vector<int> V_E_z(V.size(), -1);

    int n_vert_x = n_cell_x + 1;
    int n_vert_y = n_cell_y + 1;
    int n_vert_xy = n_vert_x * n_vert_y;

    int v_idx = 0;
    for (int i = 0; i < n_cell_z; ++i) {
        for (int j = 0; j < n_cell_y; ++j) {
            for (int k = 0; k < n_cell_x; ++k) {
                E.emplace_back(v_idx, v_idx + 1);
                E.emplace_back(v_idx, v_idx + n_vert_x);
                E.emplace_back(v_idx, v_idx + n_vert_xy);
                V_E_x[v_idx] = E.size() -3;
                V_E_y[v_idx] = E.size() -2;
                V_E_z[v_idx] = E.size() -1;
                ++v_idx;
            }
            // k == n_cell_x
            E.emplace_back(v_idx, v_idx + n_vert_x);
            E.emplace_back(v_idx, v_idx + n_vert_xy);
            V_E_y[v_idx] = E.size() - 2;
            V_E_z[v_idx] = E.size() - 1;
            ++v_idx;
        }
        // j == n_cell_y
        for (int k = 0; k < n_cell_x; ++k) {
            E.emplace_back(v_idx, v_idx + 1);
            E.emplace_back(v_idx, v_idx + n_vert_xy);
            V_E_x[v_idx] = E.size() -2;
            V_E_z[v_idx] = E.size() -1;
            ++v_idx;
        }
        // j == n_cell_y && k == n_cell_x
        E.emplace_back(v_idx, v_idx + n_vert_xy);
        V_E_z[v_idx] = E.size() - 1;
        ++v_idx;
    }
    // i == n_cell_z
    for (int j = 0; j < n_cell_y; ++j) {
        for (int k = 0; k < n_cell_x; ++k) {
            E.emplace_back(v_idx, v_idx + 1);
            E.emplace_back(v_idx, v_idx + n_vert_x);
            V_E_x[v_idx] = E.size() -2;
            V_E_y[v_idx] = E.size() -1;
            ++v_idx;
        }
        // i == n_cell_z && k == n_cell_x
        E.emplace_back(v_idx, v_idx + n_vert_x);
        V_E_y[v_idx] = E.size() - 1;
        ++v_idx;
    }
    // i == n_cell_z && j == n_cell_y
    for (int k = 0; k < n_cell_x; ++k) {
        E.emplace_back(v_idx, v_idx + 1);
        V_E_x[v_idx] = E.size() -1;
        ++v_idx;
    }
    // i == n_cell_z && j == n_cell_y && k == n_cell_x
    // nothing here

    /// step 3: create faces
    std::vector<int> V_F_x(V.size(), -1);
    std::vector<int> V_F_y(V.size(), -1);
    std::vector<int> V_F_z(V.size(), -1);

    v_idx = 0;
    for (int i = 0; i < n_cell_z; ++i) {
        for (int j = 0; j < n_cell_y; ++j) {
            for (int k = 0; k < n_cell_x; ++k) {
                // f_x
                F.emplace_back(std::vector<int>{V_E_z[v_idx], V_E_y[v_idx], V_E_z[v_idx + n_vert_x], V_E_y[v_idx + n_vert_xy]});
                // f_y
                F.emplace_back(std::vector<int>{V_E_z[v_idx], V_E_x[v_idx], V_E_z[v_idx +1], V_E_x[v_idx + n_vert_xy]});
                // f_z
                F.emplace_back(std::vector<int>{V_E_y[v_idx], V_E_x[v_idx], V_E_y[v_idx+1], V_E_x[v_idx + n_vert_x]});
                V_F_x[v_idx] = F.size() - 3;
                V_F_y[v_idx] = F.size() - 2;
                V_F_z[v_idx] = F.size() - 1;
                ++v_idx;
            }
            // k == n_cell_x
            // f_x
            F.emplace_back(std::vector<int>{V_E_z[v_idx], V_E_y[v_idx], V_E_z[v_idx + n_vert_x], V_E_y[v_idx + n_vert_xy]});
            V_F_x[v_idx] = F.size() - 1;
            ++v_idx;
        }
        // j == n_cell_y
        for (int k = 0; k < n_cell_x; ++k) {
            // f_y
            F.emplace_back(std::vector<int>{V_E_z[v_idx], V_E_x[v_idx], V_E_z[v_idx +1], V_E_x[v_idx + n_vert_xy]});
            V_F_y[v_idx] = F.size() - 1;
            ++v_idx;
        }
        // j == n_cell_y && k == n_cell_x
        ++v_idx;
    }
    // i == n_cell_z
    for (int j = 0; j < n_cell_y; ++j) {
        for (int k = 0; k < n_cell_x; ++k) {
            // f_z
            F.emplace_back(std::vector<int>{V_E_y[v_idx], V_E_x[v_idx], V_E_y[v_idx+1], V_E_x[v_idx + n_vert_x]});
            V_F_z[v_idx] = F.size() - 1;
            ++v_idx;
        }
        // i == n_cell_z && k == n_cell_x
        ++v_idx;
    }

    /// step 4: create cells
    C.reserve(n_cell_x * n_cell_y * n_cell_z);

    v_idx = 0;
    for (int i = 0; i < n_cell_z; ++i) {
        for (int j = 0; j < n_cell_y; ++j) {
            for (int k = 0; k < n_cell_x; ++k) {
                C.emplace_back(std::vector<int>{V_F_x[v_idx], V_F_y[v_idx], V_F_z[v_idx],
                                                V_F_x[v_idx +1], V_F_y[v_idx + n_vert_x], V_F_z[v_idx + n_vert_xy]});
                ++v_idx;
            }
            // k == n_cell_x
            ++v_idx;
        }
        // j == n_cell_y
        v_idx += n_vert_x;
    }

    // step 5: initialize implicits
    Impl.clear();
    V_Impl.clear(); E_Impl.clear(); F_Impl.clear();
    V_Impl.resize(V.size(), std::vector<int>{-1});
    E_Impl.resize(E.size(), std::vector<int>{-1});
    F_Impl.resize(F.size(), -1); // -1: faces of the regular grid
    ///

}
