//
// Created by Charles Du on 8/26/20.
//

#include "Grid.h"

#include <map>
#include <queue>
#include <set>

void Grid::compute_arrangement(const Sampled_Implicit &sImplicit)
{
    // initialize
    std::vector<Point> Vertices = V;
    std::vector<Edge>  Edges = E;
    std::vector<std::vector<int>> Faces = F;
    std::vector<std::vector<int>> Cells = C;

    std::vector<int> E_vert(E.size(), -1); // -1 for null
    std::vector<std::pair<int,int>> E_new(E.size());
    std::vector<bool> is_split_edge(E.size(), false);

    std::vector<std::vector<int>> F_edge(F.size());
    std::vector<std::vector<int>> F_new(F.size());

    std::vector<bool> C_split(C.size(), false);

    // step 1: get implicit function value at each vertex
    std::vector<double> V_value(V.size());
    for (int i = 0; i < V.size(); ++i) {
        V_value[i] = sImplicit.function_at(V[i]);
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
            // split edge
//            Edges.push_back(Edge(vi, v_index));
            Edges.emplace_back(vi, v_index);
//            Edges.push_back(Edge(v_index, vj));
            Edges.emplace_back(v_index, vj);
            is_split_edge.push_back(false);
            is_split_edge.push_back(false);
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
//                    Edges.push_back(Edge(v_start, v_end));
                    Edges.emplace_back(v_start, v_end);
                    is_split_edge.push_back(true);
                    F_edge[f].push_back(e_split);
                    split_face.push_back(e_split);
                    //
                    int f_split = Faces.size();
                    Faces.push_back(split_face);
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
                    Q.push(f);
                    while (! Q.empty()) {
                        int f_front = Q.front();
                        Q.pop();
                        new_cell.push_back(f_front);
                        visited[f_front] = true;
                        for (auto e : Faces[f_front]) {
                            if (is_split_edge[e]) {
                                boundary_loops.insert(split_edge_to_loop[e]);
                            }
                            else {
                                const auto &neighbor_faces = edge_to_face[e];
                                int f_opposite = (neighbor_faces[0] == f_front ? neighbor_faces[1] : neighbor_faces[0]);
                                if (! visited[f_opposite]) {
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
    int cur_index = 0;
    for (int e = 0; e < n_E; ++e) {
        if (E_vert[e] == -1) { // edge e is not split
            E.emplace_back(Edges[e]);
            E_new_index[e] = cur_index;
            ++cur_index;
        }
    }
    auto eit = Edges.begin() + n_E;
    E.insert(E.end(), eit, Edges.end());
    for (int e = n_E; e < Edges.size(); ++e) {
        E_new_index[e] = cur_index;
        ++cur_index;
    }

    // update faces F
    std::vector<int> F_new_indices(Faces.size());
    int n_F = F.size();
    F.clear();
    cur_index = 0;
    for (int f = 0; f < n_F; ++f) {
        if (F_edge[f].empty()) { // face f is not split
            F.emplace_back(Faces[f]);
            F_new_indices[f] = cur_index;
            ++cur_index;
        }
    }
    auto fit = Faces.begin() + n_F;
    F.insert(F.end(), fit, Faces.end());
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
