//
// Created by Charles Du on 8/26/20.
//

#include "Grid.h"

#include <map>
#include <queue>
#include <set>

#include <fstream>
#include <iostream>

Grid::Grid(const Point &p_min, const Point &p_max, int n_cell_x, int n_cell_y, int n_cell_z) {
    init_grid(p_min, p_max, n_cell_x, n_cell_y, n_cell_z);
}

Grid::Grid(const Point &p_min, const Point &p_max) {
    init_grid(p_min, p_max, 1, 1, 1);
}




void Grid::compute_arrangement(const Sampled_Implicit &sImplicit) {
    // record the implicit
    Impl.push_back(&sImplicit);
    int cur_Impl = Impl.size() -1;

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
            V_Impl.emplace_back(std::vector<int>{cur_Impl});
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
                    Edge_Implicits.emplace_back(std::vector<int>{F_Impl[f],cur_Impl});
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
    fout << F_Impl.size() << std::endl;
    for (auto &impl : F_Impl) {
        fout << impl << " ";
    }
    fout << std::endl;


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
    ///

}
