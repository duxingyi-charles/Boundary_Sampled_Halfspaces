//
// Created by Charles Du on 8/14/20.
//

#ifndef PSI_GRID_H
#define PSI_GRID_H

#include "Sampled_Implicit.h"
#include <vector>
#include <string>

typedef std::pair<int,int> Edge;

class Grid {
public:

    Grid() = default;
    // construct a single cuboid grid with lower corner p_min and upper corner p_max
    Grid(const Point &p_min, const Point &p_max);
    // construct of uniform n_cell_x * n_cell_y * n_cell_z grid with lower corner p_min and upper corner p_max
    Grid(const Point &p_min, const Point &p_max, int n_cell_x, int n_cell_y, int n_cell_z);

    ~Grid() = default;

    void compute_arrangement(const Sampled_Implicit &);

    void prepare_graph_data();

    void graph_cut();

    bool export_grid(const std::string &filename) const;

private:
    // vertices
    std::vector<Point> V;
    // edges
    std::vector<Edge>  E;
    // faces
    std::vector<std::vector<int>> F;
    // cells
    std::vector<std::vector<int>> C;

    // implicits
    std::vector<const Sampled_Implicit*> Impl;
    // indices of implicits passing each vertex
    std::vector<std::vector<int>> V_Impl;
    // indices of implicits passing each edge
    std::vector<std::vector<int>> E_Impl;
    // the index of the implicit surface passing each face
    std::vector<int> F_Impl;

    // --- data for graph-cut ---
    // (unordered) list of faces on each patch
    std::vector<std::vector<int>> P;
    // index of implicit for each patch
    std::vector<int> P_Impl;
    // indices of sample points on each patch
    std::vector<std::vector<int>> P_samples;
    // distance weighted area of each patch
    std::vector<double> P_dist;
    // (unordered) list of boundary patches of each block
    std::vector<std::vector<int>> B_patch;
    // (unordered) list of cells of each block
    std::vector<std::vector<int>> B_cell;
    // blocks incident to each patch
    std::vector<std::vector<int>> P_block;
    // sign of implicit on blocks incident to each patch
    std::vector<std::vector<int>> P_sign;

    // --- result of graph-cut ---
    // block labels: object -> true, background -> false
    std::vector<bool> B_label;
    // patch labels: surface -> true, not surface -> false
    std::vector<bool> P_label;



    // initialize V,E,F,G as a uniform n_cell_x * n_cell_y * n_cell_z grid with lower corner p_min and upper corner p_max
    void init_grid(const Point &p_min, const Point &p_max,
                   int n_cell_x, int n_cell_y, int n_cell_z);

};

#endif //PSI_GRID_H
