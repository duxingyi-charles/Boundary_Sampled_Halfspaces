//
// Created by Charles Du on 11/5/20.
//

#ifndef PSI_TOPO_PSI_H
#define PSI_TOPO_PSI_H

#include "PSI.h"

class Topo_PSI : public PSI {

public:
    Topo_PSI() : PSI() {};

    ~Topo_PSI() override = default;


private:
    void compute_arrangement_for_graph_cut(
            const GridSpec &grid,
            const std::vector<std::unique_ptr<Sampled_Implicit>> &implicits) override;

    // initialize V,E,F,G as a uniform n_cell_x * n_cell_y * n_cell_z grid with lower corner p_min and upper corner p_max
    void init_grid(const Point &p_min, const Point &p_max,
                   int n_cell_x, int n_cell_y, int n_cell_z);

    void compute_arrangement(const Sampled_Implicit &sImplicit, int cur_Impl);

    void prepare_graph_data();

    // the algorithm represent faces as list of boundary edge indices
    // this function convert them to faces represented by a list of boundary vertex indices
    // and save them in F
    void update_F();
private:
    // edges
    std::vector<Edge>  E;
    // cells
    std::vector<std::vector<int>> C;
    // face: each face is a list of boundary edge indices
    std::vector<std::vector<int>> Face_edges;

    // implicits
    // indices of implicits passing each vertex
    std::vector<std::vector<int>> V_Impl;
    // indices of implicits passing each edge
    std::vector<std::vector<int>> E_Impl;
    // the index of the implicit surface passing each face

    // --- data for graph-cut ---
    // (unordered) list of cells of each block
    std::vector<std::vector<int>> B_cell;

};


#endif //PSI_TOPO_PSI_H
