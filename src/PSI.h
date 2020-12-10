//
// Created by Charles Du on 11/5/20.
//

#ifndef PSI_PSI_H
#define PSI_PSI_H

#include <string>
#include <Eigen/Core>
#include "Sampled_Implicit.h"

#include "GridSpec.h"


typedef std::pair<int,int> Edge;

class PSI {

public:
    PSI() : Impl_ptr(), arrangement_ready(false),ready_for_graph_cut(false), graph_cut_finished(false) {};
    virtual ~PSI() = default;

    void run(const GridSpec &grid,
             std::vector<std::unique_ptr<Sampled_Implicit>> &implicits);

    bool export_data(const std::string &filename) const;

protected:
    virtual void compute_arrangement_for_graph_cut(
            const GridSpec &grid,
            const std::vector<std::unique_ptr<Sampled_Implicit>> &implicits) = 0;

    void process_samples();
    static void process_samples(
            //input
            const std::vector<Point> &V,
            const std::vector<std::vector<int>> &F,
            const std::vector<std::vector<int>> &P,
            const std::vector<int> &P_Impl,
            const std::vector<std::unique_ptr<Sampled_Implicit>> *Impl_ptr,
            //output
            std::vector<std::vector<int>> &P_samples,
            std::vector<double> &P_dist
    );

    void graph_cut();
    static void graph_cut(
            //input
            const std::vector<double> &P_dist,
            const std::vector<std::vector<int>> &P_samples,
            const std::vector<std::vector<int>> &P_block,
            const std::vector<std::vector<int>> &P_sign,
            const std::vector<std::vector<int>> &B_patch,
            // output
            std::vector<bool> &B_label,
            std::vector<bool> &P_label
    );

public:
    const std::vector<Point>& get_vertices() const { return V; }
    const std::vector<std::vector<int>>& get_faces() const { return F; }
    const std::vector<std::vector<int>>& get_patches() const { return P; }
    const std::vector<std::vector<int>>& get_cells() const { return B_patch; }
    const std::vector<bool>& get_cell_labels() const { return B_label; }

protected:
    // vertices
    std::vector<Point> V;
    // faces
    std::vector<std::vector<int>> F;

    // implicits
    // todo: is this needed?
    std::vector<std::unique_ptr<Sampled_Implicit>> *Impl_ptr;
    // the index of the implicit surface passing each face
    std::vector<int> F_Impl;

    // --- data for graph-cut ---
    bool arrangement_ready;
    bool ready_for_graph_cut;
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
    // blocks incident to each patch
    std::vector<std::vector<int>> P_block;
    // sign of implicit on blocks incident to each patch
    std::vector<std::vector<int>> P_sign;

    // --- result of graph-cut ---
    bool graph_cut_finished;
    // block labels: object -> true, background -> false
    std::vector<bool> B_label;
    // patch labels: surface -> true, not surface -> false
    std::vector<bool> P_label;


};


#endif //PSI_PSI_H
