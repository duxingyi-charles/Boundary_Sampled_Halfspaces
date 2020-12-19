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

struct PSI_Search_State {
    std::vector<int> prohibited_patches;
    std::vector<bool> P_label;
    std::vector<bool> B_label;
    double cost;
};

class PSI {

public:
    PSI() : Impl_ptr(), arrangement_ready(false),ready_for_graph_cut(false),use_distance_weighted_area(false),graph_cut_finished(false) {};
    virtual ~PSI() = default;

    void run(const GridSpec &grid,
             std::vector<std::unique_ptr<Sampled_Implicit>> &implicits);

    bool export_data(const std::string &filename) const;

    bool export_state(const std::string &filename,
                             const std::vector<bool> &P_label,
                             const std::vector<std::vector<int>> &components) const;

    // Algorithms for reverse engineering

    // after running PSI on a dense sampling,
    // find a subset of samples that produce the same result as the original dense sampling
    void reduce_samples();

    void reduce_samples(const std::vector<std::unique_ptr<Sampled_Implicit>> *Impl_ptr_sparse);

    void search_for_connected_result(int topK);

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
            bool  use_distance_weighted_area,
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

    static void graph_cut(
            //input
            const std::vector<double> &P_dist,
            const std::vector<std::vector<int>> &P_samples,
            const std::vector<std::vector<int>> &P_block,
            const std::vector<std::vector<int>> &P_sign,
            const std::vector<std::vector<int>> &B_patch,
            const std::vector<int> &prohibited_patches,
            // output
            std::vector<bool> &B_label,
            std::vector<bool> &P_label,
            double &cut_cost
    );

    void compute_patch_adjacency(std::vector<std::vector<int>> &P_Adj_same, std::vector<std::vector<int>> &P_Adj_diff);

    static void compute_patch_adjacency(
            //input
            const std::vector<std::vector<int>> &P,
            const std::vector<std::vector<int>> &F,
            const std::vector<Point> &V,
            const std::vector<int> &P_Impl,
            //output
            std::vector<std::vector<int>> &P_Adj_same,
            std::vector<std::vector<int>> &P_Adj_diff
            );

    static void get_unsampled_patch_components(
            //input
            const std::vector<std::vector<int>> &P_Adj_same,
            const std::vector<std::vector<int>> &P_Adj_diff,
            const std::vector<std::vector<int>> &P_samples,
            const std::vector<bool> P_label,
            //output
            std::vector<std::vector<int>> &components
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
    // flag: use distance weighted area or not
    bool use_distance_weighted_area;
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
