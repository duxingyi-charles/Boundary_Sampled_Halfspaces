//
// Created by Charles Du on 11/5/20.
//

#ifndef PSI_PSI_H
#define PSI_PSI_H

#include <string>
#include <Eigen/Core>
#include "Sampled_Implicit.h"

#include "GridSpec.h"
#include "PSI_Param.h"


typedef std::pair<int,int> Edge;

struct PSI_Search_State {
    std::vector<int> prohibited_patches;
    std::vector<bool> P_label;
    std::vector<bool> B_label;
    double cost;
};


class PSI {

public:
    PSI() : Impl_ptr(),
    arrangement_ready(false),ready_for_graph_cut(false),ready_for_connected_graph_cut(false),graph_cut_finished(false),
    use_distance_weighted_area(false),use_state_space_graph_cut(false),topK(1),consider_adj_diff(true) {};



    virtual ~PSI() = default;

    void run(const GridSpec &grid,
             std::vector<std::unique_ptr<Sampled_Implicit>> &implicits);
    virtual void update_implicit(const GridSpec &grid,
            const std::unique_ptr<Sampled_Implicit>& fn, size_t i) {}
    void process_samples();
    void graph_cut();

    bool export_data(const std::string &filename) const;

    void set_parameters(const PSI_Param &param_spec);

    // Algorithms for reverse engineering

    // after running PSI on a dense sampling,
    // find a subset of samples that produce the same result as the original dense sampling
    void reduce_samples(const std::vector<std::unique_ptr<Sampled_Implicit>> *Impl_ptr_sparse);

protected:
    virtual void compute_arrangement_for_graph_cut(
            const GridSpec &grid,
            const std::vector<std::unique_ptr<Sampled_Implicit>> &implicits) = 0;

    static void process_samples(
            //input
            const std::vector<Point> &V,
            const std::vector<std::vector<int>> &F,
            const std::vector<std::vector<int>> &P,
            const std::vector<int> &P_Impl,
            const std::vector<bool> &P_touch_bbox,
            const double &bbox_area,
            const std::vector<std::unique_ptr<Sampled_Implicit>> *Impl_ptr,
            bool  use_distance_weighted_area,
            //output
            std::vector<std::vector<int>> &P_samples,
            std::vector<double> &P_dist
    );

    static void simple_graph_cut(
            //input
            const std::vector<double> &P_dist,
            const std::vector<std::vector<int>> &P_samples,
            const std::vector<std::vector<int>> &P_block,
            const std::vector<std::vector<int>> &P_sign,
            const std::vector<std::vector<int>> &B_patch,
            // output
            std::vector<bool> &B_label,
            std::vector<bool> &P_label,
            double &cut_cost
    );

    static void simple_graph_cut(
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

    static void simple_graph_cut(
            //input
            const std::vector<double> &P_dist,
            const std::vector<std::vector<int>> &P_samples,
            const std::vector<std::vector<int>> &P_block,
            const std::vector<std::vector<int>> &P_sign,
            const std::vector<std::vector<int>> &B_patch,
            double Delta,
            // output
            std::vector<bool> &B_label,
            std::vector<bool> &P_label,
            double &cut_cost
    );

    void connected_graph_cut();

    static void connected_graph_cut(
            //input
            const std::vector<double> &P_dist,
            const std::vector<std::vector<int>> &P_samples,
            const std::vector<std::vector<int>> &P_block,
            const std::vector<std::vector<int>> &P_sign,
            const std::vector<std::vector<int>> &B_patch,
            const std::vector<std::vector<int>> &P_Adj_same,
            const std::vector<std::vector<int>> &P_Adj_diff,
            int topK, bool consider_adj_diff,
            //output
            std::vector<bool> &B_label,
            std::vector<bool> &P_label
            );

    // export intermediate state for state-space search
    static bool export_state(const std::string &filename,
                             const std::vector<bool> &P_label,
                             const std::vector<std::vector<int>> &components);

    void compute_patch_adjacency();

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
            const std::vector<bool> &P_label,
            bool consider_adj_diff,
            //output
            std::vector<std::vector<int>> &components
            );

    // helpers
    static double point_triangle_distance(const Point &p1, const Point &p2, const Point& p3, const Point& q);

public:
    const std::vector<Point>& get_vertices() const { return V; }
    const std::vector<std::vector<int>>& get_faces() const { return F; }
    const std::vector<std::vector<int>>& get_patches() const { return P; }
    const std::vector<std::vector<int>>& get_cells() const { return B_patch; }
    const std::vector<bool>& get_cell_labels() const { return B_label; }
    const std::vector<std::vector<int>>& get_cell_patch_adjacency() const { return P_block; }
    const std::vector<std::vector<int>>& get_cell_patch_sign() const { return P_sign; }
    const std::vector<bool>& get_patch_labels() const { return P_label; }
    const std::vector<int>& get_patch_implicit() const { return P_Impl; }
    size_t get_num_implicits() const { return Impl_ptr->size(); }
    size_t get_num_patches() const { return P.size(); }
    size_t get_num_cells() const { return B_patch.size(); }

    bool export_sampled_implicits(const std::string &output_dir) const;

protected:
    // bounding box
    double bbox_area;

    // vertices
    std::vector<Point> V;
    // faces
    std::vector<std::vector<int>> F;

    // implicits
    std::vector<std::unique_ptr<Sampled_Implicit>> *Impl_ptr;
    // the index of the implicit surface passing each face
    std::vector<int> F_Impl;

    // --- data for graph-cut ---
    bool arrangement_ready;
    bool ready_for_graph_cut;
    bool ready_for_connected_graph_cut;
    bool graph_cut_finished;
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
    // flag: whether the patch touches the bounding box
    std::vector<bool> P_touch_bbox;

    // --- parameters for connected state search ---
    // flag: use state-space graph-cut or not
    bool use_state_space_graph_cut;
    // in state space search, expand a state into k children states with least graph-cut cost.
    // if topK=0, explore all children states.
    int topK;
    // consider adjacent patches from other implicits or not
    bool consider_adj_diff;
    // adjacency list of patches from the same implicit
    std::vector<std::vector<int>> P_Adj_same;
    // adjacency list of patches from different implicits
    std::vector<std::vector<int>> P_Adj_diff;

    // --- result of graph-cut ---
    // block labels: object -> true, background -> false
    std::vector<bool> B_label;
    // patch labels: surface -> true, not surface -> false
    std::vector<bool> P_label;


};


#endif //PSI_PSI_H
