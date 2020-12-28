//
// Created by Charles Du on 12/28/20.
//

#ifndef PSI_PSI_PARAM_H
#define PSI_PSI_PARAM_H

#include <nlohmann/json.hpp>
#include <fstream>

class PSI_Param {

public:
    PSI_Param() :use_distance_weighted_area(false),use_state_space_graph_cut(false),topK(1),consider_adj_diff(true) {};
    ~PSI_Param() = default;

    bool use_distance_weighted_area;
    bool use_state_space_graph_cut;
    int topK;
    bool consider_adj_diff;

    static PSI_Param parse_psi_param(const std::string& param_spec);

};


#endif //PSI_PSI_PARAM_H
