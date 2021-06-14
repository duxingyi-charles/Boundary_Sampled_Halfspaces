//
// Created by Charles Du on 12/28/20.
//

#include "BSH_Param.h"

BSH_Param BSH_Param::parse_bsh_param(const std::string &param_spec) {
    using json = nlohmann::json;
    std::ifstream fin(param_spec.c_str());
    if (!fin) {
        throw std::runtime_error("Parameter Config file does not exist!");
    }
    json config;
    fin >> config;

    BSH_Param spec;

    assert(config.contains("use_distance_weighted_area"));
    spec.use_distance_weighted_area = config["use_distance_weighted_area"].get<bool>();

    if (config.contains("state_space_search")) {
        spec.use_state_space_graph_cut = true;
        auto entry = config["state_space_search"];
        assert(config.contains("topK"));
        spec.topK = entry["topK"].get<int>();
        assert(config.contains("consider_adj_diff"));
        spec.consider_adj_diff = entry["consider_adj_diff"].get<bool>();
    }

    return spec;
}