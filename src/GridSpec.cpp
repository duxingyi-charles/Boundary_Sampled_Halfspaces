//
// Created by Charles Du on 11/5/20.
//

#include "GridSpec.h"


GridSpec GridSpec::parse_grid_spec(const std::string& grid_spec) {
    using json = nlohmann::json;
    std::ifstream fin(grid_spec.c_str());
    if (!fin) {
        throw std::runtime_error("Config file does not exist!");
    }
    json config;
    fin >> config;

    GridSpec spec;
    spec.resolution <<
                    config["resolution"][0],
            config["resolution"][1],
            config["resolution"][2];
    spec.bbox_min <<
                  config["bbox_min"][0],
            config["bbox_min"][1],
            config["bbox_min"][2];
    spec.bbox_max <<
                  config["bbox_max"][0],
            config["bbox_max"][1],
            config["bbox_max"][2];

    return spec;
}

