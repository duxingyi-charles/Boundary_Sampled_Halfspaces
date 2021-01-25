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


double GridSpec::get_bbox_area() const {
    double dx = fabs(bbox_max.x() - bbox_min.x());
    double dy = fabs(bbox_max.y() - bbox_min.y());
    double dz = fabs(bbox_max.z() - bbox_min.z());

    return 2 * (dx*dy + dx*dz + dy*dz);

}

