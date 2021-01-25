//
// Created by Charles Du on 11/5/20.
//

#ifndef PSI_GRIDSPEC_H
#define PSI_GRIDSPEC_H

#include <Eigen/Core>

#include <nlohmann/json.hpp>
#include <fstream>

//struct GridSpec {
//    Eigen::Vector3i resolution;
//    Eigen::Vector3d bbox_min;
//    Eigen::Vector3d bbox_max;
//};

class GridSpec {

public:
    GridSpec() = default;
    ~GridSpec() = default;

    Eigen::Vector3i resolution;
    Eigen::Vector3d bbox_min;
    Eigen::Vector3d bbox_max;

    static GridSpec parse_grid_spec(const std::string& grid_spec);

    double get_bbox_area() const;
};






#endif //PSI_GRIDSPEC_H
