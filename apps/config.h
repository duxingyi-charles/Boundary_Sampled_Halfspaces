#pragma once

#include <Sampled_Implicit.h>

#include <nlohmann/json.hpp>
#include <string>

std::unique_ptr<Sampled_Implicit> initialize_RBF(const nlohmann::json& entry,
                                                 const std::string& path_name);

std::unique_ptr<Sampled_Implicit> initialize_Sphere(
    const nlohmann::json& entry, const std::string& path_name);

std::unique_ptr<Sampled_Implicit> initialize_Cylinder(
    const nlohmann::json& entry, const std::string& path_name);

std::vector<std::unique_ptr<Sampled_Implicit>>
initialize_sampled_implicit_functions(const std::string& config_file);

struct GridSpec {
    Eigen::Vector3i resolution;
    Eigen::Vector3d bbox_min;
    Eigen::Vector3d bbox_max;
};

GridSpec parse_grid_spec(const std::string& grid_spec);
