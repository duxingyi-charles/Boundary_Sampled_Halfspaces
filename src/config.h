#pragma once

#include <Sampled_Implicit.h>

#include <memory>
#include <nlohmann/json.hpp>
#include <string>

std::unique_ptr<Sampled_Implicit> initialize_RBF(const nlohmann::json& entry,
                                                 const std::string& path_name);

std::unique_ptr<Sampled_Implicit> initialize_Plane(
    const nlohmann::json& entry, const std::string& path_name);

std::unique_ptr<Sampled_Implicit> initialize_Quadric(
        const nlohmann::json& entry, const std::string& path_name);

std::unique_ptr<Sampled_Implicit> initialize_Sphere(
    const nlohmann::json& entry, const std::string& path_name);

std::unique_ptr<Sampled_Implicit> initialize_Cylinder(
    const nlohmann::json& entry, const std::string& path_name);

std::unique_ptr<Sampled_Implicit> initialize_Cone(const nlohmann::json& entry,
                                                  const std::string& path_name);

std::vector<std::unique_ptr<Sampled_Implicit>>
initialize_sampled_implicit_functions(const std::string& config_file);
