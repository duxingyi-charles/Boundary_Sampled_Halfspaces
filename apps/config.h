#pragma once

#include <Hermite_RBF_sImplicit.h>
#include <Plane_sImplicit.h>
#include <Sampled_Implicit.h>

#include <cassert>
#include <exception>
#include <nlohmann/json.hpp>
#include <string>

std::unique_ptr<Sampled_Implicit> initialize_RBF(const nlohmann::json& entry,
                                                 const std::string& path_name) {
    std::string point_file, coeff_file, sample_file;

    assert(entry.contains("points"));
    point_file = path_name + entry["points"].get<std::string>();
    assert(entry.contains("rbf_coeffs"));
    coeff_file = path_name + entry["rbf_coeffs"].get<std::string>();

    auto fn = std::make_unique<Hermite_RBF_sImplicit>();

    if (entry.contains("samples")) {
        sample_file = path_name + entry["samples"].get<std::string>();
        fn->import_sampled_Hermite_RBF(point_file, coeff_file, sample_file);
    } else {
        fn->import_Hermite_RBF(point_file, coeff_file);
    }

    return fn;
}

std::vector<std::unique_ptr<Sampled_Implicit>>
initialize_sampled_implicit_functions(const std::string& config_file) {
    using json = nlohmann::json;
    std::ifstream fin(config_file.c_str());
    if (!fin) {
        throw std::runtime_error("Config file does not exist!");
    }
    json config;
    fin >> config;

    if (!config.contains("input")) {
        throw std::runtime_error("Invalid config file");
    }

    auto pos = config_file.find_last_of("/\\");
    auto path_name = config_file.substr(0, pos + 1);

    std::vector<std::unique_ptr<Sampled_Implicit>> implicit_functions;

    for (auto entry : config["input"]) {
        if (entry["type"] == "rbf") {
            implicit_functions.push_back(initialize_RBF(entry, path_name));
        } else {
            throw std::runtime_error(
                "Unsupported implicit function type detected");
        }
    }

    return implicit_functions;
}

struct GridSpec {
    Eigen::Vector3i resolution;
    Eigen::Vector3d bbox_min;
    Eigen::Vector3d bbox_max;
};

GridSpec parse_grid_spec(const std::string& grid_spec) {
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
