#pragma once

#include <Hermite_RBF_sImplicit.h>
#include <Plane_sImplicit.h>
#include <Sphere_sImplicit.h>
#include <Cylinder_sImplicit.h>
#include <Cone_sImplicit.h>
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

std::unique_ptr<Sampled_Implicit> initialize_Plane(const nlohmann::json& entry,
                                                    const std::string& path_name) {
    int dimension = 3;

    assert(entry.contains("point"));
    assert(entry["point"].size()==dimension);
    Point p(0,0,0);
    for (int i = 0; i < dimension; ++i) {
        p[i] = entry["point"][i].get<double>();
    }

    assert(entry.contains("normal"));
    assert(entry["normal"].size()==dimension);
    Eigen::Vector3d n(0,0,0);
    for (int i = 0; i < dimension; ++i) {
        n[i] = entry["normal"][i].get<double>();
    }

    auto fn = std::make_unique<Plane_sImplicit>(p,n);

    if (entry.contains("samples")) {
        std::string sample_file = path_name + entry["samples"].get<std::string>();
        std::vector<Point> pts;
        fn->import_xyz(sample_file,pts);
        fn->set_sample_points(pts);
    }

    return fn;
}

std::unique_ptr<Sampled_Implicit> initialize_Sphere(const nlohmann::json& entry,
                                                 const std::string& path_name) {
    assert(entry.contains("center"));
    int dimension = 3;
    assert(entry["center"].size()==dimension);
    Point center(0,0,0);
    for (int i = 0; i < dimension; ++i) {
        center[i] = entry["center"][i].get<double>();
    }
    assert(entry.contains("radius"));
    double radius = entry["radius"].get<double>();
    assert(entry.contains("is_flipped"));
    bool is_flipped = entry["is_flipped"].get<bool>();

    auto fn = std::make_unique<Sphere_sImplicit>(center,radius,is_flipped);

    if (entry.contains("samples")) {
        std::string sample_file = path_name + entry["samples"].get<std::string>();
        std::vector<Point> pts;
        fn->import_xyz(sample_file,pts);
        fn->set_sample_points(pts);
    }

    return fn;
}

std::unique_ptr<Sampled_Implicit> initialize_Cylinder(const nlohmann::json& entry,
                                                    const std::string& path_name) {
    int dimension = 3;

    assert(entry.contains("axis_point1"));
    assert(entry["axis_point1"].size()==dimension);
    Point p1(0,0,0);
    for (int i = 0; i < dimension; ++i) {
        p1[i] = entry["axis_point1"][i].get<double>();
    }
    assert(entry.contains("axis_point2"));
    assert(entry["axis_point2"].size()==dimension);
    Point p2(0,0,0);
    for (int i = 0; i < dimension; ++i) {
        p2[i] = entry["axis_point2"][i].get<double>();
    }
    assert(entry.contains("radius"));
    double radius = entry["radius"].get<double>();
    assert(entry.contains("is_flipped"));
    bool is_flipped = entry["is_flipped"].get<bool>();

    auto fn = std::make_unique<Cylinder_sImplicit>(std::make_pair(p1,p2),radius,is_flipped);

    if (entry.contains("samples")) {
        std::string sample_file = path_name + entry["samples"].get<std::string>();
        std::vector<Point> pts;
        fn->import_xyz(sample_file,pts);
        fn->set_sample_points(pts);
    }

    return fn;
}

std::unique_ptr<Sampled_Implicit> initialize_Cone(const nlohmann::json& entry,
                                                      const std::string& path_name) {
    int dimension = 3;

    assert(entry.contains("apex"));
    assert(entry["apex"].size()==dimension);
    Point p(0,0,0);
    for (int i = 0; i < dimension; ++i) {
        p[i] = entry["apex"][i].get<double>();
    }
    assert(entry.contains("axis_vector"));
    assert(entry["axis_vector"].size()==dimension);
    Eigen::Vector3d v(0,0,0);
    for (int i = 0; i < dimension; ++i) {
        v[i] = entry["axis_vector"][i].get<double>();
    }
    assert(entry.contains("apex_angle"));
    double a = entry["apex_angle"].get<double>();
    assert(entry.contains("is_flipped"));
    bool is_flipped = entry["is_flipped"].get<bool>();

    auto fn = std::make_unique<Cone_sImplicit>(p,v,a,is_flipped);

    if (entry.contains("samples")) {
        std::string sample_file = path_name + entry["samples"].get<std::string>();
        std::vector<Point> pts;
        fn->import_xyz(sample_file,pts);
        fn->set_sample_points(pts);
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
        }
        else if (entry["type"] == "plane") {
            implicit_functions.push_back(initialize_Plane(entry, path_name));
        }
        else if (entry["type"] == "sphere") {
            implicit_functions.push_back(initialize_Sphere(entry, path_name));
        }
        else if (entry["type"] == "cylinder") {
            implicit_functions.push_back(initialize_Cylinder(entry, path_name));
        }
        else if (entry["type"] == "cone") {
            implicit_functions.push_back(initialize_Cone(entry, path_name));
        }
        else {
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
