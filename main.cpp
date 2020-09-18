#include <CLI/CLI.hpp>
#include <cassert>
#include <exception>
#include <iostream>
#include <nlohmann/json.hpp>
#include <string>

#include "src/Grid.h"
#include "src/Hermite_RBF_sImplicit.h"
#include "src/Plane_sImplicit.h"
#include "src/Sampled_Implicit.h"

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

int main(int argc, char** argv) {
    struct {
        std::string config_file;
        std::string output_grid_file;
        int grid_size;
    } args;

    CLI::App app{"Piecewise implicit surface demo"};
    app.add_option("-s,--grid-size", args.grid_size, "Grid size")
        ->default_val(32);
    app.add_option("config_file", args.config_file, "Configuration file")
        ->required();
    app.add_option("output_grid_file", args.output_grid_file,
                   "Output grid file")
        ->required();
    CLI11_PARSE(app, argc, argv);

    auto implicit_functions =
        initialize_sampled_implicit_functions(args.config_file);
    Grid grid(Point(-2, -2, -2), Point(2, 2, 2), args.grid_size, args.grid_size,
              args.grid_size);

    // before
    grid.export_grid("init.grid");

    // compute arrangement
    for (const auto& rbf : implicit_functions) {
        grid.compute_arrangement(*rbf);
    }

    // prepare for graph-cut
    grid.prepare_graph_data();

    // after
    grid.export_grid(args.output_grid_file);

    return 0;
}

