#include <CLI/CLI.hpp>
#include <string>

#include <Grid.h>

#include "config.h"

int main(int argc, char** argv) {
    struct {
        std::string grid_file;
        std::string config_file;
        std::string output_grid_file;
        int grid_size;
    } args;

    CLI::App app{"Piecewise implicit surface demo"};
    app.add_option("-G,--grid-file", args.grid_file, "Grid spec file")
        ->required();
    app.add_option("config_file", args.config_file, "Configuration file")
        ->required();
    app.add_option("output_grid_file", args.output_grid_file,
                   "Output grid file")
        ->required();
    CLI11_PARSE(app, argc, argv);

    auto implicit_functions =
        initialize_sampled_implicit_functions(args.config_file);

    auto grid_spec = parse_grid_spec(args.grid_file);
    Grid grid(
        {grid_spec.bbox_min[0], grid_spec.bbox_min[1], grid_spec.bbox_min[2]},
        {grid_spec.bbox_max[0], grid_spec.bbox_max[1], grid_spec.bbox_max[2]},
        grid_spec.resolution[0], grid_spec.resolution[1], grid_spec.resolution[2]);

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

