#include <CLI/CLI.hpp>
#include <string>

#include <Grid.h>

#include "config.h"

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

