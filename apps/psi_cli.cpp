#include <CLI/CLI.hpp>
#include <string>

//#include <Grid.h>
//#include "PSI.h"
#include "Mesh_PSI.h"

#include "config.h"
#include "ScopedTimer.h"

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


    // grid specification and implicit functions
    auto implicit_functions =
        initialize_sampled_implicit_functions(args.config_file);

    auto grid_spec = GridSpec::parse_grid_spec(args.grid_file);

    // PSI
    Mesh_PSI psi;
    psi.run(grid_spec, implicit_functions);

    // export result
    psi.export_data(args.output_grid_file);

//--------------------- old ----------------
//    Grid grid(
//        {grid_spec.bbox_min[0], grid_spec.bbox_min[1], grid_spec.bbox_min[2]},
//        {grid_spec.bbox_max[0], grid_spec.bbox_max[1], grid_spec.bbox_max[2]},
//        grid_spec.resolution[0], grid_spec.resolution[1], grid_spec.resolution[2]);
//
//    // compute arrangement
//    {
//        ScopedTimer<> timer("topological arrangement");
//        for (const auto& rbf : implicit_functions) {
//            grid.compute_arrangement(*rbf);
//        }
//    }
//
//    // graph-cut
//    {
//        ScopedTimer<> timer("graph cut");
//        grid.prepare_graph_data();
//        grid.graph_cut();
//    }
//
//    // after
//    grid.export_grid(args.output_grid_file);
// -------------------------------------------------

    return 0;
}

