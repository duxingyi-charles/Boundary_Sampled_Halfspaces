//
// Created by Charles Du on 11/25/20.
//

#include <CLI/CLI.hpp>
#include <string>

//#include "BSH.h"
#include "Mesh_BSH.h"
#include "Topo_BSH.h"

#include "config.h"
//#include "ScopedTimer.h"

#include "Sphere_sImplicit.h"

int main(int argc, char** argv) {
    struct {
        std::string grid_file;
        std::string config_file;
        std::string output_grid_file;
        int grid_size;
        std::string arrangement_algorithm;
    } args;

    CLI::App app{"BSH demo"};
    app.add_option("-G,--grid-file", args.grid_file, "Grid spec file")
            ->required();
    app.add_option("-A,--arr-algo", args.arrangement_algorithm, "Arrangement algorithm " )
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

    // BSH
    Topo_BSH topo_bsh;
    Mesh_BSH mesh_bsh;

    BSH *bsh;
    if (args.arrangement_algorithm[0] == 't') {
        bsh = &topo_bsh;
    } else {
        bsh = &mesh_bsh;
    }

    bsh->run(grid_spec, implicit_functions);

    // add implicit
    Point center(1,1,0);
//    Point center(0,0,0);
    double radius = 1.0;
    std::vector<Point> samples;
    samples.push_back(Point(1,0,0));
    implicit_functions.push_back(std::make_unique<Sphere_sImplicit>(center, radius));
    implicit_functions.back()->set_sample_points(samples);
    bsh->run(grid_spec, implicit_functions);

    // export result
    bsh->export_data(args.output_grid_file);


    return 0;
}

