//
// Created by Charles Du on 12/10/20.
//

#include <CLI/CLI.hpp>
#include <string>

#include "Mesh_BSH.h"
#include "Topo_BSH.h"

#include "config.h"
#include "ScopedTimer.h"



int main(int argc, char** argv) {
    struct {
        std::string grid_file;
        std::string config_file;
        std::string output_grid_file;
        std::string config_file_sparse;
        std::string param_file;
        std::string arrangement_algorithm;
    } args;

    CLI::App app{"BSH demo"};
    app.add_option("-G,--grid-file", args.grid_file, "Grid spec file")
            ->required();
    app.add_option("-P,--param-file", args.param_file, "Parameter spec file")
            ->required();
    app.add_option("-A,--arr-algo", args.arrangement_algorithm, "Arrangement algorithm " )
            ->required();
    app.add_option("config_file", args.config_file, "Configuration file")
            ->required();
    app.add_option("config_file_sparse", args.config_file_sparse, "Configuration file for sparse sampling")
            ->required();
    app.add_option("output_grid_file", args.output_grid_file,
                   "Output grid file")
            ->required();
    CLI11_PARSE(app, argc, argv);


    // grid specification and implicit functions
    auto implicit_functions =
            initialize_sampled_implicit_functions(args.config_file);

    auto grid_spec = GridSpec::parse_grid_spec(args.grid_file);

    auto param_spec = BSH_Param::parse_bsh_param(args.param_file);

    // BSH
    Topo_BSH topo_bsh;
    Mesh_BSH mesh_bsh;

    BSH *bsh;
    if (args.arrangement_algorithm[0] == 't') {
        bsh = &topo_bsh;
    } else {
        bsh = &mesh_bsh;
    }
    bsh->set_parameters(param_spec);

    // run BSH on dense samples
    bsh->run(grid_spec, implicit_functions);

    // sparse samples
    auto implicit_functions_sparse =
            initialize_sampled_implicit_functions(args.config_file_sparse);
    //
    {
        ScopedTimer<> timer("sample reduction");
//        bsh->reduce_samples(nullptr);
        bsh->reduce_samples(&implicit_functions_sparse);
    }


    // export result
    bsh->export_data(args.output_grid_file);


    return 0;
}

