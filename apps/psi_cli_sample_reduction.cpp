//
// Created by Charles Du on 12/10/20.
//

#include <CLI/CLI.hpp>
#include <string>

#include "Mesh_PSI.h"
#include "Topo_PSI.h"

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

    CLI::App app{"Piecewise implicit surface demo"};
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

    auto param_spec = PSI_Param::parse_psi_param(args.param_file);

    // PSI
    Topo_PSI topo_psi;
    Mesh_PSI mesh_psi;

    PSI *psi;
    if (args.arrangement_algorithm[0] == 't') {
        psi = &topo_psi;
    } else {
        psi = &mesh_psi;
    }
    psi->set_parameters(param_spec);

    // run PSI on dense samples
    psi->run(grid_spec, implicit_functions);

    // sparse samples
    auto implicit_functions_sparse =
            initialize_sampled_implicit_functions(args.config_file_sparse);
    //
    {
        ScopedTimer<> timer("sample reduction");
//        psi->reduce_samples(nullptr);
        psi->reduce_samples(&implicit_functions_sparse);
    }


    // export result
    psi->export_data(args.output_grid_file);


    return 0;
}

