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
        int grid_size;
        std::string arrangement_algorithm;
        int topK;
        bool consider_adj_diff;
    } args;

    CLI::App app{"Piecewise implicit surface demo"};
    app.add_option("-G,--grid-file", args.grid_file, "Grid spec file")
            ->required();
    app.add_option("-A,--arr-algo", args.arrangement_algorithm, "Arrangement algorithm " )
            ->required();
    app.add_option("-k,--top-k", args.topK, "topK " )
            ->required();
    app.add_option("-a,--consider-adj-diff", args.consider_adj_diff,
                   "whether to consider patches adjacent to disconnected components." )
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
    Topo_PSI topo_psi;
    Mesh_PSI mesh_psi;

    PSI *psi;
    if (args.arrangement_algorithm[0] == 't') {
        psi = &topo_psi;
    } else {
        psi = &mesh_psi;
    }

    psi->run(grid_spec, implicit_functions);
    psi->export_data("/Users/charlesdu/Downloads/research/implicit_modeling/code/piecewise_sampled_implicits/data/state_search_test/init_res.grid");
    {
        ScopedTimer<> timer("search for connected result");
//        int topK = 1;
        psi->search_for_connected_result(args.topK, args.consider_adj_diff);
    }


    // export result
    psi->export_data(args.output_grid_file);


    return 0;
}

