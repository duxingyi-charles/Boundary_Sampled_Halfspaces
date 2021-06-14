//
// Created by Charles Du on 1/23/21.
//

#include <CLI/CLI.hpp>
#include <string>

//#include "BSH.h"
#include "Mesh_BSH.h"
#include "Topo_BSH.h"

#include "config.h"
//#include "ScopedTimer.h"

int main(int argc, char** argv) {
    struct {
        std::string grid_file;
        std::string config_file;
        std::string output_dir;
    } args;

    CLI::App app{"Piecewise implicit surface demo"};
    app.add_option("-G,--grid-file", args.grid_file, "Grid spec file")
            ->required();
    app.add_option("config_file", args.config_file, "Configuration file")
            ->required();
    app.add_option("output_dir", args.output_dir,
                   "Output directory")
            ->required();
    CLI11_PARSE(app, argc, argv);


    // grid specification and implicit functions
    auto implicit_functions =
            initialize_sampled_implicit_functions(args.config_file);

    auto grid_spec = GridSpec::parse_grid_spec(args.grid_file);


    // BSH
    Mesh_BSH bsh;

    bsh.run(grid_spec, implicit_functions);

    // export sampled implicits
    std::cout << "export sampled implcits..." << std::endl;
    bsh.export_sampled_implicits(args.output_dir);

    return 0;
}

