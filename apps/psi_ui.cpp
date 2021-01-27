#include <config.h>
#include <CLI/CLI.hpp>

#include <igl/opengl/glfw/Viewer.h>

#include <Eigen/Core>
#include <string>

#include "PSIStates.h"
#include "UI.h"

int main(int argc, char** argv)
{
    struct
    {
        std::string config_file;
        std::string grid_file;
    } args;

    CLI::App app{"Piecewise implicit surface demo"};
    app.add_option("-G,--grid", args.grid_file, "Grid spec file")->required();
    app.add_option("config_file", args.config_file, "Configuration file")->required();
    CLI11_PARSE(app, argc, argv);

    // grid specification and implicit functions
    auto implicit_functions = initialize_sampled_implicit_functions(args.config_file);

    auto grid_spec = GridSpec::parse_grid_spec(args.grid_file);

    PSIStates states(grid_spec, implicit_functions);

    igl::opengl::glfw::Viewer viewer;
    viewer.core().background_color.setOnes();

    UI ui;
    ui.set_states(&states);
    ui.initialize(viewer);

    viewer.launch_init(true, false, "BSH Modeling");
    auto bbox = states.get_bbox();
    viewer.core().align_camera_center(bbox);
    viewer.launch_rendering(true);
    viewer.launch_shut();

    return 0;
}
