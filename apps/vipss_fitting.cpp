//
// Created by Charles Du on 1/12/21.
//

#include <CLI/CLI.hpp>
#include <string>

#include "Hermite_RBF_sImplicit.h"

//#include "ScopedTimer.h"

int main(int argc, char** argv) {
    struct {
        std::string sample_file;
        double error_bound;
        std::string output_coef_file;
        std::string output_pts_file;
    } args;

    CLI::App app{"vipss implicit fitting"};
    app.add_option("-S,--sample-file", args.sample_file, "Input sample points file")
            ->required();
    app.add_option("-E,--error-bound", args.error_bound, "error bound")
            ->required();
    app.add_option("output_rbf_coef_file", args.output_coef_file, "Output RBF coef file")
            ->required();
    app.add_option("output_control_point_file", args.output_pts_file,
                   "Output control points file")
            ->required();
    CLI11_PARSE(app, argc, argv);


    // input points
    std::vector<Point> sample_pts;
    Sampled_Implicit::import_xyz(args.sample_file, sample_pts);

    // vipss fitting
    Hermite_RBF_sImplicit rbf;
    rbf.fit_RBF(sample_pts, args.error_bound);

    // export result
    rbf.export_RBF_coeff(args.output_coef_file);
    Sampled_Implicit::export_xyz(args.output_pts_file,rbf.get_control_points());

    return 0;
}

