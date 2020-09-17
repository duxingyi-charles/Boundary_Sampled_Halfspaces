
#include "src/Hermite_RBF_sImplicit.h"
#include "src/Plane_sImplicit.h"
#include <string>
#include <iostream>

#include "src/Grid.h"


// sampled implicits
void test1() {
    Hermite_RBF_sImplicit rbf;
    // import
    std::string pts_file = "/Users/charlesdu/Downloads/research/implicit_modeling/code/VIPSS/data/MMA_tmp/input.xyz";
    std::string coeff_file = "/Users/charlesdu/Downloads/research/implicit_modeling/code/VIPSS/data/MMA_tmp/input_rbfCoeff";
    rbf.import_Hermite_RBF(pts_file,coeff_file);

    // test: compute function and gradient
    Point p(0.5,1,-2);
    double f = rbf.function_at(p);
    std::cout << "f(" << p(0) << "," << p(1) << "," << p(2) << ") = " << f << std::endl;

    Eigen::Vector3d g = rbf.gradient_at(p);
    std::cout << "g(" << p(0) << "," << p(1) << "," << p(2) << ") = " << std::endl;
    std::cout << g << std::endl;
}

// arrangement in a single-cell grid
void test2() {
    Plane_sImplicit plane1(Point(0,-1,1), Point(-1,0,1), Point(1,-1,0));
    Plane_sImplicit plane2(Point(-1,-1,1), Point(1,1,-1), Point(-1,1,-1));
    Plane_sImplicit plane3(Point(0,0,0), Point(1,0,0), Point(0,1,0));

    Grid grid(Point(-1,-1,-1), Point(1,1,1));
    // before
    grid.export_grid("/Users/charlesdu/Downloads/research/implicit_modeling/code/VIPSS/data/MMA_tmp/init.grid");
    // compute arrangement
    grid.compute_arrangement(plane1);
    grid.compute_arrangement(plane2);
    grid.compute_arrangement(plane3);
    // after
    grid.export_grid("/Users/charlesdu/Downloads/research/implicit_modeling/code/VIPSS/data/MMA_tmp/result.grid");
}

// uniform grid initialization
void test3() {
    Grid grid(Point(-1,-1,-1), Point(1,1,1), 1, 2, 3);
    grid.export_grid("/Users/charlesdu/Downloads/research/implicit_modeling/code/VIPSS/data/MMA_tmp/init.grid");
}

// arrangement in uniform grid (planar implicits)
// 2020/08/29 known issues: implicit surfaces passing grid vertex/edge/face will produce near degenerate cells
void test4() {
    Plane_sImplicit plane1(Point(0,-1,1), Point(-1,0,1), Point(1,-1,0));
    Plane_sImplicit plane2(Point(-1,-1,1), Point(1,1,-1), Point(-1,1,-1));
    Plane_sImplicit plane3(Point(0,0,0), Point(1,0,0), Point(0,1,0));

    Grid grid(Point(-1,-1,-1), Point(1,1,1), 2, 2, 2);
    // before
    grid.export_grid("/Users/charlesdu/Downloads/research/implicit_modeling/code/VIPSS/data/MMA_tmp/init.grid");
    // compute arrangement
    grid.compute_arrangement(plane1);
    grid.compute_arrangement(plane2);
    grid.compute_arrangement(plane3);
    // after
    grid.export_grid("/Users/charlesdu/Downloads/research/implicit_modeling/code/VIPSS/data/MMA_tmp/result.grid");
}

// arrangement in uniform grid (RBF implicits)
void test5() {
    // import RBF
    Hermite_RBF_sImplicit rbf;
    std::string pts_file = "/Users/charlesdu/Downloads/research/implicit_modeling/code/VIPSS/data/MMA_tmp/input.xyz";
    std::string coeff_file = "/Users/charlesdu/Downloads/research/implicit_modeling/code/VIPSS/data/MMA_tmp/input_rbfCoeff";
    rbf.import_Hermite_RBF(pts_file,coeff_file);

    // init grid
    Grid grid(Point(-2,-2,-2), Point(2,2,2), 64, 64, 64);
    // before
    grid.export_grid("/Users/charlesdu/Downloads/research/implicit_modeling/code/VIPSS/data/MMA_tmp/init.grid");
    // compute arrangement
    grid.compute_arrangement(rbf);
    // after
    grid.export_grid("/Users/charlesdu/Downloads/research/implicit_modeling/code/VIPSS/data/MMA_tmp/result.grid");
}

// arrangement in uniform grid (multiple implicits)
void test6() {
    // import RBF
    Hermite_RBF_sImplicit rbf;
    std::string pts_file = "/Users/charlesdu/Downloads/research/implicit_modeling/code/VIPSS/data/MMA_tmp/input.xyz";
    std::string coeff_file = "/Users/charlesdu/Downloads/research/implicit_modeling/code/VIPSS/data/MMA_tmp/input_rbfCoeff";
    rbf.import_Hermite_RBF(pts_file,coeff_file);

    // plane
    Plane_sImplicit plane1(Point(0,-1,1), Point(-1,0,1), Point(1,-1,0));

    // init grid
    Grid grid(Point(-2,-2,-2), Point(2,2,2), 8, 8, 8);
    // before
    grid.export_grid("/Users/charlesdu/Downloads/research/implicit_modeling/code/VIPSS/data/MMA_tmp/init.grid");
    // compute arrangement
    grid.compute_arrangement(rbf);
    grid.compute_arrangement(plane1);
    // after
    grid.export_grid("/Users/charlesdu/Downloads/research/implicit_modeling/code/VIPSS/data/MMA_tmp/result.grid");
}

// get ready for graph-cut
void test7() {
    // import RBF
    Hermite_RBF_sImplicit rbf;
    std::string pts_file = "/Users/charlesdu/Downloads/research/implicit_modeling/code/VIPSS/data/MMA_tmp/input.xyz";
    std::string coeff_file = "/Users/charlesdu/Downloads/research/implicit_modeling/code/VIPSS/data/MMA_tmp/input_rbfCoeff";
    rbf.import_Hermite_RBF(pts_file,coeff_file);

    // plane
    Plane_sImplicit plane1(Point(0,-1,1), Point(-1,0,1), Point(1,-1,0));

    // init grid
    Grid grid(Point(-2,-2,-2), Point(2,2,2), 15, 15, 15);
    // before
    grid.export_grid("/Users/charlesdu/Downloads/research/implicit_modeling/code/VIPSS/data/MMA_tmp/init.grid");
    // compute arrangement
    grid.compute_arrangement(rbf);
    grid.compute_arrangement(plane1);
    // prepare for graph-cut
    grid.prepare_graph_data();
    // after
    grid.export_grid("/Users/charlesdu/Downloads/research/implicit_modeling/code/VIPSS/data/MMA_tmp/result.grid");
}

// within a pipeline
void test8(int num_rbf, int grid_resolution) {
    assert(num_rbf >0 && grid_resolution > 0);

    std::vector<Hermite_RBF_sImplicit> rbf_list;
    for (int i = 1; i <= num_rbf; ++i) {
        Hermite_RBF_sImplicit rbf;
        std::string pts_file = "/Users/charlesdu/Downloads/research/implicit_modeling/code/VIPSS/data/MMA_tmp/input_"+std::to_string(i)+".xyz";
        std::string coeff_file = "/Users/charlesdu/Downloads/research/implicit_modeling/code/VIPSS/data/MMA_tmp/input_"+std::to_string(i)+"_rbfCoeff";
        std::string sample_file = "/Users/charlesdu/Downloads/research/implicit_modeling/code/VIPSS/data/MMA_tmp/input_"+std::to_string(i)+"_sample.xyz";
        rbf.import_sampled_Hermite_RBF(pts_file,coeff_file,sample_file);
        rbf_list.push_back(rbf);
    }

    // init grid
    Grid grid(Point(-2,-2,-2), Point(2,2,2), grid_resolution, grid_resolution, grid_resolution);
    // before
    grid.export_grid("/Users/charlesdu/Downloads/research/implicit_modeling/code/VIPSS/data/MMA_tmp/init.grid");
    // compute arrangement
    for (const auto &rbf : rbf_list) {
        grid.compute_arrangement(rbf);
    }
    // prepare for graph-cut
    grid.prepare_graph_data();
    // after
    grid.export_grid("/Users/charlesdu/Downloads/research/implicit_modeling/code/VIPSS/data/MMA_tmp/result.grid");
};

int main(int argc, char** argv)
{
    test8(atoi(argv[1]),atoi(argv[2]));
    return 0;
}

