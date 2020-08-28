
#include "src/Hermite_RBF_sImplicit.h"
#include "src/Plane_sImplicit.h"
#include <string>
#include <iostream>

#include "src/Grid.h"

int main()
{
    // ------- sampled implicits -------------
//    Hermite_RBF_sImplicit rbf;
//    // import
//    std::string pts_file = "/Users/charlesdu/Downloads/research/implicit_modeling/code/VIPSS/data/MMA_tmp/input.xyz";
//    std::string coeff_file = "/Users/charlesdu/Downloads/research/implicit_modeling/code/VIPSS/data/MMA_tmp/input_rbfCoeff";
//    rbf.import_Hermite_RBF(pts_file,coeff_file);
//
//    // test
//    Point p(0.5,1,-2);
//    double f = rbf.function_at(p);
//    std::cout << "f(" << p(0) << "," << p(1) << "," << p(2) << ") = " << f << std::endl;
//
//    Eigen::Vector3d g = rbf.gradient_at(p);
//    std::cout << "g(" << p(0) << "," << p(1) << "," << p(2) << ") = " << std::endl;
//    std::cout << g << std::endl;


    // -------- cell arrangement ----------
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

    return 0;
}

