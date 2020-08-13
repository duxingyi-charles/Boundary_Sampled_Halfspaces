
#include "src/Hermite_RBF_sImplicit.h"
#include <string>
#include <iostream>

int main()
{
    Hermite_RBF_sImplicit rbf;
    // import
    std::string pts_file = "/Users/charlesdu/Downloads/research/implicit_modeling/code/VIPSS/data/MMA_tmp/input.xyz";
    std::string coeff_file = "/Users/charlesdu/Downloads/research/implicit_modeling/code/VIPSS/data/MMA_tmp/input_rbfCoeff";
    rbf.import_Hermite_RBF(pts_file,coeff_file);

    // test
    Point p(0.5,1,-2);
    double f = rbf.function_at(p);
    std::cout << "f(" << p(0) << "," << p(1) << "," << p(2) << ") = " << f << std::endl;

    Eigen::Vector3d g = rbf.gradient_at(p);
    std::cout << "g(" << p(0) << "," << p(1) << "," << p(2) << ") = " << std::endl;
    std::cout << g << std::endl;

    return 0;
}