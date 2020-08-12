//
// Created by Charles Du on 8/12/20.
//

#include <iostream>
#include "Hermite_RBF_sImplicit.h"



bool Hermite_RBF_sImplicit::import_Hermite_RBF(const std::string &pts_file, const std::string &coeff_file)
{
    // import sample points
    std::vector<Point> pts;
    bool succeed = import_sample_points(pts_file, pts);
    if (!succeed) {
        std::cout << "Fail to import sample points." << std::endl;
        return false;
    }

    // import RBF coefficients
    Eigen::VectorXd a;
    Eigen::Vector4d b;
    succeed = import_RBF_coeff(coeff_file, a, b);
    if (!succeed) {
        std::cout << "Fail to import RBF coefficients." << std::endl;
        return false;
    }

    // succeed == true
    sample_points = pts;
    coeff_a = a;
    coeff_b = b;
    return true;
}

bool Hermite_RBF_sImplicit::import_sample_points(const std::string &filename, std::vector<Point> &pts) {

}

bool Hermite_RBF_sImplicit::import_RBF_coeff(const std::string &filename, Eigen::VectorXd &a, Eigen::Vector4d &b) {

}