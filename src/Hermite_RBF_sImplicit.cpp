//
// Created by Charles Du on 8/12/20.
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
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

bool Hermite_RBF_sImplicit::import_sample_points(const std::string &filename, std::vector<Point> &pts) const {

    std::ifstream reader(filename.data(), std::ofstream::in);

    if (!reader.good()) {
        std::cout << "Can not open the file " << filename << std::endl;
        return false;
    }else {
        std::cout << "Reading: "<<filename<< std::endl;
    }

    // first number: dimension (2 or 3)
    int dim;
    reader >> dim;
    // read point coordinates
    if (dim != 3) {
        std::cout << "Can't handle non-3D points." << std::endl;
        reader.close();
        return false;
    }
    pts.clear();
    double x,y,z;
    while(!reader.eof()){
        reader >> x >> y >> z;
        pts.emplace_back(Point(x,y,z));
    }

    reader.close();
    return true;
}

bool Hermite_RBF_sImplicit::import_RBF_coeff(const std::string &filename, Eigen::VectorXd &a, Eigen::Vector4d &b) const {
    std::ifstream reader(filename.data(), std::ofstream::in);
    if (!reader.good()) {
        std::cout << "Can not open the file " << filename << std::endl;
        return false;
    }else {
        std::cout << "Reading: "<<filename<< std::endl;
    }

    // first line: coefficient a
    std::string line;
    std::getline(reader, line);
    std::istringstream iss(line);
    double val;
    std::vector<double> tmp_a;
    while (iss >> val) tmp_a.push_back(val);
    a.resize(tmp_a.size());
    for (size_t i = 0; i < tmp_a.size(); ++i) {
        a(i) = tmp_a[i];
    }

    // second line: coefficient b
    std::getline(reader, line);
    iss.str(line);
    double d,c0,c1,c2;
    if (!(iss >> d >> c0 >> c1 >> c2)) {
        std::cout << "coeff_b should have 4 elements (in 3D)." << std::endl;
        reader.close();
        return false;
    }
    b << d, c0, c1, c2;


    reader.close();
    return true;
}

double Hermite_RBF_sImplicit::kernel_function(const Point &p1, const Point &p2) const {
    return pow((p1-p2).norm(), 3);
}

Eigen::Vector3d Hermite_RBF_sImplicit::kernel_gradient(const Point &p1, const Point &p2) const {
    return 3 * (p1-p2).norm() * (p1-p2);
}

double Hermite_RBF_sImplicit::function_at(const Point &p) const {
    size_t npt = sample_points.size();
    int dim = 3;

    Eigen::VectorXd kern(npt*(dim+1));
    for (size_t i = 0; i < npt; ++i) {
        kern(i) = kernel_function(p,sample_points[i]);
    }
    Eigen::Vector3d G;
    for (size_t i = 0; i < npt; ++i) {
        G = kernel_gradient(p,sample_points[i]);
        for (int j = 0; j < 3; ++j) {
            kern(npt+i+j*npt) = G(j);
        }
    }
    double loc_part = kern.dot(coeff_a);

    Eigen::Vector4d kb(1,p(0),p(1),p(2));
    double poly_part = kb.dot(coeff_b);

    double re = loc_part + poly_part;
    return re;
}

Eigen::Vector3d Hermite_RBF_sImplicit::gradient_at(const Point &) const {

}
