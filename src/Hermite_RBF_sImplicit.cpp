//
// Created by Charles Du on 8/12/20.
//

#include <iostream>
#include <fstream>
#include <cmath>
#include "Hermite_RBF_sImplicit.h"



bool Hermite_RBF_sImplicit::import_Hermite_RBF(const std::string &pts_file, const std::string &coeff_file)
{
    // import sample points
    std::vector<Point> pts;
    bool succeed = import_xyz(pts_file, pts);
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

bool Hermite_RBF_sImplicit::import_xyz(const std::string &filename, std::vector<Point> &pts) {

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

bool Hermite_RBF_sImplicit::import_RBF_coeff(const std::string &filename, Eigen::VectorXd &a, Eigen::Vector4d &b) {
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

    // second line: coefficient b (d,c0,c1,c2)
    std::getline(reader, line);
    std::istringstream iss2(line);
    double d,c0,c1,c2;
    if (!(iss2 >> d >> c0 >> c1 >> c2)) {
        std::cout << "coeff_b should have 4 elements (in 3D)." << std::endl;
        reader.close();
        return false;
    }
    b << d, c0, c1, c2;

    reader.close();
    return true;
}

double Hermite_RBF_sImplicit::kernel_function(const Point &p1, const Point &p2) {
    return pow((p1-p2).norm(), 3);
}

Eigen::Vector3d Hermite_RBF_sImplicit::kernel_gradient(const Point &p1, const Point &p2) {
    return 3 * (p1-p2).norm() * (p1-p2);
}

Eigen::Matrix3d Hermite_RBF_sImplicit::kernel_Hessian(const Point &p1, const Point &p2) {
    Eigen::Vector3d diff = p1 - p2;
    double len = diff.norm();
    if (len < 1e-8) {
        return Eigen::Matrix3d::Zero();
    }

    Eigen::Matrix3d hess = diff * (diff.transpose()/len);
    hess(0,0) += len;
    hess(1,1) += len;
    hess(2,2) += len;
    hess *= 3;

    return hess;
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

Eigen::Vector3d Hermite_RBF_sImplicit::gradient_at(const Point &p) const {
    size_t npt = sample_points.size();

    Eigen::Vector3d grad;
    grad.setZero();
    // sum(ai * fi)
    for (size_t i = 0; i < npt; ++i) {
        grad += kernel_gradient(p,sample_points[i]) * coeff_a[i];
    }
    // sum(hi * bi)
    Eigen::Matrix3d H;
    for (size_t i = 0; i < npt; ++i) {
        H = kernel_Hessian(p,sample_points[i]);
        grad += H.col(0) * coeff_a[npt+i];
        grad += H.col(1) * coeff_a[2*npt+i];
        grad += H.col(2) * coeff_a[3*npt+i];
    }
    // c
    grad(0) += coeff_b(1);
    grad(1) += coeff_b(2);
    grad(2) += coeff_b(3);

    return grad;
}
