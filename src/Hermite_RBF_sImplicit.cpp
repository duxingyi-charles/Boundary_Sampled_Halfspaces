//
// Created by Charles Du on 8/12/20.
//

#include <iostream>
#include <fstream>
#include <cmath>
#include <limits>
#include "Hermite_RBF_sImplicit.h"

#include "vipss/rbfcore.h"


void Hermite_RBF_sImplicit::compute_RBF_coeff(const std::vector<Point> &points, Eigen::VectorXd &a, Eigen::Vector4d &b) {
    //
    vector<double> Vs;
    RBF_Core rbf_core;
    RBF_Paras para = RBF_Paras::Set_RBF_PARA();

    // read points into Vs
    para.point_dimension = points[0].size();
    for (const auto &p : points) {
        for (int i=0; i<p.size(); ++i) {
            Vs.push_back(p(i));
        }
    }
    // vipss computation
    rbf_core.InjectData(Vs,para);
    rbf_core.BuildK(para);
    rbf_core.InitNormal(para);
    rbf_core.OptNormal(0);

    // read result of vipss
    if (rbf_core.b.size() != 4) {
        std::cout << "coeff_b should have 4 elements (in 3D)." << std::endl;
        return;
    }

    // coefficient a
    a.resize(rbf_core.a.size());
    for (int i = 0; i < rbf_core.a.size(); ++i) {
        a(i) = rbf_core.a(i);
    }
    // coefficient b
    for (int i = 0; i < rbf_core.b.size(); ++i) {
        b(i) = rbf_core.b(i);
    }

}


void Hermite_RBF_sImplicit::fit_RBF(const std::vector<Point> &points, double error_bound)
{
    if (error_bound <= 0 || points.size() <= 3) {   // interpolation
        update_RBF_coeff(points);
        return;
    }

    // find mass center of points
    Point mass_center(0,0,0);
    for (const auto &p : points) {
        mass_center += p;
    }
    mass_center /= points.size();

    // first control point: the point closest to mass center
    std::vector<bool> selected(points.size(),false);
    std::vector<Point> control_pts;

    double min_dist = std::numeric_limits<double>::max();
    int selected_id = -1;
    for (int i = 0; i < points.size(); ++i) {
        double dist = (points[i] - mass_center).norm();
        if (dist < min_dist) {
            min_dist = dist;
            selected_id = i;
        }
    }
    selected[selected_id] = true;
    Point p1 = points[selected_id];
    control_pts.push_back(p1);

    // second control point: the point farthest from the first control point
    double max_dist = 0;
    selected_id = -1;
    for (int i = 0; i < points.size(); ++i) {
        double dist = (points[i] - p1).norm();
        if (dist > max_dist) {
            max_dist = dist;
            selected_id = i;
        }
    }
    if (max_dist == 0) {  // all points coincide
        update_RBF_coeff(control_pts);
        return;
    }
    selected[selected_id] = true;
    Point p2 = points[selected_id];
    control_pts.push_back(p2);

    // third control point: the point farthest from the first two control points
    max_dist = 0;
    selected_id = -1;
    for (int i = 0; i < points.size(); ++i) {
        double dist = min((points[i]-p1).norm(), (points[i]-p2).norm());
        if (dist > max_dist) {
            max_dist = dist;
            selected_id = i;
        }
    }
    if (max_dist == 0) {  // there are only two distinct points
        update_RBF_coeff(control_pts);
        return;
    }
    selected[selected_id] = true;
    Point p3 = points[selected_id];
    control_pts.push_back(p3);

    // initialize rbf with three points
    update_RBF_coeff(control_pts);

    // add control point until error falls below error bound
    max_dist = 0;
    selected_id = -1;
    for (int i = 0; i < points.size(); ++i) {
        if (!selected[i]) {
            double f_i = function_at(points[i]);
            auto   g_i = gradient_at(points[i]);
            double dist = (g_i.norm() > 0) ? fabs(f_i)/(g_i.norm()) : fabs(f_i);
            if (dist > max_dist) {
                max_dist = dist;
                selected_id = i;
            }
        }
    }

    while ((max_dist > error_bound) && (selected_id != -1)) {
        selected[selected_id] = true;
        control_pts.push_back(points[selected_id]);
        update_RBF_coeff(control_pts);

        max_dist = 0;
        selected_id = -1;
        for (int i = 0; i < points.size(); ++i) {
            if (!selected[i]) {
                double f_i = function_at(points[i]);
                auto   g_i = gradient_at(points[i]);
                double dist = (g_i.norm() > 0) ? fabs(f_i)/(g_i.norm()) : fabs(f_i);
                if (dist > max_dist) {
                    max_dist = dist;
                    selected_id = i;
                }
            }
        }
    }

}


bool Hermite_RBF_sImplicit::import_Hermite_RBF(const std::string &pts_file, const std::string &coeff_file)
{
    // import control points
    std::vector<Point> pts;
    bool succeed = import_xyz(pts_file, pts);
    if (!succeed) {
        std::cout << "Fail to import RBF control points." << std::endl;
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
    control_points = pts;
    coeff_a = a;
    coeff_b = b;
    return true;
}

bool Hermite_RBF_sImplicit::import_sampled_Hermite_RBF(const std::string &pts_file, const std::string &coeff_file,
                                                       const std::string &sample_file) {
    bool succeed = import_Hermite_RBF(pts_file, coeff_file);
    if (!succeed) {
        std::cout << "Fail to import Hermite RBF control points or coefficients." << std::endl;
        return false;
    }

    //debug
//    print_control_points();
//    print_coeff();
//    std::cout << "function value at (0,0,0) = " << function_at(Point(0,0,0)) << std::endl;

    // import sample points
    std::vector<Point> pts;
    succeed = import_xyz(sample_file, pts);
    if (!succeed) {
        std::cout << "Fail to import sample points." << std::endl;
        return false;
    }

    // succeed == true
    sample_points = pts;
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

bool Hermite_RBF_sImplicit::export_RBF_coeff(const std::string &filename) const {
    std::ofstream fout(filename, std::ofstream::out);
    if (!fout.good()) {
        std::cout << "Can not create output file " << filename << std::endl;
        return false;
    }

    //precision of output
    fout.precision(std::numeric_limits<double>::max_digits10);

    // coef_a
    for (int i = 0; i < coeff_a.size(); ++i) {
        fout << coeff_a(i) << " ";
    }
    fout << std::endl;

    // coef_b
    for (int i = 0; i < coeff_b.size(); ++i) {
        fout << coeff_b(i) << " ";
    }
    fout << std::endl;

    fout.close();
    std::cout << "export_RBF_coeff finish: " << filename << std::endl;
    return true;
}

void Hermite_RBF_sImplicit::consistent_update_RBF_coeff(const std::vector<Point> &points) {
    int n_total = sample_points.size() + control_points.size();
    if (n_total == 0) {  // nothing to keep consistent
        update_RBF_coeff(points);
        return;
    }
    //before update, compute gradient at samples and control points
    std::vector<Point> prev_control_points = control_points;
    std::vector<Eigen::Vector3d> prev_grads;
    prev_grads.reserve(n_total);
    for (const auto &p : sample_points) {
        prev_grads.emplace_back(gradient_at(p));
    }
    for (const auto &p : prev_control_points) {
        prev_grads.emplace_back(gradient_at(p));
    }

    update_RBF_coeff(points);

    // after update, compute gradients
    std::vector<Eigen::Vector3d> new_grads;
    new_grads.reserve(n_total);
    for (const auto &p : sample_points) {
        new_grads.emplace_back(gradient_at(p));
    }
    for (const auto &p : prev_control_points) {
        new_grads.emplace_back(gradient_at(p));
    }

    // orientation consistency check
    int n_consistent = 0;
    for (int i = 0; i < n_total; ++i) {
        if (prev_grads[i].dot(new_grads[i]) > 0) {
            ++n_consistent;
        }
    }

    if (2 * n_consistent < n_total) {
        flip_sign();
    }
}

void Hermite_RBF_sImplicit::print_coeff() const {
    // coef_a
    std::cout << "coef_a: " << std::endl;
    for (int i = 0; i < coeff_a.size(); ++i) {
        std::cout << coeff_a(i) << " ";
    }
    std::cout << std::endl;

    // coef_b
    std::cout << "coef_b: " << std::endl;
    for (int i = 0; i < coeff_b.size(); ++i) {
        std::cout << coeff_b(i) << " ";
    }
    std::cout << std::endl;
}

void Hermite_RBF_sImplicit::print_control_points() const {
    std::cout << "control points: " << std::endl;
    for (const auto &p : control_points) {
        std::cout << p.x() << " " << p.y() << " " << p.z() << std::endl;
    }
    std::cout << std::endl;
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
    size_t npt = control_points.size();
    int dim = 3;

    Eigen::VectorXd kern(npt*(dim+1));
    for (size_t i = 0; i < npt; ++i) {
        kern(i) = kernel_function(p,control_points[i]);
    }
    Eigen::Vector3d G;
    for (size_t i = 0; i < npt; ++i) {
        G = kernel_gradient(p,control_points[i]);
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
    size_t npt = control_points.size();

    Eigen::Vector3d grad;
    grad.setZero();
    // sum(ai * fi)
    for (size_t i = 0; i < npt; ++i) {
        grad += kernel_gradient(p,control_points[i]) * coeff_a[i];
    }
    // sum(hi * bi)
    Eigen::Matrix3d H;
    for (size_t i = 0; i < npt; ++i) {
        H = kernel_Hessian(p,control_points[i]);
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


bool Hermite_RBF_sImplicit::save(const std::string &dir, const std::string &name, nlohmann::json &json_obj) const
{
    json_obj.clear();
    json_obj["type"] = "rbf";
    json_obj["name"] = name;

    // control points
    std::string point_filename = name + ".xyz";
    json_obj["points"] = point_filename;

    std::string point_file = dir + point_filename;
    if (!export_xyz(point_file, control_points)) {
        return false;
    }

    // rbf coef
    std::string coef_filename = name + "_rbfCoeff";
    json_obj["rbf_coeffs"] = coef_filename;

    std::string coef_file = dir + coef_filename;
    if(!export_RBF_coeff(coef_file)) {
        return false;
    }

    // sample points
    std::string sample_filename = name + "_sample.xyz";
    json_obj["samples"] = sample_filename;

    std::string sample_file = dir + sample_filename;
    if (!export_xyz(sample_file,sample_points)) {
        return false;
    }

    //
    return true;
}