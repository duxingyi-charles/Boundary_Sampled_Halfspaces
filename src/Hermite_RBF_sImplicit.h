//
// Created by Charles Du on 8/12/20.
//

#ifndef PSI_HERMITE_RBF_SIMPLICIT_H
#define PSI_HERMITE_RBF_SIMPLICIT_H

#include "Sampled_Implicit.h"


class Hermite_RBF_sImplicit : public Sampled_Implicit {
public:
    // constant 0 function
    Hermite_RBF_sImplicit() : Sampled_Implicit(), coeff_a(), coeff_b(0,0,0,0) {};

    ~Hermite_RBF_sImplicit() override = default;

    // initialize with control points
    Hermite_RBF_sImplicit(const std::vector<Point> &control_pts,
                          const std::vector<Point> &sample_pts)
                          : Sampled_Implicit(sample_pts) { update_RBF_coeff(control_pts);}

    bool import_Hermite_RBF(const std::string &pts_file, const std::string &coeff_file);
    bool import_sampled_Hermite_RBF(const std::string &pts_file,
                                    const std::string &coeff_file,
                                    const std::string &sample_file);


    double function_at(const Point &) const override;
    Eigen::Vector3d gradient_at(const Point &) const override;

    static void compute_RBF_coeff(const std::vector<Point> &points, Eigen::VectorXd &a, Eigen::Vector4d &b);

    void fit_RBF(const std::vector<Point> &points, double error_bound);

    void update_RBF_coeff(const std::vector<Point> &points);

    void flip_sign() {
        coeff_a *= -1;
        coeff_b *= -1;
    }

protected:


private:
    Eigen::VectorXd coeff_a;
    Eigen::Vector4d coeff_b;

    std::vector<Point> control_points;

    static bool import_RBF_coeff(const std::string &filename, Eigen::VectorXd &a, Eigen::Vector4d &b);

    // |p1-p2|^3
    static double kernel_function(const Point &p1, const Point &p2);
    // 3 |p1-p2| (p1-p2)
    static Eigen::Vector3d kernel_gradient(const Point &p1, const Point &p2);
    // 3 [ |p1-p2|I + (p1-p2)*(p1-p2)^T/|p1-p1| ]
    static Eigen::Matrix3d kernel_Hessian(const Point &p1, const Point &p2);


};


#endif //PSI_HERMITE_RBF_SIMPLICIT_H
