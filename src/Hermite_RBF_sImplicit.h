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

    ~Hermite_RBF_sImplicit() = default;

    bool import_Hermite_RBF(const std::string& pts_file, const std::string& coeff_file);


    double function_at(const Point &) const override;
    Eigen::Vector3d gradient_at(const Point &) const override;

protected:


private:
    Eigen::VectorXd coeff_a;
    Eigen::Vector4d coeff_b;

    bool import_sample_points(const std::string &filename, std::vector<Point> &pts);
    bool import_RBF_coeff(const std::string &filename, Eigen::VectorXd &a, Eigen::Vector4d &b);

};


#endif //PSI_HERMITE_RBF_SIMPLICIT_H
