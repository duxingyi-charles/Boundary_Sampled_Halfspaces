//
// Created by Charles Du on 11/25/20.
//

#ifndef PSI_QUADRIC_SIMPLICIT_H
#define PSI_QUADRIC_SIMPLICIT_H


#include "Sampled_Implicit.h"
#include <Eigen/Geometry>

typedef Eigen::Matrix<double,10,1> Vector10d;

class Quadric_sImplicit : public Sampled_Implicit
{
public:
    // default: plane x==0
    Quadric_sImplicit() : Sampled_Implicit(), coef() {
        coef = Vector10d::Zero(10,1);
        coef(1) = 1;
    };
    // initialize using quadratic coefficients
    Quadric_sImplicit(const Vector10d &c) : coef(c) {};

    ~Quadric_sImplicit() override = default;

    double function_at(const Point &p) const override;
    Eigen::Vector3d   gradient_at(const Point &p) const override;
    void flip() override { coef *= -1; }

private:
    // 10 coefficients: 1, x, y, z, x^2, x*y, x*z, y^2, y*z, z^2
    Vector10d coef;
};



#endif //PSI_QUADRIC_SIMPLICIT_H
