//
// Created by Charles Du on 1/5/21.
//

#ifndef PSI_TORUS_SIMPLICIT_H
#define PSI_TORUS_SIMPLICIT_H

#include "Sampled_Implicit.h"


class Torus_sImplicit : public Sampled_Implicit {
public:
    // standard torus
    Torus_sImplicit()
            : Sampled_Implicit(),
              center(0, 0, 0),
              axis_unit_vector(0, 0, 1),
              major_radius(1), minor_radius(0.5),
              is_flipped(false){};

    Torus_sImplicit(const Point &p, const Eigen::Vector3d &v, double R, double r,
                       bool flip)
            : Sampled_Implicit(),
              center(p),
              axis_unit_vector(v.normalized()),
              major_radius(fabs(R)),
              minor_radius(fmin(fabs(R), fabs(r))), // minor radius can't be larger than major radius
              is_flipped(flip){};


    ~Torus_sImplicit() override = default;

    double function_at(const Point &x) const override;
    Eigen::Vector3d gradient_at(const Point &x) const override;

private:
    // p: center point of torus
    Point center;
    // v: unit vector along axis
    Eigen::Vector3d axis_unit_vector;
    // R: major radius
    double major_radius;
    // r: minor radius
    double minor_radius;
    // is_flipped = false: function value is negative inside the torus
    // is_flipped = true : function value is positive inside the torus
    bool is_flipped;

};


#endif //PSI_TORUS_SIMPLICIT_H
