//
// Created by Charles Du on 10/9/20.
//

#ifndef PSI_CONE_SIMPLICIT_H
#define PSI_CONE_SIMPLICIT_H

#include "Sampled_Implicit.h"


class Cone_sImplicit : public Sampled_Implicit
{
public:
    // standard cone
    Cone_sImplicit() : Sampled_Implicit(),
    apex(0,0,0), axis_unit_vector(0,0,-1), apex_angle(M_PI_4), is_flipped(false)
    {};

    Cone_sImplicit(const Point &p, const Eigen::Vector3d &v, double a, bool flip) : Sampled_Implicit(),
    apex(p), axis_unit_vector(v.normalized()), apex_angle(fmod(a,M_PI)), is_flipped(flip)
    {};

    ~Cone_sImplicit() override = default;

    double function_at(const Point &x) const override;
    Eigen::Vector3d   gradient_at(const Point &x) const override;

private:
    // p: cone apex
    Point apex;
    // v: unit vector along cone axis
    Eigen::Vector3d axis_unit_vector;
    // a: angle between generatrix and axis vector (range [0,pi))
    double apex_angle;
    // is_flipped = false: f(x) = cos(a) ||x-p|| - dot(x-p,v)
    // is_flipped = true : f(x) = dot(x-p,v) - cos(a) ||x-p||
    bool is_flipped;

};

#endif //PSI_CONE_SIMPLICIT_H
