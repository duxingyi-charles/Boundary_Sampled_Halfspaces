//
// Created by Charles Du on 10/8/20.
//

#ifndef PSI_SPHERE_SIMPLICIT_H
#define PSI_SPHERE_SIMPLICIT_H

#include "Sampled_Implicit.h"

class Sphere_sImplicit : public Sampled_Implicit
{
public:
    // unit sphere
    Sphere_sImplicit() : Sampled_Implicit(), center(0,0,0), radius(1), is_flipped(false) {};

    Sphere_sImplicit(const Point &o, double r, bool flip)
    : center(o), radius(fabs(r)), is_flipped(flip) {};

    Sphere_sImplicit(const Point &o, double r)
    : center(o), radius(fabs(r)), is_flipped(false) {};

    Sphere_sImplicit(const Point &o, const Point &s, bool flip)
    : center(o), radius((o-s).norm()), is_flipped(flip) {};

    Sphere_sImplicit(const Point &o, const Point &s)
    : center(o), radius((o-s).norm()), is_flipped(false) {};

    ~Sphere_sImplicit() override = default;

    double function_at(const Point &p) const override;
    Eigen::Vector3d   gradient_at(const Point &p) const override;

private:
    Point center;
    double radius;
    // is_flipped = false: f(x) = (||x-center||^2 - r^2)
    // is_flipped = true : f(x) = (r^2 - ||x-center||^2)
    bool is_flipped;

};

double Sphere_sImplicit::function_at(const Point &p) const {
    if (is_flipped) {
        return (radius*radius - (p-center).squaredNorm());
    }
    else {
        return ((p-center).squaredNorm() - radius*radius);
    }
}

Eigen::Vector3d Sphere_sImplicit::gradient_at(const Point &p) const {
    if (is_flipped) {
        return 2*(center - p);
    }
    else {
        return 2*(p - center);
    }
}




#endif //PSI_SPHERE_SIMPLICIT_H
