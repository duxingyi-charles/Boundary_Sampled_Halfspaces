//
// Created by Charles Du on 10/9/20.
//

#ifndef PSI_CYLINDER_SIMPLICIT_H
#define PSI_CYLINDER_SIMPLICIT_H


#include "Sampled_Implicit.h"

class Cylinder_sImplicit : public Sampled_Implicit
{
public:
    // standard cylinder
    Cylinder_sImplicit() : Sampled_Implicit(),
        axis_point(0,0,0), axis_unit_vector(0,0,1), radius(1), is_flipped(false) {};

    Cylinder_sImplicit(const Point &p, const Eigen::Vector3d &v, double r, bool flip) : Sampled_Implicit(),
        axis_point(p), axis_unit_vector(v.normalized()), radius(fabs(r)), is_flipped(flip) {};

    Cylinder_sImplicit(std::pair<Point,Point> p1p2, double r, bool flip) : Sampled_Implicit(),
        axis_point(p1p2.first), axis_unit_vector((p1p2.second-p1p2.first).normalized()), radius(fabs(r)), is_flipped(flip) {};

    ~Cylinder_sImplicit() override = default;

    double function_at(const Point &x) const override;
    Eigen::Vector3d   gradient_at(const Point &x) const override;

private:
    // p: point on cylinder axis
    Point axis_point;
    // v: unit vector along axis
    Eigen::Vector3d axis_unit_vector;
    double radius;
    // is_flipped = false: f(x) = (||x - (p + dot(x-p,v)v)||^2 - r^2)
    // is_flipped = true : f(x) = (r^2 - ||x - (p + dot(x-p,v)v)||^2)
    bool is_flipped;

};

double Cylinder_sImplicit::function_at(const Point &x) const {
    if (is_flipped) {
        return radius*radius - (x - axis_point - (axis_unit_vector.dot(x-axis_point))*axis_unit_vector).squaredNorm();
    } else {
        return (x - axis_point - (axis_unit_vector.dot(x-axis_point))*axis_unit_vector).squaredNorm() - radius*radius;
    }
}

Eigen::Vector3d Cylinder_sImplicit::gradient_at(const Point &x) const {
    if (is_flipped) {
        return (axis_point + (axis_unit_vector.dot(x-axis_point))*axis_unit_vector - x);
    } else {
       return (x - axis_point - (axis_unit_vector.dot(x-axis_point))*axis_unit_vector);
    }
}


#endif //PSI_CYLINDER_SIMPLICIT_H