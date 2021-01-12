//
// Created by Charles Du on 10/9/20.
//

#ifndef PSI_CYLINDER_SIMPLICIT_H
#define PSI_CYLINDER_SIMPLICIT_H

#include "Sampled_Implicit.h"

class Cylinder_sImplicit : public Sampled_Implicit {
   public:
    // standard cylinder
    Cylinder_sImplicit()
        : Sampled_Implicit(),
          axis_point(0, 0, 0),
          axis_unit_vector(0, 0, 1),
          radius(1),
          is_flipped(false){};

    Cylinder_sImplicit(const Point &p, const Eigen::Vector3d &v, double r,
                       bool flip)
        : Sampled_Implicit(),
          axis_point(p),
          axis_unit_vector(v.normalized()),
          radius(fabs(r)),
          is_flipped(flip){};

    Cylinder_sImplicit(std::pair<Point, Point> p1p2, double r, bool flip)
        : Sampled_Implicit(),
          axis_point(p1p2.first),
          axis_unit_vector((p1p2.second - p1p2.first).normalized()),
          radius(fabs(r)),
          is_flipped(flip){};

    ~Cylinder_sImplicit() override = default;

    double function_at(const Point &x) const override;
    Eigen::Vector3d gradient_at(const Point &x) const override;

    std::string get_type() const override { return "cylinder"; }

    void get_axis_point(Point &p) const override { p = axis_point; };
    void get_axis_unit_vector(Eigen::Vector3d &vec) const override { vec = axis_unit_vector; };
    void get_radius(double &r) const override { r = radius; };
    void get_is_flipped(bool &flip) const override { flip = is_flipped; };

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

#endif  // PSI_CYLINDER_SIMPLICIT_H
