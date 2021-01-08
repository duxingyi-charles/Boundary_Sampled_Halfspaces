//
// Created by Charles Du on 1/5/21.
//

#ifndef PSI_TORUS_SIMPLICIT_H
#define PSI_TORUS_SIMPLICIT_H

#include "Sampled_Implicit.h"

#include <iostream>

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

    bool has_control_points() const override { return true; }
    const std::vector<Point>& get_control_points() const override {
        if (m_control_pts.empty()) {
            m_control_pts.push_back(center);
        }
        return m_control_pts;
    }
    void set_control_points(const std::vector<Point>& pts) override {
        if (pts.size() < 1) {
            std::cerr << "Torus primitive expects at least 1 control points";
            return;
        }
        m_control_pts = pts;
        center = pts[0];
    }

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

    mutable std::vector<Point> m_control_pts;
};


#endif //PSI_TORUS_SIMPLICIT_H
