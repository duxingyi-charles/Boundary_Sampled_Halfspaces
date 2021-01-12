//
// Created by Charles Du on 10/8/20.
//

#ifndef PSI_SPHERE_SIMPLICIT_H
#define PSI_SPHERE_SIMPLICIT_H

#include "Sampled_Implicit.h"

#include <iostream>

class Sphere_sImplicit : public Sampled_Implicit {
   public:
    // unit sphere
    Sphere_sImplicit()
        : Sampled_Implicit(), center(0, 0, 0), radius(1), is_flipped(false){};

    Sphere_sImplicit(const Point &o, double r, bool flip)
        : Sampled_Implicit(), center(o), radius(fabs(r)), is_flipped(flip){};

    Sphere_sImplicit(const Point &o, double r)
        : Sampled_Implicit(), center(o), radius(fabs(r)), is_flipped(false){};

    Sphere_sImplicit(const Point &o, const Point &s, bool flip)
        : Sampled_Implicit(),
          center(o),
          radius((o - s).norm()),
          is_flipped(flip){};

    Sphere_sImplicit(const Point &o, const Point &s)
        : Sampled_Implicit(),
          center(o),
          radius((o - s).norm()),
          is_flipped(false){};

    ~Sphere_sImplicit() override = default;

    double function_at(const Point &p) const override;
    Eigen::Vector3d gradient_at(const Point &p) const override;

    bool has_control_points() const override { return true; }
    const std::vector<Point>& get_control_points() const override {
        if (m_control_pts.empty()) {
            m_control_pts.push_back(center);
        }
        return m_control_pts;
    }
    void set_control_points(const std::vector<Point>& pts) override {
        if (pts.size() < 1) {
            std::cerr << "Sphere primitive expects at least 1 control points";
            return;
        }
        m_control_pts = pts;
        center = pts[0];
    }

    double get_radius() const {
        return radius;
    }

    void set_radius(double r) {
        radius = r;
    }

   private:
    Point center;
    double radius;
    // is_flipped = false: f(x) = (||x-center||^2 - r^2)
    // is_flipped = true : f(x) = (r^2 - ||x-center||^2)
    bool is_flipped;

    mutable std::vector<Point> m_control_pts;
};

#endif  // PSI_SPHERE_SIMPLICIT_H
