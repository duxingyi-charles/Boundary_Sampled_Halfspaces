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
    void flip() override { is_flipped = !is_flipped; }

    bool has_control_points() const override { return true; }
    const std::vector<Point>& get_control_points() const override {
        if (m_control_pts.empty()) {
            m_control_pts.push_back(center);
            m_control_pts.push_back(center + Point(radius, 0, 0));
        }
        return m_control_pts;
    }
    void set_control_points(const std::vector<Point>& pts) override {
        if (pts.size() < 2) {
            std::cerr << "Sphere primitive expects at least 1 control points";
            return;
        }
        m_control_pts = pts;
        if ((center - pts[0]).squaredNorm() < 1e-12) {
            // Update radius
            radius = (pts[1] - pts[0]).norm();
        } else {
            // Update center;
            center = pts[0];
        }
        m_control_pts[1] = center + Point(radius, 0, 0);
    }

    std::string get_type() const override { return "sphere"; }

    void get_radius(double &r) const override { r = radius; };
    void get_is_flipped(bool &flip) const override { flip = is_flipped; };
    void get_center(Point &p) const override { p = center; };

    double get_radius() const { return radius; }
    void set_radius(double& r) { radius = r; }

    const Point& get_center() const { return center; }
    void set_center(const Point& p) { center = p; }

    bool save(const std::string &dir, const std::string &name, nlohmann::json &json_obj) const override;

    void translate(const Point& t) override {
        Sampled_Implicit::translate(t);
        center += t;
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
