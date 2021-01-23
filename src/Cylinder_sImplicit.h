//
// Created by Charles Du on 10/9/20.
//

#ifndef PSI_CYLINDER_SIMPLICIT_H
#define PSI_CYLINDER_SIMPLICIT_H

#include "Sampled_Implicit.h"

#include <iostream>
#include <Eigen/Dense>

class Cylinder_sImplicit : public Sampled_Implicit
{
public:
    // standard cylinder
    Cylinder_sImplicit()
        : Sampled_Implicit()
        , axis_point(0, 0, 0)
        , axis_unit_vector(0, 0, 1)
        , radius(1)
        , is_flipped(false){};

    Cylinder_sImplicit(const Point &p, const Eigen::Vector3d &v, double r, bool flip)
        : Sampled_Implicit()
        , axis_point(p)
        , axis_unit_vector(v.normalized())
        , radius(fabs(r))
        , is_flipped(flip){};

    Cylinder_sImplicit(std::pair<Point, Point> p1p2, double r, bool flip)
        : Sampled_Implicit()
        , axis_point(p1p2.first)
        , axis_unit_vector((p1p2.second - p1p2.first).normalized())
        , radius(fabs(r))
        , is_flipped(flip){};

    ~Cylinder_sImplicit() override = default;

    double function_at(const Point &x) const override;
    Eigen::Vector3d gradient_at(const Point &x) const override;

    void flip() override { is_flipped = !is_flipped; }

    bool has_control_points() const override { return true; }
    const std::vector<Point> &get_control_points() const override
    {
        if (m_control_pts.empty()) {
            m_control_pts.push_back(axis_point);
            m_control_pts.push_back(axis_point + axis_unit_vector);
            if (std::abs(axis_unit_vector[1]) < 0.9) {
                Point dir = Point(0, 1, 0).cross(axis_unit_vector);
                m_control_pts.push_back(axis_point + radius * dir);
            } else {
                Point dir = Point(0, 0, 1).cross(axis_unit_vector);
                m_control_pts.push_back(axis_point + radius * dir);
            }
        }
        return m_control_pts;
    }
    void set_control_points(const std::vector<Point> &pts) override
    {
        if (pts.size() != 3) {
            std::cerr << "Cylinder primitive expects 3 control points";
            return;
        }
        m_control_pts = pts;
        if ((axis_point - pts[0]).norm() < 1e-6) {
            const Point dir = (pts[1] - pts[0]).normalized();
            if ((dir - axis_unit_vector).norm() < 1e-6) {
                radius = (pts[2] - pts[0]).norm();
            } else {
                axis_unit_vector = dir;
            }
        } else {
            axis_point = pts[0];
        }

        m_control_pts[1] = axis_point + axis_unit_vector;
        if (std::abs(axis_unit_vector[1]) < 0.9) {
            Point dir = Point(0, 1, 0).cross(axis_unit_vector);
            m_control_pts[2] = axis_point + radius * dir;
        } else {
            Point dir = Point(0, 0, 1).cross(axis_unit_vector);
            m_control_pts[2] = axis_point + radius * dir;
        }
    }

    std::string get_type() const override { return "cylinder"; }

    void get_axis_point(Point &p) const override { p = axis_point; };
    void get_axis_unit_vector(Eigen::Vector3d &vec) const override { vec = axis_unit_vector; };
    void get_radius(double &r) const override { r = radius; };
    void get_is_flipped(bool &flip) const override { flip = is_flipped; };

    bool save(const std::string &dir, const std::string &name, nlohmann::json &json_obj) const override;

private:
    // p: point on cylinder axis
    Point axis_point;
    // v: unit vector along axis
    Eigen::Vector3d axis_unit_vector;
    double radius;
    // is_flipped = false: f(x) = (||x - (p + dot(x-p,v)v)||^2 - r^2)
    // is_flipped = true : f(x) = (r^2 - ||x - (p + dot(x-p,v)v)||^2)
    bool is_flipped;

    mutable std::vector<Point> m_control_pts;
};

#endif // PSI_CYLINDER_SIMPLICIT_H
