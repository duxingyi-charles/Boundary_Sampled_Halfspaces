//
// Created by Charles Du on 1/5/21.
//

#ifndef PSI_TORUS_SIMPLICIT_H
#define PSI_TORUS_SIMPLICIT_H

#include "Sampled_Implicit.h"

#include <Eigen/Dense>
#include <iostream>

class Torus_sImplicit : public Sampled_Implicit
{
public:
    // standard torus
    Torus_sImplicit()
        : Sampled_Implicit()
        , center(0, 0, 0)
        , axis_unit_vector(0, 0, 1)
        , major_radius(1)
        , minor_radius(0.5)
        , is_flipped(false){};

    Torus_sImplicit(const Point &p, const Eigen::Vector3d &v, double R, double r, bool flip)
        : Sampled_Implicit()
        , center(p)
        , axis_unit_vector(v.normalized())
        , major_radius(fabs(R))
        , minor_radius(fmin(fabs(R), fabs(r)))
        , // minor radius can't be larger than major radius
        is_flipped(flip){};


    ~Torus_sImplicit() override = default;

    double function_at(const Point &x) const override;
    Eigen::Vector3d gradient_at(const Point &x) const override;
    void flip() override { is_flipped = !is_flipped; }

    bool has_control_points() const override { return true; }
    const std::vector<Point> &get_control_points() const override
    {
        if (m_control_pts.empty()) {
            m_control_pts.push_back(center);
            m_control_pts.push_back(center + axis_unit_vector * major_radius);
            Point dir;
            if (std::abs(axis_unit_vector[1]) < 0.9) {
                dir = Point(0, 1, 0).cross(axis_unit_vector).normalized();
            } else {
                dir = Point(0, 0, 1).cross(axis_unit_vector).normalized();
            }
            m_control_pts.push_back(center + dir * major_radius);
            m_control_pts.push_back(center + dir * (major_radius + minor_radius));
        }
        return m_control_pts;
    }
    void set_control_points(const std::vector<Point> &pts) override
    {
        if (pts.size() < 1) {
            std::cerr << "Torus primitive expects at least 1 control points";
            return;
        }
        Eigen::Transform<double, 3, Eigen::AffineCompact> transform;
        transform.setIdentity();
        if ((center - pts[0]).norm() < 1e-6) {
            Point dir = (pts[1] - pts[0]).normalized();
            if ((dir - axis_unit_vector).norm() < 1e-6) {
                auto R = (pts[2] - pts[0]).norm();
                auto r = (pts[3] - pts[2]).norm();
                if (std::abs(R - major_radius) < 1e-6) {
                    minor_radius = r;
                } else {
                    major_radius = R;
                }
            } else {
                Point axis = axis_unit_vector.cross(dir);
                double theta = std::atan2(axis.norm(), axis_unit_vector.dot(dir));
                transform.rotate(Eigen::AngleAxis<double>(theta, axis));
                axis_unit_vector = dir;
            }
        } else {
            transform.translate(pts[0] - center);
            center = pts[0];
        }
        m_control_pts = pts;
        m_control_pts[1] = center + axis_unit_vector * major_radius;

        // Reproject sample points.
        for (auto &p : Sampled_Implicit::sample_points) {
            p = transform * p;
            auto d = p - center;
            auto h = d.dot(axis_unit_vector);
            auto d2 = d - axis_unit_vector * h;
            auto q = center + d2.normalized() * major_radius;
            p = q + (p-q).normalized() * minor_radius;
        }


        Point dir;
        if (std::abs(axis_unit_vector[1]) < 0.9) {
            dir = Point(0, 1, 0).cross(axis_unit_vector).normalized();
        } else {
            dir = Point(0, 0, 1).cross(axis_unit_vector).normalized();
        }
        m_control_pts[2] = center + dir * major_radius;
        m_control_pts[3] = center + dir * (major_radius + minor_radius);
    }

    std::string get_type() const override { return "torus"; }

    void get_center(Point &p) const override { p = center; };
    void get_axis_unit_vector(Eigen::Vector3d &vec) const override { vec = axis_unit_vector; };
    void get_major_radius(double &R) const override { R = major_radius; };
    void get_minor_radius(double &r) const override { r = minor_radius; };
    void get_is_flipped(bool &flip) const override { flip = is_flipped; };

    bool save(
        const std::string &dir, const std::string &name, nlohmann::json &json_obj) const override;

    void translate(const Point &t) override
    {
        Sampled_Implicit::translate(t);
        center += t;
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


#endif // PSI_TORUS_SIMPLICIT_H
