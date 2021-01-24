//
// Created by Charles Du on 10/9/20.
//

#ifndef PSI_CONE_SIMPLICIT_H
#define PSI_CONE_SIMPLICIT_H

#include "Sampled_Implicit.h"

#include <iostream>
#include <Eigen/Dense>


class Cone_sImplicit : public Sampled_Implicit
{
public:
    // standard cone
    Cone_sImplicit() : Sampled_Implicit(),
    apex(0,0,0), axis_unit_vector(0,0,-1), apex_angle(M_PI/4), is_flipped(false)
    {};

    Cone_sImplicit(const Point &p, const Eigen::Vector3d &v, double a, bool flip) : Sampled_Implicit(),
    apex(p), axis_unit_vector(v.normalized()), apex_angle(fmod(a,M_PI)), is_flipped(flip)
    {};

    ~Cone_sImplicit() override = default;

    std::string get_type() const override { return "cone"; }


    double function_at(const Point &x) const override;
    Eigen::Vector3d   gradient_at(const Point &x) const override;
    void flip() override { is_flipped = !is_flipped; }

    bool has_control_points() const override { return true; }
    const std::vector<Point>& get_control_points() const override {
        if (m_control_pts.empty()) {
            m_control_pts.push_back(apex);
            m_control_pts.push_back(apex + axis_unit_vector);
            if (std::abs(axis_unit_vector[1]) < 0.9) {
                Point dir = (Eigen::Vector3d(0, 1, 0).cross(axis_unit_vector.transpose())).transpose();
                m_control_pts.push_back(apex + axis_unit_vector + dir * std::sin(apex_angle));
            } else {
                Point dir = (Eigen::Vector3d(0, 0, 1).cross(axis_unit_vector.transpose())).transpose();
                m_control_pts.push_back(apex + axis_unit_vector + dir * std::sin(apex_angle));
            }
        }
        return m_control_pts;
    }
    void set_control_points(const std::vector<Point>& pts) override {
        if (pts.size() != 3) {
            std::cerr << "Cone primitive expects 3 control points";
            return;
        }
        m_control_pts = pts;
        if ((pts[0] - apex).norm() < 1e-6) {
            Point dir = pts[1] - pts[0];
            if ((dir - axis_unit_vector).norm() < 1e-6) {
                dir = pts[2] - pts[0];
                apex_angle = std::atan2(dir.cross(axis_unit_vector).norm(), dir.dot(axis_unit_vector));
            } else {
                axis_unit_vector = dir;
            }
        } else {
            apex = pts[0];
        }

        m_control_pts[1] = apex + axis_unit_vector;
        if (std::abs(axis_unit_vector[1]) < 0.9) {
            Point dir = (Eigen::Vector3d(0, 1, 0).cross(axis_unit_vector.transpose())).transpose();
            m_control_pts[2] = apex + axis_unit_vector + dir * std::sin(apex_angle);
        } else {
            Point dir = (Eigen::Vector3d(0, 0, 1).cross(axis_unit_vector.transpose())).transpose();
            m_control_pts[2] = apex + axis_unit_vector + dir * std::sin(apex_angle);
        }
    }

    void get_apex(Point &p) const override { p = apex; };
    void get_axis_unit_vector(Eigen::Vector3d &vec) const override { vec = axis_unit_vector; };
    void get_apex_angle(double &a) const override { a = apex_angle; };
    void get_is_flipped(bool &flip) const override { flip = is_flipped; };

    bool save(const std::string &dir, const std::string &name, nlohmann::json &json_obj) const override;

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

    mutable std::vector<Point> m_control_pts;
};

#endif //PSI_CONE_SIMPLICIT_H
