#ifndef PLANE_S_IMPLICIT_H
#define PLANE_S_IMPLICIT_H

#include "Sampled_Implicit.h"
#include <Eigen/Geometry>

#include <iostream>

class Plane_sImplicit : public Sampled_Implicit
{
public:
    // standard plane
    Plane_sImplicit() : Sampled_Implicit(), p(0,0,0), normal(0,0,1) {};
    // plane from a point and a normal vector
    Plane_sImplicit(const Point &pt, const Eigen::Vector3d &n) : Sampled_Implicit(), p(pt), normal(n.normalized()) {};
    // plane from three points
	Plane_sImplicit(const Point &p1, const Point &p2, const Point &p3);
	~Plane_sImplicit() override = default;

	double function_at(const Point &x) const override { return normal.dot(x-p); }
	Eigen::Vector3d   gradient_at(const Point &x) const override { return normal; }

    bool has_control_points() const override { return true; }
    const std::vector<Point>& get_control_points() const override {
        if (m_control_pts.empty()) {
            m_control_pts.push_back(p);
        }
        return m_control_pts;
    }

    void set_control_points(const std::vector<Point>& pts) override {
        if (pts.size() < 1) {
            std::cerr << "Plane primitive expects at least 1 control points";
            return;
        }
        m_control_pts = pts;
        p = pts[0];
    }

private:
    Point p;
    Eigen::Vector3d normal;
    mutable std::vector<Point> m_control_pts;
};

#endif
