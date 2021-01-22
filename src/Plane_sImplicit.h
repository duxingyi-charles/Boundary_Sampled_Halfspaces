#ifndef PLANE_S_IMPLICIT_H
#define PLANE_S_IMPLICIT_H

#include <Eigen/Geometry>
#include "Sampled_Implicit.h"

#include <iostream>

class Plane_sImplicit : public Sampled_Implicit
{
public:
    // standard plane
    Plane_sImplicit()
        : Sampled_Implicit()
        , p(0, 0, 0)
        , normal(0, 0, 1){};
    // plane from a point and a normal vector
    Plane_sImplicit(const Point &pt, const Eigen::Vector3d &n)
        : Sampled_Implicit()
        , p(pt)
        , normal(n.normalized()){};
    // plane from three points
    Plane_sImplicit(const Point &p1, const Point &p2, const Point &p3);
    ~Plane_sImplicit() override = default;

    double function_at(const Point &x) const override { return normal.dot(x - p); }
    Eigen::Vector3d gradient_at(const Point &x) const override { return normal; }
    void flip() override { normal *= -1; }

    bool has_control_points() const override { return true; }
    const std::vector<Point> &get_control_points() const override
    {
        if (m_control_pts.empty()) {
            m_control_pts.push_back(p);
            m_control_pts.push_back(p + normal);
        }
        return m_control_pts;
    }

    void set_control_points(const std::vector<Point> &pts) override
    {
        if (pts.size() < 2) {
            std::cerr << "Plane primitive expects 2 control points";
            return;
        }
        if ((p - pts[0]).norm() < 1e-6) {
            normal = (pts[1] - pts[0]).normalized();
        } else {
            p = pts[0];
        }
        m_control_pts = pts;
        m_control_pts[1] = p + normal;
    }

    // getter
    virtual void get_point(Point &pt) const override { pt = p; };
    virtual void get_normal(Eigen::Vector3d &n) const override { n = normal; };

    std::string get_type() const override { return "plane"; }

private:
    Point p;
    Eigen::Vector3d normal;
    mutable std::vector<Point> m_control_pts;
};

#endif
