#ifndef PLANE_S_IMPLICIT_H
#define PLANE_S_IMPLICIT_H

#include "Sampled_Implicit.h"
#include <Eigen/Geometry>

class Plane_sImplicit : public Sampled_Implicit
{
public:
    // standard plane
    Plane_sImplicit() : Sampled_Implicit(), p(0,0,0), normal(0,0,1) {};
    // plane from a point and a normal vector
    Plane_sImplicit(const Point &pt, const Eigen::Vector3d &n): Sampled_Implicit(), p(pt), normal(n.normalized()) {};
    // plane from three points
	Plane_sImplicit(const Point &p1, const Point &p2, const Point &p3);
	~Plane_sImplicit() override = default;

	double function_at(const Point &x) const override { return normal.dot(x-p); }
	Eigen::Vector3d   gradient_at(const Point &x) const override { return normal; }

	// getter
    virtual void get_point(Point& pt) const override { pt = p;};
    virtual void get_normal(Eigen::Vector3d& n) const override { n = normal; };

    std::string get_type() const override { return "plane"; }

private:
    Point p;
    Eigen::Vector3d normal;
};


Plane_sImplicit::Plane_sImplicit(const Point &p1, const Point &p2, const Point &p3)
	: Sampled_Implicit(), p(p1)
{
    normal = (p2-p1).cross(p3-p1);
	normal.normalize();
}


#endif
