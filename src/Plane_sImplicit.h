#ifndef PLANE_S_IMPLICIT_H
#define PLANE_S_IMPLICIT_H

#include "Sampled_Implicit.h"

class Plane_sImplicit : public Sampled_Implicit
{
public:
	Plane_sImplicit(const Point &p1, const Point &p2, const Point &p3);
	~Plane_sImplicit() = default;

	double function_at(const Point &p) const override { return normal.dot(p-sample_points[0]); }
	Eigen::Vector3d   gradient_at(const Point &p) const override { return normal; }

protected:
	Eigen::Vector3d normal;
	
};


Plane_sImplicit::Plane_sImplicit(const Point &p1, const Point &p2, const Point &p3)
	: Sampled_Implicit(std::vector<Point>{p1,p2,p3})
{
    normal = (p2-p1).cross(p3-p1);
	normal.normalize();
}


#endif