#ifndef PLANE_SIMPLICIT_H
#define PLANE_SIMPLICIT_H

#include "Sampled_Implicit.h"
#include "Vec3.h"

class Plane_sImplicit : public Sampled_Implicit
{
public:
	Plane_sImplicit(const Point &p1, const Point &p2, const Point &p3);
	~Plane_sImplicit() = default;

	double function_at(const Point &p) const override;
	Vec3   gradient_at(const Point &p) const override { return normal; }

protected:
	Vec3 normal;
	
};


Plane_sImplicit::Plane_sImplicit(const Point &p1, const Point &p2, const Point &p3)
	: Sampled_Implicit()
{
	std::vector<Point> pts = {p1,p2,p3};
	sample_points = pts;

	normal = cross(p2-p1,p3-p1);
	normal.normalize();
}


double Plane_sImplicit::function_at(const Point &p) const {
	return dot(normal, p-sample_points[0]);
}

#endif