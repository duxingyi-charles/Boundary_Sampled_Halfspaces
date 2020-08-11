#ifndef PLANE_SIMPLICIT_H
#define PLANE_SIMPLICIT_H

#include "Sampled_Implicit.h"
#include "Vec3.h"

class Plane_sImplicit : public Sampled_Implicit
{
public:
	Plane_sImplicit(const Point &p1, const Point &p2, const Point &p3);
	~Plane_sImplicit() = default;

	double evaluate_at(const Point &p) const override;

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


double Plane_sImplicit::evaluate_at(const Point &p) {
	return dot(normal, p-p1);
}

#endif