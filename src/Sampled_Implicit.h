#ifndef SAMPLED_INPLICIT_H
#define SAMPLED_INPLICIT_H

#include <vector>
#include "Vec3.h"

typedef Vec3 Point;

class Sampled_Implicit
{
public:
	Sampled_Implicit() = default;
	Sampled_Implicit(const std::vector<Point> &pts)
	: sample_points(pts) {};

	virtual ~Sampled_Implicit() = default;

	std::vector<Point> get_sample_points() const { return sample_points; }

	virtual double evaluate_at(const Point &p) const = 0;

protected:

	// void set_sample_points(const std::vector<Point> &pts) { sample_points = pts; } 

	std::vector<Point> sample_points;

};

#endif