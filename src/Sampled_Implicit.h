#ifndef SAMPLED_INPLICIT_H
#define SAMPLED_INPLICIT_H

#include <vector>
#include "Point.h"

class Sampled_Implicit
{
public:
	Sampled_Implicit(const std::vector<Point> &pts)
	: sample_points(pts) {}
	virtual ~Sampled_Implicit() = default;

	std::vector<Point> get_sample_points() const { return sample_points; }

	virtual double evaluate_at(const Point &p) const = 0;

protected:

	std::vector<Point> sample_points;

};

#endif