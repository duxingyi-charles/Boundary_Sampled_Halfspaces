#ifndef SAMPLED_IMPLICIT_H
#define SAMPLED_IMPLICIT_H

#include <vector>
#include <Eigen/Core>


typedef Eigen::Vector3d Point;

class Sampled_Implicit
{
public:
	Sampled_Implicit() = default;
	explicit Sampled_Implicit(const std::vector<Point> &pts)
	: sample_points(pts) {};

	virtual ~Sampled_Implicit() = default;

	std::vector<Point> get_sample_points() const { return sample_points; }

	virtual double function_at(const Point &) const = 0;
	virtual Eigen::Vector3d gradient_at(const Point &) const = 0;

protected:

	std::vector<Point> sample_points;

};

#endif
