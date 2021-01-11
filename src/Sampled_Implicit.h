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
	void set_sample_points(const std::vector<Point> &samples) { sample_points = samples; }

    virtual double function_at(const Point &) const = 0;
	virtual Eigen::Vector3d gradient_at(const Point &) const = 0;

    static bool import_xyz(const std::string &filename, std::vector<Point> &pts);

    virtual std::string get_type() const { return "unknown"; }

    // public getter
    virtual void get_point(Point&) const {};    //plane
    virtual void get_normal(Eigen::Vector3d&) const {};  //plane

    virtual void get_axis_point(Point &) const {};   //cylinder
    virtual void get_axis_unit_vector(Eigen::Vector3d &) const {}; //cylinder, cone, torus
    virtual void get_radius(double &) const {};     //cylinder
    virtual void get_is_flipped(bool &) const {};   //cylinder, cone, sphere, torus

    virtual void get_apex(Point&) const {};  //cone
    virtual void get_apex_angle(double&) const {}; //cone

    virtual void get_center(Point&) const {}; //torus
    virtual void get_major_radius(double&) const {}; //torus
    virtual void get_minor_radius(double&) const {}; //torus


protected:

	std::vector<Point> sample_points;

};

#endif
