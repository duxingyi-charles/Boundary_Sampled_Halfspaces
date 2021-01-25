#ifndef SAMPLED_IMPLICIT_H
#define SAMPLED_IMPLICIT_H

#include <vector>
#include <Eigen/Core>

#include <nlohmann/json.hpp>


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
    void add_sample_point(const Point& p) { sample_points.push_back(p); }

    virtual bool has_control_points() const { return false; }
    virtual const std::vector<Point>& get_control_points() const { throw "Not supported yet!"; }
    virtual void set_control_points(const std::vector<Point>& pts) { }

    virtual double function_at(const Point &) const = 0;
	virtual Eigen::Vector3d gradient_at(const Point &) const = 0;
    virtual void flip() =0;

    static bool import_xyz(const std::string &filename, std::vector<Point> &pts);
    static bool export_xyz(const std::string &filename, const std::vector<Point> &pts);

    virtual bool save(const std::string &dir, const std::string &name, nlohmann::json &json_obj) const = 0;

    virtual std::string get_type() const { return "unknown"; }

    // public getter
    virtual void get_point(Point&) const {};    //plane
    virtual void get_normal(Eigen::Vector3d&) const {};  //plane

    virtual void get_axis_point(Point &) const {};   //cylinder
    virtual void get_axis_unit_vector(Eigen::Vector3d &) const {}; //cylinder, cone, torus
    virtual void get_radius(double &) const {};     //cylinder, sphere
    virtual void get_is_flipped(bool &) const {};   //cylinder, cone, sphere, torus

    virtual void get_apex(Point&) const {};  //cone
    virtual void get_apex_angle(double&) const {}; //cone

    virtual void get_center(Point&) const {}; //torus, sphere
    virtual void get_major_radius(double&) const {}; //torus
    virtual void get_minor_radius(double&) const {}; //torus

    virtual void translate(const Point& t) {
        for (auto& p : sample_points) {
            p += t;
        }
    }


protected:

	std::vector<Point> sample_points;

};

#endif
