#include "Cylinder_sImplicit.h"

double Cylinder_sImplicit::function_at(const Point &x) const {
    if (is_flipped) {
        return radius * radius -
               (x - axis_point -
                (axis_unit_vector.dot(x - axis_point)) * axis_unit_vector)
                   .squaredNorm();
    } else {
        return (x - axis_point -
                (axis_unit_vector.dot(x - axis_point)) * axis_unit_vector)
                   .squaredNorm() -
               radius * radius;
    }
}

Eigen::Vector3d Cylinder_sImplicit::gradient_at(const Point &x) const {
    if (is_flipped) {
        return (axis_point +
                (axis_unit_vector.dot(x - axis_point)) * axis_unit_vector - x);
    } else {
        return (x - axis_point -
                (axis_unit_vector.dot(x - axis_point)) * axis_unit_vector);
    }
}


bool Cylinder_sImplicit::save(const std::string &dir, const std::string &name, nlohmann::json &json_obj) const
{
    json_obj.clear();
    json_obj["type"] = "cylinder";
    json_obj["axis_point1"] = {axis_point.x(), axis_point.y(), axis_point.z()};
    Point point2 = axis_point + axis_unit_vector;
    json_obj["axis_point2"] = {point2.x(), point2.y(), point2.z()};
    json_obj["radius"] = radius;
    json_obj["is_flipped"] = is_flipped;
    json_obj["name"] = name;
    json_obj["color"] = {
        Sampled_Implicit::m_color[0],
        Sampled_Implicit::m_color[1],
        Sampled_Implicit::m_color[2],
        Sampled_Implicit::m_color[3]
    };

    std::string sample_filename = name + "_sample.xyz";
    json_obj["samples"] = sample_filename;

    std::string sample_file = dir + sample_filename;
    if (export_xyz(sample_file,sample_points)) {
        return true;
    } else {
        return false;
    }
}
