#include "Cone_sImplicit.h"

double Cone_sImplicit::function_at(const Point &x) const {
    if (is_flipped) {
        return (x-apex).dot(axis_unit_vector) - cos(apex_angle) * (x-apex).norm();
    } else {
        return cos(apex_angle) * (x-apex).norm() - (x-apex).dot(axis_unit_vector);
    }
}

Eigen::Vector3d Cone_sImplicit::gradient_at(const Point &x) const {
    if (is_flipped) {
        return axis_unit_vector - cos(apex_angle) * (x-apex).normalized();
    } else {
        return cos(apex_angle) * (x-apex).normalized() - axis_unit_vector;
    }
}


bool Cone_sImplicit::save(const std::string &dir, const std::string &name, nlohmann::json &json_obj) const
{
    json_obj.clear();
    json_obj["type"] = "cone";
    json_obj["apex"] = {apex.x(), apex.y(), apex.z()};
    json_obj["axis_vector"] = {axis_unit_vector.x(), axis_unit_vector.y(), axis_unit_vector.z()};
    json_obj["apex_angle"] = apex_angle;
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
