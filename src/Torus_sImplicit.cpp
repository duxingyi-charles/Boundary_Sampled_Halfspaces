//
// Created by Charles Du on 1/5/21.
//

#include "Torus_sImplicit.h"


double Torus_sImplicit::function_at(const Point &x) const {
    double val;
    Eigen::Vector3d vec = x - center;
    Eigen::Vector3d vec_para = vec.dot(axis_unit_vector) * axis_unit_vector;
    Eigen::Vector3d vec_perp = vec - vec_para;
    if (vec_perp.norm() == 0) { // point x lies on torus axis
        val = vec_para.squaredNorm() + major_radius * major_radius - minor_radius * minor_radius;
    } else {
        vec_perp.normalize();
        auto q = center + major_radius * vec_perp;
        val = (x - q).squaredNorm() - minor_radius * minor_radius;
    }
    //
    if (is_flipped) {
        return -val;
    } else {
        return val;
    }
}

Eigen::Vector3d Torus_sImplicit::gradient_at(const Point &x) const {
    Eigen::Vector3d grad;
    Eigen::Vector3d vec = x - center;
    Eigen::Vector3d vec_para = vec.dot(axis_unit_vector) * axis_unit_vector;
    Eigen::Vector3d vec_perp = vec - vec_para;
    if (vec_perp.norm() == 0) { // point x lies on torus axis
        grad = 2 * vec;
    } else {
        vec_perp.normalize();
        auto q = center + major_radius * vec_perp;
        grad = 2 * (x - q);
    }
    //
    if (is_flipped) {
        return -grad;
    } else {
        return grad;
    }

}


bool Torus_sImplicit::save(const std::string &dir, const std::string &name, nlohmann::json &json_obj) const
{
    json_obj.clear();
    json_obj["type"] = "torus";
    json_obj["center"] = {center.x(), center.y(), center.z()};
    json_obj["axis_vector"] = {axis_unit_vector.x(), axis_unit_vector.y(), axis_unit_vector.z()};
    json_obj["major_radius"] = major_radius;
    json_obj["minor_radius"] = minor_radius;
    json_obj["is_flipped"] = is_flipped;
    json_obj["name"] = name;

    std::string sample_filename = name + "_sample.xyz";
    json_obj["samples"] = sample_filename;

    std::string sample_file = dir + sample_filename;
    if (export_xyz(sample_file,sample_points)) {
        return true;
    } else {
        return false;
    }
}