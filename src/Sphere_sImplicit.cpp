#include "Sphere_sImplicit.h"

double Sphere_sImplicit::function_at(const Point &p) const {
    if (is_flipped) {
        return (radius * radius - (p - center).squaredNorm());
    } else {
        return ((p - center).squaredNorm() - radius * radius);
    }
}

Eigen::Vector3d Sphere_sImplicit::gradient_at(const Point &p) const {
    if (is_flipped) {
        return 2 * (center - p);
    } else {
        return 2 * (p - center);
    }
}

bool Sphere_sImplicit::save(const std::string &dir, const std::string &name, nlohmann::json &json_obj) const
{
    json_obj.clear();
    json_obj["type"] = "sphere";
    json_obj["center"] = {center.x(), center.y(), center.z()};
    json_obj["radius"] = radius;
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