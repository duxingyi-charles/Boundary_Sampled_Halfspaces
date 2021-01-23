//
// Created by Charles Du on 11/25/20.
//

#include "Quadric_sImplicit.h"

double Quadric_sImplicit::function_at(const Point &p) const {
    double x = p(0);
    double y = p(1);
    double z = p(2);
    // coef: coefficents for {1, x, y, z, x^2, x*y, x*z, y^2, y*z, z^2}
    return coef(0) + coef(1)*x + coef(2)*y + coef(3)*z
            + coef(4)*x*x + coef(5)*x*y + coef(6)*x*z
            + coef(7)*y*y + coef(8)*y*z + coef(9)*z*z;
}

Eigen::Vector3d  Quadric_sImplicit::gradient_at(const Point &p) const {
    double x = p(0);
    double y = p(1);
    double z = p(2);
    return Eigen::Vector3d(coef(1) + coef(4)*2*x + coef(5)*y + coef(6)*z,
                           coef(2) + coef(5)*x + coef(7)*2*y + coef(8)*z,
                           coef(3) + coef(6)*x + coef(8)*y + coef(9)*2*z);
}


bool Quadric_sImplicit::save(const std::string &dir, const std::string &name, nlohmann::json &json_obj) const
{
    json_obj.clear();
    json_obj["type"] = "quadric";
    json_obj["coef"] =  {
            coef(0),coef(1),coef(2),coef(3),coef(4),
            coef(5),coef(6),coef(7),coef(8),coef(9)
    };
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
