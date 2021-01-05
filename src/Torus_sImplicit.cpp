//
// Created by Charles Du on 1/5/21.
//

#include "Torus_sImplicit.h"


double Torus_sImplicit::function_at(const Point &x) const {
    double val;
    auto vec = x - center;
    auto vec_para = vec.dot(axis_unit_vector) * axis_unit_vector;
    auto vec_perp = vec - vec_para;
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
    auto vec = x - center;
    auto vec_para = vec.dot(axis_unit_vector) * axis_unit_vector;
    auto vec_perp = vec - vec_para;
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