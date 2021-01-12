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

