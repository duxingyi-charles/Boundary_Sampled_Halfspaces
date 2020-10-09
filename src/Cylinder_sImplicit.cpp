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
