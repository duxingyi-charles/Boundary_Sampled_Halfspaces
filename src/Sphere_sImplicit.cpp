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

