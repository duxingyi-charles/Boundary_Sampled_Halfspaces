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