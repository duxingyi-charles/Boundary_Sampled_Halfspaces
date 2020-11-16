//
// Created by Charles Du on 10/8/20.
//

#ifndef PSI_SPHERE_SIMPLICIT_H
#define PSI_SPHERE_SIMPLICIT_H

#include "Sampled_Implicit.h"

class Sphere_sImplicit : public Sampled_Implicit {
   public:
    // unit sphere
    Sphere_sImplicit()
        : Sampled_Implicit(), center(0, 0, 0), radius(1), is_flipped(false){};

    Sphere_sImplicit(const Point &o, double r, bool flip)
        : Sampled_Implicit(), center(o), radius(fabs(r)), is_flipped(flip){};

    Sphere_sImplicit(const Point &o, double r)
        : Sampled_Implicit(), center(o), radius(fabs(r)), is_flipped(false){};

    Sphere_sImplicit(const Point &o, const Point &s, bool flip)
        : Sampled_Implicit(),
          center(o),
          radius((o - s).norm()),
          is_flipped(flip){};

    Sphere_sImplicit(const Point &o, const Point &s)
        : Sampled_Implicit(),
          center(o),
          radius((o - s).norm()),
          is_flipped(false){};

    ~Sphere_sImplicit() override = default;

    double function_at(const Point &p) const override;
    Eigen::Vector3d gradient_at(const Point &p) const override;

   private:
    Point center;
    double radius;
    // is_flipped = false: f(x) = (||x-center||^2 - r^2)
    // is_flipped = true : f(x) = (r^2 - ||x-center||^2)
    bool is_flipped;
};

#endif  // PSI_SPHERE_SIMPLICIT_H
