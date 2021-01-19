#include "Plane_sImplicit.h"

Plane_sImplicit::Plane_sImplicit(const Point &p1, const Point &p2, const Point &p3)
    : Sampled_Implicit()
    , p(p1)
{
    normal = (p2 - p1).cross(p3 - p1);
    normal.normalize();
}
