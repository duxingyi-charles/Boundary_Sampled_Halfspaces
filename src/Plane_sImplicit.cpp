#include "Plane_sImplicit.h"

Plane_sImplicit::Plane_sImplicit(const Point &p1, const Point &p2, const Point &p3)
    : Sampled_Implicit()
    , p(p1)
{
    normal = (p2 - p1).cross(p3 - p1);
    normal.normalize();
}


bool Plane_sImplicit::save(const std::string &dir, const std::string &name, nlohmann::json &json_obj) const
{
    json_obj.clear();
    json_obj["type"] = "plane";
    json_obj["point"] = {p.x(), p.y(), p.z()};
    json_obj["normal"] = {normal.x(), normal.y(), normal.z()};
    json_obj["name"] = name;
    json_obj["color"] = {
        Sampled_Implicit::m_color[0],
        Sampled_Implicit::m_color[1],
        Sampled_Implicit::m_color[2],
        Sampled_Implicit::m_color[3]
    };

    std::string sample_filename = name + "_sample.xyz";
    json_obj["samples"] = sample_filename;

    std::string sample_file = dir + sample_filename;
    if (export_xyz(sample_file,sample_points)) {
        return true;
    } else {
        return false;
    }
}
