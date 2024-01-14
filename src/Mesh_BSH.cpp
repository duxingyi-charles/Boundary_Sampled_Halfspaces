//
// Created by Charles Du on 11/5/20.
//

#include "Mesh_BSH.h"

#include <igl/copyleft/marching_cubes.h>

#include <limits>

#include <numeric> //std::partial_sum

#include "PyMesh/Arrangement.h"
#include "ScopedTimer.h"

#include <Eigen/Geometry> //cross product

IGL_Mesh Mesh_BSH::generate_cube(const GridSpec& grid_spec) {
    ScopedTimer<> timer("Generate cube");
    constexpr double eps = 0;
    auto bbox_min = grid_spec.bbox_min.array() + eps;
    auto bbox_max = grid_spec.bbox_max.array() - eps;
    IGL_Mesh cube;
    cube.vertices.resize(8, 3);
    cube.vertices <<
                  bbox_min[0], bbox_min[1], bbox_max[2],
            bbox_min[0], bbox_min[1], bbox_min[2],
            bbox_max[0], bbox_min[1], bbox_min[2],
            bbox_max[0], bbox_min[1], bbox_max[2],
            bbox_min[0], bbox_max[1], bbox_max[2],
            bbox_max[0], bbox_max[1], bbox_max[2],
            bbox_max[0], bbox_max[1], bbox_min[2],
            bbox_min[0], bbox_max[1], bbox_min[2];

    cube.faces.resize(12, 3);
    cube.faces << 0, 1, 2, 0, 2, 3, 4, 5, 6, 4, 6, 7, 0, 3, 5, 0, 5, 4, 1, 7, 6,
            1, 6, 2, 2, 6, 5, 2, 5, 3, 0, 4, 7, 0, 7, 1;

    return cube;
}

IGL_Mesh Mesh_BSH::marching_cubes(Sampled_Implicit& fn, const GridSpec& grid_spec) {
    ScopedTimer<> timer("Marching cube");
    size_t num_grid_pts = (grid_spec.resolution[0] + 1) *
                          (grid_spec.resolution[1] + 1) *
                          (grid_spec.resolution[2] + 1);

    Eigen::Matrix<double, Eigen::Dynamic, 1> values(num_grid_pts);
    Eigen::Matrix<double, Eigen::Dynamic, 3> grid(num_grid_pts, 3);
    for (int i = 0; i <= grid_spec.resolution[2]; i++) {
        for (int j = 0; j <= grid_spec.resolution[1]; j++) {
            for (int k = 0; k <= grid_spec.resolution[0]; k++) {
                int idx = i * (grid_spec.resolution[0] + 1) *
                          (grid_spec.resolution[1] + 1) +
                          j * (grid_spec.resolution[0] + 1) + k;
                double x = double(k) / double(grid_spec.resolution[0]) *
                           (grid_spec.bbox_max[0] - grid_spec.bbox_min[0]) +
                           grid_spec.bbox_min[0];
                double y = double(j) / double(grid_spec.resolution[1]) *
                           (grid_spec.bbox_max[1] - grid_spec.bbox_min[1]) +
                           grid_spec.bbox_min[1];
                double z = double(i) / double(grid_spec.resolution[2]) *
                           (grid_spec.bbox_max[2] - grid_spec.bbox_min[2]) +
                           grid_spec.bbox_min[2];
                grid.row(idx) << x, y, z;
                values[idx] = fn.function_at({x, y, z});
            }
        }
    }

    IGL_Mesh mesh;
    igl::copyleft::marching_cubes(
            values, grid, grid_spec.resolution[0]+1, grid_spec.resolution[1]+1,
            grid_spec.resolution[2]+1, mesh.vertices, mesh.faces);
    return mesh;
}

void Mesh_BSH::merge_meshes(const std::vector<IGL_Mesh>& meshes,
                      // output
                      IGL_Mesh &merged_mesh,
                      Eigen::VectorXi &face_to_mesh) {
    ScopedTimer<> timer("Merge meshes");
    size_t num_meshes = meshes.size();
    std::vector<size_t> num_vertices(num_meshes + 1, 0);
    std::vector<size_t> num_faces(num_meshes + 1, 0);

    for (size_t i = 0; i < num_meshes; i++) {
        num_vertices[i + 1] = meshes[i].vertices.rows();
        num_faces[i + 1] = meshes[i].faces.rows();
    }

    std::partial_sum(num_vertices.begin(), num_vertices.end(),
                     num_vertices.begin());
    std::partial_sum(num_faces.begin(), num_faces.end(), num_faces.begin());

    merged_mesh.vertices.resize(num_vertices.back(), 3);
    merged_mesh.faces.resize(num_faces.back(), 3);

    // face id to mesh id
    face_to_mesh = Eigen::VectorXi::Ones(num_faces.back());


    for (size_t i = 0; i < num_meshes; i++) {
        merged_mesh.vertices.block(num_vertices[i], 0,
                                   num_vertices[i + 1] - num_vertices[i], 3) =
                meshes[i].vertices;
        merged_mesh.faces.block(num_faces[i], 0,
                                num_faces[i + 1] - num_faces[i], 3) =
                meshes[i].faces.array() + num_vertices[i];
        //
        face_to_mesh.segment(num_faces[i],num_faces[i+1]-num_faces[i]) *= i;
    }

}

void Mesh_BSH::add_implicit(const GridSpec &grid, const std::unique_ptr<Sampled_Implicit>& fn) {
    m_implicit_meshes.emplace_back();
    update_implicit(grid, fn, m_implicit_meshes.size()-2);
}

void Mesh_BSH::update_implicit(
        const GridSpec &grid,
        const std::unique_ptr<Sampled_Implicit>& fn, size_t i) {
    // Index of implicit goes from 1 to num_implicits.  Index 0 is for bbox.
    auto& meshes = m_implicit_meshes;
    assert(i+1 < meshes.size());

    double grid_diag = (grid.bbox_max - grid.bbox_min).norm();
    Eigen::Vector3d delta_diag((grid.bbox_max.x()-grid.bbox_min.x())/grid.resolution.x(),
                               (grid.bbox_max.y()-grid.bbox_min.y())/grid.resolution.y(),
                               (grid.bbox_max.z()-grid.bbox_min.z())/grid.resolution.z());
    double error_bound = delta_diag.norm();

    if (fn->get_type() == "plane") {
        Point p;
        Eigen::Vector3d normal;
        fn->get_point(p);
        fn->get_normal(normal);
        meshes[i+1] = generate_plane(grid, p, normal);
    } else if (fn->get_type() == "cylinder") {
        Point axis_point;
        Eigen::Vector3d axis_unit_vector;
        double radius;
        bool is_flipped;
        fn->get_axis_point(axis_point);
        fn->get_axis_unit_vector(axis_unit_vector);
        fn->get_radius(radius);
        fn->get_is_flipped(is_flipped);
        int n = (2 * radius < error_bound) ? 15: (int)round(abs(M_PI/asin(0.5*error_bound/radius)));
//            std::cout << "cylinder n: " << n << std::endl;
        n = (n < 15) ? 15 : n;
        meshes[i+1] = generate_cylinder(grid,n,axis_point,axis_unit_vector,radius,is_flipped);
    } else if (fn->get_type() == "cone") {
        Point apex;
        Eigen::Vector3d axis_unit_vector;
        double apex_angle;
        bool is_flipped;
        fn->get_apex(apex);
        fn->get_axis_unit_vector(axis_unit_vector);
        fn->get_apex_angle(apex_angle);
        fn->get_is_flipped(is_flipped);
//        int n = (2 * (grid_diag*abs(tan(apex_angle))) < error_bound) ? 15 : (int)round(M_PI/abs(asin(0.5*error_bound/(grid_diag*abs(tan(apex_angle))))));
        int n = (2 * (grid_diag*abs(sin(apex_angle))) < error_bound) ? 15 : (int)round(M_PI/abs(asin(0.5*error_bound/(grid_diag*abs(sin(apex_angle))))));
//            std::cout << "cone n: " << n << std::endl;
        n = (n < 15) ? 15 : n;
        meshes[i+1] = generate_cone(grid,n,apex,axis_unit_vector,apex_angle,is_flipped);
    } else if (fn->get_type() == "sphere") {
        Point center;
        double radius;
        bool is_flipped;
        fn->get_center(center);
        fn->get_radius(radius);
        fn->get_is_flipped(is_flipped);
        int n_latitude = (2 * radius < error_bound) ? 15 : (int)round(abs(M_PI/asin(0.5*error_bound/radius)));
        int n_longitude = (2 * radius < error_bound) ? 15 : (int)round(abs(0.25*M_PI/asin(0.5*error_bound/radius)));
        n_latitude = (n_latitude < 15) ? 15 : n_latitude;
        n_longitude = (n_longitude < 15) ? 15 : n_longitude;
        meshes[i+1] = generate_sphere(grid,n_latitude, n_longitude, center,radius,is_flipped);
    } else if (fn->get_type() == "torus") {
        Point center;
        Eigen::Vector3d axis_unit_vector;
        double major_radius;
        double minor_radius;
        bool is_flipped;
        fn->get_center(center);
        fn->get_axis_unit_vector(axis_unit_vector);
        fn->get_major_radius(major_radius);
        fn->get_minor_radius(minor_radius);
        fn->get_is_flipped(is_flipped);
        int n_major = (2 * major_radius < error_bound) ? 15: (int)round(abs(M_PI/asin(0.5*error_bound/major_radius)));
        int n_minor = (2 * minor_radius < error_bound) ? 15: (int)round(abs(M_PI/asin(0.5*error_bound/minor_radius)));
        n_major = (n_major < 15) ? 15 : n_major;
        n_minor = (n_minor < 15) ? 15 : n_minor;
        meshes[i+1] = generate_torus(grid,n_major,n_minor,center,axis_unit_vector,major_radius,minor_radius,is_flipped);
    }
    else {
        meshes[i+1] = marching_cubes(*fn, grid);
    }
}

void Mesh_BSH::compute_arrangement_for_graph_cut(
        const GridSpec &grid,
        const std::vector<std::unique_ptr<Sampled_Implicit>> &implicits) {
    // arrangement computation independent of samples
    generate_meshes(grid, implicits);
    compute_arrangement(grid, implicits);
    arrangement_ready = true;
    //
    process_samples();
}

IGL_Mesh Mesh_BSH::generate_plane(const GridSpec &grid, const Point &p, const Eigen::Vector3d &normal)
{
    auto pmin = grid.bbox_min;
    auto pmax = grid.bbox_max;
    auto center = (pmin + pmax) /2;

    Eigen::Vector3d xAxi(pmax.x()-center.x(),0,0);
    Eigen::Vector3d yAxi(0,pmax.y()-center.y(),0);
    Eigen::Vector3d zAxi(0,0,pmax.z()-center.z());

    // find the corner farthest to point p
    double max_dist = 0;
    for (int i = 0; i < 2; ++i) {
        int sx = (i == 0 ? -1 : 1);
        for (int j = 0; j < 2; ++j) {
            int sy = (j == 0 ? -1 : 1);
            for (int k = 0; k < 2; ++k) {
                int sz = (k == 0 ? -1 : 1);
                auto q_ijk = center + sx * xAxi + sy * yAxi + sz * zAxi;
                double dist = (p-q_ijk).norm();
                if (dist > max_dist) {
                    max_dist = dist;
                }
            }
        }
    }
    max_dist *= 1.001;

//    Eigen::Vector3d axis1 = p.cross(p+normal);
    Eigen::Vector3d axis1 = normal.cross(xAxi);
    if (axis1.norm() == 0) {
        axis1 = normal.cross(yAxi);
    }
    axis1.normalize();
    Eigen::Vector3d axis2 = normal.cross(axis1);

//    std::cout << "plane center: " << p << std::endl;
//    std::cout << "max_dist: " << max_dist << std::endl;

    Point p0 = p + max_dist * (-axis1-axis2);
    Point p1 = p + max_dist * (axis1 -axis2);
    Point p2 = p + max_dist * (axis1 +axis2);
    Point p3 = p + max_dist * (-axis1+axis2);

    IGL_Mesh plane;
    plane.vertices.resize(4, 3);
    plane.vertices.row(0) = p0;
    plane.vertices.row(1) = p1;
    plane.vertices.row(2) = p2;
    plane.vertices.row(3) = p3;

//    std::cout << "point0: " << p0 << std::endl;


    plane.faces.resize(2, 3);
    plane.faces << 0, 1, 2, 0, 2, 3;

    return plane;
}


//IGL_Mesh Mesh_BSH::generate_cylinder(const GridSpec &grid, int n,
//                           const Point& axis_point, const Eigen::Vector3d &axis_unit_vector,
//                           double radius, bool is_flipped)
//{
//    auto pmin = grid.bbox_min;
//    auto pmax = grid.bbox_max;
//    double diag_dist = (pmax - pmin).norm();
//
//    // find intersection of the cylinder axis and the bounding box
//    double px = axis_point.x();
//    double py = axis_point.y();
//    double pz = axis_point.z();
//
//    double vx = axis_unit_vector.x();
//    double vy = axis_unit_vector.y();
//    double vz = axis_unit_vector.z();
//
//    double pmin_x = pmin.x();
//    double pmin_y = pmin.y();
//    double pmin_z = pmin.z();
//
//    double pmax_x = pmax.x();
//    double pmax_y = pmax.y();
//    double pmax_z = pmax.z();
//
//    //
//    bool find_intersection = false;
//    Point q(0,0,0);
//    while (true) {
//        if (vx != 0) {
//            // large x
//            double t = (pmax_x - px)/vx;
//            q = axis_point + t * axis_unit_vector;
//            if (pmin_y <= q.y() <= pmax_y && pmin_z <= q.z() <= pmax_z) {
//                find_intersection = true;
//                break;
//            }
//            // small x
//            t = (pmin_x - px)/vx;
//            q = axis_point + t * axis_unit_vector;
//            if (pmin_y <= q.y() <= pmax_y && pmin_z <= q.z() <= pmax_z) {
//                find_intersection = true;
//                break;
//            }
//        }
//
//        if (vy != 0) {
//            // large y
//            double t = (pmax_y - py)/vy;
//            q = axis_point + t * axis_unit_vector;
//            if (pmin_x <= q.x() <= pmax_x && pmin_z <= q.z() <= pmax_z) {
//                find_intersection = true;
//                break;
//            }
//            // small y
//            t = (pmin_y - py)/vy;
//            q = axis_point + t * axis_unit_vector;
//            if (pmin_x <= q.x() <= pmax_x && pmin_z <= q.z() <= pmax_z) {
//                find_intersection = true;
//                break;
//            }
//        }
//
//        if (vz != 0) {
//            // large z
//            double t = (pmax_z - pz)/vz;
//            q = axis_point + t * axis_unit_vector;
//            if (pmin_x <= q.x() <= pmax_x && pmin_y <= q.y() <= pmax_y) {
//                find_intersection = true;
//                break;
//            }
//            // small z
//            t = (pmin_z - pz)/vz;
//            q = axis_point + t * axis_unit_vector;
//            if (pmin_x <= q.x() <= pmax_x && pmin_y <= q.y() <= pmax_y) {
//                find_intersection = true;
//                break;
//            }
//        }
//        // find_intersection = false
//        break;
//    }
//    if (!find_intersection) {
//        if (vx != 0) {
//            double t = (pmax_x - px)/vx;
//            q = axis_point + t * axis_unit_vector;
//        } else if (vy != 0) {
//            double t = (pmax_y - py)/vy;
//            q = axis_point + t * axis_unit_vector;
//        } else { // vz != 0
//            double t = (pmax_z - pz)/vz;
//            q = axis_point + t * axis_unit_vector;
//        }
//    }
//
//    std::cout << "q:" << q << std::endl;
//
//
//
//    //
//    Point p1 = q - 2 * diag_dist * axis_unit_vector;
////    Point p2 = q + 2 * diag_dist * axis_unit_vector;
//    std::cout << "p1: " << p1 << std::endl;
//
//    //
//    Eigen::Vector3d axis1 = axis_unit_vector.cross(Eigen::Vector3d(1,0,0));
//    if (axis1.norm() == 0) {
//        axis1 = axis_unit_vector.cross(Eigen::Vector3d(0,1,0));
//    }
//    axis1.normalize();
//    Eigen::Vector3d axis2 = axis_unit_vector.cross(axis1);
//
//    IGL_Mesh cylinder;
//    cylinder.vertices.resize(2 * n, 3);
//    for (int i = 0; i < n; ++i) {
//        cylinder.vertices.row(i) = p1 + radius * cos(1.0 * i * 2 * M_PI / n) * axis1
//                + radius * sin(1.0 * i * 2 * M_PI / n) * axis2;
//    }
//    for (int i = 0; i < n; ++i) {
////        cylinder.vertices.row(n+i) = cylinder.vertices.row(i) + 2 * diag_dist * axis_unit_vector;
//        cylinder.vertices.row(n+i) = cylinder.vertices.row(i);
//        cylinder.vertices(n+i,0) += 2 * diag_dist * axis_unit_vector.x();
//        cylinder.vertices(n+i,1) += 2 * diag_dist * axis_unit_vector.y();
//        cylinder.vertices(n+i,2) += 2 * diag_dist * axis_unit_vector.z();
//    }
//
//    // flip or not
//    Point q0 = cylinder.vertices.row(0);
//    Point q1 = cylinder.vertices.row(1);
//    Point q2 = cylinder.vertices.row(n);
//    Eigen::Vector3d face_normal = (q1-q0).cross(q2-q0);
//    Eigen::Vector3d radius_vector = q0 - p1;
//    bool need_flip = false;
//
//    if (radius_vector.dot(face_normal) > 0) {
//        if (is_flipped) need_flip = true;
//    } else {
//        if (!is_flipped) need_flip = true;
//    }
//
//    cylinder.faces.resize(2 * n, 3);
//    int num_vert = 2*n;
//    for (int i = 0; i < n; ++i) {
//        if (need_flip) {
//            cylinder.faces <<  (i+n % num_vert), (i+1 % num_vert), (i % num_vert);
//            cylinder.faces <<  ((i+n+1) % num_vert), ((i+1) % num_vert), ((i+n) % num_vert);
//        } else {
//            cylinder.faces << (i % num_vert), (i+1 % num_vert), (i+n % num_vert);
//            cylinder.faces << ((i+n) % num_vert), ((i+1) % num_vert), ((i+n+1) % num_vert);
//        }
//
//    }
//
//    return cylinder;
//}


IGL_Mesh Mesh_BSH::generate_cylinder(const GridSpec &grid, int n,
                                     const Point& axis_point, const Eigen::Vector3d &axis_unit_vector,
                                     double radius, bool is_flipped)
{
    n = (n < 3) ? 3 : n;
    auto pmin = grid.bbox_min;
    auto pmax = grid.bbox_max;
    auto center = (pmin + pmax) /2;

    Eigen::Vector3d xAxi(pmax.x()-center.x(),0,0);
    Eigen::Vector3d yAxi(0,pmax.y()-center.y(),0);
    Eigen::Vector3d zAxi(0,0,pmax.z()-center.z());

    // find the corner farthest to axis_point
    double max_dist = 0;
    for (int i = 0; i < 2; ++i) {
        int sx = (i == 0 ? -1 : 1);
        for (int j = 0; j < 2; ++j) {
            int sy = (j == 0 ? -1 : 1);
            for (int k = 0; k < 2; ++k) {
                int sz = (k == 0 ? -1 : 1);
                auto q_ijk = center + sx * xAxi + sy * yAxi + sz * zAxi;
                double dist = (axis_point-q_ijk).norm();
                if (dist > max_dist) {
                    max_dist = dist;
                }
            }
        }
    }
    max_dist *= 1.001;


    // center of the base disk of cylinder
    Point p1 = axis_point - max_dist * axis_unit_vector;
//    std::cout << "p1: " << p1 << std::endl;

    //
    Eigen::Vector3d axis1 = axis_unit_vector.cross(Eigen::Vector3d(1,0,0));
    if (axis1.norm() == 0) {
        axis1 = axis_unit_vector.cross(Eigen::Vector3d(0,1,0));
    }
    axis1.normalize();
    Eigen::Vector3d axis2 = axis_unit_vector.cross(axis1);

    IGL_Mesh cylinder;
    cylinder.vertices.resize(2 * n, 3);
    for (int i = 0; i < n; ++i) { // base circle
        cylinder.vertices.row(i) = p1 + radius * cos(i * 2 * M_PI / n) * axis1
                                   + radius * sin(i * 2 * M_PI / n) * axis2;
    }
    for (int i = 0; i < n; ++i) { // top circle
//        cylinder.vertices.row(n+i) = cylinder.vertices.row(i) + 2 * diag_dist * axis_unit_vector;
        cylinder.vertices.row(n+i) = cylinder.vertices.row(i);
        cylinder.vertices(n+i,0) += 2 * max_dist * axis_unit_vector.x();
        cylinder.vertices(n+i,1) += 2 * max_dist * axis_unit_vector.y();
        cylinder.vertices(n+i,2) += 2 * max_dist * axis_unit_vector.z();
    }

    // flip orientation or not
    Point q0 = cylinder.vertices.row(0);
    Point q1 = cylinder.vertices.row(1);
    Point q2 = cylinder.vertices.row(n);
    Eigen::Vector3d face_normal = (q1-q0).cross(q2-q0);
    Eigen::Vector3d radius_vector = q0 - p1;
    bool need_flip = false;

    if (radius_vector.dot(face_normal) > 0) {
        if (is_flipped) need_flip = true;
    } else {
        if (!is_flipped) need_flip = true;
    }

    cylinder.faces.resize(2 * n, 3);
    for (int i = 0; i < n; ++i) {
        if (need_flip) {
            cylinder.faces.row(i) << (i+n), (i+1)%n, (i);
            cylinder.faces.row(i+n) << (i+1)%n + n, (i+1)%n, (i+n);
        } else {
            cylinder.faces.row(i) << (i), (i+1)%n, (i+n);
            cylinder.faces.row(i+n) << (i+n), (i+1)%n, (i+1)%n + n;
        }
    }
//    std::cout << cylinder.faces << std::endl;

    return cylinder;
}

IGL_Mesh Mesh_BSH::generate_cone(const GridSpec &grid, int n,
                       const Point& apex, const Eigen::Vector3d& axis_unit_vector,
                       double apex_angle, bool is_flipped)
{
    if (cos(apex_angle) == 0) {
        if (is_flipped) {
            return generate_plane(grid,apex,axis_unit_vector);
        } else {
            return generate_plane(grid,apex,-axis_unit_vector);
        }
    }

    n = (n < 3) ? 3 : n;
    auto pmin = grid.bbox_min;
    auto pmax = grid.bbox_max;
    auto center = (pmin + pmax) /2;

    Eigen::Vector3d xAxi(pmax.x()-center.x(),0,0);
    Eigen::Vector3d yAxi(0,pmax.y()-center.y(),0);
    Eigen::Vector3d zAxi(0,0,pmax.z()-center.z());

    // find the corner farthest to apex
    double max_dist = 0;
    for (int i = 0; i < 2; ++i) {
        int sx = (i == 0 ? -1 : 1);
        for (int j = 0; j < 2; ++j) {
            int sy = (j == 0 ? -1 : 1);
            for (int k = 0; k < 2; ++k) {
                int sz = (k == 0 ? -1 : 1);
                auto q_ijk = center + sx * xAxi + sy * yAxi + sz * zAxi;
                double dist = (apex-q_ijk).norm();
                if (dist > max_dist) {
                    max_dist = dist;
                }
            }
        }
    }
    max_dist *= 1.001;

    Point base_center(0,0,0);
    if (cos(apex_angle) > 0) {
//        base_center = apex + max_dist * axis_unit_vector;
        base_center = apex + max_dist * cos(apex_angle) * axis_unit_vector;
    } else {
//        base_center = apex - max_dist * axis_unit_vector;
        base_center = apex - max_dist * cos(apex_angle) * axis_unit_vector;
    }

    //
    Eigen::Vector3d axis1 = axis_unit_vector.cross(Eigen::Vector3d(1,0,0));
    if (axis1.norm() == 0) {
        axis1 = axis_unit_vector.cross(Eigen::Vector3d(0,1,0));
    }
    axis1.normalize();
    Eigen::Vector3d axis2 = axis_unit_vector.cross(axis1);

//    double radius = fabs(max_dist * tan(apex_angle));
    double radius = fabs(max_dist * sin(apex_angle));

    IGL_Mesh cone;
    cone.vertices.resize(n+1, 3);
    if (is_flipped) {
        for (int i = 0; i < n; ++i) {
            cone.vertices.row(i) = base_center + radius * cos(i*2*M_PI/n) * axis1
                                    + radius * sin(i*2*M_PI/n) * axis2;
        }
    } else {
        for (int i = 0; i < n; ++i) {
            cone.vertices.row(i) = base_center + radius * cos(-i*2*M_PI/n) * axis1
                                   + radius * sin(-i*2*M_PI/n) * axis2;
        }
    }
    cone.vertices.row(n) = apex;

    cone.faces.resize(n,3);
    for (int i = 0; i < n; ++i) {
        cone.faces.row(i) << n, i, (i+1)%n;
    }

    return cone;
}

IGL_Mesh Mesh_BSH::generate_torus(const GridSpec &grid, int n_major, int n_minor,
                            const Point& center, const Eigen::Vector3d& axis_unit_vector,
                            double major_radius, double minor_radius,
                            bool is_flipped)
{
    Eigen::Vector3d axis1 = axis_unit_vector.cross(Eigen::Vector3d(1,0,0));
    if (axis1.norm() == 0) {
        axis1 = axis_unit_vector.cross(Eigen::Vector3d(0,1,0));
    }
    axis1.normalize();
    Eigen::Vector3d axis2 = axis_unit_vector.cross(axis1);

    IGL_Mesh torus;
    torus.vertices.resize(n_major*n_minor, 3);

    for (int i = 0; i < n_major; ++i) {
        Point m_center = center + major_radius * cos(i*2*M_PI/n_major) * axis1
                         + major_radius * sin(i*2*M_PI/n_major) * axis2;
        Eigen::Vector3d m_axis = (m_center - center);
        m_axis.normalize();

        for (int j = 0; j < n_minor; ++j) {
            torus.vertices.row(i*n_minor + j) = m_center
                                                    + minor_radius * cos(j*2*M_PI/n_minor) * m_axis
                                                    + minor_radius * sin(j*2*M_PI/n_minor) * axis_unit_vector;
        }
    }

    torus.faces.resize(n_major*n_minor*2, 3);
    if (is_flipped) {
        for (int i = 0; i < n_major; ++i) {
            for (int j = 0; j < n_minor; ++j) {
                torus.faces.row(i*n_minor+j) << ((i+1)%n_major)*n_minor+j,
                        i*n_minor+j,
                        i*n_minor+((j+1)%n_minor);
                torus.faces.row(i*n_minor+j+n_major*n_minor) << ((i+1)%n_major)*n_minor+j,
                        i*n_minor+((j+1)%n_minor),
                        ((i+1)%n_major)*n_minor+((j+1)%n_minor);
            }
        }
    } else {
        for (int i = 0; i < n_major; ++i) {
            for (int j = 0; j < n_minor; ++j) {
                torus.faces.row(i*n_minor+j) << i*n_minor+j,
                                                ((i+1)%n_major)*n_minor+j,
                                                i*n_minor+((j+1)%n_minor);
                torus.faces.row(i*n_minor+j+n_major*n_minor) << i*n_minor+((j+1)%n_minor),
                                                                ((i+1)%n_major)*n_minor+j,
                                                                ((i+1)%n_major)*n_minor+((j+1)%n_minor);
            }
        }
    }

    return torus;

}

IGL_Mesh Mesh_BSH::generate_random_plane(const GridSpec &grid)
{
    auto pmin = grid.bbox_min;
    auto pmax = grid.bbox_max;
    auto resolution = grid.resolution;

    int rand_x = rand() % resolution.x();
    double x = pmin.x() + (1.0 * rand_x / resolution.x()) * (pmax.x() - pmin.x());
    int rand_y = rand() % resolution.y();
    double y = pmin.y() + (1.0 * rand_y / resolution.y()) * (pmax.y() - pmin.y());
    int rand_z = rand() % resolution.z();
    double z = pmin.z() + (1.0 * rand_z / resolution.z()) * (pmax.z() - pmin.z());

    Point center(x,y,z);
    std::cout << "plane center: " << center << std::endl;

    Eigen::Vector3d normal((double)rand()/RAND_MAX, (double)rand()/RAND_MAX, (double)rand()/RAND_MAX);
    normal.normalize();

    double max_dist = (pmax - pmin).norm() * 1.001;
    std::cout << "max_dist: " << max_dist << std::endl;

    Eigen::Vector3d axis1 = center.cross(center+normal);
    axis1.normalize();
    Eigen::Vector3d axis2 = normal.cross(axis1);

    Point p0 = center + max_dist * (-axis1-axis2);
    Point p1 = center + max_dist * (axis1 -axis2);
    Point p2 = center + max_dist * (axis1 +axis2);
    Point p3 = center + max_dist * (-axis1+axis2);

    std::cout << "point0: " << p0 << std::endl;

    IGL_Mesh plane;
    plane.vertices.resize(4, 3);
    plane.vertices.row(0) = p0;
    plane.vertices.row(1) = p1;
    plane.vertices.row(2) = p2;
    plane.vertices.row(3) = p3;


    plane.faces.resize(2, 3);
    plane.faces << 0, 1, 2, 0, 2, 3;


    return plane;

}

IGL_Mesh Mesh_BSH::generate_unit_sphere(int n_latitude, int n_longitude, bool is_flipped)
{
    IGL_Mesh sphere;
    sphere.vertices.resize((2*n_longitude-1)*n_latitude+2, 3);
    
    double latitude_delta = 2 * M_PI / n_latitude;
    double longtitude_delta = M_PI / (2 * n_longitude);

    for (int i = 1; i < 2*n_longitude; ++i) {
        double longtitude_angle = -M_PI/2 + i * longtitude_delta;
        for (int j = 0; j < n_latitude; ++j) {
            double latitude_angle = j * latitude_delta;
            sphere.vertices.row((i-1)*n_latitude + j) << cos(longtitude_angle) * cos(latitude_angle),
                                                            cos(longtitude_angle) * sin(latitude_angle),
                                                            sin(longtitude_angle);
        }
    }
    int south_pole = (2*n_longitude-1)*n_latitude;
    sphere.vertices.row(south_pole) << 0, 0, -1;
    int north_pole = (2*n_longitude-1)*n_latitude + 1;
    sphere.vertices.row(north_pole) << 0, 0, 1;
    

    sphere.faces.resize((2*n_longitude-1)*2*n_latitude,3);
    if (is_flipped) {
        for (int i = 1; i < 2*n_longitude-1; ++i) {
            for (int j = 0; j < n_latitude; ++j) {
                sphere.faces.row((i-1)*2*n_latitude + j) << i*n_latitude+j,
                        (i-1)*n_latitude + (j+1)%n_latitude,
                        (i-1)*n_latitude+j;

                sphere.faces.row((i-1)*2*n_latitude + n_latitude + j) << i*n_latitude + (j+1)%n_latitude,
                        (i-1)*n_latitude + (j+1)%n_latitude,
                        i*n_latitude+j;
            }
        }

        for (int j = 0; j < n_latitude; ++j) {
            sphere.faces.row((2*n_longitude-2)*2*n_latitude + j) << (j+1)%n_latitude, south_pole, j;
        }

        for (int j = 0; j < n_latitude; ++j) {
            sphere.faces.row((2*n_longitude-2)*2*n_latitude + n_latitude + j) << (2*n_longitude-2)*n_latitude + (j+1)%n_latitude,
                    (2*n_longitude-2)*n_latitude + j,
                    north_pole;
        }
    } else {
        for (int i = 1; i < 2*n_longitude-1; ++i) {
            for (int j = 0; j < n_latitude; ++j) {
                sphere.faces.row((i-1)*2*n_latitude + j) << (i-1)*n_latitude+j,
                                                            (i-1)*n_latitude + (j+1)%n_latitude,
                                                            i*n_latitude+j;

                sphere.faces.row((i-1)*2*n_latitude + n_latitude + j) << i*n_latitude+j,
                        (i-1)*n_latitude + (j+1)%n_latitude,
                        i*n_latitude + (j+1)%n_latitude;
            }
        }

        for (int j = 0; j < n_latitude; ++j) {
            sphere.faces.row((2*n_longitude-2)*2*n_latitude + j) << j, south_pole, (j+1)%n_latitude;
        }

        for (int j = 0; j < n_latitude; ++j) {
            sphere.faces.row((2*n_longitude-2)*2*n_latitude + n_latitude + j) << north_pole,
                    (2*n_longitude-2)*n_latitude + j,
                    (2*n_longitude-2)*n_latitude + (j+1)%n_latitude;
        }
    }

    return sphere;
}

IGL_Mesh Mesh_BSH::generate_sphere(const GridSpec &grid, int n_latitude, int n_longitude,
                         const Point& center, double radius, bool is_flipped)
{
    IGL_Mesh sphere = generate_unit_sphere(n_latitude, n_longitude, is_flipped);

    sphere.vertices *= radius;
    for (int i = 0; i < sphere.vertices.rows(); ++i) {
        sphere.vertices.row(i) += center;
    }

    return sphere;
}

void Mesh_BSH::generate_meshes(
        const GridSpec &grid,
        const std::vector<std::unique_ptr<Sampled_Implicit>> &implicits) {
    auto& meshes = m_implicit_meshes;
    meshes.clear();
    meshes.resize(implicits.size()+1);

    meshes[0] = generate_cube(grid);
    size_t count = 0;
    for (const auto& fn : implicits) {
        update_implicit(grid, fn, count);
        count++;
    }
}

void Mesh_BSH::compute_arrangement(
        const GridSpec &grid,
        const std::vector<std::unique_ptr<Sampled_Implicit>> &implicits)
{
    // generate meshes
    // merge meshes
    IGL_Mesh merged_mesh;
    Eigen::VectorXi face_to_mesh;
    merge_meshes(m_implicit_meshes,merged_mesh,face_to_mesh);


    // compute arrangement
    ScopedTimer<> timer("mesh arrangement for graph-cut");
    auto engine = PyMesh::Arrangement::create_fast_arrangement(
        merged_mesh.vertices, merged_mesh.faces, face_to_mesh);
    {
        ScopedTimer<> timer("mesh arrangement");
        engine->run();
    }

    auto vertices = engine->get_vertices();
    auto faces = engine->get_faces();

    //copy vertices to V
    V.clear();
    for (int i = 0; i < vertices.rows(); ++i) {
        V.emplace_back(vertices(i,0),vertices(i,1),vertices(i,2));
    }
    //copy faces to F
    F.clear();
    F.resize(faces.rows());
    for (int i = 0; i < faces.rows(); ++i) {
        for (int j = 0; j < faces.cols(); ++j) {
            F[i].push_back(faces(i,j));
        }
    }

    //find vertices/faces outside bounding box
    const Point& bbox_min = grid.bbox_min;
    const Point& bbox_max = grid.bbox_max;
    std::vector<bool> V_outside(V.size(), false);
    for (int i = 0; i < V.size(); ++i) {
        double x = V[i].x();
        double y = V[i].y();
        double z = V[i].z();
        if (x < bbox_min.x() || x > bbox_max.x()
            || y < bbox_min.y() || y > bbox_max.y()
            || z < bbox_min.z() || z > bbox_max.z()){
            V_outside[i] = true;
        }
    }
    std::vector<bool> F_outside(F.size(), false);
    for (int i = 0; i < F.size(); ++i) {
        for (auto vi : F[i]) {
            if (V_outside[vi]) {
                F_outside[i] = true;
                break;
            }
        }
    }

    // compute face to implicit mapping.  -1 means mapped to bounding box.
    Eigen::VectorXi result_face_to_implicit = engine->get_out_face_labels().array() -1;

    auto cells = engine->get_cells();
    int num_patch = cells.rows();

    std::vector<int> patch_impl(num_patch, -1);
    auto face_to_patch = engine->get_patches();

    std::vector<bool> vertex_on_bbox(V.size(), false);
    double eps_x = (bbox_max.x()-bbox_min.x()) / 1000;
    double eps_y = (bbox_max.y()-bbox_min.y()) / 1000;
    double eps_z = (bbox_max.z()-bbox_min.z()) / 1000;
    for (int i = 0; i < V.size(); ++i) {
        double x = V[i].x();
        double y = V[i].y();
        double z = V[i].z();
        if (fabs(x - bbox_min.x())< eps_x || fabs(x - bbox_max.x()) < eps_x
            || fabs(y - bbox_min.y())< eps_y || fabs(y - bbox_max.y()) < eps_y
            || fabs(z - bbox_min.z())< eps_z || fabs(z - bbox_max.z()) < eps_z){
            vertex_on_bbox[i] = true;
        }
    }


    for (int i=0; i<F.size(); i++) {
        int implicit_id = result_face_to_implicit(i);
        // the current face (belongs to OR resides outside) grid cube
//        if (implicit_id == -1 || F_outside[i]) continue;
        if (implicit_id == -1) {
            for (auto vi : F[i]) {
                vertex_on_bbox[vi] = true;
            }
            continue;
        }
        if (F_outside[i]) continue;
        int patch_id = face_to_patch(i);
        patch_impl[patch_id] = implicit_id;
    }

    // label faces that touch the bounding box
    std::vector<bool> face_touch_bbox(F.size(), false);
    for (int i = 0; i < F.size(); ++i) {
        for (auto vi : F[i]) {
            if (vertex_on_bbox[vi]) {
                face_touch_bbox[i] = true;
                break;
            }
        }
    }

    //blocks incident to each patch
    std::vector<std::vector<int>> patch_block;
    //sign of implicit on blocks incident to each patch
    std::vector<std::vector<int>> patch_sign;

    for (int i=0; i < num_patch; ++i) {
        patch_block.emplace_back(std::vector<int> {cells(i, 0), cells(i, 1)});
        patch_sign.emplace_back(std::vector<int> {1, -1});
    }

    // remove patches from the grid cube OR outside the grid cube
    std::vector<int> P_old_to_new_index(patch_impl.size());

    P_block.clear();
    P_sign.clear();
    P_Impl.clear();

    int patch_count = 0;
    for (int i=0; i < patch_impl.size(); ++i) {
        if (patch_impl[i] != -1) {
            P_old_to_new_index[i] = patch_count;
            P_block.emplace_back(patch_block[i]);
            P_sign.emplace_back(patch_sign[i]);
            P_Impl.push_back(patch_impl[i]);
            patch_count++;
        }
        else {
            P_old_to_new_index[i] = -1; //patch deleted
        }
    }
    num_patch = P_Impl.size();

    // P_block: make sure the minimal index is 0
    int min_block_index = std::numeric_limits<int>::max();
    int max_block_index = -1;
    for (auto &blocks : P_block) {
        for (auto b : blocks) {
            min_block_index = (b < min_block_index) ? b : min_block_index;
            max_block_index = (b > max_block_index) ? b : max_block_index;
        }
    }
    if (min_block_index > 0) {
        for (auto &blocks : P_block) {
            for (auto &b : blocks) {
                b -= min_block_index;
            }
        }
    }
    int num_block = max_block_index - min_block_index + 1;


    // collect faces belonging to each patches
    P.clear();
    P.resize(num_patch);
    P_touch_bbox.clear();
    P_touch_bbox.resize(num_patch, false);
    F_Impl.clear();
    F_Impl.resize(F.size(), -1);  // -1 for faces on/outside grid cube
    for (int i=0; i<F.size(); ++i) {
        int patch_id = P_old_to_new_index[face_to_patch(i)];
        if (patch_id != -1) {  // current patch not deleted
            P[patch_id].push_back(i);
            F_Impl[i] = P_Impl[patch_id];
            if (face_touch_bbox[i]) {
                P_touch_bbox[patch_id] = true;
            }
        }
    }

    // collect patches belonging to each block
    B_patch.clear();
    B_patch.resize(num_block);
    for (int i=0; i < P_block.size(); ++i) {
        for (int j=0; j < P_block[i].size(); ++j) {
            B_patch[P_block[i][j]].push_back(i);
        }
    }

}
