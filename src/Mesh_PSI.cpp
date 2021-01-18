//
// Created by Charles Du on 11/5/20.
//

#include "Mesh_PSI.h"

#include <igl/copyleft/marching_cubes.h>

#include <limits>

#include <numeric> //std::partial_sum

#include "PyMesh/CellPartition.h"
#include "ScopedTimer.h"

#include <Eigen/Geometry> //cross product

IGL_Mesh Mesh_PSI::generate_cube(const GridSpec& grid_spec) {
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

IGL_Mesh Mesh_PSI::marching_cubes(Sampled_Implicit& fn, const GridSpec& grid_spec) {
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

void Mesh_PSI::merge_meshes(const std::vector<IGL_Mesh>& meshes,
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


void Mesh_PSI::compute_arrangement_for_graph_cut(
        const GridSpec &grid,
        const std::vector<std::unique_ptr<Sampled_Implicit>> &implicits) {
    // arrangement computation independent of samples
    compute_arrangement(grid, implicits);
    arrangement_ready = true;
    //
    process_samples();
}

IGL_Mesh Mesh_PSI::generate_plane(const GridSpec &grid, const Point &p, const Eigen::Vector3d &normal)
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


//IGL_Mesh Mesh_PSI::generate_cylinder(const GridSpec &grid, int n,
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


IGL_Mesh Mesh_PSI::generate_cylinder(const GridSpec &grid, int n,
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

IGL_Mesh Mesh_PSI::generate_cone(const GridSpec &grid, int n,
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
        base_center = apex + max_dist * axis_unit_vector;
    } else {
        base_center = apex - max_dist * axis_unit_vector;
    }

    //
    Eigen::Vector3d axis1 = axis_unit_vector.cross(Eigen::Vector3d(1,0,0));
    if (axis1.norm() == 0) {
        axis1 = axis_unit_vector.cross(Eigen::Vector3d(0,1,0));
    }
    axis1.normalize();
    Eigen::Vector3d axis2 = axis_unit_vector.cross(axis1);

    double radius = fabs(max_dist * tan(apex_angle));

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

IGL_Mesh Mesh_PSI::generate_torus(const GridSpec &grid, int n_major, int n_minor,
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

IGL_Mesh Mesh_PSI::generate_random_plane(const GridSpec &grid)
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

IGL_Mesh Mesh_PSI::generate_unit_sphere(bool is_flipped)
{
    IGL_Mesh sphere;
    sphere.vertices.resize(162, 3);
    sphere.vertices << 0.,0.,-1.,0.,0.,1.,-0.894427,0.,-0.447214,0.894427,0.,0.447214,0.723607,-0.525731,-0.447214,0.723607,
    0.525731,-0.447214,-0.723607,-0.525731,0.447214,-0.723607,0.525731,0.447214,-0.276393,-0.850651,-0.447214,-0.276393,0.850651,
    -0.447214,0.276393,-0.850651,0.447214,0.276393,0.850651,0.447214,0.951057,0.309017,0.,0.951057,-0.309017,0.,0.850651,0.,
    -0.525731,0.587785,0.809017,0.,0.688191,0.5,0.525731,0.,-1.,0.,0.262866,-0.809017,-0.525731,0.587785,-0.809017,0.,
    -0.262866,-0.809017,0.525731,-0.587785,-0.809017,0.,0.16246,-0.5,0.850651,0.688191,-0.5,0.525731,0.525731,0.,0.850651,
    0.16246,0.5,0.850651,0.262866,0.809017,-0.525731,0.425325,0.309017,-0.850651,-0.16246,0.5,-0.850651,0.425325,-0.309017,
    -0.850651,-0.525731,0.,-0.850651,-0.688191,0.5,-0.525731,-0.16246,-0.5,-0.850651,-0.688191,-0.5,-0.525731,-0.951057,
    0.309017,0.,-0.587785,0.809017,0.,-0.951057,-0.309017,0.,-0.850651,0.,0.525731,-0.262866,0.809017,0.525731,-0.425325,
    0.309017,0.850651,-0.425325,-0.309017,0.850651,0.,1.,0.,0.870463,0.433889,-0.232454,0.947214,0.16246,-0.276393,0.818274,
    0.273267,-0.505721,0.959253,0.160622,0.232454,0.959253,-0.160622,0.232454,1.,0.,0.,0.947214,-0.16246,-0.276393,0.870463,
    -0.433889,-0.232454,0.818274,-0.273267,-0.505721,0.861803,0.425325,0.276393,0.822619,0.259892,0.505721,0.68164,0.69378,
    -0.232454,0.809017,0.587785,0.,0.67082,0.688191,0.276393,0.449186,0.862668,0.232454,0.501375,0.702046,0.505721,0.143665,
    -0.961938,0.232454,0.309017,-0.951057,0.,0.449186,-0.862668,0.232454,-0.143665,-0.961938,-0.232454,-0.00703145,-0.862668,
    -0.505721,0.138197,-0.951057,-0.276393,0.447214,-0.850651,-0.276393,0.512752,-0.69378,-0.505721,0.68164,-0.69378,-0.232454,
    -0.309017,-0.951057,0.,-0.449186,-0.862668,-0.232454,0.00703145,-0.862668,0.505721,-0.138197,-0.951057,0.276393,-0.447214,
    -0.850651,0.276393,-0.512752,-0.69378,0.505721,-0.68164,-0.69378,0.232454,0.084444,-0.259892,0.961938,0.361803,-0.262866,
    0.894427,0.273267,0.,0.961938,0.228109,-0.702046,0.674609,0.501375,-0.702046,0.505721,0.447214,-0.525731,0.723607,0.638197,
    -0.262866,0.723607,0.822619,-0.259892,0.505721,0.738175,0.,0.674609,0.361803,0.262866,0.894427,0.084444,0.259892,0.961938,
    0.638197,0.262866,0.723607,0.447214,0.525731,0.723607,0.228109,0.702046,0.674609,-0.00703145,0.862668,-0.505721,0.0527864,
    0.688191,-0.723607,-0.228109,0.702046,-0.674609,0.512752,0.69378,-0.505721,0.597196,0.433889,-0.674609,0.361803,0.587785,
    -0.723607,0.138197,0.425325,-0.894427,0.221077,0.160622,-0.961938,-0.084444,0.259892,-0.961938,0.67082,0.16246,-0.723607,
    0.597196,-0.433889,-0.674609,0.67082,-0.16246,-0.723607,0.447214,0.,-0.894427,0.221077,-0.160622,-0.961938,-0.447214,
    0.525731,-0.723607,-0.501375,0.702046,-0.505721,-0.273267,0.,-0.961938,-0.361803,0.262866,-0.894427,-0.638197,0.262866,
    -0.723607,-0.738175,0.,-0.674609,-0.822619,0.259892,-0.505721,-0.084444,-0.259892,-0.961938,-0.361803,-0.262866,-0.894427,
    -0.228109,-0.702046,-0.674609,-0.501375,-0.702046,-0.505721,-0.447214,-0.525731,-0.723607,-0.638197,-0.262866,-0.723607,
    -0.822619,-0.259892,-0.505721,-0.959253,0.160622,-0.232454,-0.861803,0.425325,-0.276393,-0.870463,0.433889,0.232454,
    -0.68164,0.69378,0.232454,-0.809017,0.587785,0.,-0.67082,0.688191,-0.276393,-0.449186,0.862668,-0.232454,-0.947214,
    0.16246,0.276393,-0.818274,0.273267,0.505721,-0.959253,-0.160622,-0.232454,-1.,0.,0.,-0.947214,-0.16246,0.276393,
    -0.870463,-0.433889,0.232454,-0.818274,-0.273267,0.505721,0.00703145,0.862668,0.505721,-0.0527864,0.688191,0.723607,
    -0.512752,0.69378,0.505721,-0.597196,0.433889,0.674609,-0.361803,0.587785,0.723607,-0.138197,0.425325,0.894427,
    -0.221077,0.160622,0.961938,-0.447214,0.,0.894427,-0.221077,-0.160622,0.961938,-0.67082,0.16246,0.723607,-0.67082,
    -0.16246,0.723607,-0.597196,-0.433889,0.674609,0.138197,-0.425325,-0.894427,0.361803,-0.587785,-0.723607,0.0527864,
    -0.688191,-0.723607,-0.861803,-0.425325,-0.276393,-0.67082,-0.688191,-0.276393,-0.809017,-0.587785,0.,0.861803,
    -0.425325,0.276393,0.67082,-0.688191,0.276393,0.809017,-0.587785,0.,-0.0527864,-0.688191,0.723607,-0.138197,
    -0.425325,0.894427,-0.361803,-0.587785,0.723607,-0.143665,0.961938,-0.232454,0.138197,0.951057,-0.276393,0.143665,
    0.961938,0.232454,0.309017,0.951057,0.,0.447214,0.850651,-0.276393,-0.138197,0.951057,0.276393,-0.309017,0.951057,0.,
    -0.447214,0.850651,0.276393;

    sphere.faces.resize(320,3);
    if (is_flipped) {
        sphere.faces << 44,42,5,43,12,42,14,43,44,44,43,42,47,45,12,46,3,45,13,46,47,47,46,45,50,48,14,49,13,48,4,49,50,50,49,48,43,47,12,48,13,47,14,48,43,43,48,47,52,45,3,51,12,45,16,51,52,52,51,45,54,42,12,53,5,42,15,53,54,54,53,42,57,55,16,56,15,55,11,56,57,57,56,55,51,54,12,55,15,54,16,55,51,51,55,54,60,58,10,59,17,58,19,59,60,60,59,58,63,61,17,62,8,61,18,62,63,63,62,61,66,64,19,65,18,64,4,65,66,66,65,64,59,63,17,64,18,63,19,64,59,59,64,63,68,61,8,67,17,61,21,67,68,68,67,61,70,58,17,69,10,58,20,69,70,70,69,58,73,71,21,72,20,71,6,72,73,73,72,71,67,70,17,71,20,70,21,71,67,67,71,70,76,74,1,75,22,74,24,75,76,76,75,74,79,77,22,78,10,77,23,78,79,79,78,77,82,80,24,81,23,80,3,81,82,82,81,80,75,79,22,80,23,79,24,80,75,75,80,79,84,76,1,83,24,76,25,83,84,84,83,76,85,82,24,52,3,82,16,52,85,85,52,82,87,86,25,57,16,86,11,57,87,87,57,86,83,85,24,86,16,85,25,86,83,83,86,85,90,88,9,89,26,88,28,89,90,90,89,88,93,91,26,92,5,91,27,92,93,93,92,91,96,94,28,95,27,94,0,95,96,96,95,94,89,93,26,94,27,93,28,94,89,89,94,93,92,44,5,97,14,44,27,97,92,92,97,44,99,50,14,98,4,50,29,98,99,99,98,50,95,100,27,101,29,100,0,101,95,95,101,100,97,99,14,100,29,99,27,100,97,97,100,99,103,90,9,102,28,90,31,102,103,103,102,90,105,96,28,104,0,96,30,104,105,105,104,96,108,106,31,107,30,106,2,107,108,108,107,106,102,105,28,106,30,105,31,106,102,102,106,105,104,109,0,110,32,109,30,110,104,104,110,109,113,111,32,112,8,111,33,112,113,113,112,111,107,114,30,115,33,114,2,115,107,107,115,114,110,113,32,114,33,113,30,114,110,110,114,113,108,116,2,117,34,116,31,117,108,108,117,116,120,118,34,119,7,118,35,119,120,120,119,118,103,121,31,122,35,121,9,122,103,103,122,121,117,120,34,121,35,120,31,121,117,117,121,120,124,118,7,123,34,118,37,123,124,124,123,118,126,116,34,125,2,116,36,125,126,126,125,116,129,127,37,128,36,127,6,128,129,129,128,127,123,126,34,127,36,126,37,127,123,123,127,126,87,130,11,131,38,130,25,131,87,87,131,130,134,132,38,133,7,132,39,133,134,134,133,132,84,135,25,136,39,135,1,136,84,84,136,135,131,134,38,135,39,134,25,135,131,131,135,134,138,136,1,137,39,136,40,137,138,138,137,136,139,133,39,124,7,133,37,124,139,139,124,133,141,140,40,129,37,140,6,129,141,141,129,140,137,139,39,140,37,139,40,140,137,137,140,139,109,101,0,142,29,101,32,142,109,109,142,101,143,98,29,65,4,98,18,65,143,143,65,98,111,144,32,62,18,144,8,62,111,111,62,144,142,143,29,144,18,143,32,144,142,142,144,143,125,115,2,145,33,115,36,145,125,125,145,115,146,112,33,68,8,112,21,68,146,146,68,112,128,147,36,73,21,147,6,73,128,128,73,147,145,146,33,147,21,146,36,147,145,145,147,146,46,81,3,148,23,81,13,148,46,46,148,81,149,78,23,60,10,78,19,60,149,149,60,78,49,150,13,66,19,150,4,66,49,49,66,150,148,149,23,150,19,149,13,150,148,148,150,149,69,77,10,151,22,77,20,151,69,69,151,77,152,74,22,138,1,74,40,138,152,152,138,74,72,153,20,141,40,153,6,141,72,72,141,153,151,152,22,153,40,152,20,153,151,151,153,152,88,154,9,155,41,154,26,155,88,88,155,154,157,156,41,56,11,156,15,56,157,157,56,156,91,158,26,53,15,158,5,53,91,91,53,158,155,157,41,158,15,157,26,158,155,155,158,157,130,156,11,159,41,156,38,159,130,130,159,156,160,154,41,122,9,154,35,122,160,160,122,154,132,161,38,119,35,161,7,119,132,132,119,161,159,160,41,161,35,160,38,161,159,159,161,160;
    } else {
        sphere.faces << 5,42,44,42,12,43,44,43,14,42,43,44,12,45,47,45,3,46,47,46,13,45,46,47,14,48,50,48,13,49,50,49,4,48,49,50,12,47,43,47,13,48,43,48,14,47,48,43,3,45,52,45,12,51,52,51,16,45,51,52,12,42,54,42,5,53,54,53,15,42,53,54,16,55,57,55,15,56,57,56,11,55,56,57,12,54,51,54,15,55,51,55,16,54,55,51,10,58,60,58,17,59,60,59,19,58,59,60,17,61,63,61,8,62,63,62,18,61,62,63,19,64,66,64,18,65,66,65,4,64,65,66,17,63,59,63,18,64,59,64,19,63,64,59,8,61,68,61,17,67,68,67,21,61,67,68,17,58,70,58,10,69,70,69,20,58,69,70,21,71,73,71,20,72,73,72,6,71,72,73,17,70,67,70,20,71,67,71,21,70,71,67,1,74,76,74,22,75,76,75,24,74,75,76,22,77,79,77,10,78,79,78,23,77,78,79,24,80,82,80,23,81,82,81,3,80,81,82,22,79,75,79,23,80,75,80,24,79,80,75,1,76,84,76,24,83,84,83,25,76,83,84,24,82,85,82,3,52,85,52,16,82,52,85,25,86,87,86,16,57,87,57,11,86,57,87,24,85,83,85,16,86,83,86,25,85,86,83,9,88,90,88,26,89,90,89,28,88,89,90,26,91,93,91,5,92,93,92,27,91,92,93,28,94,96,94,27,95,96,95,0,94,95,96,26,93,89,93,27,94,89,94,28,93,94,89,5,44,92,44,14,97,92,97,27,44,97,92,14,50,99,50,4,98,99,98,29,50,98,99,27,100,95,100,29,101,95,101,0,100,101,95,14,99,97,99,29,100,97,100,27,99,100,97,9,90,103,90,28,102,103,102,31,90,102,103,28,96,105,96,0,104,105,104,30,96,104,105,31,106,108,106,30,107,108,107,2,106,107,108,28,105,102,105,30,106,102,106,31,105,106,102,0,109,104,109,32,110,104,110,30,109,110,104,32,111,113,111,8,112,113,112,33,111,112,113,30,114,107,114,33,115,107,115,2,114,115,107,32,113,110,113,33,114,110,114,30,113,114,110,2,116,108,116,34,117,108,117,31,116,117,108,34,118,120,118,7,119,120,119,35,118,119,120,31,121,103,121,35,122,103,122,9,121,122,103,34,120,117,120,35,121,117,121,31,120,121,117,7,118,124,118,34,123,124,123,37,118,123,124,34,116,126,116,2,125,126,125,36,116,125,126,37,127,129,127,36,128,129,128,6,127,128,129,34,126,123,126,36,127,123,127,37,126,127,123,11,130,87,130,38,131,87,131,25,130,131,87,38,132,134,132,7,133,134,133,39,132,133,134,25,135,84,135,39,136,84,136,1,135,136,84,38,134,131,134,39,135,131,135,25,134,135,131,1,136,138,136,39,137,138,137,40,136,137,138,39,133,139,133,7,124,139,124,37,133,124,139,40,140,141,140,37,129,141,129,6,140,129,141,39,139,137,139,37,140,137,140,40,139,140,137,0,101,109,101,29,142,109,142,32,101,142,109,29,98,143,98,4,65,143,65,18,98,65,143,32,144,111,144,18,62,111,62,8,144,62,111,29,143,142,143,18,144,142,144,32,143,144,142,2,115,125,115,33,145,125,145,36,115,145,125,33,112,146,112,8,68,146,68,21,112,68,146,36,147,128,147,21,73,128,73,6,147,73,128,33,146,145,146,21,147,145,147,36,146,147,145,3,81,46,81,23,148,46,148,13,81,148,46,23,78,149,78,10,60,149,60,19,78,60,149,13,150,49,150,19,66,49,66,4,150,66,49,23,149,148,149,19,150,148,150,13,149,150,148,10,77,69,77,22,151,69,151,20,77,151,69,22,74,152,74,1,138,152,138,40,74,138,152,20,153,72,153,40,141,72,141,6,153,141,72,22,152,151,152,40,153,151,153,20,152,153,151,9,154,88,154,41,155,88,155,26,154,155,88,41,156,157,156,11,56,157,56,15,156,56,157,26,158,91,158,15,53,91,53,5,158,53,91,41,157,155,157,15,158,155,158,26,157,158,155,11,156,130,156,41,159,130,159,38,156,159,130,41,154,160,154,9,122,160,122,35,154,122,160,38,161,132,161,35,119,132,119,7,161,119,132,41,160,159,160,35,161,159,161,38,160,161,159;
    }

    return sphere;
}

IGL_Mesh Mesh_PSI::generate_sphere(const GridSpec &grid,
                         const Point& center, double radius, bool is_flipped)
{
    IGL_Mesh sphere = generate_unit_sphere(is_flipped);

    sphere.vertices *= radius;
    for (int i = 0; i < sphere.vertices.rows(); ++i) {
        sphere.vertices.row(i) += center;
    }

    return sphere;
}

void Mesh_PSI::compute_arrangement(
        const GridSpec &grid,
        const std::vector<std::unique_ptr<Sampled_Implicit>> &implicits)
{
    // generate meshes
    std::vector<IGL_Mesh> meshes;
    meshes.reserve(implicits.size()+1);

    meshes.push_back(generate_cube(grid));
    double grid_diag = (grid.bbox_max - grid.bbox_min).norm();
    Eigen::Vector3d delta_diag((grid.bbox_max.x()-grid.bbox_min.x())/grid.resolution.x(),
                               (grid.bbox_max.y()-grid.bbox_min.y())/grid.resolution.y(),
                               (grid.bbox_max.z()-grid.bbox_min.z())/grid.resolution.z());
    double error_bound = delta_diag.norm();
    for (const auto& fn : implicits) {
        //
        if (fn->get_type() == "plane") {
            Point p;
            Eigen::Vector3d normal;
            fn->get_point(p);
            fn->get_normal(normal);
            meshes.push_back(generate_plane(grid, p, normal));
        } else if (fn->get_type() == "cylinder") {
            Point axis_point;
            Eigen::Vector3d axis_unit_vector;
            double radius;
            bool is_flipped;
            fn->get_axis_point(axis_point);
            fn->get_axis_unit_vector(axis_unit_vector);
            fn->get_radius(radius);
            fn->get_is_flipped(is_flipped);
            int n = (radius < error_bound) ? 15: (int)round(abs(M_PI/asin(0.5*error_bound/radius)));
//            std::cout << "cylinder n: " << n << std::endl;
            n = (n < 15) ? 15 : n;
            meshes.push_back(generate_cylinder(grid,n,axis_point,axis_unit_vector,radius,is_flipped));
        } else if (fn->get_type() == "cone") {
            Point apex;
            Eigen::Vector3d axis_unit_vector;
            double apex_angle;
            bool is_flipped;
            fn->get_apex(apex);
            fn->get_axis_unit_vector(axis_unit_vector);
            fn->get_apex_angle(apex_angle);
            fn->get_is_flipped(is_flipped);
            int n = (int)round(M_PI/abs(asin(0.5*error_bound/(grid_diag*abs(tan(apex_angle))))));
//            std::cout << "cone n: " << n << std::endl;
            n = (n < 15) ? 15 : n;
            meshes.push_back(generate_cone(grid,n,apex,axis_unit_vector,apex_angle,is_flipped));
        } else if (fn->get_type() == "sphere") {
            Point center;
            double radius;
            bool is_flipped;
            fn->get_center(center);
            fn->get_radius(radius);
            fn->get_is_flipped(is_flipped);
            meshes.push_back(generate_sphere(grid,center,radius,is_flipped));
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
            int n_major = (major_radius < error_bound) ? 15: (int)round(abs(M_PI/asin(0.5*error_bound/major_radius)));
            int n_minor = (minor_radius < error_bound) ? 15: (int)round(abs(M_PI/asin(0.5*error_bound/minor_radius)));
//            std::cout << "n_major: " << n_major << std::endl;
//            std::cout << "n_minor: " << n_minor << std::endl;
            n_major = (n_major < 15) ? 15 : n_major;
            n_minor = (n_minor < 15) ? 15 : n_minor;
            meshes.push_back(generate_torus(grid,n_major,n_minor,center,axis_unit_vector,major_radius,minor_radius,is_flipped));
        }
        else {
            meshes.push_back(marching_cubes(*fn, grid));
        }
//
//        meshes.push_back(marching_cubes(*fn, grid));
        // test: random planes
//        meshes.push_back(generate_random_plane(grid));
    }


    // merge meshes
    IGL_Mesh merged_mesh;
    Eigen::VectorXi face_to_mesh;
    merge_meshes(meshes,merged_mesh,face_to_mesh);


    // compute arrangement
    ScopedTimer<> timer("mesh arrangement for graph-cut");
    auto engine = PyMesh::CellPartition::create_raw(merged_mesh.vertices, merged_mesh.faces);
    engine->run();

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

    //find vertices/faces outside of bounding box
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

    // compute map: result face id -> input implicit id
    auto source_faces = engine->get_source_faces();
    Eigen::VectorXi result_face_to_implicit(source_faces.size());
    for (int i=0; i<result_face_to_implicit.size(); ++i) {
        // -1: grid bounding cube
        result_face_to_implicit(i) = face_to_mesh(source_faces(i)) - 1;
    }

    auto cells = engine->get_cells();
    int num_patch = cells.rows();

    std::vector<int> patch_impl(num_patch, -1);
    auto face_to_patch = engine->get_patches();

    for (int i=0; i<F.size(); i++) {
        int implicit_id = result_face_to_implicit(i);
        // the current face (belongs to OR resides outside) grid cube
        if (implicit_id == -1 || F_outside[i]) continue;
        int patch_id = face_to_patch(i);
        patch_impl[patch_id] = implicit_id;
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
    F_Impl.clear();
    F_Impl.resize(F.size(), -1);  // -1 for faces on/outside grid cube
    for (int i=0; i<F.size(); ++i) {
        int patch_id = P_old_to_new_index[face_to_patch(i)];
        if (patch_id != -1) {  // current patch not deleted
            P[patch_id].push_back(i);
            F_Impl[i] = P_Impl[patch_id];
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
