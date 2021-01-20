#pragma once
#include <Mesh_PSI.h>
#include <PSI.h>
#include <Sphere_sImplicit.h>
#include <Cylinder_sImplicit.h>
#include <Plane_sImplicit.h>
#include <Cone_sImplicit.h>

#include <Eigen/Core>

#include <memory>
#include <vector>

class PSIStates
{
public:
    using VertexArray = Eigen::Matrix<double, Eigen::Dynamic, 3>;
    using FaceArray = Eigen::Matrix<int, Eigen::Dynamic, 3>;
    using MapArray = Eigen::Matrix<int, Eigen::Dynamic, 1>;

public:
    PSIStates(const GridSpec& grid, std::vector<std::unique_ptr<Sampled_Implicit>>& implicits)
        : m_grid(grid)
        , m_implicits(implicits)
    {
        m_psi = std::make_unique<Mesh_PSI>();
        m_psi->run(grid, implicits);
        initialize_states();
        initialize_colors();
    }

    const VertexArray& get_vertices() const { return m_vertices; }
    const FaceArray& get_faces() const { return m_faces; }
    const auto& get_patches() const { return m_psi->get_patches(); }
    const auto& get_cells() const { return m_psi->get_cells(); }
    const auto& get_cell_labels() const { return m_psi->get_cell_labels(); }
    const auto& get_patch_labels() const { return m_psi->get_patch_labels(); }
    const auto& get_patch_implicit() const { return m_psi->get_patch_implicit(); }
    const auto& get_cell_patch_adjacency() const { return m_psi->get_cell_patch_adjacency(); }
    const size_t get_num_implicits() const { return m_psi->get_num_implicits(); }
    const size_t get_num_patches() const { return m_psi->get_num_patches(); }
    const size_t get_num_cells() const { return m_psi->get_num_cells(); }
    const auto& get_bbox() const { return m_bbox; }
    const auto get_grid_diag() const { return (m_grid.bbox_max - m_grid.bbox_min).norm(); }

    const int get_implicit_from_patch(int patch_id) const { return get_patch_implicit()[patch_id]; }
    const Eigen::RowVector4d get_implicit_color(int implicit_id) const
    {
        return m_colors.row(implicit_id);
    }

    void set_resolution(int res) {
        m_grid.resolution << res, res, res;
    }
    int get_resolution() const {
        return m_grid.resolution.maxCoeff();
    }

    void add_plane(const Point& p, const Point& n)
    {
        m_implicits.push_back(std::make_unique<Plane_sImplicit>(p, n));
        auto& fn = m_implicits.back();
        fn->set_sample_points({p});
        m_psi->run(m_grid, m_implicits);
        initialize_states();
        initialize_colors();
    }

    void add_sphere(const Point& center, double radius)
    {
        m_implicits.push_back(std::make_unique<Sphere_sImplicit>(center, radius));
        auto& fn = m_implicits.back();
        fn->set_sample_points({Point(center + Point(radius, 0, 0))});
        m_psi->run(m_grid, m_implicits);
        initialize_states();
        initialize_colors();
    }

    void add_cylinder(const Point& center, const Point& axis, double radius)
    {
        m_implicits.push_back(std::make_unique<Cylinder_sImplicit>(center, axis, radius, false));
        auto& fn = m_implicits.back();
        Point d = Point(0, 1, 0).cross(axis);
        if (d.norm() < 1e-3) {
            d = Point(0, 0, 1).cross(axis);
        }
        fn->set_sample_points({Point(center + d * radius)});
        m_psi->run(m_grid, m_implicits);
        initialize_states();
        initialize_colors();
    }

    void add_cone(const Point& center, const Point& axis, double angle)
    {
        m_implicits.push_back(std::make_unique<Cone_sImplicit>(center, axis, angle, false));
        auto& fn = m_implicits.back();
        fn->set_sample_points({center});
        m_psi->run(m_grid, m_implicits);
        initialize_states();
        initialize_colors();
    }

    Sampled_Implicit& get_implicit_function(size_t i) { return *m_implicits[i]; }

    void update_control_points(size_t implicit_id, const std::vector<Eigen::Vector3d>& pts)
    {
        m_implicits[implicit_id]->set_control_points(pts);
    }

    void update_sample_points(size_t implicit_id, const std::vector<Eigen::Vector3d>& pts)
    {
        m_implicits[implicit_id]->set_sample_points(pts);
        m_psi->process_samples();
        m_psi->graph_cut();
    }

    void refresh()
    {
        m_psi->run(m_grid, m_implicits);
        initialize_states();
        initialize_colors();
    }

private:
    void initialize_states()
    {
        const auto& V = m_psi->get_vertices();
        const auto& F = m_psi->get_faces();

        const size_t num_vertices = V.size();
        const size_t num_faces = F.size();

        m_vertices.resize(num_vertices, 3);
        for (size_t i = 0; i < num_vertices; i++) {
            m_vertices.row(i) = V[i];
        }
        m_bbox.row(0) = m_vertices.colwise().minCoeff();
        m_bbox.row(1) = m_vertices.colwise().maxCoeff();

        size_t num_triangles = 0;
        for (size_t i = 0; i < num_faces; i++) {
            if (F[i].size() == 3) {
                num_triangles++;
            } else {
                throw "Non-triangle face is not yet supported";
            }
        }

        m_faces.resize(num_triangles, 3);
        for (size_t i = 0; i < num_faces; i++) {
            m_faces(i, 0) = F[i][0];
            m_faces(i, 1) = F[i][1];
            m_faces(i, 2) = F[i][2];
        }
    }

    Eigen::Matrix<double, 1, 3> map_color(double t)
    {
        t = std::max<double>(std::min<double>(t, 1), 0);
        Eigen::Matrix<double, 8, 3> pallette;
        //pallette << 0x8d, 0xd3, 0xc7, 0xff, 0xff, 0xb3, 0xbe, 0xba, 0xda, 0xfb, 0x80, 0x72, 0x80,
        //    0xb1, 0xd3, 0xfd, 0xb4, 0x62, 0xb3, 0xde, 0x69, 0xfc, 0xcd, 0xe5;

        pallette <<
                0xe4, 0x1a, 0x1c,
                0x37, 0x7e, 0xb8,
                0x4d, 0xaf, 0x4a,
                0x98, 0x4e, 0xa3,
                0xff, 0x7f, 0x00,
                0xff, 0xff, 0x33,
                0xa6, 0x56, 0x28,
                0xf7, 0x81, 0xbf;
        size_t i0 = static_cast<size_t>(std::floor(t * 8));
        size_t i1 = static_cast<size_t>(std::ceil(t * 8));
        i0 = std::min<size_t>(i0, 7);
        i1 = std::min<size_t>(i1, 7);

        const double d = t * 8 - i0;
        return (pallette.row(i0) * (1 - d) + pallette.row(i1) * d) / 255.0;
    }

    void initialize_colors()
    {
        const auto num_implicits = m_implicits.size();
        m_colors.resize(num_implicits, 4);
        for (int i = 0; i < num_implicits; i++) {
            double t = double(i) / double(num_implicits - 1);
            m_colors.row(i) << map_color(t), 1;
        }

        // m_colors.setRandom(num_implicits, 4);
        // m_colors = m_colors.array() * 0.5 + 0.5;
        // m_colors.col(3).setConstant(0.5);
    }

private:
    // Input states.
    GridSpec m_grid;
    std::vector<std::unique_ptr<Sampled_Implicit>>& m_implicits;
    std::unique_ptr<PSI> m_psi;

    // Derived states.
    VertexArray m_vertices;
    FaceArray m_faces;
    Eigen::Matrix<double, 2, 3> m_bbox;
    Eigen::Matrix<double, Eigen::Dynamic, 4> m_colors;
};

