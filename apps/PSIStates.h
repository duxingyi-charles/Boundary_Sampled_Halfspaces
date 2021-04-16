#pragma once
#include <Cone_sImplicit.h>
#include <Cylinder_sImplicit.h>
#include <Hermite_RBF_sImplicit.h>
#include <Mesh_PSI.h>
#include <PSI.h>
#include <Plane_sImplicit.h>
#include <Sphere_sImplicit.h>
#include <Torus_sImplicit.h>

#include <Eigen/Core>

#include <igl/write_triangle_mesh.h>

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
        m_bbox.resize(2, 3);
        m_bbox.row(0) = grid.bbox_min;
        m_bbox.row(1) = grid.bbox_max;

        m_params.use_distance_weighted_area = false;
        m_params.use_state_space_graph_cut = true;
        m_params.topK = 1;
        m_params.consider_adj_diff = true;
        m_psi->set_parameters(m_params);

        refresh();
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
    const auto& get_implicit_meshes() const { return m_psi->get_input_meshes(); }

    const int get_implicit_from_patch(int patch_id) const { return get_patch_implicit()[patch_id]; }
    const Eigen::RowVector4d get_implicit_color(int implicit_id) const
    {
        return m_implicits[implicit_id]->get_color();
    }

    void set_resolution(int res) { m_grid.resolution << res, res, res; }
    int get_resolution() const { return m_grid.resolution.maxCoeff(); }
    int get_k() const { return m_params.topK; }
    void set_k(int k) { m_params.topK = k; }
    int get_max_search_count() const {
        constexpr int MAX_INT = std::numeric_limits<int>::max();
        if (m_params.max_search_count > 0 && m_params.max_search_count < MAX_INT) {
            return m_params.max_search_count;
        } else {
            return 0;
        }
    }
    void set_max_search_count(int c) {
        if (c <= 0) {
            m_params.max_search_count = std::numeric_limits<int>::max();
        } else {
            m_params.max_search_count = c;
        }
    }

    void add_plane(const Point& p, const Point& n)
    {
        m_implicits.push_back(std::make_unique<Plane_sImplicit>(p, n));
        auto& fn = m_implicits.back();
        fn->set_sample_points({p});
        fn->set_reference_length(get_grid_diag() / 10);
        fn->set_color(get_new_color());
        m_psi->add_implicit(m_grid, fn);
    }

    void add_sphere(const Point& center, double radius)
    {
        m_implicits.push_back(std::make_unique<Sphere_sImplicit>(center, radius));
        auto& fn = m_implicits.back();
        fn->set_sample_points({Point(center + Point(radius, 0, 0))});
        fn->set_reference_length(get_grid_diag() / 10);
        fn->set_color(get_new_color());
        m_psi->add_implicit(m_grid, fn);
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
        fn->set_reference_length(get_grid_diag() / 10);
        fn->set_color(get_new_color());
        m_psi->add_implicit(m_grid, fn);
    }

    void add_cone(const Point& center, const Point& axis, double angle)
    {
        m_implicits.push_back(std::make_unique<Cone_sImplicit>(center, axis, angle, false));
        auto& fn = m_implicits.back();
        fn->set_sample_points({center});
        fn->set_reference_length(get_grid_diag() / 10);
        fn->set_color(get_new_color());
        m_psi->add_implicit(m_grid, fn);
    }

    void add_torus(const Point& center, const Point& axis, double R, double r)
    {
        m_implicits.push_back(std::make_unique<Torus_sImplicit>(center, axis, R, r, false));
        auto& fn = m_implicits.back();
        Point d = Point(0, 1, 0).cross(axis);
        if (d.norm() < 1e-3) {
            d = Point(0, 0, 1).cross(axis);
        }
        fn->set_sample_points({center + d * (R + r)});
        fn->set_reference_length(get_grid_diag() / 10);
        fn->set_color(get_new_color());
        m_psi->add_implicit(m_grid, fn);
    }

    void add_vipss(const std::vector<Point>& ctrl_pts, const std::vector<Point>& sample_pts)
    {
        m_implicits.push_back(std::make_unique<Hermite_RBF_sImplicit>(ctrl_pts, sample_pts));
        auto& fn = m_implicits.back();
        fn->set_reference_length(get_grid_diag() / 10);
        fn->set_color(get_new_color());
        m_psi->add_implicit(m_grid, fn);
    }

    void remove_implicit(size_t i) { m_implicits.erase(m_implicits.begin() + i); }

    Sampled_Implicit& get_implicit_function(size_t i) { return *m_implicits[i]; }
    const Sampled_Implicit& get_implicit_function(size_t i) const { return *m_implicits[i]; }

    void update_implicit(size_t implicit_id)
    {
        m_psi->update_implicit(m_grid, m_implicits[implicit_id], implicit_id);
    }

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
        m_psi->set_parameters(m_params);
        m_psi->run(m_grid, m_implicits);
        initialize_states();
        initialize_reference_length();
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
        // m_bbox.row(0) = m_vertices.colwise().minCoeff();
        // m_bbox.row(1) = m_vertices.colwise().maxCoeff();

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

    void initialize_reference_length()
    {
        // TODO: reference length is for UI only.  Should not live inside the
        // implicit classes.
        for (auto& fn : m_implicits) {
            fn->set_reference_length(get_grid_diag() / 10);
        }
    }

    Eigen::Matrix<double, Eigen::Dynamic, 3> get_old_pallette() const
    {
        Eigen::Matrix<double, 8, 3> pallette;
        // pallette << 0x8d, 0xd3, 0xc7, 0xff, 0xff, 0xb3, 0xbe, 0xba, 0xda, 0xfb, 0x80, 0x72, 0x80,
        //    0xb1, 0xd3, 0xfd, 0xb4, 0x62, 0xb3, 0xde, 0x69, 0xfc, 0xcd, 0xe5;

        pallette << 0xe4, 0x1a, 0x1c, 0x37, 0x7e, 0xb8, 0x4d, 0xaf, 0x4a, 0x98, 0x4e, 0xa3, 0xff,
            0x7f, 0x00, 0xff, 0xff, 0x33, 0xa6, 0x56, 0x28, 0xf7, 0x81, 0xbf;
        return pallette;
    }

    Eigen::Matrix<double, Eigen::Dynamic, 3> get_rainbow_pallette() const
    {
        Eigen::Matrix<double, Eigen::Dynamic, 3> pallette(11, 3);
        pallette << 126, 41, 135, 82, 70, 156, 71, 98, 173, 45, 149, 196, 69, 174, 153, 113, 189,
            103, 164, 192, 60, 210, 183, 43, 239, 151, 33, 241, 91, 34, 237, 28, 36;
        return pallette;
    }

    Eigen::Matrix<double, 1, 3> map_color(double t) const
    {
        auto pallette = get_rainbow_pallette();
        const auto n = pallette.rows();
        t = std::max<double>(std::min<double>(t, 1), 0);
        size_t i0 = static_cast<size_t>(std::floor(t * n));
        size_t i1 = static_cast<size_t>(std::ceil(t * n));
        i0 = std::min<size_t>(i0, n - 1);
        i1 = std::min<size_t>(i1, n - 1);

        const double d = t * n - i0;
        return (pallette.row(i0) * (1 - d) + pallette.row(i1) * d) / 255.0;
    }

    void initialize_colors()
    {
        const auto num_implicits = m_implicits.size();
        for (size_t i = 0; i < num_implicits; i++) {
            auto& color = m_implicits[i]->get_color();
            if (color.template segment<3>(0).norm() == 0) {
                color = get_new_color(m_color_count);
                m_color_count++;
            }
        }
    }

    Eigen::RowVector4d get_new_color() const { return get_new_color(m_color_count++); }

    Eigen::RowVector4d get_new_color(size_t i) const
    {
        auto halton_sequence = [](int i, int b) -> double {
            double f = 1;
            double r = 0;
            while (i > 0) {
                f = f / b;
                r = r + f * (i % b);
                i = i / b;
            }
            return r;
        };
        Eigen::Matrix<double, 1, 4> color;
        color << map_color(halton_sequence(i, 3)), 1;
        return color;
    }

public:
    void save_all(const std::string& filename) const
    {
        save_output(filename);
        save_sample_points(filename);
        save_implicit_meshes(filename);
        save_implicits("psi_cache");
    }

    void save_sample_points(const std::string& filename) const
    {
        const auto pos = filename.rfind('.');
        const auto basename = filename.substr(0, pos);
        const auto ext = filename.substr(pos);

        const double l = get_grid_diag();
        const size_t num_implicits = get_num_implicits();
        for (size_t i = 0; i < num_implicits; i++) {
            const auto& pts = get_implicit_function(i).get_sample_points();
            std::vector<IGL_Mesh> point_meshes;
            point_meshes.reserve(pts.size());
            for (const auto& p : pts) {
                point_meshes.push_back(Mesh_PSI::generate_unit_sphere(8, 16, false));
                auto& mesh = point_meshes.back();
                mesh.vertices = mesh.vertices.array() * (l / 100);
                mesh.vertices.rowwise() += p.transpose();
            }

            IGL_Mesh point_mesh;
            Eigen::VectorXi face_to_mesh;
            Mesh_PSI::merge_meshes(point_meshes, point_mesh, face_to_mesh);

            igl::write_triangle_mesh(basename + "_sample_points_" + std::to_string(i) + ext,
                point_mesh.vertices,
                point_mesh.faces);
        }
    }

    void save_implicit_meshes(const std::string& filename) const
    {
        const auto pos = filename.rfind('.');
        const auto basename = filename.substr(0, pos);
        const auto ext = filename.substr(pos);
        const size_t num_implicits = get_num_implicits();

        for (size_t i = 0; i < num_implicits; i++) {
            const auto& mesh = get_implicit_meshes()[i + 1];
            igl::write_triangle_mesh(
                basename + "_implicit_" + std::to_string(i) + ext, mesh.vertices, mesh.faces);
        }
    }

    void save_output(const std::string& filename) const
    {
        const auto pos = filename.rfind('.');
        const auto basename = filename.substr(0, pos);
        const auto ext = filename.substr(pos);
        const size_t num_patches = get_num_patches();
        const size_t num_implicits = get_num_implicits();
        const auto& patches = get_patches();
        auto visible = get_patch_labels();
        assert(num_patches == patches.size());
        assert(num_patches == visible.size());

        const size_t num_faces = m_faces.rows();
        FaceArray visible_faces(num_faces, 3);
        FaceArray patch_faces(num_faces, 3);
        std::vector<FaceArray> visible_implicit_faces(num_implicits);
        std::vector<size_t> face_counts(num_implicits, 0);
        for (size_t i = 0; i < num_implicits; i++) {
            visible_implicit_faces[i].resize(num_faces, 3);
        }

        size_t count = 0;
        for (size_t i = 0; i < num_patches; i++) {
            if (!visible[i]) continue;

            size_t local_count = 0;
            const auto& patch = patches[i];
            const auto& implicit_id = get_implicit_from_patch(i);
            patch_faces.resize(patch.size(), 3);
            for (auto j : patch) {
                visible_faces.row(count) = m_faces.row(j);
                patch_faces.row(local_count) = m_faces.row(j);
                visible_implicit_faces[implicit_id].row(face_counts[implicit_id] + local_count) =
                    m_faces.row(j);

                count++;
                local_count++;
            }
            face_counts[implicit_id] += local_count;
            // igl::write_triangle_mesh(basename + "_patch_" + std::to_string(i) + ext, m_vertices,
            // patch_faces);
        }
        visible_faces.conservativeResize(count, 3);
        igl::write_triangle_mesh(filename, m_vertices, visible_faces);

        for (size_t i = 0; i < num_implicits; i++) {
            visible_implicit_faces[i].conservativeResize(face_counts[i], 3);
            igl::write_triangle_mesh(basename + "_visible_implicit_" + std::to_string(i) + ext,
                m_vertices,
                visible_implicit_faces[i]);
        }
    }

    void save_implicits(const std::string& output_dir) const
    {
        m_psi->export_sampled_implicits(output_dir);
    }

private:
    // Input states.
    GridSpec m_grid;
    PSI_Param m_params;
    std::vector<std::unique_ptr<Sampled_Implicit>>& m_implicits;
    std::unique_ptr<Mesh_PSI> m_psi;
    mutable int m_color_count = 0;

    // Derived states.
    VertexArray m_vertices;
    FaceArray m_faces;
    Eigen::Matrix<double, 2, 3> m_bbox;
};

