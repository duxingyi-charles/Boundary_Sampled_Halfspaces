#include <PSI.h>
#include <config.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/project.h>
#include <igl/unproject_on_plane.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/write_triangle_mesh.h>

#include <CLI/CLI.hpp>
#include <Eigen/Core>
#include <array>
#include <limits>
#include <map>
#include <string>

#include "PSIStates.h"
#include "ScopedTimer.h"

class MeshArrangementMenu
{
public:
    MeshArrangementMenu()
        : m_states(nullptr)
    {}

    void set_states(PSIStates* states)
    {
        assert(states != nullptr);
        m_states = states;
    }

    void initialize(igl::opengl::glfw::Viewer& viewer)
    {
        initialize_data(viewer); // Data must be initialized first.
        initialize_manu(viewer);
        initialize_picking(viewer);
    }

private:
    void reset_patch_visibility(igl::opengl::glfw::Viewer& viewer) {
        static auto set_patch_visible = [&](auto& viewer, size_t i, bool val) {
            assert(i < m_data_ids.size());
            m_patch_visible[i] = val;
            auto pid = m_data_ids[i];
            viewer.data(pid).set_visible(val);
            assert(viewer.core().is_set(viewer.data(pid).is_visible) == val);
        };

        const auto& patches = m_states->get_patches();
        const size_t num_patches = patches.size();

        m_patch_visible = m_states->get_patch_labels();

        for (size_t i = 0; i < num_patches; i++) {
            set_patch_visible(viewer, i, m_patch_visible[i]);
        }
    }

    void initialize_manu(igl::opengl::glfw::Viewer& viewer)
    {
        //static auto is_patch_visible = [&](auto& viewer, size_t i) -> bool {
        //    assert(i < m_data_ids.size());
        //    auto pid = m_data_ids[i];
        //    return viewer.core().is_set(viewer.data(pid).is_visible);
        //};

        static auto set_patch_visible = [&](auto& viewer, size_t i, bool val) {
            assert(i < m_data_ids.size());
            m_patch_visible[i] = val;
            auto pid = m_data_ids[i];
            viewer.data(pid).set_visible(val);
            assert(viewer.core().is_set(viewer.data(pid).is_visible) == val);
        };

        static auto update_cell_visibility = [&](auto& viewer, size_t i) {
            const auto& cell = m_states->get_cells()[i];
            const bool curr_cell_visible = m_cell_visible[i];
            const auto& adj_list = m_states->get_cell_patch_adjacency();
            for (auto patch_index : cell) {
                const auto& adj_cells = adj_list[patch_index];
                bool val = false;
                for (auto cell_index : adj_cells) {
                    if (m_cell_visible[cell_index] != curr_cell_visible) {
                        val = true;
                    }
                }

                set_patch_visible(viewer, patch_index, val);
            }
        };

        static auto update_implicit_visibility = [&](auto& viewer, size_t i) {
            const auto& patch_implicit = m_states->get_patch_implicit();
            const auto num_patches = patch_implicit.size();
            for (size_t patch_index = 0; patch_index < num_patches; patch_index++) {
                if (patch_implicit[patch_index] == i) {
                    set_patch_visible(viewer, patch_index, m_implicit_visible[i]);
                }
            }

            const auto control_id = m_control_view_ids[i];
            viewer.data(control_id).set_visible(m_implicit_visible[i]);
            const auto sample_id = m_sample_view_ids[i];
            viewer.data(sample_id).set_visible(m_implicit_visible[i]);
        };

        static auto hide_control_pts = [&](auto& viewer) {
            const auto num_implicits = m_states->get_num_implicits();
            for (size_t i = 0; i < num_implicits; i++) {
                const int control_id = m_control_view_ids[i];
                viewer.data(control_id).set_visible(false);
                const int sample_id = m_sample_view_ids[i];
                viewer.data(sample_id).set_visible(true);
                m_implicit_visible[i] = false;
            }
        };

        viewer.plugins.push_back(&m_menu);
        m_menu.callback_draw_viewer_menu = [&]() {
            // m_menu.draw_viewer_menu();
            if (ImGui::Button("Center object", ImVec2(-1, 0))) {
                auto bbox = m_states->get_bbox();
                viewer.core().align_camera_center(bbox);
            }

            static int view_type = 1;
            ImGui::RadioButton("Cells", &view_type, 0);
            ImGui::SameLine();
            ImGui::RadioButton("Patches", &view_type, 1);
            ImGui::SameLine();
            ImGui::RadioButton("Implicits", &view_type, 2);

            auto add_cell_menu = [&]() {
                hide_control_pts(viewer);
                if (ImGui::Button("Refresh", ImVec2(-1, 0))) {
                    const auto& cells = m_states->get_cells();
                    const size_t num_cells = cells.size();

                    m_cell_visible = m_states->get_cell_labels();
                    for (size_t i = 0; i < num_cells; i++) {
                        update_cell_visibility(viewer, i);
                    }
                }

                if (ImGui::CollapsingHeader("Cells", ImGuiTreeNodeFlags_DefaultOpen)) {
                    const auto& cells = m_states->get_cells();
                    const size_t num_cells = cells.size();

                    for (size_t i = 0; i < num_cells; i++) {
                        std::string name = "Cells " + std::to_string(i);
                        if (ImGui::Checkbox(
                                name.c_str(),
                                [&]() { return m_cell_visible[i]; },
                                [&](bool val) { m_cell_visible[i] = val; })) {
                            update_cell_visibility(viewer, i);
                        }
                    }
                }
            };

            auto add_patch_menu = [&]() {
                hide_control_pts(viewer);
                if (ImGui::Button("Refresh", ImVec2(-1, 0))) {
                    reset_patch_visibility(viewer);
                }

                if (ImGui::CollapsingHeader("Patches", ImGuiTreeNodeFlags_DefaultOpen)) {
                    const auto& patches = m_states->get_patches();
                    const size_t num_patches = patches.size();

                    for (size_t i = 0; i < num_patches; i++) {
                        std::string name = "Patch " + std::to_string(i);
                        if (ImGui::Checkbox(
                                name.c_str(),
                                [&]() { return m_patch_visible[i]; },
                                [&](bool val) { m_patch_visible[i] = val; })) {
                            set_patch_visible(viewer, i, m_patch_visible[i]);
                        }
                    }
                }
            };

            auto add_implicit_menu = [&]() {
                if (ImGui::Button("Refresh", ImVec2(-1, 0))) {
                    const size_t num_implicits = m_states->get_num_implicits();
                    m_implicit_visible.assign(num_implicits, true);
                    for (size_t i = 0; i < num_implicits; i++) {
                        update_implicit_visibility(viewer, i);
                    }
                }

                if (ImGui::CollapsingHeader("Implicits", ImGuiTreeNodeFlags_DefaultOpen)) {
                    const size_t num_implicits = m_states->get_num_implicits();
                    for (size_t i = 0; i < num_implicits; i++) {
                        std::string name = "Implicit " + std::to_string(i);
                        if (ImGui::Checkbox(
                                name.c_str(),
                                [&]() { return m_implicit_visible[i]; },
                                [&](bool val) { m_implicit_visible[i] = val; })) {
                            update_implicit_visibility(viewer, i);
                        }
                    }
                }
            };

            if (view_type == 0) {
                add_cell_menu();
            } else if (view_type == 1) {
                add_patch_menu();
            } else {
                add_implicit_menu();
            }

            // Primitive menu
            if (ImGui::Button("Add sphere", ImVec2(-1, 0))) {
                auto bbox = m_states->get_bbox();
                double diag = (bbox.row(0) - bbox.row(1)).norm();
                m_states->add_sphere(0.5 * (bbox.row(0) + bbox.row(1)), diag / 5);
                initialize_data(viewer);
            }
        };
    }

    void clear_data(igl::opengl::glfw::Viewer& viewer)
    {
        // Clear viewer data.
        viewer.selected_data_index = viewer.data_list.size() - 1;
        while (viewer.erase_mesh(viewer.selected_data_index)) {
        };
        viewer.data().clear();

        m_data_ids.clear();
        m_control_view_ids.clear();
        m_sample_view_ids.clear();
    }

    void initialize_data(igl::opengl::glfw::Viewer& viewer)
    {
        clear_data(viewer);
        initialize_visibility();
        initialize_patch_data(viewer);
        initialize_implicit_data(viewer);
    }

    void initialize_visibility()
    {
        const auto num_implicits = m_states->get_num_implicits();
        m_patch_visible = m_states->get_patch_labels();
        m_cell_visible = m_states->get_cell_labels();
        m_implicit_visible.assign(num_implicits, false);
    }

    void initialize_patch_data(igl::opengl::glfw::Viewer& viewer)
    {
        const auto& vertices = m_states->get_vertices();
        const auto& faces = m_states->get_faces();

        const auto& patches = m_states->get_patches();
        const auto num_patches = patches.size();

        m_data_ids.reserve(num_patches);
        assert(m_patch_visible.size() == num_patches);

        for (size_t i = 0; i < num_patches; i++) {
            const auto& patch = patches[i];
            const auto patch_size = patch.size();
            PSIStates::FaceArray patch_faces(patch_size, 3);
            for (size_t j = 0; j < patch_size; j++) {
                patch_faces.row(j) = faces.row(patch[j]);
            }

            int id = viewer.append_mesh();
            m_data_ids.push_back(id);

            const int implicit_id = m_states->get_implicit_from_patch(i);

            // Change the color intensity for different patches within the same
            // implicit.
            Eigen::RowVector4d color = m_states->get_implicit_color(implicit_id);
            color.template segment<3>(0) +=
                (Eigen::RowVector3d::Ones() - color.template segment<3>(0)) * (double)(i) /
                (double)(num_patches * 2 + 1);

            viewer.data(id).set_mesh(vertices, patch_faces);
            viewer.data(id).set_colors(color);
            viewer.data(id).set_visible(m_patch_visible[i]);
            viewer.data(id).double_sided = true;
        }
    }

    void initialize_implicit_data(igl::opengl::glfw::Viewer& viewer)
    {
        const auto num_implicits = m_states->get_num_implicits();
        m_control_view_ids.reserve(num_implicits);

        for (size_t i = 0; i < num_implicits; i++) {
            int id = viewer.append_mesh();
            m_control_view_ids.push_back(id);

            const auto& fn = m_states->get_implicit_function(i);
            if (!fn.has_control_points()) continue;
            const auto& pts = fn.get_control_points();

            for (const auto& p : pts) {
                viewer.data(id).add_points(p.transpose(), get_control_pt_color(i));
            }
        }

        m_sample_view_ids.reserve(num_implicits);
        for (size_t i = 0; i < num_implicits; i++) {
            int id = viewer.append_mesh();
            m_sample_view_ids.push_back(id);

            const auto& fn = m_states->get_implicit_function(i);
            const auto& pts = fn.get_sample_points();

            for (const auto& p : pts) {
                viewer.data(id).add_points(p.transpose(), get_sample_pt_color(i));
            }
        }
    }

    void initialize_picking(igl::opengl::glfw::Viewer& viewer)
    {
        viewer.callback_mouse_down = [&](igl::opengl::glfw::Viewer& viewer, int, int) -> bool {
            int fid;
            Eigen::Vector3f bc;
            // Cast a ray in the view direction starting from the mouse position
            double x = viewer.current_mouse_x;
            double y = viewer.core().viewport(3) - viewer.current_mouse_y;
            m_down_x = x;
            m_down_y = y;
            m_hit = false;
            m_hit_implicit = -1;
            std::cout << "x: " << x << " y:" << y << std::endl;

            const auto num_patches = m_states->get_patches().size();
            const auto num_implicits = m_states->get_num_implicits();

            // Check of control points.
            m_active_control_point = -1;
            for (size_t i = 0; i < num_implicits; i++) {
                if (!m_implicit_visible[i]) continue;
                auto view_id = m_control_view_ids[i];
                const auto& points = viewer.data(view_id).points;
                const size_t num_pts = points.rows();
                for (size_t j = 0; j < num_pts; j++) {
                    const Eigen::RowVector3d p = points.row(j).template segment<3>(0);
                    auto screen_p = igl::project(p.template cast<float>().transpose().eval(),
                        viewer.core().view,
                        viewer.core().proj,
                        viewer.core().viewport);
                    if (std::abs(screen_p[0] - x) < 5 && std::abs(screen_p[1] - y) < 5) {
                        m_hit = true;
                        m_active_control_point = j;
                        m_hit_implicit = i;
                        return false;
                    }
                }
            }

            // Check of sample points.
            m_active_sample_point = -1;
            for (size_t i = 0; i < num_implicits; i++) {
                auto view_id = m_sample_view_ids[i];
                const auto& points = viewer.data(view_id).points;
                const size_t num_pts = points.rows();
                for (size_t j = 0; j < num_pts; j++) {
                    const Eigen::RowVector3d p = points.row(j).template segment<3>(0);
                    auto screen_p = igl::project(p.template cast<float>().transpose().eval(),
                        viewer.core().view,
                        viewer.core().proj,
                        viewer.core().viewport);
                    if (std::abs(screen_p[0] - x) < 5 && std::abs(screen_p[1] - y) < 5) {
                        m_hit = true;
                        m_active_sample_point = j;
                        m_hit_implicit = i;
                        return false;
                    }
                }
            }

            std::vector<Eigen::RowVector3d> hits;
            std::vector<int> hit_patches;
            hits.reserve(num_patches);
            hit_patches.reserve(num_patches);
            for (size_t i = 0; i < num_patches; i++) {
                if (!m_patch_visible[i]) continue;
                if (!m_implicit_visible[m_states->get_implicit_from_patch(i)]) continue;
                const auto pid = m_data_ids[i];
                const auto& V = viewer.data(pid).V;
                const auto& F = viewer.data(pid).F;
                if (igl::unproject_onto_mesh(Eigen::Vector2f(x, y),
                        viewer.core().view,
                        viewer.core().proj,
                        viewer.core().viewport,
                        V,
                        F,
                        fid,
                        bc)) {
                    const int v0 = F(fid, 0);
                    const int v1 = F(fid, 1);
                    const int v2 = F(fid, 2);
                    Eigen::RowVector3d p =
                        V.row(v0) * bc[0] + V.row(v1) * bc[1] + V.row(v2) * bc[2];

                    hits.push_back(p);
                    hit_patches.push_back(i);
                }
            }

            if (hits.size() > 0) {
                Eigen::RowVector3d best_hit;
                double best_z = std::numeric_limits<double>::max();
                for (size_t i = 0; i < hits.size(); i++) {
                    const auto& p = hits[i];
                    const auto patch_index = hit_patches[i];
                    Eigen::RowVector4d q;
                    q << p, 1;
                    double z = -(viewer.core().view.template cast<double>() * q.transpose())[2];
                    if (z < best_z) {
                        best_hit = p;
                        best_z = z;
                        m_hit_implicit = m_states->get_implicit_from_patch(patch_index);
                    }
                }

                m_hit_pt = best_hit;
                m_hit = true;
                return false;
            }
            m_active_control_point = -1;
            return false;
        };

        viewer.callback_mouse_move = [&](igl::opengl::glfw::Viewer& viewer, int, int) -> bool {
            auto update_point = [&](int view_id, int point_id) {
                Eigen::RowVector3d p =
                    viewer.data(view_id).points.row(point_id).template segment<3>(0);
                Eigen::RowVector3d n(0, 0, 1);
                n = (viewer.core().view.block(0, 0, 3, 3).template cast<double>().inverse() *
                     n.transpose())
                        .transpose();
                Eigen::RowVector4d plane;
                plane << n, -n.dot(p);

                double x = viewer.current_mouse_x;
                double y = viewer.core().viewport(3) - viewer.current_mouse_y;

                Eigen::RowVector3d q;
                igl::unproject_on_plane(Eigen::Vector2f(x, y),
                    viewer.core().proj * viewer.core().view,
                    viewer.core().viewport,
                    plane,
                    q);
                viewer.data(view_id).points.row(point_id).template segment<3>(0) = q;
                viewer.data(view_id).dirty |= igl::opengl::MeshGL::DIRTY_OVERLAY_POINTS;
            };

            if (m_hit) {
                if (m_active_control_point >= 0) {
                    const auto view_id = m_control_view_ids[m_hit_implicit];
                    update_point(view_id, m_active_control_point);
                    return true;
                } else if (m_active_sample_point >= 0) {
                    const auto view_id = m_sample_view_ids[m_hit_implicit];
                    update_point(view_id, m_active_sample_point);
                    return true;
                }
            }
            return false;
        };

        viewer.callback_mouse_up = [&](igl::opengl::glfw::Viewer& viewer, int, int) -> bool {
            double x = viewer.current_mouse_x;
            double y = viewer.core().viewport(3) - viewer.current_mouse_y;
            if (x == m_down_x && y == m_down_y && m_hit) {
                if (m_active_control_point < 0 && m_active_sample_point < 0) {
                    const auto view_id = m_sample_view_ids[m_hit_implicit];
                    viewer.data(view_id).add_points(m_hit_pt, get_sample_pt_color(m_hit_implicit));

                    const auto& fn = m_states->get_implicit_function(m_hit_implicit);
                    auto pts = fn.get_sample_points();
                    pts.push_back(m_hit_pt);

                    m_states->update_sample_points(m_hit_implicit, pts);
                    initialize_data(viewer);

                    m_hit = false;
                    m_hit_implicit = -1;
                    m_active_control_point = -1;
                    m_active_sample_point = -1;
                    return true;
                }
            } else if (m_hit) {
                assert(m_hit_implicit >= 0);
                assert(m_active_control_point >= 0 || m_active_sample_point >= 0);

                const auto& fn = m_states->get_implicit_function(m_hit_implicit);
                if (m_active_control_point >= 0) {
                    // Control point moved.
                    assert(fn.has_control_points());
                    auto view_id = m_control_view_ids[m_hit_implicit];
                    auto pts = fn.get_control_points();
                    pts[m_active_control_point] =
                        viewer.data(view_id).points.row(m_active_control_point).template segment<3>(0);
                    m_states->update_control_points(m_hit_implicit, pts);
                    initialize_data(viewer);
                } else {
                    // Sample point moved.
                    auto view_id = m_sample_view_ids[m_hit_implicit];
                    auto pts = fn.get_sample_points();
                    pts[m_active_sample_point] =
                        viewer.data(view_id).points.row(m_active_sample_point).template segment<3>(0);
                    m_states->update_sample_points(m_hit_implicit, pts);
                    reset_patch_visibility(viewer);
                }
            }
            m_hit = false;
            m_hit_implicit = -1;
            m_active_control_point = -1;
            m_active_sample_point = -1;
            return false;
        };
    }

    Eigen::RowVector3d get_control_pt_color(int implicit_id) const
    {
        Eigen::RowVector4d implicit_color = m_states->get_implicit_color(implicit_id);
        Eigen::Vector3d pt_color =
            implicit_color.template segment<3>(0) +
            (Eigen::RowVector3d::Ones() - implicit_color.template segment<3>(0)) * 0.5;
        return pt_color;
    }

    Eigen::RowVector3d get_sample_pt_color(int implicit_id) const
    {
        Eigen::Vector3d pt_color(1, 1, 0);
        return pt_color;
    }

private:
    std::vector<int> m_data_ids; // data per patch.
    std::vector<int> m_control_view_ids; // control point views.
    std::vector<int> m_sample_view_ids; // sample point views.

    std::vector<bool> m_patch_visible;
    std::vector<bool> m_cell_visible;
    std::vector<bool> m_implicit_visible;

    igl::opengl::glfw::imgui::ImGuiMenu m_menu;
    PSIStates* m_states;
    int m_active_control_point = -1;
    int m_active_sample_point = -1;
    bool m_hit;
    int m_hit_implicit = -1;
    Eigen::RowVector3d m_hit_pt;
    double m_down_x, m_down_y;
};

int main(int argc, char** argv)
{
    struct
    {
        std::string config_file;
        std::string grid_file;
    } args;

    CLI::App app{"Piecewise implicit surface demo"};
    app.add_option("-G,--grid", args.grid_file, "Grid spec file")->required();
    app.add_option("config_file", args.config_file, "Configuration file")->required();
    CLI11_PARSE(app, argc, argv);

    // grid specification and implicit functions
    auto implicit_functions = initialize_sampled_implicit_functions(args.config_file);

    auto grid_spec = GridSpec::parse_grid_spec(args.grid_file);

    PSIStates states(grid_spec, implicit_functions);

    igl::opengl::glfw::Viewer viewer;

    MeshArrangementMenu menu;
    menu.set_states(&states);
    menu.initialize(viewer);

    viewer.launch();

    return 0;
}
