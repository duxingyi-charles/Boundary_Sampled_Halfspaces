#pragma once
#include <PSI.h>

#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/project.h>
#include <igl/unproject_on_plane.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/write_triangle_mesh.h>

#include <array>
#include <limits>
#include <map>
#include <string>

#include "PSIStates.h"

struct ActiveState
{
    int active_implicit_id = -1;
    int active_point_id = -1;
    void reset()
    {
        active_implicit_id = -1;
        active_point_id = -1;
    }
};

struct PickState
{
    bool hit = false;
    int hit_implicit_id = -1;
    int hit_point_id = -1;
    Eigen::RowVector3d hit_point;
    void reset()
    {
        hit = false;
        hit_implicit_id = -1;
        hit_point_id = -1;
    }
};

class UI
{
public:
    UI()
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
        initialize_menu(viewer);
        initialize_mouse_behaviors(viewer);
        initialize_hotkeys(viewer);
    }

private:
    void reset_patch_visibility(igl::opengl::glfw::Viewer& viewer)
    {
        auto set_patch_visible = [&](auto& viewer, size_t i, bool val) {
            assert(i < m_data_ids.size());
            m_patch_visible[i] = val;
            const auto pid = m_data_ids[i];

            const int implicit_id = m_states->get_implicit_from_patch(i);
            Eigen::RowVector4d color = m_states->get_implicit_color(implicit_id);

            if (val) {
                viewer.data(pid).set_visible(true);
                viewer.data(pid).show_lines = false;
                viewer.data(pid).show_faces = true;
            } else {
                viewer.data(pid).set_visible(false);
                viewer.data(pid).show_lines = false;
                viewer.data(pid).show_faces = false;
                viewer.data(pid).line_color = color.transpose().template cast<float>();
            }
        };

        auto set_implicit_visible = [&](auto& viewer, size_t implicit_id) {
            const auto pid = m_implicit_data_ids[implicit_id];
            Eigen::RowVector4d color = m_states->get_implicit_color(implicit_id);

            // Show hit implicit patches.
            if (implicit_id == m_active_state.active_implicit_id) {
                viewer.data(pid).set_visible(true);
                viewer.data(pid).show_lines = true;
                viewer.data(pid).show_faces = true;
                viewer.data(pid).line_color = color.transpose().template cast<float>();
            } else {
                viewer.data(pid).set_visible(true);
                viewer.data(pid).show_lines = m_show_wire_frame;
                viewer.data(pid).show_faces = false;
                viewer.data(pid).line_color = color.transpose().template cast<float>();
            }
        };

        const size_t num_patches = m_states->get_num_patches();
        const size_t num_implicits = m_states->get_num_implicits();

        m_patch_visible = m_states->get_patch_labels();

        for (size_t i = 0; i < num_patches; i++) {
            set_patch_visible(viewer, i, m_patch_visible[i]);
        }
        for (size_t i = 0; i < num_implicits; i++) {
            set_implicit_visible(viewer, i);
        }

        update_mode(viewer);
    }

    void update_mode(igl::opengl::glfw::Viewer& viewer)
    {
        const auto num_implicits = m_states->get_num_implicits();
        for (size_t i = 0; i < num_implicits; i++) {
            const int control_id = m_control_view_ids[i];
            const bool is_active = m_active_state.active_implicit_id == i;
            viewer.data(control_id).set_visible(m_ui_mode == 1 && is_active);

            const int sample_id = m_sample_view_ids[i];
            viewer.data(sample_id).set_visible(m_ui_mode == 0 && is_active);
        }
        m_pick_state.reset();
    }

    void initialize_menu(igl::opengl::glfw::Viewer& viewer)
    {
        update_mode(viewer);

        static bool imgui_demo = false;
        viewer.plugins.push_back(&m_menu);
        m_menu.callback_draw_viewer_menu = [&]() {
            auto post_update_geometry = [&]() {
                initialize_data(viewer);
                m_active_state.reset();
                reset_patch_visibility(viewer);
            };

            ImGui::StyleColorsLight();
            const auto width = ImGui::GetContentRegionAvailWidth();
            // ImGui::Checkbox("imgui demo", &imgui_demo);
            if (imgui_demo) ImGui::ShowDemoWindow(&imgui_demo);
            if (ImGui::RadioButton("Sample Points", &m_ui_mode, 0)) {
                m_active_state.active_point_id = -1;
                reset_patch_visibility(viewer);
            }
            if (ImGui::RadioButton("Control Points", &m_ui_mode, 1)) {
                m_active_state.active_point_id = -1;
                reset_patch_visibility(viewer);
            }
            if (ImGui::Checkbox("Wire frame", &m_show_wire_frame)) {
                reset_patch_visibility(viewer);
            }
            if (ImGui::Checkbox("Interactive", &m_interactive)) {
                reset_patch_visibility(viewer);
            }
            static int res = m_states->get_resolution();
            ImGui::SetNextItemWidth(-1);
            if (ImGui::SliderInt("", &res, 16, 128, "grid: %d")) {
                m_states->set_resolution(res);
            }
            if (ImGui::Button("Plane", ImVec2(width / 2.1, 0.0f))) {
                const auto& bbox = m_states->get_bbox();
                m_states->add_plane(bbox.colwise().mean(), Point(0, 0, 1));
                post_update_geometry();
            }
            ImGui::SameLine();
            if (ImGui::Button("Sphere", ImVec2(width / 2.1, 0.0f))) {
                const auto& bbox = m_states->get_bbox();
                const auto l = (bbox.row(1) - bbox.row(0)).minCoeff();
                m_states->add_sphere(bbox.colwise().mean(), l / 3);
                post_update_geometry();
            }
            if (ImGui::Button("Cylinder", ImVec2(width / 2.1, 0.0f))) {
                const auto& bbox = m_states->get_bbox();
                const auto l = (bbox.row(1) - bbox.row(0)).minCoeff();
                m_states->add_cylinder(bbox.colwise().mean(), Point(0, 1, 0), l / 10);
                post_update_geometry();
            }
            ImGui::SameLine();
            if (ImGui::Button("Cone", ImVec2(width / 2.1, 0.0f))) {
                const auto& bbox = m_states->get_bbox();
                m_states->add_cone(bbox.colwise().mean(), Point(0, 0, -1), M_PI / 4);
                post_update_geometry();
            }
            if (ImGui::Button("Torus", ImVec2(width / 2.1, 0.0f))) {
                const auto& bbox = m_states->get_bbox();
                const auto l = (bbox.row(1) - bbox.row(0)).minCoeff();
                m_states->add_torus(bbox.colwise().mean(), Point(0, 0, 1), l/3, l/20);
                post_update_geometry();
            }
            ImGui::SameLine();
            if (ImGui::Button("Implicit", ImVec2(width / 2.1, 0.0f))) {
                const auto& bbox = m_states->get_bbox();
                const auto l = (bbox.row(1) - bbox.row(0)).minCoeff();
                std::vector<Point> pts(3);
                pts[0] = bbox.colwise().mean();
                pts[1] = pts[0] + Point(l/20, 0, 0);
                pts[2] = pts[0] + Point(0, l/20, 0);
                m_states->add_vipss(pts, pts);
                post_update_geometry();
            }
            ImGui::PushItemWidth(width);
            if (ImGui::SliderInt("Implicit id",
                    &m_active_state.active_implicit_id,
                    -1,
                    m_states->get_num_implicits() - 1,
                    "Implicit #%d")) {
                const int num_implicits = m_states->get_num_implicits();
                m_active_state.active_implicit_id =
                    std::min(m_active_state.active_implicit_id, num_implicits - 1);
                m_active_state.active_implicit_id = std::max(m_active_state.active_implicit_id, -1);
                reset_patch_visibility(viewer);
            }
            ImGui::PopItemWidth();
            if (ImGui::Button("Flip", ImVec2(width, 0.0f))) {
                if (m_active_state.active_implicit_id >= 0) {
                    auto& fn = m_states->get_implicit_function(m_active_state.active_implicit_id);
                    fn.flip();
                    m_states->refresh();
                    post_update_geometry();
                }
            }
            if (ImGui::Button("Update", ImVec2(width, 0.0f))) {
                m_states->refresh();
                post_update_geometry();
            }
            if (ImGui::Button("Save mesh", ImVec2(width, 0.0f))) {
                m_states->save_all("psi_output.obj");
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
        m_implicit_data_ids.clear();
        m_control_view_ids.clear();
        m_sample_view_ids.clear();
    }

    void initialize_data(igl::opengl::glfw::Viewer& viewer)
    {
        clear_data(viewer);
        initialize_visibility();
        initialize_patch_data(viewer);
        initialize_implicit_mesh_data(viewer);
        initialize_implicit_data(viewer);
    }

    void initialize_visibility() { m_patch_visible = m_states->get_patch_labels(); }

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
            viewer.data(id).show_lines = false;
        }
    }

    void initialize_implicit_mesh_data(igl::opengl::glfw::Viewer& viewer)
    {
        const auto num_implicits = m_states->get_num_implicits();
        m_implicit_data_ids.reserve(num_implicits);
        const auto& implicit_meshes = m_states->get_implicit_meshes();
        assert(num_implicits + 1 == implicit_meshes.size());

        for (size_t i = 0; i < num_implicits; i++) {
            int id = viewer.append_mesh();
            m_implicit_data_ids.push_back(id);

            const auto& mesh = implicit_meshes[i + 1]; // +1 to skip 1st mesh, which is bbox.
            Eigen::RowVector4d color = m_states->get_implicit_color(i);
            color[3] = 0.2;

            viewer.data(id).set_mesh(mesh.vertices, mesh.faces);
            viewer.data(id).set_colors(color);
            viewer.data(id).set_visible(false);
            viewer.data(id).double_sided = true;
            viewer.data(id).show_lines = true;
            viewer.data(id).show_faces = true;
        }
    }

    void initialize_implicit_data(igl::opengl::glfw::Viewer& viewer)
    {
        const auto num_implicits = m_states->get_num_implicits();
        m_control_view_ids.reserve(num_implicits);
        const auto l = m_states->get_grid_diag();

        for (size_t i = 0; i < num_implicits; i++) {
            int id = viewer.append_mesh();
            m_control_view_ids.push_back(id);

            const auto& fn = m_states->get_implicit_function(i);
            if (!fn.has_control_points()) continue;
            const auto& pts = fn.get_control_points();

            if (fn.get_type() == "sphere") {
                viewer.data(id).add_points(pts[0].transpose(), get_control_pt_color(i));
                viewer.data(id).add_points(
                    pts[1].transpose(), Eigen::Matrix<double, 1, 3>(1, 1, 0));
            } else if (fn.get_type() == "cone") {
                viewer.data(id).add_points(pts[0].transpose(), get_control_pt_color(i));
                viewer.data(id).add_points(
                    pts[1].transpose(), Eigen::Matrix<double, 1, 3>(1, 1, 0));
                viewer.data(id).add_points(
                    pts[2].transpose(), Eigen::Matrix<double, 1, 3>(0, 1, 0));
            } else if (fn.get_type() == "cylinder") {
                viewer.data(id).add_points(pts[0].transpose(), get_control_pt_color(i));
                viewer.data(id).add_points(
                    pts[1].transpose(), Eigen::Matrix<double, 1, 3>(1, 1, 0));
                viewer.data(id).add_points(
                    pts[2].transpose(), Eigen::Matrix<double, 1, 3>(0, 1, 0));
            } else if (fn.get_type() == "plane") {
                viewer.data(id).add_points(pts[0].transpose(), get_control_pt_color(i));
                viewer.data(id).add_points(
                    pts[1].transpose(), Eigen::Matrix<double, 1, 3>(1, 1, 0));
            } else if (fn.get_type() == "torus") {
                viewer.data(id).add_points(pts[0].transpose(), get_control_pt_color(i));
                viewer.data(id).add_points(
                    pts[1].transpose(), Eigen::Matrix<double, 1, 3>(1, 1, 0));
                viewer.data(id).add_points(
                    pts[2].transpose(), Eigen::Matrix<double, 1, 3>(0, 1, 0));
                viewer.data(id).add_points(
                    pts[3].transpose(), Eigen::Matrix<double, 1, 3>(0.25, 0.75, 0));
            } else {
                for (const auto& p : pts) {
                    viewer.data(id).add_points(p.transpose(), get_control_pt_color(i));
                }
            }

            // add_secondary_implicit_data(viewer, id, fn);
            viewer.data(id).show_overlay_depth = 0;
        }

        m_sample_view_ids.reserve(num_implicits);
        for (size_t i = 0; i < num_implicits; i++) {
            int id = viewer.append_mesh();
            m_sample_view_ids.push_back(id);

            const auto& fn = m_states->get_implicit_function(i);
            const auto& pts = fn.get_sample_points();

            for (const auto& p : pts) {
                viewer.data(id).add_points(p.transpose(), get_sample_pt_color(i));
                auto n = fn.gradient_at(p);
                viewer.data(id).add_edges(p.transpose(),
                    (p + n.normalized() * l / 20).transpose(),
                    Eigen::RowVector3d(0, 0, 0));
            }
            viewer.data(id).show_overlay_depth = 0;
            viewer.data(id).show_lines = true;
        }
    }

    void add_secondary_implicit_data(
        igl::opengl::glfw::Viewer& viewer, int id, const Sampled_Implicit& fn)
    {
        assert(id >= 0);
        if (const Sphere_sImplicit* sphere = dynamic_cast<const Sphere_sImplicit*>(&fn)) {
            std::cout << "Sphere!" << std::endl;
            double r = sphere->get_radius();
            Eigen::Vector3d dir(1, 0, 0);
            Eigen::MatrixXd p(1, 3);
            p = (sphere->get_center() + dir.normalized() * r).transpose();
            Eigen::Matrix<double, 1, 3> c(1, 1, 0);
            viewer.data(id).add_points(p, c);
        } else if (const Cylinder_sImplicit* cylinder =
                       dynamic_cast<const Cylinder_sImplicit*>(&fn)) {
            std::cout << "Cylinder!" << std::endl;
        } else if (const Plane_sImplicit* plane = dynamic_cast<const Plane_sImplicit*>(&fn)) {
            std::cout << "Plane!" << std::endl;
        } else if (const Cone_sImplicit* cone = dynamic_cast<const Cone_sImplicit*>(&fn)) {
            std::cout << "Cone!" << std::endl;
        } else {
            std::cout << "Other!" << std::endl;
        }
    }

    /**
     * Shoot a ray from camera through the screen point (x,y) and see which
     * patch it hits.  The result will be stored in m_pick_state.
     */
    void pick(igl::opengl::glfw::Viewer& viewer, double x, double y)
    {
        m_pick_state.reset();
        const size_t num_patches = m_states->get_num_patches();

        std::vector<Eigen::RowVector3d> hits;
        std::vector<int> hit_patches;
        hits.reserve(num_patches);
        hit_patches.reserve(num_patches);

        int fid;
        Eigen::Vector3f bc;

        // Check if any visible patch is picked.
        for (size_t i = 0; i < num_patches; i++) {
            if (!m_patch_visible[i]) continue;
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
                Eigen::RowVector3d p = V.row(v0) * bc[0] + V.row(v1) * bc[1] + V.row(v2) * bc[2];

                hits.push_back(p);
                hit_patches.push_back(i);
            }
        }

        if (hits.size() > 0) {
            Eigen::RowVector3d best_hit;
            double best_z = std::numeric_limits<double>::max();
            int implicit_id;
            for (size_t i = 0; i < hits.size(); i++) {
                const auto& p = hits[i];
                const auto patch_index = hit_patches[i];
                Eigen::RowVector4d q;
                q << p, 1;
                double z = -(viewer.core().view.template cast<double>() * q.transpose())[2];
                if (z < best_z) {
                    best_hit = p;
                    best_z = z;
                    implicit_id = m_states->get_implicit_from_patch(patch_index);
                }
            }

            m_pick_state.hit = true;
            m_pick_state.hit_implicit_id = implicit_id;
            m_pick_state.hit_point = best_hit;
        } else if (m_active_state.active_implicit_id >= 0) {
            // No visible patch is picked.
            // Check if active impliciti surface it is picked.
            const auto pid = m_implicit_data_ids[m_active_state.active_implicit_id];
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
                Eigen::RowVector3d p = V.row(v0) * bc[0] + V.row(v1) * bc[1] + V.row(v2) * bc[2];

                m_pick_state.hit = true;
                m_pick_state.hit_implicit_id = m_active_state.active_implicit_id;
                m_pick_state.hit_point = p;
            }
        }
    }

    void initialize_mouse_down_behaviors(igl::opengl::glfw::Viewer& viewer)
    {
        static auto has_hit_point = [&](double x, double y) -> int {
            if (m_active_state.active_implicit_id < 0) return -1;

            constexpr float point_radius = 10;
            const auto& ids = (m_ui_mode == 0) ? m_sample_view_ids : m_control_view_ids;
            const auto view_id = ids[m_active_state.active_implicit_id];
            const auto& points = viewer.data(view_id).points;

            const size_t num_pts = points.rows();
            for (size_t i = 0; i < num_pts; i++) {
                const Eigen::RowVector3d p = points.row(i).template segment<3>(0);
                auto screen_p = igl::project(p.template cast<float>().transpose().eval(),
                    viewer.core().view,
                    viewer.core().proj,
                    viewer.core().viewport);
                if (std::abs(screen_p[0] - x) < point_radius &&
                    std::abs(screen_p[1] - y) < point_radius) {
                    return i;
                }
            }
            return -1;
        };

        viewer.callback_mouse_down = [&](igl::opengl::glfw::Viewer& viewer, int key, int modifier) -> bool {
            m_mouse_down = true;
            m_pick_state.reset();

            double x = viewer.current_mouse_x;
            double y = viewer.core().viewport(3) - viewer.current_mouse_y;
            m_down_x = x;
            m_down_y = y;
            std::cout << "Mouse down: " << x << ", " << y;
            std::cout << " Key: " << key << ", " << modifier << std::endl;
            m_shift_down = modifier == 1;

            auto hit_point_id = has_hit_point(x, y);
            if (hit_point_id >= 0) {
                // Hit a control/sample point in the active implicit.
                m_active_state.active_point_id = hit_point_id;
            } else {
                // Pick to see if we hit an implicit surface.
                m_active_state.active_point_id = -1;
                pick(viewer, x, y);
            }
            return false;
        };
    }

    void initialize_mouse_move_behaviors(igl::opengl::glfw::Viewer& viewer)
    {
        viewer.callback_mouse_move = [&](igl::opengl::glfw::Viewer& viewer, int, int) -> bool {
            double x = viewer.current_mouse_x;
            double y = viewer.core().viewport(3) - viewer.current_mouse_y;

            auto update_active_point = [&]() {
                int view_id = -1;
                int point_id = m_active_state.active_point_id;
                if (m_ui_mode == 0) {
                    view_id = m_sample_view_ids[m_active_state.active_implicit_id];
                } else {
                    view_id = m_control_view_ids[m_active_state.active_implicit_id];
                }
                assert(view_id >= 0);

                if (m_ui_mode == 1) {
                    // Update control point.
                    Eigen::RowVector3d p =
                        viewer.data(view_id).points.row(point_id).template segment<3>(0);
                    Eigen::RowVector3d n(0, 0, 1);
                    n = (viewer.core().view.block(0, 0, 3, 3).template cast<double>().inverse() *
                         n.transpose())
                            .transpose();
                    Eigen::RowVector4d plane;
                    plane << n, -n.dot(p);

                    Eigen::RowVector3d q;
                    igl::unproject_on_plane(Eigen::Vector2f(x, y),
                        viewer.core().proj * viewer.core().view,
                        viewer.core().viewport,
                        plane,
                        q);
                    viewer.data(view_id).points.row(point_id).template segment<3>(0) = q;
                    viewer.data(view_id).dirty |= igl::opengl::MeshGL::DIRTY_OVERLAY_POINTS;
                } else {
                    // Update sample point.  Keep sample point on the implicit
                    // surface.
                    auto pid = m_implicit_data_ids[m_active_state.active_implicit_id];
                    const auto& V = viewer.data(pid).V;
                    const auto& F = viewer.data(pid).F;
                    int fid;
                    Eigen::Vector3f bc;
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

                        viewer.data(view_id).points.row(point_id).template segment<3>(0) = p;
                        const auto& fn =
                            m_states->get_implicit_function(m_active_state.active_implicit_id);
                        auto n = fn.gradient_at(p);
                        const auto l = m_states->get_grid_diag();
                        viewer.data(view_id).lines.row(point_id).template segment<3>(0) = p;
                        viewer.data(view_id).lines.row(point_id).template segment<3>(3) =
                            p + n.transpose().normalized() * l / 20;
                        viewer.data(view_id).dirty |= igl::opengl::MeshGL::DIRTY_OVERLAY_POINTS;
                        viewer.data(view_id).dirty |= igl::opengl::MeshGL::DIRTY_OVERLAY_LINES;
                    }
                }
            };

            if (m_mouse_down) {
                const bool is_implicit_active = m_active_state.active_implicit_id >= 0;
                const bool is_point_active = m_active_state.active_point_id >= 0;
                if (is_implicit_active && !is_point_active) {
                    if (m_shift_down && m_pick_state.hit) {
                        return true;
                    }
                } else if (is_implicit_active && is_point_active) {
                    update_active_point();
                    return true;
                }
            } else {
                // Hover. Do nothing for now.
            }
            return false;
        };
    }

    void initialize_mouse_up_behaviors(igl::opengl::glfw::Viewer& viewer)
    {
        viewer.callback_mouse_up = [&](igl::opengl::glfw::Viewer& viewer, int key, int modifier) -> bool {
            if (!m_mouse_down) {
                // Clicked on menu item will not trigger mouse down callback.
                return false;
            }
            double x = viewer.current_mouse_x;
            double y = viewer.core().viewport(3) - viewer.current_mouse_y;
            std::cout << "Mouse up: " << x << ", " << y;
            std::cout << " Key: " << key << ", " << modifier << std::endl;
            m_shift_down = modifier == 1;
            const bool mouse_moved = (x != m_down_x) || (y != m_down_y);
            const bool is_patch_active = m_active_state.active_implicit_id >= 0;
            const bool is_point_active = is_patch_active && m_active_state.active_point_id >= 0;
            const auto point_id = m_active_state.active_point_id;
            const auto implicit_id = m_active_state.active_implicit_id;

            auto translate_implicit = [&]() {
                auto& fn = m_states->get_implicit_function(m_active_state.active_implicit_id);
                const auto& p = m_pick_state.hit_point;
                Eigen::RowVector3d n(0, 0, 1);
                n = (viewer.core().view.block(0, 0, 3, 3).template cast<double>().inverse() *
                        n.transpose())
                    .transpose();
                Eigen::RowVector4d plane;
                plane << n, -n.dot(p);

                Eigen::RowVector3d q;
                igl::unproject_on_plane(Eigen::Vector2f(x, y),
                        viewer.core().proj * viewer.core().view,
                        viewer.core().viewport,
                        plane,
                        q);

                fn.translate(q-p);
            };

            if (is_point_active && mouse_moved) {
                // Active point dragged.
                const auto& fn = m_states->get_implicit_function(implicit_id);
                if (m_ui_mode == 0) {
                    // Sample mode.
                    const auto view_id = m_sample_view_ids[implicit_id];
                    auto pts = fn.get_sample_points();
                    pts[point_id] =
                        viewer.data(view_id).points.row(point_id).template segment<3>(0);
                    m_states->update_sample_points(implicit_id, pts);
                } else {
                    // Control mode.
                    const auto view_id = m_control_view_ids[implicit_id];
                    auto pts = fn.get_control_points();
                    pts[point_id] =
                        viewer.data(view_id).points.row(point_id).template segment<3>(0);
                    m_states->update_control_points(implicit_id, pts);
                    if (m_interactive) {
                        // Todo: find ways to only update changed implicit.
                        m_states->refresh();
                        initialize_data(viewer);
                    } else {
                        m_states->update_implicit(implicit_id);
                        initialize_data(viewer);
                    }
                }
            } else if (is_point_active) {
                // A point is clicked.  Do nothing.
            } else if (is_patch_active && !mouse_moved && m_pick_state.hit &&
                       m_pick_state.hit_implicit_id == m_active_state.active_implicit_id) {
                // Insert new point to the active patch.
                const auto& hit_point = m_pick_state.hit_point;
                const auto& fn = m_states->get_implicit_function(implicit_id);
                const auto l = m_states->get_grid_diag();

                if (m_ui_mode == 0) {
                    const auto view_id = m_sample_view_ids[implicit_id];
                    viewer.data(view_id).add_points(hit_point, get_sample_pt_color(implicit_id));

                    auto n = fn.gradient_at(hit_point);
                    viewer.data(view_id).add_edges(hit_point,
                        hit_point + n.transpose().normalized() * l / 20,
                        Eigen::RowVector3d(0, 0, 0));

                    auto pts = fn.get_sample_points();
                    pts.push_back(hit_point);
                    m_states->update_sample_points(implicit_id, pts);
                } else {
                    const auto view_id = m_control_view_ids[implicit_id];
                    viewer.data(view_id).add_points(hit_point, get_control_pt_color(implicit_id));

                    const auto& fn = m_states->get_implicit_function(implicit_id);
                    auto pts = fn.get_control_points();
                    pts.push_back(hit_point);
                    m_states->update_control_points(implicit_id, pts);
                    if (m_interactive) {
                        m_states->refresh();
                        initialize_data(viewer);
                    }
                }
            } else if (is_patch_active && m_pick_state.hit && mouse_moved && modifier == 1) {
                // Translated.
                translate_implicit();
                if (m_interactive) {
                    m_states->refresh();
                    initialize_data(viewer);
                } else {
                    m_states->update_implicit(implicit_id);
                    initialize_data(viewer);
                }
            } else if (!mouse_moved) {
                // Nothing was active, or another patch is being activated.
                if (m_pick_state.hit) {
                    // Clicked a non-active patch, just activate it.
                    m_active_state.active_implicit_id = m_pick_state.hit_implicit_id;
                    m_active_state.active_point_id = -1;
                } else {
                    // Clicked on empty space.  Deactivate all.
                    m_active_state.reset();
                }
            }

            reset_patch_visibility(viewer);
            m_mouse_down = false;
            return true;
        };
    }

    void initialize_mouse_behaviors(igl::opengl::glfw::Viewer& viewer)
    {
        initialize_mouse_down_behaviors(viewer);
        initialize_mouse_move_behaviors(viewer);
        initialize_mouse_up_behaviors(viewer);
    }

    void initialize_hotkeys(igl::opengl::glfw::Viewer& viewer)
    {
        viewer.callback_key_down =
            [&](igl::opengl::glfw::Viewer& viewer, unsigned int key, int modifier) -> bool {
            std::cout << "Key down -- Key: " << key << " Modifier: " << modifier << std::endl;
            if (key == 32) {
                m_show_wire_frame = !m_show_wire_frame;
                reset_patch_visibility(viewer);
                return true;
            }
            return false;
        };
        viewer.callback_key_up =
            [&](igl::opengl::glfw::Viewer& viewer, unsigned int key, int modifier) -> bool {
            std::cout << "Key up   -- Key: " << key << " Modifier: " << modifier << std::endl;
            if (key == 32) {
                m_show_wire_frame = !m_show_wire_frame;
                reset_patch_visibility(viewer);
                return true;
            }
            if (key == 259 || key == 88) {
                remove_active_point(viewer);
                reset_patch_visibility(viewer);
                return true;
            }
            return false;
        };
    }

    void remove_active_point(igl::opengl::glfw::Viewer& viewer)
    {
        if (m_active_state.active_implicit_id < 0) return;
        if (m_active_state.active_point_id < 0) return;

        const auto implicit_id = m_active_state.active_implicit_id;
        const auto point_id = m_active_state.active_point_id;
        const auto& fn = m_states->get_implicit_function(implicit_id);
        const auto l = m_states->get_grid_diag();

        if (m_ui_mode == 0) {
            // Working with sample points.
            auto view_id = m_sample_view_ids[implicit_id];
            auto pts = fn.get_sample_points();
            pts.erase(pts.begin() + point_id);
            m_states->update_sample_points(implicit_id, pts);

            viewer.data(view_id).clear_points();
            viewer.data(view_id).clear_edges();
            auto c = get_sample_pt_color(implicit_id);
            for (auto& p : pts) {
                viewer.data(view_id).add_points(p.transpose(), c);

                auto n = fn.gradient_at(p);
                viewer.data(view_id).add_edges(p.transpose(),
                    (p + n.normalized() * l / 20).transpose(),
                    Eigen::RowVector3d(0, 0, 0));
            }
            m_active_state.active_point_id = -1;
        } else {
            // Working with control points.
            auto view_id = m_control_view_ids[implicit_id];
            auto pts = fn.get_control_points();
            pts.erase(pts.begin() + point_id);
            m_states->update_control_points(implicit_id, pts);
            m_states->refresh();

            viewer.data(view_id).clear_points();
            auto c = get_control_pt_color(implicit_id);
            for (auto& p : pts) {
                viewer.data(view_id).add_points(p.transpose(), c);
            }
            m_active_state.active_point_id = -1;
        }
    }

    Eigen::RowVector3d get_control_pt_color(int implicit_id) const
    {
        Eigen::RowVector4d implicit_color = m_states->get_implicit_color(implicit_id);
        Eigen::Vector3d pt_color =
            implicit_color.template segment<3>(0) +
            (Eigen::RowVector3d::Ones() - implicit_color.template segment<3>(0)) * 0.2;
        return pt_color;
    }

    Eigen::RowVector3d get_sample_pt_color(int implicit_id) const
    {
        Eigen::RowVector4d implicit_color = m_states->get_implicit_color(implicit_id);
        Eigen::Vector3d pt_color =
            implicit_color.template segment<3>(0) +
            (Eigen::RowVector3d::Ones() - implicit_color.template segment<3>(0)) * 0.1;
        return pt_color;
    }

private:
    std::vector<int> m_data_ids; // data per patch.
    std::vector<int> m_implicit_data_ids; // data per implicit surface.
    std::vector<int> m_control_view_ids; // control point views.
    std::vector<int> m_sample_view_ids; // sample point views.

    std::vector<bool> m_patch_visible;

    igl::opengl::glfw::imgui::ImGuiMenu m_menu;
    PSIStates* m_states;
    double m_down_x, m_down_y;
    int m_ui_mode = 0;
    bool m_show_wire_frame = false;
    bool m_interactive = true;

    bool m_mouse_down = false;
    bool m_shift_down = false;
    PickState m_pick_state;
    ActiveState m_active_state;
};

