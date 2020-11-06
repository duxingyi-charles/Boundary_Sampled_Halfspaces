#include <PSI.h>
#include <config.h>
#include <igl/copyleft/marching_cubes.h>
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
        viewer.plugins.push_back(&m_menu);
        m_menu.callback_draw_viewer_menu = [&]() {
            // m_menu.draw_viewer_menu();
            if (ImGui::Button("Center object", ImVec2(-1, 0))) {
                auto bbox = m_states->get_bbox();
                viewer.core().align_camera_center(bbox);
            }
            if (ImGui::CollapsingHeader("Cells", ImGuiTreeNodeFlags_DefaultOpen)) {
                const auto& cells = m_states->get_cells();
                const size_t num_cells = cells.size();
                for (size_t i = 0; i < num_cells; i++) {
                    std::string name = "Cell " + std::to_string(i);
                    if (ImGui::Checkbox(
                            name.c_str(),
                            [&]() { return m_visible[i]; },
                            [&](bool val) { m_visible[i] = val; })) {
                        viewer.data(i).set_visible(m_visible[i]);
                    }
                }
            }
            if (ImGui::Button("Add sphere", ImVec2(-1, 0))) {
                auto bbox = m_states->get_bbox();
                double diag = (bbox.row(0) - bbox.row(1)).norm();
                m_states->add_sphere(0.5 * (bbox.row(0) + bbox.row(1)), diag / 5);
                initialize_data(viewer);
            }
        };

        initialize_data(viewer);
        initialize_picking(viewer);
    }

private:
    void initialize_data(igl::opengl::glfw::Viewer& viewer)
    {
        // Clear viewer data.
        viewer.selected_data_index = viewer.data_list.size()-1;
        while(viewer.erase_mesh(viewer.selected_data_index)){};
        viewer.data().clear();
        m_colors.clear();

        const auto& vertices = m_states->get_vertices();
        const auto& faces = m_states->get_faces();

        const auto& cells = m_states->get_cells();
        m_visible = m_states->get_cell_labels();
        const auto& patches = m_states->get_patches();
        const auto num_cells = cells.size();
        assert(m_visible.size() == num_cells);
        for (size_t i = 0; i < num_cells; i++) {
            const auto cell_size = cells[i].size();

            size_t num_faces = 0;
            for (size_t j = 0; j < cell_size; j++) {
                const auto& patch = patches[cells[i][j]];
                const auto patch_size = patch.size();
                num_faces += patch_size;
            }

            PSIStates::FaceArray cell_faces(num_faces, 3);
            size_t fid = 0;
            for (size_t j = 0; j < cell_size; j++) {
                const auto& patch = patches[cells[i][j]];
                const auto patch_size = patch.size();
                for (size_t k = 0; k < patch_size; k++) {
                    cell_faces.row(fid) = faces.row(patch[k]);
                    fid++;
                }
            }

            std::string filename = "cell_" + std::to_string(i) + ".obj";
            igl::write_triangle_mesh(filename, vertices, cell_faces);

            viewer.append_mesh();
            viewer.data(i).V.resize(0, 3);
            viewer.data(i).F.resize(0, 3);
            viewer.data(i).set_mesh(vertices, cell_faces);
            int id = viewer.data(i).id;
            std::cout << i << " " << id << std::endl;
            m_colors.emplace(id, 0.5 * Eigen::RowVector3d::Random().array() + 0.5);
            viewer.data(i).set_colors(m_colors[id]);
            viewer.data(i).set_visible(m_visible[i]);
            viewer.data(i).double_sided = true;
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
            m_hit_cell = -1;
            std::cout << "x: " << x << " y:" << y << std::endl;
            const auto num_cells = m_states->get_cells().size();

            // Check of points.
            m_active_point = -1;
            for (size_t i = 0; i < num_cells; i++) {
                const auto& points = viewer.data(i).points;
                const size_t num_pts = points.rows();
                for (size_t j = 0; j < num_pts; j++) {
                    const Eigen::RowVector3d p = points.row(j).template segment<3>(0);
                    auto screen_p = igl::project(p.template cast<float>().transpose().eval(),
                        viewer.core().view,
                        viewer.core().proj,
                        viewer.core().viewport);
                    if (std::abs(screen_p[0] - x) < 5 && std::abs(screen_p[1] - y) < 5) {
                        m_hit = true;
                        m_active_point = j;
                        m_hit_cell = i;
                        return false;
                    }
                }
            }

            std::vector<Eigen::RowVector3d> hits;
            std::vector<int> hit_cells;
            hits.reserve(num_cells);
            hit_cells.reserve(num_cells);
            for (size_t i = 0; i < num_cells; i++) {
                if (!m_visible[i]) continue;
                const auto& V = viewer.data(i).V;
                const auto& F = viewer.data(i).F;
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
                    hit_cells.push_back(i);
                }
            }

            if (hits.size() > 0) {
                Eigen::RowVector3d best_hit;
                double best_z = std::numeric_limits<double>::max();
                for (size_t i = 0; i < hits.size(); i++) {
                    const auto& p = hits[i];
                    const auto cell_id = hit_cells[i];
                    Eigen::RowVector4d q;
                    q << p, 1;
                    double z = -(viewer.core().view.template cast<double>() * q.transpose())[2];
                    if (z < best_z) {
                        best_hit = p;
                        best_z = z;
                        m_hit_cell = cell_id;
                    }
                }

                m_hit_pt = best_hit;
                m_hit = true;
                return false;
            }
            m_active_point = -1;
            return false;
        };

        viewer.callback_mouse_move = [&](igl::opengl::glfw::Viewer& viewer, int, int) -> bool {
            if (m_hit && m_active_point >= 0) {
                assert(m_hit_cell >= 0);
                Eigen::RowVector3d p =
                    viewer.data(m_hit_cell).points.row(m_active_point).template segment<3>(0);
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
                m_points[m_active_point] = q;
                viewer.data(m_hit_cell).points.row(m_active_point).template segment<3>(0) = q;
                viewer.data(m_hit_cell).dirty |= igl::opengl::MeshGL::DIRTY_OVERLAY_POINTS;
                return true;
            }
            return false;
        };

        viewer.callback_mouse_up = [&](igl::opengl::glfw::Viewer& viewer, int, int) -> bool {
            double x = viewer.current_mouse_x;
            double y = viewer.core().viewport(3) - viewer.current_mouse_y;
            if (x == m_down_x && y == m_down_y && m_hit) {
                if (m_active_point < 0) {
                    m_points.push_back(m_hit_pt);
                    Eigen::Matrix<double, 1, 3> color(1, 1, 0);
                    viewer.data(m_hit_cell).add_points(m_points.back(), color);
                    m_hit = false;
                    m_hit_cell = -1;
                    m_active_point = -1;
                    return true;
                }
            }
            m_hit = false;
            m_hit_cell = -1;
            m_active_point = -1;
            return false;
        };
    }

private:
    std::vector<bool> m_visible;
    igl::opengl::glfw::imgui::ImGuiMenu m_menu;
    PSIStates* m_states;
    std::map<int, Eigen::RowVector3d> m_colors;
    std::vector<Eigen::RowVector3d> m_points;
    int m_active_point = -1;
    bool m_hit;
    int m_hit_cell = -1;
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
