#include <igl/copyleft/marching_cubes.h>
#include <igl/write_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/unproject_on_plane.h>
#include <igl/project.h>

#include <CLI/CLI.hpp>
#include <Eigen/Core>
#include <array>
#include <string>
#include <map>
#include <limits>

#include "PyMesh/CellPartition.h"
#include "ScopedTimer.h"
#include "config.h"
#include "mesh_arrangement_utils.h"

class MeshArrangementStates {
    public:
        MeshArrangementStates(const std::string& grid_file, const std::string& config_file) {
            auto implicit_functions =
                initialize_sampled_implicit_functions(config_file);
            auto grid_spec = parse_grid_spec(grid_file);

            m_input_meshes.reserve(implicit_functions.size() + 1);

            m_input_meshes.push_back(generate_cube(grid_spec));
            for (const auto& fn : implicit_functions) {
                m_input_meshes.push_back(marching_cubes(*fn, grid_spec));
            }

            m_cells = compute_arrangement(m_input_meshes);
        }

        size_t get_num_cells() const { return m_cells.size(); }
        std::vector<IGL_Mesh>& get_cells() { return m_cells; }
        const std::vector<IGL_Mesh>& get_cells() const { return m_cells; }

        Eigen::Matrix<double, 2, 3> get_bbox() const {
            const auto num_cells = get_num_cells();
            if (num_cells == 0) {
                return Eigen::Matrix<double, 2, 3>::Zero();
            }
            Eigen::Matrix<double, 2, 3> bbox;
            bbox.row(0) = m_cells[0].vertices.colwise().minCoeff();
            bbox.row(1) = m_cells[0].vertices.colwise().maxCoeff();
            for (size_t i=1; i<num_cells; i++) {
                bbox.row(0) = bbox.row(0).cwiseMin(m_cells[i].vertices.colwise().minCoeff());
                bbox.row(1) = bbox.row(1).cwiseMax(m_cells[i].vertices.colwise().maxCoeff());
            }
            return bbox;
        }

    private:
        std::vector<IGL_Mesh> m_input_meshes;
        std::vector<IGL_Mesh> m_cells;
};

class MeshArrangementMenu {
    public:
        MeshArrangementMenu() : m_states(nullptr) {}

        void set_states(MeshArrangementStates* states) {
            assert(states != nullptr);
            m_states = states;
        }

        void initialize(igl::opengl::glfw::Viewer& viewer) {
            const auto& cells = m_states->get_cells();
            const size_t num_cells = m_states->get_num_cells();
            m_visible.assign(num_cells, true);

            for (size_t i=1; i<num_cells; i++) {
                viewer.append_mesh();
                viewer.data(i).set_mesh(cells[i].vertices, cells[i].faces);
                int id = viewer.data(i).id;
                m_colors.emplace(id, 0.5*Eigen::RowVector3d::Random().array() + 0.5);
                viewer.data(i).set_colors(m_colors[id]);
                viewer.data(i).set_visible(true);
            }

            viewer.plugins.push_back(&m_menu);
            m_menu.callback_draw_viewer_menu = [&]() {
                //m_menu.draw_viewer_menu();
                if (ImGui::Button("Center object", ImVec2(-1, 0))) {
                    auto bbox = m_states->get_bbox();
                    viewer.core().align_camera_center(bbox);
                }
                if (ImGui::CollapsingHeader("Cells", ImGuiTreeNodeFlags_DefaultOpen)) {
                    for (size_t i=1; i<m_states->get_num_cells(); i++) {
                        std::string name = "Cell " + std::to_string(i);
                        if (ImGui::Checkbox(name.c_str(),
                                    [&](){ return m_visible[i];},
                                    [&](bool val){ m_visible[i] = val;}
                                    )) {
                            viewer.data(i).set_visible(m_visible[i]);
                        }
                    }
                }
            };

            initialize_picking(viewer);
        }

    private:
        void initialize_picking(igl::opengl::glfw::Viewer& viewer) {
            viewer.callback_mouse_down =
                [&](igl::opengl::glfw::Viewer& viewer, int, int)->bool {
                    int fid;
                    Eigen::Vector3f bc;
                    // Cast a ray in the view direction starting from the mouse position
                    double x = viewer.current_mouse_x;
                    double y = viewer.core().viewport(3) - viewer.current_mouse_y;
                    m_down_x = x;
                    m_down_y = y;
                    m_hit = false;
                    std::cout << "x: " << x << " y:" << y << std::endl;

                    // Check of points.
                    m_active_point = -1;
                    for (size_t i=0; i<m_points.size(); i++) {
                        const auto& p = m_points[i];
                        auto screen_p = igl::project(p.template cast<float>().transpose().eval(), viewer.core().view, viewer.core().proj, viewer.core().viewport);
                        if (std::abs(screen_p[0] - x) < 5 && std::abs(screen_p[1] - y) < 5) {
                            m_hit = true;
                            m_active_point = i;
                            return false;
                        }
                    }

                    const auto& cells = m_states->get_cells();
                    std::vector<Eigen::RowVector3d> hits;
                    hits.reserve(cells.size());
                    for (size_t i=1; i<cells.size(); i++) {
                        if (!m_visible[i]) continue;
                        const auto& cell = cells[i];
                        const auto& V = cell.vertices;
                        const auto& F = cell.faces;
                        if(igl::unproject_onto_mesh(Eigen::Vector2f(x,y),
                                    viewer.core().view, viewer.core().proj, viewer.core().viewport, V, F, fid, bc))
                        {
                            const int v0 = F(fid, 0);
                            const int v1 = F(fid, 1);
                            const int v2 = F(fid, 2);
                            Eigen::RowVector3d p =
                                    V.row(v0) * bc[0] + V.row(v1) * bc[1] + V.row(v2) * bc[2];

                            hits.push_back(p);
                        }
                    }


                    if (hits.size() > 0) {
                        Eigen::RowVector3d best_hit;
                        double best_z = std::numeric_limits<double>::max();
                        for (const auto& p : hits) {
                            Eigen::RowVector4d q;
                            q << p, 1;
                            double z = -(viewer.core().view.template cast<double>() * q.transpose())[2];
                            if (z < best_z) {
                                best_hit = p;
                                best_z = z;
                            }
                        }

                        m_hit_pt = best_hit;
                        m_hit = true;
                        return false;
                    }
                    m_active_point = -1;
                    return false;
                };

            viewer.callback_mouse_move =
                [&](igl::opengl::glfw::Viewer& viewer, int, int)->bool {
                    if (m_hit && m_active_point >= 0) {
                        auto& p = m_points[m_active_point];
                        Eigen::RowVector3d n(0, 0, 1);
                        n = (viewer.core().view.block(0, 0, 3, 3).template cast<double>().inverse() * n.transpose()).transpose();
                        Eigen::RowVector4d plane;
                        plane << n, -n.dot(p);

                        double x = viewer.current_mouse_x;
                        double y = viewer.core().viewport(3) - viewer.current_mouse_y;

                        Eigen::RowVector3d q;
                        igl::unproject_on_plane(Eigen::Vector2f(x,y),
                                    viewer.core().proj * viewer.core().view,
                                    viewer.core().viewport, plane, q);
                        m_points[m_active_point] = q;
                        viewer.data(0).points.row(m_active_point).template segment<3>(0) = q;
                        viewer.data(0).dirty |= igl::opengl::MeshGL::DIRTY_OVERLAY_POINTS;
                        return true;
                    }
                    return false;
                };

            viewer.callback_mouse_up =
                [&](igl::opengl::glfw::Viewer& viewer, int, int)->bool {
                    double x = viewer.current_mouse_x;
                    double y = viewer.core().viewport(3) - viewer.current_mouse_y;
                    if (x == m_down_x && y == m_down_y && m_hit) {
                        if (m_active_point < 0) {
                            m_points.push_back(m_hit_pt);
                            Eigen::Matrix<double, 1, 3> color(1, 1, 0);
                            viewer.data(0).add_points(m_points.back(), color);
                            m_hit = false;
                            m_active_point = -1;
                            return true;
                        }
                    }
                    m_hit = false;
                    m_active_point = -1;
                    return false;
                };
        }

    private:
        std::vector<bool> m_visible;
        igl::opengl::glfw::imgui::ImGuiMenu m_menu;
        MeshArrangementStates* m_states;
        std::map<int, Eigen::RowVector3d> m_colors;
        std::vector<Eigen::RowVector3d> m_points;
        int m_active_point = -1;
        bool m_hit;
        Eigen::RowVector3d m_hit_pt;
        double m_down_x, m_down_y;
};

int main(int argc, char** argv) {
    struct {
        std::string config_file;
        std::string grid_file;
    } args;

    CLI::App app{"Piecewise implicit surface demo"};
    app.add_option("-G,--grid", args.grid_file, "Grid spec file")->required();
    app.add_option("config_file", args.config_file, "Configuration file")
        ->required();
    CLI11_PARSE(app, argc, argv);

    MeshArrangementStates states(args.grid_file, args.config_file);

    igl::opengl::glfw::Viewer viewer;

    MeshArrangementMenu menu;
    menu.set_states(&states);
    menu.initialize(viewer);

    viewer.launch();

    return 0;
}
