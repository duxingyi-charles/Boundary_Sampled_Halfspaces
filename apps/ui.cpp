#include <igl/copyleft/marching_cubes.h>
#include <igl/write_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>

#include <CLI/CLI.hpp>
#include <Eigen/Core>
#include <array>
#include <string>
#include <map>

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
        }

    private:
        std::vector<bool> m_visible;
        igl::opengl::glfw::imgui::ImGuiMenu m_menu;
        MeshArrangementStates* m_states;
        std::map<int, Eigen::RowVector3d> m_colors;
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
