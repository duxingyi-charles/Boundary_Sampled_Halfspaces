#include <Mesh_PSI.h>
#include <PSI.h>
#include <Sphere_sImplicit.h>

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
        : m_grid(grid), m_implicits(implicits)
    {
        m_psi = std::make_unique<Mesh_PSI>();
        m_psi->run(grid, implicits);
        initialize_states();
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
    const auto& get_bbox() const { return m_bbox; }

    void add_sphere(const Point& center, double radius) {
        m_implicits.push_back(std::make_unique<Sphere_sImplicit>(center, radius));
        m_psi->run(m_grid, m_implicits);
        initialize_states();
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

private:
    std::unique_ptr<PSI> m_psi;
    VertexArray m_vertices;
    FaceArray m_faces;
    Eigen::Matrix<double, 2, 3> m_bbox;
    const GridSpec& m_grid;
    std::vector<std::unique_ptr<Sampled_Implicit>>& m_implicits;
};

