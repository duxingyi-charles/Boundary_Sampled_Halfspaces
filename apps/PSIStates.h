#include <Mesh_PSI.h>
#include <PSI.h>

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
    PSIStates(
        const GridSpec& grid, std::vector<std::unique_ptr<Sampled_Implicit>>& implicit_fns)
    {
        m_psi = std::make_unique<Mesh_PSI>();
        m_psi->run(grid, implicit_fns);
        initialize_states();
    }

    const VertexArray& get_vertices() const { return m_vertices; }
    const FaceArray& get_faces() const { return m_faces; }
    const MapArray& get_face_mapping() const { return m_face_mapping; }
    const auto& get_patches() const { return m_psi->get_patches(); }
    const auto& get_cells() const { return m_psi->get_cells(); }
    const auto& get_bbox() const { return m_bbox; }

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
        m_face_mapping.resize(num_triangles, 1);
        for (size_t i = 0; i < num_faces; i++) {
            m_faces(i, 0) = F[i][0];
            m_faces(i, 1) = F[i][1];
            m_faces(i, 2) = F[i][2];
            m_face_mapping[i] = i;
        }
    }

private:
    std::unique_ptr<PSI> m_psi;
    VertexArray m_vertices;
    FaceArray m_faces;
    MapArray m_face_mapping;
    Eigen::Matrix<double, 2, 3> m_bbox;
};
