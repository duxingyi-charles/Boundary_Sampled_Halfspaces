/* This file is part of PyMesh. Copyright (c) 2016 by Qingnan Zhou */
#pragma once
#include <memory>

#include "EigenTypedef.h"

namespace PyMesh {

class Arrangement
{
public:
    typedef std::shared_ptr<Arrangement> Ptr;
    static Ptr create_mesh_arrangement(const MatrixFr& vertices, const MatrixIr& faces);
    static Ptr create_fast_arrangement(const MatrixFr& vertices, const MatrixIr& faces);

public:
    Arrangement(const MatrixFr& vertices, const MatrixIr& faces)
        : m_vertices(vertices)
        , m_faces(faces)
    {}
    virtual ~Arrangement() = default;

    virtual void run() =0;

    MatrixFr get_vertices() const { return m_vertices; }
    MatrixIr get_faces() const { return m_faces; }
    VectorI get_source_faces() const { return m_source_faces; }

    size_t get_num_cells() const
    {
        if (m_cells.rows() > 0)
            return m_cells.maxCoeff();
        else
            return 0;
    }
    MatrixIr get_cell_faces(const size_t cell_id) const
    {
        const size_t num_faces = m_faces.rows();
        MatrixIr faces(num_faces, 3);
        size_t face_count = 0;
        for (size_t i = 0; i < num_faces; i++) {
            const size_t patch_id = m_patches[i];
            if (m_cells(patch_id, 0) == cell_id) {
                faces.row(face_count) = m_faces.row(i).reverse();
                face_count++;
            } else if (m_cells(patch_id, 1) == cell_id) {
                faces.row(face_count) = m_faces.row(i);
                face_count++;
            }
        }
        faces.conservativeResize(face_count, 3);
        return faces;
    }
    MatrixIr get_cells() const { return m_cells; }

    size_t get_num_patches() const
    {
        if (m_patches.size() > 0)
            return m_patches.maxCoeff();
        else
            return 0;
    }
    VectorI get_patches() const { return m_patches; }

    MatrixIr get_winding_number() const { return m_winding_number; }

protected:
    MatrixFr m_vertices;
    MatrixIr m_faces;
    VectorI m_source_faces;
    MatrixIr m_cells;
    VectorI m_patches;
    MatrixIr m_winding_number;
};

} // namespace PyMesh

