#pragma once

#include "Arrangement.h"

namespace PyMesh {

class MeshArrangement final : public Arrangement
{
public:
    using Base = Arrangement;

public:
    MeshArrangement(const MatrixFr& vertices, const MatrixIr& faces)
        : Base(vertices, faces)
    {}
    ~MeshArrangement() = default;
    void run() override;

private:
    using Base::m_cells;
    using Base::m_faces;
    using Base::m_patches;
    using Base::m_source_faces;
    using Base::m_vertices;
    using Base::m_winding_number;
};

} // namespace PyMesh
