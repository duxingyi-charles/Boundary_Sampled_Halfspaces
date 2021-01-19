#pragma once

#include "Arrangement.h"

namespace PyMesh {

class FastArrangement final : public Arrangement
{
public:
    using Base = Arrangement;

public:
    FastArrangement(const MatrixFr& vertices, const MatrixIr& faces)
        : Base(vertices, faces)
    {}
    ~FastArrangement() = default;

    [[clang::optnone]]
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
