#include "Arrangement.h"
#include "MeshArrangement.h"
#include "FastArrangement.h"

using namespace PyMesh;

Arrangement::Ptr Arrangement::create_mesh_arrangement(const MatrixFr& vertices, const MatrixIr& faces)
{
    return std::make_shared<MeshArrangement>(vertices, faces);
}

Arrangement::Ptr Arrangement::create_fast_arrangement(const MatrixFr& vertices, const MatrixIr& faces)
{
    return std::make_shared<FastArrangement>(vertices, faces);
}
