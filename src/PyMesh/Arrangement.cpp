#include "Arrangement.h"
#include "MeshArrangement.h"

using namespace PyMesh;

Arrangement::Ptr Arrangement::create_raw(const MatrixFr& vertices, const MatrixIr& faces)
{
    return std::make_shared<MeshArrangement>(vertices, faces);
}

