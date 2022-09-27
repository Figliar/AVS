/**
 * @file    tree_mesh_builder.cpp
 *
 * @author  FULL NAME <xreset00@stud.fit.vutbr.cz>
 *
 * @brief   Parallel Marching Cubes implementation using OpenMP tasks + octree early elimination
 *
 * @date    DATE
 **/

#include <iostream>
#include <math.h>
#include <limits>

#include "tree_mesh_builder.h"

TreeMeshBuilder::TreeMeshBuilder(unsigned gridEdgeSize)
    : BaseMeshBuilder(gridEdgeSize, "Octree")
{

}

unsigned int TreeMeshBuilder::divide(float size, const ParametricScalarField &field, Vec3_t<float> cubeOffset) {
    float Length = ((float) size * mGridResolution);
    float half_Length = Length / 2.0f;
    const Vec3_t<float> middle(
            cubeOffset.x * mGridResolution + half_Length,
            cubeOffset.y * mGridResolution + half_Length,
            cubeOffset.z * mGridResolution + half_Length
            );

    float is_empty = mIsoLevel + sqrtf(3.F) / 2.F * Length;
    if (is_empty < evaluateFieldAt(middle, field)){
        return 0;
    }
    else if (size <= 1) {
        return buildCube(cubeOffset, field);
    }

    unsigned total_triangleCount = 0;
    const float new_size = ((float) size / 2.0f);

    for (const Vec3_t<float> vertex : sc_vertexNormPos) {
        #pragma omp task default(none) firstprivate(vertex) shared(cubeOffset, new_size, field, total_triangleCount)
        {
            const Vec3_t<float> newCubeOffset(
                    vertex.x * new_size + cubeOffset.x,
                    vertex.y * new_size + cubeOffset.y,
                    vertex.z * new_size + cubeOffset.z
            );
            const unsigned count = divide(new_size, field, newCubeOffset);
            #pragma omp atomic update
            total_triangleCount += count;
        }
    }
    #pragma omp taskwait
    return total_triangleCount;
}

unsigned TreeMeshBuilder::marchCubes(const ParametricScalarField &field)
{
    unsigned return_value;
    // Suggested approach to tackle this problem is to add new method to
    // this class. This method will call itself to process the children.
    // It is also strongly suggested to first implement Octree as sequential
    // code and only when that works add OpenMP tasks to achieve parallelism.
    // 1. Compute total number of cubes in the grid.
    #pragma omp parallel default(none) shared(return_value, field)
    #pragma omp single nowait
    return_value = divide((float) mGridSize, field, Vec3_t<float>());

    return return_value;
}

float TreeMeshBuilder::evaluateFieldAt(const Vec3_t<float> &pos, const ParametricScalarField &field)
{
    const Vec3_t<float> *pPoints = field.getPoints().data();
    const auto count = unsigned(field.getPoints().size());

    float value = std::numeric_limits<float>::max();

    // 2. Find minimum square distance from points "pos" to any point in the
    //    field.
    for(unsigned i = 0; i < count; ++i)
    {
        float distanceSquared  = (pos.x - pPoints[i].x) * (pos.x - pPoints[i].x);
        distanceSquared       += (pos.y - pPoints[i].y) * (pos.y - pPoints[i].y);
        distanceSquared       += (pos.z - pPoints[i].z) * (pos.z - pPoints[i].z);

        // Comparing squares instead of real distance to avoid unnecessary
        // "sqrt"s in the loop.
        value = std::min(value, distanceSquared);
    }

    // 3. Finally take square root of the minimal square distance to get the real distance
    return sqrt(value);
}

void TreeMeshBuilder::emitTriangle(const BaseMeshBuilder::Triangle_t &triangle)
{
    #pragma omp critical(tree_emitTriangle)
    mTriangles.push_back(triangle);
}
