#pragma once

#include <Eigen/Core>

using namespace Eigen;

class Light
{
    public:
        Vector3f origin;
        Vector3f dir;
        Vector3f emittedLight;
        float area;
};