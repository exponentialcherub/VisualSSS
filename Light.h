#pragma once

#include <Eigen/Core>

using namespace Eigen;

class Light
{
    public:
        Vector3f dir;
        float ambientI;
        float localI;
};