#pragma once

#include <Eigen/Core>
#include <stdlib.h>

using namespace Eigen;

class Light
{
    public:
        Vector3f origin;
        Vector3f emittedLight;
        float radius;

        Vector3f randomPoint()
        {
            Vector3f point;
            point[0] = rand() / RAND_MAX;
            point[1] = rand() / RAND_MAX;
            point[2] = rand() / RAND_MAX;
            point.normalize();
            return point * radius + origin;
        }

        float getArea()
        {
            return (4/3) * EIGEN_PI * radius * radius * radius;
        }
};