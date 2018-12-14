#pragma once

#include <Eigen/Core>
#include <stdlib.h>

using namespace Eigen;

/**
 * Defines a light by a position, emitted light and radius. Making it a spherical light casting in all directions.
 **/
class Light
{
    public:
        Vector3f origin;
        Vector3f emittedLight;
        float radius;

        /**
         * Gets a random point on the surface of the sphere.
         **/
        Vector3f randomPoint()
        {
            Vector3f point;
            point[0] = rand() / RAND_MAX;
            point[1] = rand() / RAND_MAX;
            point[2] = rand() / RAND_MAX;
            point.normalize();
            return point * radius + origin;
        }

        /**
         * Area of the sphere.
         **/
        float getArea()
        {
            return (4/3) * EIGEN_PI * radius * radius * radius;
        }
};