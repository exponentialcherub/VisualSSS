#pragma once

#include "Line.h"
#include "Light.h"
#include <Eigen/Dense>

using namespace Eigen;

class Object
{
    public:
        Vector3f colour;
        Vector3f diffuse;
        Vector3f specular;
        
        Object(Vector3f c)
        {
            colour = c;
            diffuse = colour * 0.5;
            specular = {0.5, 0.5, 0.5};
        }
        virtual bool intersects(Line ray, Light light, float &t, float &ambAngle, float &specAngle, Vector3f &intersectionPoint) = 0;
        virtual bool intersects(Line ray, float &t) = 0;
};