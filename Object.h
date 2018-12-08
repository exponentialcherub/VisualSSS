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
        Vector3f sigmaS;
        Vector3f sigmaA;
        
        Object(Vector3f c)
        {
            colour = c;
            diffuse = colour * 1;
            specular = {0.5, 0.5, 0.5};
        }
        virtual bool intersects(Line ray, float &t, Vector3f &normal, Vector3f &intersectionPoint) = 0;
        virtual bool intersects(Line ray, float &t) = 0;
        virtual bool isTranslucent() = 0;

        Vector3f getSigmaT()
        {
            return sigmaS + sigmaA;
        }
};