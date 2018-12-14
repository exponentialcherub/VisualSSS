#pragma once

#include "Line.h"
#include "Light.h"
#include <Eigen/Dense>
#include <math.h>

using namespace Eigen;

class Object
{
    public:
        Vector3f colour;
        Vector3f diffuse;
        Vector3f specular;
        Vector3f sigmaS;
        Vector3f sigmaA;
        float albedo;
        
        Object(Vector3f c)
        {
            colour = c;
            diffuse = colour * 1;
            specular = {0.5, 0.5, 0.5};
        }
        virtual bool intersects(Line ray, float &t, Vector3f &normal, Vector3f &intersectionPoint, int &face) = 0;
        virtual bool intersects(Line ray, float &t) = 0;
        virtual bool isTranslucent() = 0;
        virtual Vector3f randomPoint(Vector3f &normal, int face) = 0;
        virtual float getBoundingBoxIntersect(Line ray) = 0;

        Vector3f getSigmaT()
        {
            return sigmaS + sigmaA;
        }

        Vector3f getReducedSigmaT()
        {
            // Constant phase funtion which causes isotropic scattering so g = 0; 
            // So sigmaT = reducedSigmaT
            return getSigmaT();
        }

        Vector3f getSigmaTR()
        {
            // Constant phase funtion which causes isotropic scattering so g = 0; 
            Vector3f reducedSigmaT = sigmaS + sigmaA;
            Vector3f sigmaTR;
            sigmaTR[0] = sqrt(3 * sigmaA[0] * reducedSigmaT[0]);
            sigmaTR[1] = sqrt(3 * sigmaA[1] * reducedSigmaT[1]);
            sigmaTR[2] = sqrt(3 * sigmaA[2] * reducedSigmaT[2]);
            return sigmaTR;
        }
};