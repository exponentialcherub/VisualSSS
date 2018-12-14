#pragma once

#include "Object.h"

/**
 * Defines a plane with a point and a normal.
 **/
class Plane : public Object
{
    public:

        Vector3f point;
        Vector3f normal;
        float range;

        Plane() : Object({0, 0, 0}){}
        Plane(Vector3f p, Vector3f n, Vector3f c, float range = 100);

        bool intersects(Line ray, float &t);
        bool intersects(Line ray, float &t, Vector3f &normalRet, Vector3f &intersectionPoint);
        bool isTranslucent();
        Vector3f randomPoint(Vector3f &normal);
        float getBoundingBoxIntersect(Line ray);
};