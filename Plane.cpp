#include "Plane.h"

Plane::Plane(Vector3f p, Vector3f n, Vector3f c) : Object(c)
{
    point = p;
    normal = n;
}

bool Plane::intersects(Line ray, float &t)
{
    if(ray.direction.dot(normal) == 0)
    {
        return false;
    }
            
    t = normal.dot(point - ray.origin) / (normal.dot(ray.direction));

    return true;
}

bool Plane::intersects(Line ray, float &t, Vector3f &normalRet, Vector3f &intersectionPoint)
{
    if(ray.direction.dot(normal) == 0)
    {
        return false;
    }
            
    t = normal.dot(point - ray.origin) / (normal.dot(ray.direction));
    intersectionPoint = ray.origin + t*ray.direction;
    normalRet = normal;   

    return true;
}

bool Plane::isTranslucent()
{
    return false;
}