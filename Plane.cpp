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

bool Plane::intersects(Line ray, Light light, float &t, float &ambAngle, float &specAngle, Vector3f &intersectionPoint)
{
    if(ray.direction.dot(normal) == 0)
    {
        return false;
    }
            
    t = normal.dot(point - ray.origin) / (normal.dot(ray.direction));

    intersectionPoint = ray.origin + t*ray.direction;
                        
    ambAngle = normal.dot(-light.dir);
    if(ambAngle < 0)
    {
        ambAngle = 0;
    }

    Vector3f reflection = light.dir - 2*(light.dir.dot(normal))*normal;
    specAngle = reflection.dot(normal);
    if(specAngle < 0)
    {
        specAngle = 0;
    }   

    return true;
}

bool Plane::isTranslucent()
{
    return false;
}