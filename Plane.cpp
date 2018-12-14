#include "Plane.h"

Plane::Plane(Vector3f p, Vector3f n, Vector3f c, float theRange) : Object(c)
{
    point = p;
    normal = n;
    range = theRange;
}

/**
 * Finds if line intersects with plane, returning t value.
 **/
bool Plane::intersects(Line ray, float &t)
{
    // If line is not parallel, then it intersects.
    if(ray.direction.dot(normal) == 0)
    {
        return false;
    }
            
    t = normal.dot(point - ray.origin) / (normal.dot(ray.direction));

    return true;
}

/**
 * Finds if line intersects with plane, returning t value, normal and intersection point.
 **/
bool Plane::intersects(Line ray, float &t, Vector3f &normalRet, Vector3f &intersectionPoint)
{
    // If line is not parallel, then it intersects.
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
    return true;
}

/**
 * Gets a random point on plane within a certain range.
 **/
Vector3f Plane::randomPoint(Vector3f &theNormal)
{
    Vector3f ret;
    float d = -normal[0]*point[0] - normal[1]*point[1] - normal[2]*point[2];
    float mag = sqrt(pow(normal[0], 2) + pow(normal[1], 2) + pow(normal[2], 2));

    if(normal[0] == 1 && mag == 1)
    {
        ret[0] = point[0]; 
    }
    else
    {
        ret[0] = range + rand() / (RAND_MAX / (range*2));
    }
    if(normal[1] == 1 && mag == 1)
    {
        ret[1] = point[1]; 
    }
    else
    {
        ret[1] = range + rand() / (RAND_MAX / (range*2));
    }
    if(normal[2] == 1 && mag == 1)
    {
        ret[2] = point[2]; 
    }
    else
    {
        if(normal[2] == 0)
        {
            ret[2] = range + rand() / (RAND_MAX / (range*2));
        }
        else
        {
            ret[2] = (-d - normal[0] * ret[0] - normal[1] * ret[1]) / normal[2];
        }
    }
    theNormal = normal;
    return ret;
}

float Plane::getBoundingBoxIntersect(Line ray)
{
    // No bounding box, so just return 1. This is used for edge cases so should not be an issue in the case of a plane.
    return 1;
}