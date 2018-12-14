#include "Sphere.h"

Sphere::Sphere(float r, Vector3f o, Vector3f c): Object(c) 
{
    radius = r;
    origin = o;
}

bool Sphere::intersects(Line ray, float & t) 
{
    Line rayObjectSpace = ray;
    rayObjectSpace.origin -= origin;

    float b = 2 * rayObjectSpace.direction.dot(rayObjectSpace.origin);
    float a = rayObjectSpace.direction.dot(rayObjectSpace.direction);
    float discriminant = (b * b) - 4 * a * (rayObjectSpace.origin.dot(rayObjectSpace.origin) - radius * radius);

    if (discriminant > 0) {
        float t1 = (-b + (float) sqrt(discriminant)) / (2 * a);
        float t2 = (-b - (float) sqrt(discriminant)) / (2 * a);
        if (t1 < t2 && t1 > 0) {
            t = t1;
        } else if (t2 > 0) {
            t = t2;
        }
        return true;
    }
    return false;
}

bool Sphere::intersects(Line ray, float & t, Vector3f & normal, Vector3f & intersectionPoint) 
{
    Line rayObjectSpace = ray;
    rayObjectSpace.origin -= origin;

    float b = 2 * rayObjectSpace.direction.dot(rayObjectSpace.origin);
    float a = rayObjectSpace.direction.dot(rayObjectSpace.direction);
    float discriminant = (b * b) - 4 * a * (rayObjectSpace.origin.dot(rayObjectSpace.origin) - radius * radius);

    if (discriminant > 0) {
        float t1 = (-b + (float) sqrt(discriminant)) / (2 * a);
        float t2 = (-b - (float) sqrt(discriminant)) / (2 * a);
        if (t1 < t2 && t1 > 0) {
            t = t1;
        } else if (t2 > 0) {
            t = t2;
        }

        // Intersection point
        Vector3f a = ray.origin + t * ray.direction;
        intersectionPoint = a;
        // Normal for lighting calculations
        Vector3f normal = a - origin;
        normal.normalize();

        return true;
    }
    return false;
}

bool Sphere::isTranslucent()
{
    return false;
}

Vector3f Sphere::randomPoint(Vector3f &normal)
{
    // TODO
    Vector3f point;
    return point;
}

float Sphere::getBoundingBoxIntersect(Line ray)
{
    // N/A
    return -1;
}