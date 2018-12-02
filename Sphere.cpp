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

bool Sphere::intersects(Line ray, Light light, float & t, float & ambAngle, float & specAngle, Vector3f & intersectionPoint) 
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
        Vector3f normal = a - origin;
        normal.normalize();
        ambAngle = normal.dot(-light.dir);
        if (ambAngle < 0) {
            ambAngle = 0;
        }
        Vector3f reflection = light.dir - 2 * (light.dir.dot(normal)) * normal;
        specAngle = reflection.dot(-ray.direction);
        if (specAngle < 0) {
            specAngle = 0;
        }

        return true;
    }
    return false;
}