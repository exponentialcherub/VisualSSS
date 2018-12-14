#include "BoundingBox.h"

/**Calculates if a ray instersects with bounding cube.
*  Based on code found here: https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-box-intersection.
**/ 
bool BoundingBox::intersects(Line ray)
{
    // If ray originates in the bounding box then it will always intersect with the cube.
    if(ray.origin[0] <= max[0] && ray.origin[1] <= max[1] && ray.origin[2] <= max[2] &&
       ray.origin[0] >= min[0] && ray.origin[1] >= min[1] && ray.origin[2] >= min[2])
    {
        return true;
    }

    float tMinX = (min[0] - ray.origin[0]) / ray.direction[0];
    float tMaxX = (max[0] - ray.origin[0]) / ray.direction[0];
    float tMinY = (min[1] - ray.origin[1]) / ray.direction[1];
    float tMaxY = (max[1] - ray.origin[1]) / ray.direction[1];
    if(tMinX > tMaxX)
    {
        float temp = tMinX;
        tMinX = tMaxX;
        tMaxX = temp;
    }
    if(tMinY > tMaxY)
    {
        float temp = tMinY;
        tMinY = tMaxY;
        tMaxY = temp;
    }

    if ((tMinX > tMaxY) || (tMinY > tMaxX))
    { 
        return false; 
    }

    float tMinZ = (min[2] - ray.origin[2]) / ray.direction[2];
    float tMaxZ = (max[2] - ray.origin[2]) / ray.direction[2];
    if(tMinZ > tMaxZ)
    {
        float temp = tMinZ;
        tMinZ = tMaxZ;
        tMaxZ = temp;
    }

    if (tMaxY > tMaxX)
    {
        tMaxX = tMaxY;
    } 
 
    if (tMinY < tMinX)
    { 
        tMinX = tMinY; 
    }

    if ((tMinX > tMaxZ) || (tMinZ > tMaxX))
    { 
        return false; 
    }

    return true;
}

/**
 * Gets the t intersect of a ray with the bounding box. All plane faces are parallel with world axis so we know how to
 * define planes.
 **/
float BoundingBox::getIntersectT(Line ray)
{
    float tValues[6];

    // Intersection with x plane.
    Plane planeMinX(min, {1, 0, 0}, {0, 0, 0});
    Plane planeMaxX(max, {1, 0, 0}, {0, 0, 0});
    planeMinX.intersects(ray, tValues[0]);
    planeMaxX.intersects(ray, tValues[1]);

    // Intersection with y plane.
    Plane planeMinY(min, {0, 1, 0}, {0, 0, 0});
    Plane planeMaxY(max, {0, 1, 0}, {0, 0, 0});
    planeMinY.intersects(ray, tValues[2]);
    planeMaxY.intersects(ray, tValues[3]);

    // Intersection with z plane.
    Plane planeMinZ(min, {0, 0, 1}, {0, 0, 0});
    Plane planeMaxZ(max, {0, 0, 1}, {0, 0, 0});
    planeMinZ.intersects(ray, tValues[4]);
    planeMaxZ.intersects(ray, tValues[5]);

    float t = std::numeric_limits<float>::max();
    for(int i=0; i<6; i++)
    {
        // Get smallest value (closest intersection).
        if(tValues[i] < t && tValues[i] > 0)
        {
            t = tValues[i];
        }
    }

    return t;
}

void BoundingBox::setValues(Vector3f minimum, Vector3f maximum)
{
    min = minimum;
    max = maximum;
}