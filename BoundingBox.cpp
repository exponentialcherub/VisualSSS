#include "BoundingBox.h"

// Calculates if a ray instersects with bounding cube.
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

void BoundingBox::setValues(Vector3f minimum, Vector3f maximum)
{
    min = minimum;
    max = maximum;
}