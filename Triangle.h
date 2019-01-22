#pragma once

#include "Plane.h"
#include "BoundingBox.h"
#include <Eigen/Dense>
#include <math.h>
#include <iostream>

using namespace Eigen;

/**
 * Defines a triangle, used as a face on a mesh.
 **/
class Triangle
{
    public:
        Plane plane;
        Vector3f point1;
        Vector3f point2;
        Vector3f point3;
        Vector3f normal1;
        Vector3f normal2;
        Vector3f normal3;
        BoundingBox boundingBox;

        void setValues(Plane thePlane, Vector3f p1, Vector3f p2, Vector3f p3, Vector3f n1, Vector3f n2, Vector3f n3)
        {
            plane = thePlane;
            point1 = p1;
            point2 = p2;
            point3 = p3;

            boundingBox = BoundingBox();
            float minX = point1[0];float maxX = point1[0];
            float minY = point1[1]; float maxY = point1[1];
            float minZ = point1[2]; float maxZ = point1[2];
            if(point2[0] < minX)
            {
                minX = point2[0];
            }
            if(point2[1] < minY)
            {
                minY = point2[1];
            }
            if(point2[2] < minZ)
            {
                minZ = point2[2];
            }
            if(point3[0] < minX)
            {
                minX = point3[0];
            }
            if(point3[1] < minY)
            {
                minY = point3[1];
            }
            if(point3[2] < minZ)
            {
                minZ = point3[2];
            }

            if(point2[0] > maxX)
            {
                maxX = point2[0];
            }
            if(point2[1] > maxY)
            {
                maxY = point2[1];
            }
            if(point2[2] > maxZ)
            {
                maxZ = point2[2];
            }
            if(point3[0] > maxX)
            {
                maxX = point3[0];
            }
            if(point3[1] > maxY)
            {
                maxY = point3[1];
            }
            if(point3[2] > maxZ)
            {
                maxZ = point3[2];
            }
            Vector3f min = {minX, minY, minZ};
            Vector3f max = {maxX, maxY, maxZ};
            boundingBox.setValues(min, max);
        }

        Vector3f midPoint()
        {
            Vector3f midpt;
            midpt[0] = (point1[0] + point2[0] + point3[0]) / 3;
            midpt[1] = (point1[1] + point2[1] + point3[1]) / 3;
            midpt[2] = (point1[2] + point2[2] + point3[2]) / 3;
            return midpt;
        }

        bool intersect(Line ray, Vector3f &intersectionPoint, float &t, Vector3f &normal)
        {
            bool intersects = false;
            float newT = plane.normal.dot(plane.point - ray.origin) /
            (plane.normal.dot(ray.direction));

        // If intersection is on or behind origin then ignore. 
        if (newT < 0.00001) {
            return false;
        }

        // Intersection point
        Vector3f a = ray.origin + newT * ray.direction;

        Vector3f normal1 = (a - point1).cross(point2 - point1);
        Vector3f normal2 = (a - point2).cross(point3 - point2);
        Vector3f normal3 = (a - point3).cross(point1 - point3);

        if (normal1.dot(normal2) >= 0) {
            if (normal2.dot(normal3) >= 0) {
                
                t = newT;
                intersectionPoint = a;
                intersects = true;

                    normal = plane.normal;
            
                
                normal.normalize();
            }
        }
        return intersects;
        }
};