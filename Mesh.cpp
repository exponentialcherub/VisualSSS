#include "Mesh.h"
#include <iostream>
Mesh::Mesh(Eigen::MatrixXf v, Eigen::MatrixXi f, Vector3f c, bool isTranslucent): Object(c)
{
    calculateTriangles(v, f);

    translucent = isTranslucent;

    float maxLimit = numeric_limits<float>::max();
    float minLimit = numeric_limits<float>::min();
    Vector3f min = {maxLimit, maxLimit, maxLimit};
    Vector3f max = {minLimit, minLimit, minLimit};
    for(int i=0; i< v.rows(); i++)
    {
        // Need to change this, shouldn't be making points smaller by 50, get updated model.
        Vector3f vertice = v.block < 1, 3 > (i, 0);// 50;
        Vector3f translate = {0, 0, -0.5};
        vertice = (vertice / 40) + translate;
        if(vertice[0] < min[0])
            min[0] = vertice[0];
        if(vertice[1] < min[1])
            min[1] = vertice[1];
        if(vertice[2] < min[2])
            min[2] = vertice[2];
        if(vertice[0] > max[0])
            max[0] = vertice[0];
        if(vertice[1] > max[1])
            max[1] = vertice[1];
        if(vertice[2] > max[2])
            max[2] = vertice[2];
    }

    boundingBox.setValues(min, max);
}

bool Mesh::intersects(Line ray, float & t)
{
    if(!boundingBox.intersects(ray))
    {
        return false;
    }

    bool intersects = false;

    for (int i = 0; i < noTriangles; i++) {
        t = triangles[i].plane.normal.dot(triangles[i].plane.point - ray.origin) /
            (triangles[i].plane.normal.dot(ray.direction));

        // Intersection point
        Vector3f a = ray.origin + t * ray.direction;

        Vector3f normal1 = (a - triangles[i].point1).cross(triangles[i].point2 - triangles[i].point1);
        Vector3f normal2 = (a - triangles[i].point2).cross(triangles[i].point3 - triangles[i].point2);
        Vector3f normal3 = (a - triangles[i].point3).cross(triangles[i].point1 - triangles[i].point3);

        if (normal1.dot(normal2) >= 0) {
            if (normal2.dot(normal3) >= 0) {
                intersects = true;
            }
        }
    }

    return intersects;
}

bool Mesh::intersects(Line ray, float & t, Vector3f & normal, Vector3f & intersectionPoint)
{
    if(!boundingBox.intersects(ray))
    {
        return false;
    }

    bool intersects = false;
    t = numeric_limits < float > ::max();

    for (int i = 0; i < noTriangles; i++) {
        float newT = triangles[i].plane.normal.dot(triangles[i].plane.point - ray.origin) /
            (triangles[i].plane.normal.dot(ray.direction));

        // If intersection is on or behind origin then ignore. 
        if (newT < 0.00001) {
            continue;
        }

        // Intersection point
        Vector3f a = ray.origin + newT * ray.direction;

        Vector3f normal1 = (a - triangles[i].point1).cross(triangles[i].point2 - triangles[i].point1);
        Vector3f normal2 = (a - triangles[i].point2).cross(triangles[i].point3 - triangles[i].point2);
        Vector3f normal3 = (a - triangles[i].point3).cross(triangles[i].point1 - triangles[i].point3);

        if (normal1.dot(normal2) >= 0) {
            if (normal2.dot(normal3) >= 0) {
                if (newT > t) {
                    continue;
                }
                t = newT;
                intersectionPoint = a;
                intersects = true;

                normal = triangles[i].plane.normal;
                normal.normalize();
            }
        }
    }

    return intersects;
}

/**
 * Adds all triangles to list, this assumes that the faces are given in sets of 3 vertices.
 **/
void Mesh::calculateTriangles(Eigen::MatrixXf vertices, Eigen::MatrixXi faces) {
    noTriangles = faces.rows();
    for (int i = 0; i < faces.rows(); i++) {
        // Assumes faces are given in sets of four vertices.
        // Triangle 1
        Vector3f translate = {0, 0, -0.5};
        Vector3f vertice1 = vertices.block < 1, 3 > (faces.coeff(i, 0), 0);// 50;
        Vector3f vertice2 = vertices.block < 1, 3 > (faces.coeff(i, 1), 0);// 50;
        Vector3f vertice3 = vertices.block < 1, 3 > (faces.coeff(i, 2), 0);// 50;
        vertice1 = (vertice1 / 40) + translate;
        vertice2 = (vertice2 / 40) + translate;
        vertice3 = (vertice3 / 40) + translate;
        // Rotate bunny towards camera
        // Ideally don't do this or resizing above
        /*Vector3f temp = {(float)(vertice1[2]-1), vertice1[1], -vertice1[0]};
        vertice1 = temp;
        temp = {(float)(vertice2[2]-1), vertice2[1], -vertice2[0]};
        vertice2 = temp;
        temp = {(float)(vertice3[2]-1), vertice3[1], -vertice3[0]};
        vertice3 = temp;*/

        Vector3f vector1 = vertice2 - vertice1;
        Vector3f vector2 = vertice3 - vertice2;
        Plane plane = Plane(vertice1, vector1.cross(vector2), colour);
        plane.normal.normalize();
        Triangle triangle;
        triangle.setValues(
            plane,
            vertice1,
            vertice2,
            vertice3
        );

        triangles.push_back(triangle);
    }
}

bool Mesh::isTranslucent()
{
    return translucent;
}

// Found formula for uniform random point in triangle here, 
// https://stackoverflow.com/questions/19654251/random-point-inside-triangle-inside-java.
Vector3f Mesh::randomPoint(Vector3f &normal)
{
    float r1 = rand() / RAND_MAX;
    float r2 = rand() / RAND_MAX;
    float r3 = rand() % noTriangles;

    Triangle tri = triangles[r3];
    Vector3f point = (1 - sqrt(r1)) * tri.point1 + (sqrt(r1)*(1-r2)) * tri.point2 + (sqrt(r1)*r2) * tri.point3;
    normal = tri.plane.normal;
    return point;
}

float Mesh::getBoundingBoxIntersect(Line ray)
{
    return boundingBox.getIntersectT(ray);
}