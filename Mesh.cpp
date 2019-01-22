#include "Mesh.h"
#include <iostream>

Mesh::Mesh(Eigen::MatrixXf v, Eigen::MatrixXi f, Eigen::MatrixXf vn, Eigen::MatrixXi fn, Vector3f c, bool isTranslucent): Object(c)
{
    calculateTriangles(v, f, vn, fn);

    translucent = isTranslucent;

    float maxLimit = numeric_limits<float>::max();
    float minLimit = numeric_limits<float>::min();
    Vector3f min = {maxLimit, maxLimit, maxLimit};
    Vector3f max = {minLimit, minLimit, minLimit};
    for(int i=0; i< v.rows(); i++)
    {
        Vector3f vertice = v.block < 1, 3 > (i, 0);
        // Positions bunny in front of camera.
        Vector3f translate = {0, -0.5, 0};
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
    //kdNode->boundingBox = boundingBox;
}

/**
 * Checks if a line intersects with the mesh and returns the t value for the closest intersection point.
 **/
bool Mesh::intersects(Line ray, float & t)
{
    if(!boundingBox.intersects(ray))
    {
        //return false;
    }
    Vector3f intersectionPoint;
    Vector3f normal;
    if(kdNode->hit(kdNode, ray, t, t, intersectionPoint, normal))
    {
        return true;
    }
    else
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

/**
 * Checks if a line intersects with the mesh and returns the t value for the closest intersection point.
 * This function also returns the normal of the face intersected with and the intersection point.
 **/
bool Mesh::intersects(Line ray, float & t, Vector3f & normal, Vector3f & intersectionPoint)
{
    if(!boundingBox.intersects(ray))
    {
        //return false;
    }
    if(kdNode->hit(kdNode, ray, t, t, intersectionPoint, normal))
    {
        return true;
    }
    else
    {
        return false;
    }

    bool interpolateNormals = false;
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

                if(interpolateNormals)
                {
                    float totalArea = getTriangleArea(triangles[i].point1, triangles[i].point2, triangles[i].point3);
                    float v = getTriangleArea(triangles[i].point3, triangles[i].point1, a) / totalArea;
                    float u = getTriangleArea(triangles[i].point2, triangles[i].point3, a) / totalArea;

                    normal = u * triangles[i].normal1 + v * triangles[i].normal2 + (1 - v - u) * triangles[i].normal3;
                }
                else
                {
                    normal = triangles[i].plane.normal;
                }
                
                normal.normalize();
            }
        }
    }

    return intersects;
}

/**
 * Adds all triangles to list, this assumes that the faces are given in sets of 3 vertices.
 **/
void Mesh::calculateTriangles(Eigen::MatrixXf vertices, Eigen::MatrixXi faces, Eigen::MatrixXf vertexNormals, Eigen::MatrixXi faceNormals) {
    noTriangles = faces.rows();
    vector<Triangle*>* tris = new vector<Triangle*>(); 
    for (int i = 0; i < noTriangles; i++) {
        // Assumes faces are given in sets of four vertices.
        // Triangle 1
        Vector3f translate = {0, 0, -1};
        Vector3f vertice1 = vertices.block < 1, 3 > (faces.coeff(i, 0), 0);
        Vector3f vertice2 = vertices.block < 1, 3 > (faces.coeff(i, 1), 0);
        Vector3f vertice3 = vertices.block < 1, 3 > (faces.coeff(i, 2), 0);
        Vector3f normal1 = vertexNormals.block < 1, 3 > (faceNormals.coeff(i, 0), 0);
        Vector3f normal2 = vertexNormals.block < 1, 3 > (faceNormals.coeff(i, 1), 0);
        Vector3f normal3 = vertexNormals.block < 1, 3 > (faceNormals.coeff(i, 2), 0);

        // Repositions bunny used so it fits in front of camera.
        vertice1 = (vertice1 / 80) + translate;
        vertice2 = (vertice2 / 80) + translate;
        vertice3 = (vertice3 / 80) + translate;

        Vector3f vector1 = vertice2 - vertice1;
        Vector3f vector2 = vertice3 - vertice2;
        Plane plane = Plane(vertice1, vector1.cross(vector2), colour);
        plane.normal.normalize();
        Triangle* triangle = new Triangle();
        triangle->setValues(
            plane,
            vertice1,
            vertice2,
            vertice3,
            normal1,
            normal2,
            normal3
        );
        tris->push_back(triangle);
        triangles.push_back(*triangle);
    }

    kdNode = kdNode->build(*tris, 0);
}

bool Mesh::isTranslucent()
{
    return translucent;
}

/** Found formula for uniform random point in triangle here, 
 * https://stackoverflow.com/questions/19654251/random-point-inside-triangle-inside-java.
**/
Vector3f Mesh::randomPoint(Vector3f &normal)
{
    float r1 = rand() / RAND_MAX;
    float r2 = rand() / RAND_MAX;
    int r3 = rand() % noTriangles;

    Triangle tri = triangles[r3];
    Vector3f point = (1 - sqrt(r1)) * tri.point1 + (sqrt(r1)*(1-r2)) * tri.point2 + (sqrt(r1)*r2) * tri.point3;
    normal = tri.plane.normal;
    return point;
}

float Mesh::getBoundingBoxIntersect(Line ray)
{
    return boundingBox.getIntersectT(ray);
}

float Mesh::getTriangleArea(Vector3f point1, Vector3f point2, Vector3f point3)
{
    Vector3f ab = point1 - point2;
    Vector3f ac = point1 - point3;
    Vector3f cross = ab.cross(ac);

    return 0.5 * sqrt(pow(cross[0], 2) + pow(cross[1], 2) + pow(cross[2], 2));
}