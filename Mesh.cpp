#include "Mesh.h"

Mesh::Mesh(Eigen::MatrixXf v, Eigen::MatrixXi f, Vector3f c): Object(c)
{
    vertices = v;
    faces = f;

    calculateTriangles();
}

bool Mesh::intersects(Line ray, float & t)
{
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

bool Mesh::intersects(Line ray, Light light, float & t, float & ambAngle, float & specAngle, Vector3f & intersectionPoint)
{
    bool intersects = false;
    t = numeric_limits < float > ::max();

    for (int i = 0; i < noTriangles; i++) {
        float newT = triangles[i].plane.normal.dot(triangles[i].plane.point - ray.origin) /
            (triangles[i].plane.normal.dot(ray.direction));

        if (newT < 0) {
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

                ambAngle = triangles[i].plane.normal.dot(-light.dir);
                if (ambAngle < 0) {
                    ambAngle = 0;
                }

                Vector3f reflection = light.dir - 2 * (light.dir.dot(triangles[i].plane.normal)) * triangles[i].plane.normal;
                specAngle = reflection.dot(triangles[i].plane.normal);
                if (specAngle < 0) {
                    specAngle = 0;
                }
            }
        }
    }

    return intersects;
}

/**
 * Crude method for extracting triangles for cube. Eventually this should be made into a generic mesh class. As
 * my cube object gives faces in sets of 4 vertices.
 **/
void Mesh::calculateTriangles() {
    noTriangles = faces.rows() * 2;
    for (int i = 0; i < faces.rows(); i++) {
        // Assumes faces are given in sets of four vertices.
        // Triangle 1
        Vector3f vertice1 = vertices.block < 1, 3 > (faces.coeff(i, 0), 0);
        Vector3f vertice2 = vertices.block < 1, 3 > (faces.coeff(i, 1), 0);
        Vector3f vertice3 = vertices.block < 1, 3 > (faces.coeff(i, 2), 0);
        Vector3f vertice4 = vertices.block < 1, 3 > (faces.coeff(i, 3), 0);

        Vector3f vector1 = vertice2 - vertice1;
        Vector3f vector2 = vertice3 - vertice2;
        Plane plane = Plane(vertice1, vector1.cross(vector2), colour);
        plane.normal.normalize();
        Triangle triangle1 = {
            plane,
            vertice1,
            vertice2,
            vertice3
        };

        // Triangle 2
        Triangle triangle2 = {
            plane,
            vertice3,
            vertice4,
            vertice1
        };

        triangles.push_back(triangle1);
        triangles.push_back(triangle2);
    }
}