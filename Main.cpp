#include <igl/readOBJ.h>
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <vector>
#include <cmath>
#include <limits>

using namespace std;
using namespace Eigen;

struct Camera
{
    Vector3f focalPoint;
    float focalLength;
    float height;
    float width;
    Vector3f upVector;
    Vector3f rightVector;
    Vector3f forwardVector;
};

struct Light
{
    Vector3f dir;
    float ambientI;
    float localI;
};

struct Line
{
    Vector3f origin;
    Vector3f direction;
};

class Object
{
    public:
        Vector3f colour;
        Vector3f diffuse;
        Vector3f specular;
        
        Object(Vector3f c)
        {
            colour = c;
            diffuse = colour * 0.5;
            specular = {0.5, 0.5, 0.5};
        }
        virtual bool intersects(Line ray, Light light, float &t, float &ambAngle, float &specAngle, Vector3f &intersectionPoint) = 0;
        virtual bool intersects(Line ray, float &t) = 0;
};

class Plane : public Object
{
    public:

        Vector3f point;
        Vector3f normal;

        Plane(Vector3f p, Vector3f n, Vector3f c) : Object(c)
        {
            point = p;
            normal = n;
        }

        bool intersects(Line ray, float &t)
        {
            if(ray.direction.dot(normal) == 0)
            {
                return false;
            }
            
            t = normal.dot(point - ray.origin) / (normal.dot(ray.direction));

            return true;
        }

        bool intersects(Line ray, Light light, float &t, float &ambAngle, float &specAngle, Vector3f &intersectionPoint)
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
};

struct Triangle
{
    Plane plane;
    Vector3f point1;
    Vector3f point2;
    Vector3f point3;
};

class Sphere : public Object
{
    float radius;
    Vector3f origin;

    public:
        Sphere(float r, Vector3f o, Vector3f c) : Object(c)
        {
            radius = r;
            origin = o;
        }

        bool intersects(Line ray, float &t)
        {
            Line rayObjectSpace = ray;
            rayObjectSpace.origin -= origin;
        
            float b = 2 * rayObjectSpace.direction.dot(rayObjectSpace.origin);
            float a = rayObjectSpace.direction.dot(rayObjectSpace.direction);
            float discriminant = (b*b) - 4*a*(rayObjectSpace.origin.dot(rayObjectSpace.origin) - radius * radius);

            if(discriminant > 0)
            {
                float t1 = (-b + (float)sqrt(discriminant)) / (2*a);
                float t2 = (-b - (float)sqrt(discriminant)) / (2*a);
                if(t1 < t2 && t1 > 0)
                {
                    t = t1;
                }
                else if(t2 > 0)
                {
                    t = t2;
                }
                return true;
            }
            return false;
        }

        bool intersects(Line ray, Light light, float &t, float &ambAngle, float &specAngle, Vector3f &intersectionPoint) override
        {
            Line rayObjectSpace = ray;
            rayObjectSpace.origin -= origin;
        
            float b = 2 * rayObjectSpace.direction.dot(rayObjectSpace.origin);
            float a = rayObjectSpace.direction.dot(rayObjectSpace.direction);
            float discriminant = (b*b) - 4*a*(rayObjectSpace.origin.dot(rayObjectSpace.origin) - radius * radius);

            if(discriminant > 0)
            {
                float t1 = (-b + (float)sqrt(discriminant)) / (2*a);
                float t2 = (-b - (float)sqrt(discriminant)) / (2*a);
                if(t1 < t2 && t1 > 0)
                {
                    t = t1;
                }
                else if(t2 > 0)
                {
                    t = t2;
                }

                // Intersection point
                Vector3f a = ray.origin + t*ray.direction;
                intersectionPoint = a;
                Vector3f normal = a - origin;
                normal.normalize();
                ambAngle = normal.dot(-light.dir);
                if(ambAngle < 0)
                {
                    ambAngle = 0;
                }
                Vector3f reflection = light.dir - 2*(light.dir.dot(normal))*normal;
                specAngle = reflection.dot(-ray.direction);
                if(specAngle < 0)
                {
                    specAngle = 0;
                }

                return true;
            }
            return false;
        }
};

class Mesh : public Object
{
    Eigen::MatrixXf vertices;
    Eigen::MatrixXi faces;
    vector<Triangle> triangles;
    int noTriangles;

    public:
        Mesh(Eigen::MatrixXf v, Eigen::MatrixXi f, Vector3f c) : Object(c)
        {
            vertices = v;
            faces = f;

            calculateTriangles();
        }
        
        bool intersects(Line ray, float &t)
        {
            bool intersects = false;

            for(int i=0; i<noTriangles; i++)
            {
                t = triangles[i].plane.normal.dot(triangles[i].plane.point - ray.origin) 
                    / (triangles[i].plane.normal.dot(ray.direction));

                // Intersection point
                Vector3f a = ray.origin + t*ray.direction;
                
                Vector3f normal1 = (a - triangles[i].point1).cross(triangles[i].point2 - triangles[i].point1);
                Vector3f normal2 = (a - triangles[i].point2).cross(triangles[i].point3 - triangles[i].point2);
                Vector3f normal3 = (a - triangles[i].point3).cross(triangles[i].point1 - triangles[i].point3);

                if(normal1.dot(normal2) >= 0)
                {
                    if(normal2.dot(normal3) >= 0)
                    {
                        intersects = true;
                    }
                }
            }

            return intersects;
        }

        bool intersects(Line ray, Light light, float &t, float &ambAngle, float &specAngle, Vector3f &intersectionPoint) override
        {
            bool intersects = false;
            t = numeric_limits<float>::max();

            for(int i=0; i<noTriangles; i++)
            {
                float newT = triangles[i].plane.normal.dot(triangles[i].plane.point - ray.origin) 
                    / (triangles[i].plane.normal.dot(ray.direction));

                if(newT < 0)
                {
                    continue;
                }

                // Intersection point
                Vector3f a = ray.origin + newT*ray.direction;
                
                Vector3f normal1 = (a - triangles[i].point1).cross(triangles[i].point2 - triangles[i].point1);
                Vector3f normal2 = (a - triangles[i].point2).cross(triangles[i].point3 - triangles[i].point2);
                Vector3f normal3 = (a - triangles[i].point3).cross(triangles[i].point1 - triangles[i].point3);

                if(normal1.dot(normal2) >= 0)
                {
                    if(normal2.dot(normal3) >= 0)
                    {
                        if(newT > t)
                        {
                            continue;
                        }
                        t = newT;
                        intersectionPoint = a;
                        intersects = true;
                        
                        ambAngle = triangles[i].plane.normal.dot(-light.dir);
                        if(ambAngle < 0)
                        {
                            ambAngle = 0;
                        }

                        Vector3f reflection = light.dir - 2*(light.dir.dot(triangles[i].plane.normal))*triangles[i].plane.normal;
                        specAngle = reflection.dot(triangles[i].plane.normal);
                        if(specAngle < 0)
                        {
                            specAngle = 0;
                        }
                    }
                }
            }

            return intersects;
        }

    private:
        /**
         * Crude method for extracting triangles for cube. Eventually this should be made into a generic mesh class. As
         * my cube object gives faces in sets of 4 vertices.
         **/
        void calculateTriangles()
        {
            noTriangles = faces.rows()*2;
            for(int i=0; i<faces.rows(); i++)
            {
                // Assumes faces are given in sets of four vertices.
                // Triangle 1
                Vector3f vertice1 = vertices.block<1, 3>(faces.coeff(i, 0), 0);
                Vector3f vertice2 = vertices.block<1, 3>(faces.coeff(i, 1), 0);
                Vector3f vertice3 = vertices.block<1, 3>(faces.coeff(i, 2), 0);
                Vector3f vertice4 = vertices.block<1, 3>(faces.coeff(i, 3), 0);
                
                Vector3f vector1 = vertice2 - vertice1;
                Vector3f vector2 = vertice3 - vertice2;
                Plane plane = Plane(vertice1, vector1.cross(vector2), colour);
                plane.normal.normalize();
                Triangle triangle1 = {plane, vertice1, vertice2, vertice3};

                // Triangle 2
                Triangle triangle2 = {plane, vertice3, vertice4, vertice1};

                triangles.push_back(triangle1);
                triangles.push_back(triangle2);
            }
        }
};

int main(int argc, char ** argv)
{
    int width = 400;
    int height = 400;
    float canvas[400][400][3] = {0};

    Vector3f sphere1Origin = {0, 1, 3};
    Vector3f sphere1Colour = {1, 0, 0};
    Sphere sphere1(1, sphere1Origin, sphere1Colour);

    Vector3f sphere2Origin = {0, -1, 5};
    Vector3f sphere2Colour = {0, 1, 0};
    Sphere sphere2(1, sphere2Origin, sphere2Colour);

    Vector3f sphere3Origin = {-1, -0.5, 4};
    Vector3f sphere3Colour = {1, 0.8, 0.8};
    Sphere sphere3(0.5, sphere3Origin, sphere3Colour);

    Vector3f sphere4Origin = {1, 0, 4};
    Vector3f sphere4Colour = {0.8, 0.6, 0.1};
    Sphere sphere4(1, sphere4Origin, sphere4Colour);
    
    /*Eigen::MatrixXf V;
    Eigen::MatrixXi F;
    igl::readOBJ("cube.obj", V, F);
    Vector3f cubeColour = {0, 1, 0};
    Mesh cube(V, F, cubeColour);*/

    Plane plane = Plane({-1, 0, 0}, {1, 0, 0}, {1, 1, 1});

    vector<Object*> objects;
    objects.push_back(&sphere1);
    objects.push_back(&sphere2);
    objects.push_back(&sphere3);
    objects.push_back(&sphere4);
    objects.push_back(&plane);

    Camera camera = {{0, 0, 0}, 1, 2, 2, {0, 1, 0}, {1, 0, 0}, {0, 0, -1}};
    Light light = {{-0.5, -0.5, 1}, 0.2, 0.8};
    light.dir.normalize();
    
    for(int i=0; i<width; i++)
    {
        for(int j=0; j<height; j++)
        {
            // Calculate point of pixel.
            float r = camera.width * (((i + 0.5) / width) - 0.5);
            float b = camera.height * (((j + 0.5) / height) - 0.5);

            float x = camera.focalPoint[0] - camera.focalLength*camera.forwardVector[0] - r*camera.rightVector[0] 
                + b *camera.upVector[0];
            float y = camera.focalPoint[1] - camera.focalLength*camera.forwardVector[1] - r*camera.rightVector[1] 
                + b *camera.upVector[1];
            float z = camera.focalPoint[2] - camera.focalLength*camera.forwardVector[2] - r*camera.rightVector[2] 
                + b *camera.upVector[2];

            Vector3f pixelPoint = {x, y, z};
            Vector3f dir = pixelPoint - camera.focalPoint;
            dir.normalize();
            
            Line ray = {pixelPoint, dir};

            float t = numeric_limits<float>::max();
            float ambAngle = 0;
            float specAngle = 0;
            int closest = -1;
            for(int k=0; k<objects.size(); k++)
            {
                float tNew = 0;
                float newAmbAngle;
                float newSpecAngle;
                Vector3f intersectionPoint;

                if(objects[k]->intersects(ray, light, tNew, newAmbAngle, newSpecAngle, intersectionPoint))
                {
                    // Check t > 0 to be sure.
                    if(tNew < t && tNew > 0)
                    {
                        t = tNew;
                        ambAngle = newAmbAngle;
                        specAngle = newSpecAngle;
                        closest = k;

                        for(int h=0; h<objects.size(); h++)
                        {
                            if(k == h)
                            {
                                continue;
                            }
                            
                            Line lightRay = {intersectionPoint + 0.01*(-light.dir), -light.dir};
                            
                            float shadowT = -1;

                            if(objects[h]->intersects(lightRay, shadowT) && shadowT > 0)
                            {
                                ambAngle = 0;
                                specAngle = 0;
                            }
                        }
                    }
                }
            }
            if(closest == -1)
            {
                continue;
            }
            
            // Calculate light here I think, getting normal from intersection method
            canvas[i][j][0] = objects[closest]->colour[0] * light.ambientI + ambAngle*objects[closest]->diffuse[0]*light.localI + pow(specAngle, 20)*objects[closest]->specular[0];
            canvas[i][j][1] = objects[closest]->colour[1] * light.ambientI + ambAngle*objects[closest]->diffuse[1]*light.localI + pow(specAngle, 20)*objects[closest]->specular[1];
            canvas[i][j][2] = objects[closest]->colour[2] * light.ambientI + ambAngle*objects[closest]->diffuse[2]*light.localI + pow(specAngle, 20)*objects[closest]->specular[2];
        }
    }

    ofstream image ("image.ppm");
    image << "P3" << endl;
    image << width << " " << height << endl;
    image << "255" << endl;
    for(int i=0; i<width; i++)
    {
        for(int j=0; j<height; j++)
        {
            image << canvas[i][j][0] * 255 << " " << canvas[i][j][1] * 255 << " " << canvas[i][j][2] * 255 << "  "; 
        }
        image << "\n";
    }

    image.close();
}