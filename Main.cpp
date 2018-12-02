#include "Plane.h"
#include "Sphere.h"
#include "Camera.h"
#include "Mesh.h"
#include <igl/readOBJ.h>
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <vector>
#include <cmath>
#include <limits>

using namespace std;
using namespace Eigen;

int main(int argc, char ** argv)
{
    int width = 100;
    int height = 100;
    float *** canvas = (float ***) malloc(width * sizeof(float **));
    for(int i = 0; i < width; i++)
    {
        canvas[i] = (float **) malloc(height * sizeof(float *));
        for(int j = 0; j < height; j++)
        {
            canvas[i][j] = (float *) malloc(3 * sizeof(float));
        }
    }

    /*Vector3f sphere1Origin = {0, 1, 3};
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
    Sphere sphere4(1, sphere4Origin, sphere4Colour);*/
    
    Eigen::MatrixXf V;
    Eigen::MatrixXi F;
    igl::readOBJ("Bunny-LowPoly.obj", V, F);
    Vector3f cubeColour = {0, 1, 0};
    Mesh bunny(V, F, cubeColour);

    Plane plane = Plane({-1, 0, 0}, {1, 0, 0}, {1, 1, 1});

    vector<Object*> objects;
    //objects.push_back(&sphere1);
    //objects.push_back(&sphere2);
    //objects.push_back(&sphere3);
    //objects.push_back(&sphere4);
    //objects.push_back(&plane);
    objects.push_back(&bunny);

    Camera camera = {{0, 0, -2}, 1, 2, 2, {0, 1, 0}, {1, 0, 0}, {0, 0, -1}};
    Light light = {{-0.5, -0.5, 1}, 0.2, 0.8};
    light.dir.normalize();
    
    for(int i=0; i<width; i++)
    {
        for(int j=0; j<height; j++)
        {
            cout << i << "   " << j << endl;
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

    for(int i = 0; i < width; i++)
    {
        for(int j = 0; j < height; j++)
        {
            free(canvas[i][j]);
        }
        free(canvas[i]);
    }
    free(canvas);
}