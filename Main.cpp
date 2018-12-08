#include "Plane.h"
#include "Sphere.h"
#include "Camera.h"
#include "Mesh.h"
#include "Scene.h"
#include <igl/readOBJ.h>
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <vector>
#include <cmath>
#include <limits>
#include <time.h>

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

    // Sub-surface scattering terms
    int singleScatteringSamples = 10;
    int multipleScatteringSamples = 5;

    // Seed random generator.
    srand(time(NULL));

    /*/ Setup uniform random generator.
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<float> dis(0, RAND_MAX);*/

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
    Vector3f bunnyColour = {0.83, 0.79, 0.75};
    Mesh bunny(V, F, bunnyColour, true);
    bunny.sigmaS = {2.19, 2.62, 3.00};
    bunny.sigmaT = {0.0021, 0.0041, 0.0071};

    Plane plane = Plane({-1, 0, 0}, {1, 0, 0}, {1, 1, 1});

    Camera camera = {{0, 0, 2}, 1, 2, 2, {0, 1, 0}, {1, 0, 0}, {0, 0, 1}};
    Light light = {{100, 100, 100}, {-0.5, -0.5, -1}, {1, 1, 1}, 40};
    light.dir.normalize();

    Scene scene = Scene(light);
    scene.addObject(bunny);
    scene.addObject(plane);
    
    int workCount = 0;

    for(int i=0; i<width; i++)
    {
        for(int j=0; j<height; j++)
        {
            //cout << i << "   " << j << endl;
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

            scene.rayTrace(ray, canvas[i][j], singleScatteringSamples, multipleScatteringSamples);

            if(++workCount % 100 == 0)
            {
                cout << (float) (workCount * 100) / (width*height) << "%" << endl;  
            }
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