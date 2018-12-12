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
#include <thread>

using namespace std;
using namespace Eigen;

int main(int argc, char ** argv)
{
    int width = 400;
    int height = 400;
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
    int singleScatteringSamples = 20;
    int multipleScatteringSamples = 200;

    // Seed random generator.
    srand(time(NULL));
    
    //These two aren't used.
    Eigen::Matrix2Xf tCo;
    Eigen::MatrixXi tF;

    Eigen::MatrixXf V; // Vertices
    Eigen::MatrixXi F; // Faces
    Eigen::MatrixXf VN; // Vertex normals
    Eigen::MatrixXi FN; // Face indices to vertex normals
    igl::readOBJ("133376.obj", V, tCo, VN, F, tF, FN);
    Vector3f bunnyColour = {0.83, 0.79, 0.75};
    Mesh bunny(V, F, VN, FN, bunnyColour, true);
    bunny.sigmaS = {2.19, 2.62, 3.00};
    bunny.sigmaA = {0.0021, 0.0041, 0.0071};
    bunny.albedo = 0.5;

    Plane plane = Plane({-1, 0, 0}, {1, 0, 0}, {1, 1, 1});

    Camera camera = {{1.5, 0.3, 0}, 1, 2, 2, {0, 0, 1}, {0, 1, 0}, {1, 0, 0}};
    Light light = {{2, 3.5, 0}, {1, 1, 1}, 1.5};

    Scene scene = Scene(light);
    scene.addObject(bunny);
    //scene.addObject(plane);
    
    int workCount = 0;

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