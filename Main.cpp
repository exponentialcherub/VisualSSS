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
    int width = 150;
    int height = 150;
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
    int multipleScatteringSamples = 50;

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

    Plane plane1 = Plane({0, -0.52, 0}, {0, 1, 0}, {1, 1, 1}, 1);
    plane1.sigmaS = {2.17, 2.62, 3.00};
    plane1.sigmaA = {0.0021, 0.0041, 0.0071};
    plane1.albedo = 0.5;

    Plane plane2 = Plane({-2, 0, 0}, {1, 0, 0}, {1, 1, 1}, 1);
    plane2.sigmaS = {2.18, 2.62, 3.00};
    plane2.sigmaA = {0.0021, 0.0041, 0.0071};
    plane2.albedo = 0.5;

    Camera camera = {{1.5, 0.3, 0}, 1, 2, 2, {0, 0, 1}, {0, 1, 0}, {1, 0, 0}};
    Light light = {{0, 1.5, 2}, {1, 1, 1}, 1.2};

    Scene scene = Scene(light, false);
    scene.addObject(bunny);
    scene.addObject(plane1);
    //scene.addObject(plane2);
    
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
            
            Line ray = {camera.focalPoint, dir};

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