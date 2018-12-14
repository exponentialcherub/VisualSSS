#include "Object.h"
#include "Light.h"
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <random>

using namespace std;

/**
 * Defines what is contained in the scene, objects, light, etc. And provides functionality to trace a ray through the scene.
 **/
class Scene
{
    vector<Object*> objects;
    Light light;
    Vector3f background;
    bool shadow;

    float FresnelTransmission(float n, float theta);
    float FresnelReflectance(float n, float theta);
    float distanceTwoPoints(Vector3f point1, Vector3f point2);
    bool shadowCheck(Line ray);

    public:
        Scene(Light theLight, bool shadows);
        void rayTrace(Line ray, float *pixel, int singleScatteringSamples, int multipleScatteringSamples);
        void addObject(Object &object);
};