#include "Face.h"
#include "Object.h"
#include "Light.h"
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <random>

using namespace std;

class Scene
{
    vector<Object*> objects;
    Light light;

    float FresnelTransmission(float n, float theta);
    float FresnelReflectance(float n, float theta);
    float distanceTwoPoints(Vector3f point1, Vector3f point2);

    public:
        Scene(Light theLight);
        void rayTrace(Line ray, float *pixel, int singleScatteringSamples, int multipleScatteringSamples);
        void addObject(Object &object);
};