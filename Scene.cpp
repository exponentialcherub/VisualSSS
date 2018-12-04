#include "Scene.h"
#include <iostream>
using namespace std;

Scene::Scene(Light theLight)
{
    light = theLight;
}

void Scene::addObject(Object &object)
{
    objects.push_back(&object);
}

void Scene::rayTrace(Line ray, float *pixel)
{
    float t = numeric_limits<float>::max();
    int closest = -1;
    Vector3f normal;
    Vector3f intersectionPoint;
    float diffAngle = 0;
    float specAngle = 0;

    for(int k=0; k<objects.size(); k++)
    {
        float tNew = 0;
        Vector3f newNormal;
        Vector3f newIntersectionPoint;

        if(objects[k]->intersects(ray, light, tNew, newNormal, newIntersectionPoint))
        {
            // Check t > 0 to be sure.
            if(tNew < t && tNew > 0)
            {
                t = tNew;
                normal = newNormal;
                intersectionPoint = newIntersectionPoint;
                closest = k;
            }
        }
    }
    if(closest == -1)
    {
        return;
    }

    bool shadow = false;

    for(int h=0; h<objects.size(); h++)
    {
        if(closest == h)
        {
            continue;
        }
        
        Line lightRay = {intersectionPoint + 0.01*(-light.dir), -light.dir};
        
        float shadowT = -1;

        if(objects[h]->intersects(lightRay, shadowT) && shadowT > 0)
        {
            shadow = true;
            break;
        }
    }

    if(!shadow)
    {
        diffAngle = normal.dot(-light.dir);
        Vector3f reflection = light.dir - 2 * (light.dir.dot(normal)) * normal;
        specAngle = reflection.dot(normal);
    }
    if(diffAngle < 0)
    {
        diffAngle = 0;
    }
    if(specAngle < 0)
    {
        specAngle = 0;
    }
    
    // Calculate light here I think, getting normal from intersection method
    pixel[0] = objects[closest]->colour[0] * light.ambientI + diffAngle*objects[closest]->diffuse[0]*light.localI + pow(specAngle, 20)*objects[closest]->specular[0];
    pixel[1] = objects[closest]->colour[1] * light.ambientI + diffAngle*objects[closest]->diffuse[1]*light.localI + pow(specAngle, 20)*objects[closest]->specular[1];
    pixel[2] = objects[closest]->colour[2] * light.ambientI + diffAngle*objects[closest]->diffuse[2]*light.localI + pow(specAngle, 20)*objects[closest]->specular[2];
}