#include "Scene.h"

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
        return;
    }
    
    // Calculate light here I think, getting normal from intersection method
    pixel[0] = objects[closest]->colour[0] * light.ambientI + ambAngle*objects[closest]->diffuse[0]*light.localI + pow(specAngle, 20)*objects[closest]->specular[0];
    pixel[1] = objects[closest]->colour[1] * light.ambientI + ambAngle*objects[closest]->diffuse[1]*light.localI + pow(specAngle, 20)*objects[closest]->specular[1];
    pixel[2] = objects[closest]->colour[2] * light.ambientI + ambAngle*objects[closest]->diffuse[2]*light.localI + pow(specAngle, 20)*objects[closest]->specular[2];
}