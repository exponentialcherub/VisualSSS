#include "Scene.h"
#include "Triangle.h"

Scene::Scene(Light theLight)
{
    light = theLight;
}

void Scene::addObject(Object &object)
{
    objects.push_back(&object);
}

void Scene::rayTrace(Line ray, float *pixel, int singleScatteringSamples, int multipleScatteringSamples)
{
    float t = numeric_limits<float>::max();
    int closest = -1;
    Vector3f normal;
    Vector3f intersectionPoint;
    float diffAngle = 0;
    float specAngle = 0;
    float refractionIndex = 1.5;
    float scale = 5;

    for(int k=0; k<objects.size(); k++)
    {
        float tNew = 0;
        Vector3f newNormal;
        Vector3f newIntersectionPoint;

        if(objects[k]->intersects(ray, tNew, newNormal, newIntersectionPoint))
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

    float s0 = distanceTwoPoints(intersectionPoint, ray.origin);//t*scale;

    if(closest == -1)
    {
        pixel[0] = 0.529;
        pixel[1] = 0.808;
        pixel[2] = 0.922;
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

    // Single scattering term
    Vector3f singleScatteringContribution = {0, 0, 0};
    if(objects[closest]->isTranslucent())
    {
        Vector3f sigmaTx0 = objects[closest]->sigmaT;
        Vector3f sigmaS = objects[closest]->sigmaS;

        float theta = acosf(normal.dot(-ray.direction));
        float fresnelTrans0 = FresnelTransmission(refractionIndex, theta);

        // Find max distance random sample can go in object
        Line rayTest = {intersectionPoint, ray.direction};
        float maxT;
        objects[closest]->intersects(ray, maxT);

        for(int i=0; i<singleScatteringSamples; i++)
        {
            // Find random point for sample
            float t = numeric_limits<float>::max();
            while(t > maxT || t < 0.0001){
                // TODO: Find out why this makes it better.
                float r = ((float) rand() / (RAND_MAX))* (1 - exp(-maxT * sigmaTx0.mean())) + exp(-maxT * sigmaTx0.mean());
                t = -log(r)/sigmaTx0.mean();
            }

            Vector3f point = rayTest.origin + t * ray.direction;

            // TODO: Currently only supports a single light.
            // Find light intersection with object
            // TODO: MAYBE WRONG, LOOK AT OTHER CODE
            Vector3f dir = light.origin - point;
            dir.normalize();
            Line lightRay = {point, dir};
            float lightT;
            Vector3f normali;
            Vector3f c;
            objects[closest]->intersects(lightRay, lightT, normali, c);
            float si = distanceTwoPoints(point, c);

            // Find distance from light to object intersection.
            float distance = distanceTwoPoints(c, light.origin);
            //float distance = ((light.origin[0] - lightRay.origin[0]) / lightRay.direction[0]) * scale;
            Vector3f outgoingRadiance = (light.emittedLight * 0.04) / (distance * distance);

            // TODO: Figure out what is right here.
            float wiDotNi = fabs(lightRay.direction.dot(normali));
            float wiDotN = fabs(lightRay.direction.dot(normali));
            float refractedDistance = (si * wiDotNi) / 
                                      sqrt(1 - (1/(refractionIndex*refractionIndex) * (1- wiDotN*wiDotN)));
            float fresnelTrans1 = FresnelTransmission(refractionIndex, acos(wiDotNi));

            // Set equal as same material.
            Vector3f sigmaTxi = sigmaTx0;
            Vector3f sigmaTc = sigmaTx0 + (fabs(normali.dot(ray.direction)) / wiDotNi) * sigmaTxi;

            singleScatteringContribution[0] += ((sigmaS[0] * fresnelTrans0 * fresnelTrans1 * (1/sigmaTc[0])) 
                                           * expf(-refractedDistance*sigmaTxi[0]) * expf(-s0 * sigmaTx0[0]) * outgoingRadiance[0])
                                            / singleScatteringSamples;
            singleScatteringContribution[1] += ((sigmaS[1] * fresnelTrans0 * fresnelTrans1 * (1/sigmaTc[1])) 
                                           * expf(-refractedDistance*sigmaTxi[1]) * expf(-s0 * sigmaTx0[1]) * outgoingRadiance[1])
                                            / singleScatteringSamples;
            singleScatteringContribution[2] += ((sigmaS[2] * fresnelTrans0 * fresnelTrans1 * (1/sigmaTc[2])) 
                                           * expf(-refractedDistance*sigmaTxi[2]) * expf(-s0 * sigmaTx0[2]) * outgoingRadiance[2])
                                            / singleScatteringSamples;
                                            singleScatteringContribution *=5;
            //cout << "Cont: " << singleScatteringContribution << endl;
        }
    }

    if(!shadow && !objects[closest]->isTranslucent())
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
    /*pixel[0] = objects[closest]->colour[0] * 0.2 + diffAngle*objects[closest]->diffuse[0]*light.localI + pow(specAngle, 20)*objects[closest]->specular[0];
    pixel[1] = objects[closest]->colour[1] * 0.2 + diffAngle*objects[closest]->diffuse[1]*light.localI + pow(specAngle, 20)*objects[closest]->specular[1];
    pixel[2] = objects[closest]->colour[2] * 0.2 + diffAngle*objects[closest]->diffuse[2]*light.localI + pow(specAngle, 20)*objects[closest]->specular[2];*/

    pixel[0] = objects[closest]->colour[0]*0.2 + diffAngle*objects[closest]->diffuse[0]*light.emittedLight[0] + 
               pow(specAngle, 20)*objects[closest]->specular[0] + singleScatteringContribution[0];
    pixel[1] = objects[closest]->colour[1]*0.2 + diffAngle*objects[closest]->diffuse[1]*light.emittedLight[1] + 
               pow(specAngle, 20)*objects[closest]->specular[1] + singleScatteringContribution[1];
    pixel[2] = objects[closest]->colour[2]*0.2 + diffAngle*objects[closest]->diffuse[2]*light.emittedLight[2] + 
               pow(specAngle, 20)*objects[closest]->specular[2] + singleScatteringContribution[2];
}

float Scene::FresnelTransmission(float n, float theta)
{
    return 1 - FresnelReflectance(n, theta);
}

float Scene::FresnelReflectance(float n, float theta)
{
    // ALSO CHANGE THIS COPIED LOLOLOLOL
    double cosi = cos(theta);
	double cost = sqrt(1 - n * n * sin(theta) * sin(theta));
	double Rs = (n * cosi - cost)/(n * cosi + cost);
	Rs *= Rs;
	double Rt = (n * cost - cosi)/(n * cost + cosi);
	Rt *= Rt;
	if(isnan(Rs)){
		Rs = 0;
	}
	if(isnan(Rt)){
		Rt = 0;
	}
    // Should this be averaged?
	double R = (Rt + Rs)/2;
	return R;
}

float Scene::distanceTwoPoints(Vector3f point1, Vector3f point2)
{
    return sqrt(pow(point1[0] - point2[0], 2) + pow(point1[1] - point2[1], 2) + pow(point1[2] - point2[2], 2));
}