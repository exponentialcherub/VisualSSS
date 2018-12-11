#include "Scene.h"
#include "Triangle.h"

Scene::Scene(Light theLight)
{
    light = theLight;

    background = {0, 0, 0};
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
    // TODO: Remove this at end, use BSSRDF or plane
    float diffAngle = 0;
    float specAngle = 0;
    float refractionIndex = 1.5;
    float singleScatterWeight = 1;

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

    if(closest == -1)
    {
        pixel[0] = background[0];
        pixel[1] = background[1];
        pixel[2] = background[2];
        return;
    }

    bool shadow = false;

    for(int h=0; h<objects.size(); h++)
    {
        if(closest == h)
        {
            continue;
        }
        
        /*Line lightRay = {intersectionPoint + 0.01*(-light.dir), -light.dir};
        
        float shadowT = -1;

        if(objects[h]->intersects(lightRay, shadowT) && shadowT > 0)
        {
            shadow = true;
            break;
        }*/
    }

    // Single scattering term
    Vector3f singleScatteringContribution = {0, 0, 0};
    if(objects[closest]->isTranslucent())
    {
        Vector3f sigmaTx0 = objects[closest]->getSigmaT();
        float sigmaTx0Mean = (sigmaTx0[0] + sigmaTx0[1] + sigmaTx0[2])/3;
        Vector3f sigmaS = objects[closest]->sigmaS;

        float theta = acosf(normal.dot(-ray.direction));
        float fresnelTrans0 = FresnelTransmission(1/refractionIndex, theta);
        if(fresnelTrans0 < 0 )
        {
            // Likely an issue with 3D model, so just assume 1.
            fresnelTrans0 = 1;
        }

        // Find max distance random sample can go in object
        Line rayTest = {intersectionPoint + 0.001 * ray.direction, ray.direction};
        float maxT;
        Vector3f maxNormal;
        Vector3f maxPoint;
        if(!objects[closest]->intersects(rayTest, maxT, maxNormal, maxPoint))
        {
            maxT = objects[closest]->getBoundingBoxIntersect(ray);
        }

        for(int i=0; i<singleScatteringSamples; i++)
        {
            // Find random point for sample
            float t = numeric_limits<float>::max();
            while(t > maxT || t < 0.0001)
            {
                if(maxT < 0.0001)
                {
                    t = maxT;
                    break;
                }
                // TODO: Find out why this makes it better. See if it's fine without this. Mean just giving first value???
                float r = ((float) rand() / (RAND_MAX)) * (1 - exp(-maxT * sigmaTx0.mean())) + exp(-maxT * sigmaTx0.mean());
                t = -log(r)/sigmaTx0.mean();
                //cout << "maxT" << maxT << "t" << t << endl;
            }

            Vector3f point = rayTest.origin + t * rayTest.direction;
            float s0 = distanceTwoPoints(intersectionPoint, point);

            // Find light intersection with object
            Vector3f lightPoint = light.randomPoint();
            Vector3f dir = lightPoint - point;
            dir.normalize();
            Line lightRay = {point, dir};
            float lightT;
            Vector3f normali;
            Vector3f c;
            float si;
            if(!objects[closest]->intersects(lightRay, lightT, normali, c))
            {
                //TODO: is this necessary, can I do this old way?
               c = point;
               si = 0;
               normali = maxNormal;
               //singleScatteringContribution = {0, 0, 1};
               //singleScatteringSamples = 1;
               //continue;
            }

            si = distanceTwoPoints(point, c);

            float wiDotNi = fabs((lightRay.direction).dot(normali));
            
            // Find distance from light to object intersection.
            float distance = distanceTwoPoints(c, lightPoint);
            Vector3f radiance = (light.emittedLight * light.getArea()) * wiDotNi / (distance * distance);
            for(int k = 0; k < 3; k++)
            {
                if(radiance[k] < 0)
                    radiance[k] = 0;
            }
            
            Vector3f averageNormal = (normali + normal) / 2;
            // TODO: Check if this is right, compare at higher res.
            float wiDotN = fabs((lightRay.direction).dot(normal));
            float refractedDistance = (si * wiDotNi) / 
                                      sqrt(1 - (1/(refractionIndex*refractionIndex) * (1- wiDotNi*wiDotNi)));
            float fresnelTrans1 = FresnelTransmission(refractionIndex, acos(wiDotNi));

            // Set equal as same object material.
            Vector3f sigmaTxi = sigmaTx0;
            Vector3f sigmaTc = sigmaTx0 + (fabs(normali.dot(ray.direction)) / wiDotNi) * sigmaTxi;
            
            singleScatteringContribution[0] += ((sigmaS[0] * fresnelTrans0 * fresnelTrans1 * (1/sigmaTc[0])) 
                                           * expf(-refractedDistance*sigmaTxi[0]) * expf(-s0 * sigmaTx0[0]) * radiance[0]);
            singleScatteringContribution[1] += ((sigmaS[0] * fresnelTrans0 * fresnelTrans1 * (1/sigmaTc[1])) 
                                           * expf(-refractedDistance*sigmaTxi[1]) * expf(-s0 * sigmaTx0[1]) * radiance[1]);
            singleScatteringContribution[2] += ((sigmaS[0] * fresnelTrans0 * fresnelTrans1 * (1/sigmaTc[2])) 
                                           * expf(-refractedDistance*sigmaTxi[2]) * expf(-s0 * sigmaTx0[2]) * radiance[2]);
            if(singleScatteringContribution[0] < 0.001)
        {
         cout << "1: " << maxNormal << endl;   
         cout << "2: " << t << endl;   
         cout << "3: " << maxT << endl;   
         cout << "4: " << refractedDistance << endl; 
         cout << "5: " << si << endl;  
        }
        }

        singleScatteringContribution =  singleScatterWeight * singleScatteringContribution / singleScatteringSamples;
    }

    // Multiple scattering
    Vector3f multipleScatteringContribution = {0, 0, 0};
    if(objects[closest]->isTranslucent())
    {
        Vector3f sigmaTx0 = objects[closest]->getSigmaT();
        Vector3f sigmaTR = objects[closest]->getSigmaTR();
        float sigmaTRMean = (sigmaTR[0] + sigmaTR[1] + sigmaTR[2]) / 3;
        Vector3f reducedSigmaT = objects[closest]->getReducedSigmaT();
        float reducedSigmaTMean = (reducedSigmaT[0] + reducedSigmaT[1] + reducedSigmaT[2]) / 3;

        // Dipole placement
        Vector3f realSourcePoint = intersectionPoint + normal * (1 / reducedSigmaTMean);

        // Fresnel approximation
        float Fdr = (1.44 * (1/pow(refractionIndex, 2))) + (0.710 / refractionIndex) + 0.668 + (0.0636 * refractionIndex);
        float rDistance = distanceTwoPoints(realSourcePoint, intersectionPoint);

        float virtualDistance = (1 / reducedSigmaTMean) + 4 * ((1 + Fdr)/(1 - Fdr)) / (3 * reducedSigmaTMean); 
        Vector3f virtualSourcePoint =  intersectionPoint - normal * virtualDistance;
        float vDistance = distanceTwoPoints(virtualSourcePoint, intersectionPoint);

        // Fresnel reflectance at x0
        float theta = acosf(normal.dot(-ray.direction));
        float fresnelTrans0 = FresnelTransmission(1/refractionIndex, theta);
        if(fresnelTrans0 < 0 )
        {
            // Likely an issue with 3D model, so just assume 1.
            fresnelTrans0 = 1;
        }

        float totalWeight = 0;
        for(int i=0; i<multipleScatteringSamples; i++)
        {
            Vector3f lightPoint = light.randomPoint();
            Vector3f normali;
            Vector3f rndPoint = objects[closest]->randomPoint(normali);
            //float weight = sigmaTRMean * expf(-distanceTwoPoints(intersectionPoint, rndPoint)*sigmaTRMean);
            float weight  = 1 - exp(-distanceTwoPoints(intersectionPoint, rndPoint) / sigmaTRMean);
            totalWeight += weight;
            // Distance between sampled point and dipole sources.
            float dr = distanceTwoPoints(rndPoint, realSourcePoint);
            float dv = distanceTwoPoints(rndPoint, virtualSourcePoint);

            // Calculate light intensity.
            float lightDistance = distanceTwoPoints(lightPoint, rndPoint);
            Vector3f radiance = (light.emittedLight * light.getArea()) * (lightPoint - rndPoint).dot(normali) / (lightDistance * lightDistance);
            for(int k = 0; k < 3; k++)
            {
                if(radiance[k] < 0)
                    radiance[k] = 0;
            }

            // Fresnel reflectance at xi
            float fresnelTrans1 = FresnelTransmission(refractionIndex, acos(fabs((lightPoint - rndPoint).dot(normali))));

            // TODO: Store this in object.
            float albedo = 0.5;
            // TODO: fresnel reflefctance
            Vector3f diffuse;
            diffuse[0] = radiance[0] * (albedo / (4 * EIGEN_PI)) * 
                         (rDistance * ((sigmaTR[0] * dr + 1) * (expf(-sigmaTR[0]*dr) / (reducedSigmaT[0]*pow(dr, 3)))) - 
                            vDistance * (sigmaTR[0]*dv + 1) * (expf(-sigmaTR[0]*dv) / (reducedSigmaT[0]*pow(dv, 3))));
            diffuse[1] = radiance[1] * (albedo / (4 * EIGEN_PI)) * 
                         (rDistance * ((sigmaTR[1] * dr + 1) * (expf(-sigmaTR[1]*dr) / (reducedSigmaT[1]*pow(dr, 3)))) - 
                            vDistance * (sigmaTR[1]*dv + 1) * (expf(-sigmaTR[1]*dv) / (reducedSigmaT[1]*pow(dv, 3))));
            diffuse[2] = radiance[2] * (albedo / (4 * EIGEN_PI)) * 
                         (rDistance * ((sigmaTR[2] * dr + 1) * (expf(-sigmaTR[2]*dr) / (reducedSigmaT[2]*pow(dr, 3)))) -
                            vDistance * (sigmaTR[2]*dv + 1) * (expf(-sigmaTR[2]*dv) / (reducedSigmaT[2]*pow(dv, 3))));
            multipleScatteringContribution += (weight * diffuse * fresnelTrans1);
        }
        multipleScatteringContribution *= 20*fresnelTrans0 / totalWeight;
        //cout << multipleScatteringContribution << endl;
    }

    if(!shadow && !objects[closest]->isTranslucent())
    {
        //diffAngle = normal.dot(-light.dir);
        //Vector3f reflection = light.dir - 2 * (light.dir.dot(normal)) * normal;
        //specAngle = reflection.dot(normal);
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

    pixel[0] = singleScatteringContribution[0] + multipleScatteringContribution[0];
    pixel[1] = singleScatteringContribution[1] + multipleScatteringContribution[1];
    pixel[2] = singleScatteringContribution[2] + multipleScatteringContribution[2];
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