#include "Scene.h"
#include "Triangle.h"

Scene::Scene(Light theLight, bool s)
{
    light = theLight;
    shadow = s;

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
    float refractionIndex = 1.5;
    float singleScatterWeight = 1;
    float multipleScatterWeight = 1;

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

    // Single scattering term
    Vector3f singleScatteringContribution = {0, 0, 0};
    if(objects[closest]->isTranslucent())
    {
        Vector3f sigmaTx0 = objects[closest]->getSigmaT();
        float sigmaTx0Mean = (sigmaTx0[0] + sigmaTx0[1] + sigmaTx0[2])/3;
        Vector3f sigmaS = objects[closest]->sigmaS;

        float theta = acosf(normal.dot(-ray.direction));
        float fresnelTrans0 = FresnelTransmission(refractionIndex, theta);
        if(fresnelTrans0 < 0 )
        {
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
            int face;
            if(!objects[closest]->intersects(lightRay, lightT, normali, c))
            {
               c = point;
               si = 0;
               normali = normal;
               cout << "WARNING: Problem calculating single scattering term, sample skipped." << endl;

               continue;
            }
            lightRay.origin = c + 0.0001 * lightRay.direction;
            if(shadowCheck(lightRay))
            {
                continue;
            }

            si = distanceTwoPoints(point, c);

            float wiDotNi = fabs((lightRay.direction).dot(normali));
            
            // Find distance from light to object intersection.
            float distance = distanceTwoPoints(c, lightPoint);
            Vector3f radiance = (light.emittedLight * light.getArea()) / (distance * distance);
            for(int k = 0; k < 3; k++)
            {
                if(radiance[k] < 0)
                    radiance[k] = 0;
            }
           
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

        // Fresnel reflectance at x0
        float theta = acosf(normal.dot(-ray.direction));
        float fresnelTrans0 = FresnelTransmission(refractionIndex, theta);
        if(fresnelTrans0 < 0 )
        {
            fresnelTrans0 = 1;
        }

        float totalWeight = 0;
        for(int i=0; i<multipleScatteringSamples; i++)
        {
            Vector3f lightPoint = light.randomPoint();
            Vector3f normali;
            Vector3f rndPoint = objects[closest]->randomPoint(normali);

            // Shadow check
            Line lightRay = {rndPoint + 0.0001 * (lightPoint - rndPoint), (lightPoint - rndPoint)};
            if(shadowCheck(lightRay))
            {
                continue;
            }

            float weight = sigmaTRMean * expf(-distanceTwoPoints(intersectionPoint, rndPoint)*sigmaTRMean);
            totalWeight += weight;
            
            // Dipole placement
            Vector3f realSourcePoint = rndPoint - normali * (1 / reducedSigmaTMean);

            // Fresnel approximation
            float Fdr = -(1.44 * (1/pow(refractionIndex, 2))) + (0.710 / refractionIndex) + 0.668 + (0.0636 * refractionIndex);
            float rDistance = distanceTwoPoints(realSourcePoint, rndPoint);
            float virtualDistance = (1 / reducedSigmaTMean) + 4 * ((1 + Fdr)/(1 - Fdr)) / (3 * reducedSigmaTMean); 
            Vector3f virtualSourcePoint =  rndPoint + normali * virtualDistance;
            float vDistance = distanceTwoPoints(virtualSourcePoint, rndPoint);

            // Distance between intersection point and dipole sources.
            float dr = distanceTwoPoints(intersectionPoint, realSourcePoint);
            float dv = distanceTwoPoints(intersectionPoint, virtualSourcePoint);

            // Calculate light intensity.
            float lightDistance = distanceTwoPoints(lightPoint, rndPoint);
            Vector3f radiance = (light.emittedLight * light.getArea()) / (lightDistance * lightDistance);
            for(int k = 0; k < 3; k++)
            {
                if(radiance[k] < 0)
                    radiance[k] = 0;
            }

            // Fresnel reflectance at xi
            float fresnelTrans1 = FresnelTransmission(refractionIndex, acos(fabs((lightPoint - rndPoint).dot(normali))));

            float albedo = objects[closest]->albedo;
            Vector3f diffuse;
            diffuse[0] = radiance[0] * (albedo / (4 * EIGEN_PI)) * 
                         (rDistance * ((sigmaTR[0] * dr + 1) * (expf(-sigmaTR[0]*dr) / (pow(dr, 3)))) -
                            vDistance * (sigmaTR[0]*dv + 1) * (expf(-sigmaTR[0]*dv) / (pow(dv, 3))));
            diffuse[1] = radiance[1] * (albedo / (4 * EIGEN_PI)) * 
                         (rDistance * ((sigmaTR[1] * dr + 1) * (expf(-sigmaTR[1]*dr) / (pow(dr, 3)))) -
                            vDistance * (sigmaTR[1]*dv + 1) * (expf(-sigmaTR[1]*dv) / (pow(dv, 3))));
            diffuse[2] = radiance[2] * (albedo / (4 * EIGEN_PI)) * 
                         (rDistance * ((sigmaTR[2] * dr + 1) * (expf(-sigmaTR[2]*dr) / (pow(dr, 3)))) -
                            vDistance * (sigmaTR[2]*dv + 1) * (expf(-sigmaTR[2]*dv) / (pow(dv, 3))));
            multipleScatteringContribution += (weight * diffuse * fresnelTrans1);
        }
        multipleScatteringContribution *= multipleScatterWeight*fresnelTrans0 / totalWeight;
    }

    pixel[0] = singleScatteringContribution[0] + multipleScatteringContribution[0];
    pixel[1] = singleScatteringContribution[1] + multipleScatteringContribution[1];
    pixel[2] = singleScatteringContribution[2] + multipleScatteringContribution[2];

    for(int i = 0; i < 3; i++)
    {
        if(pixel[i] > 1)
        {
            pixel[i] = 1;
        }
        if(pixel[i] < 0)
        {
            pixel[i] = 0;
        }
    }
}

float Scene::FresnelTransmission(float n, float theta)
{
    return 1 - FresnelReflectance(n, theta);
}

float Scene::FresnelReflectance(float n, float theta)
{
    float cosi = cos(theta);
    float sint = 1 / n * sqrtf(1 - cosi * cosi); 
    if(sint > 1)
    {
        return 1;
    }
    float cost = sqrtf(max(0.f, 1 - sint * sint));
	float Rs = (n * cosi - cost)/(n * cosi + cost);
	Rs *= Rs;
	float Rt = (n * cost - cosi)/(n * cost + cosi);
	Rt *= Rt;
	if(isnan(Rs)){
		Rs = 0;
	}
	if(isnan(Rt)){
		Rt = 0;
	}
	return (Rt + Rs)/2;
}

float Scene::distanceTwoPoints(Vector3f point1, Vector3f point2)
{
    return sqrt(pow(point1[0] - point2[0], 2) + pow(point1[1] - point2[1], 2) + pow(point1[2] - point2[2], 2));
}

bool Scene::shadowCheck(Line ray)
{
    if(!shadow)
    {
        return false;
    }

    bool shadowTrue = false;
    for(int h=0; h<objects.size(); h++)
    {   
        float shadowT = -1;

        if(objects[h]->intersects(ray, shadowT) && shadowT > 0)
        {
            shadowTrue = true;
            break;
        }
    }
    return shadowTrue;
}