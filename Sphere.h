#include "Object.h"
#include <Eigen/Core>

using namespace Eigen;

/**
 * Defines a sphere with a radius and origin.
 **/
class Sphere : public Object
{
    float radius;
    Vector3f origin;

    public:
        Sphere(float r, Vector3f o, Vector3f c);

        bool intersects(Line ray, float &t);
        bool intersects(Line ray, float &t, Vector3f &normal, Vector3f &intersectionPoint) override;
        bool isTranslucent();
        Vector3f randomPoint(Vector3f &normal);
        float getBoundingBoxIntersect(Line ray) override;
};