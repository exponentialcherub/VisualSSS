#include "Object.h"
#include <Eigen/Core>

using namespace Eigen;

class Sphere : public Object
{
    float radius;
    Vector3f origin;

    public:
        Sphere(float r, Vector3f o, Vector3f c);

        bool intersects(Line ray, float &t);

        bool intersects(Line ray, Light light, float &t, float &ambAngle, float &specAngle, Vector3f &intersectionPoint) override;
};