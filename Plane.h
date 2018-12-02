#include "Object.h"

class Plane : public Object
{
    public:

        Vector3f point;
        Vector3f normal;

        Plane(Vector3f p, Vector3f n, Vector3f c);

        bool intersects(Line ray, float &t);
        bool intersects(Line ray, Light light, float &t, float &ambAngle, float &specAngle, Vector3f &intersectionPoint);
};