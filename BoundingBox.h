
#include "Plane.h"
#include "Line.h"
#include <Eigen/Core>
#include <limits>

using namespace Eigen;

/**
 * Box which contains a mesh, used to optmise rendering.
 **/
class BoundingBox
{
    public:
        Vector3f min;
        Vector3f max;

        BoundingBox(){}

        bool intersects(Line ray);
        float getIntersectT(Line ray);
        void setValues(Vector3f min, Vector3f max);
};