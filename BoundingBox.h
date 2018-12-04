#include "Line.h"
#include <Eigen/Core>

using namespace Eigen;

class BoundingBox
{
    public:
        BoundingBox(){}

        bool intersects(Line ray);
        void setValues(Vector3f min, Vector3f max);

    private:
        Vector3f min;
        Vector3f max;
};