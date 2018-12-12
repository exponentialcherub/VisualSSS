#include "Plane.h"
#include "Face.h"
#include <Eigen/Dense>
#include <math.h>

using namespace Eigen;

class Triangle
{
    public:
        Plane plane;
        Vector3f point1;
        Vector3f point2;
        Vector3f point3;
        Vector3f normal1;
        Vector3f normal2;
        Vector3f normal3;

        void setValues(Plane thePlane, Vector3f p1, Vector3f p2, Vector3f p3, Vector3f n1, Vector3f n2, Vector3f n3)
        {
            plane = thePlane;
            point1 = p1;
            point2 = p2;
            point3 = p3;
        }
};