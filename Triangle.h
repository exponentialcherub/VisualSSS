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

        //Triangle(Plane thePlane, Vector3f p1, Vector3f p2, Vector3f p3);

        void setValues(Plane thePlane, Vector3f p1, Vector3f p2, Vector3f p3)
        {
            plane = thePlane;
            point1 = p1;
            point2 = p2;
            point3 = p3;
        }
        /*
        float getArea()
        {
            Vector3f v1 = point1 - point2;
            Vector3f v2 = point2 - point3;
            Vector3f v3 = point3 - point1;
            float a = sqrt(pow(v1[0], 2) + pow(v1[1], 2) + pow(v1[2], 2));
            float b = sqrt(pow(v2[0], 2) + pow(v2[1], 2) + pow(v2[2], 2));
            float c = sqrt(pow(v3[0], 2) + pow(v3[1], 2) + pow(v3[2], 2));

            float s = (a + b + c) / 2;
            return sqrt(s * (s - a) * (s - b) * (s - c));
        }*/
};