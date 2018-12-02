#include <Eigen/Dense>

using namespace Eigen;

class Triangle
{
    public:
        Plane plane;
        Vector3f point1;
        Vector3f point2;
        Vector3f point3;
};