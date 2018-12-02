#include "Object.h"
#include "Plane.h"
#include "Triangle.h"
#include <Eigen/Core>
#include <vector>
#include <limits>
#include <cmath>

using namespace Eigen;
using namespace std;

class Mesh : public Object
{
    vector<Triangle> triangles;
    int noTriangles;

    public:
        Mesh(Eigen::MatrixXf v, Eigen::MatrixXi f, Vector3f c);
        
        bool intersects(Line ray, float &t);
        bool intersects(Line ray, Light light, float &t, float &ambAngle, float &specAngle, Vector3f &intersectionPoint) override;

    private:
        void calculateTriangles(Eigen::MatrixXf vertices, Eigen::MatrixXi faces);
};
