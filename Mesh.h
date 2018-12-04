#include "Object.h"
#include "Plane.h"
#include "Triangle.h"
#include "BoundingBox.h"
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
    BoundingBox boundingBox;

    bool translucent;

    public:
        Mesh(Eigen::MatrixXf v, Eigen::MatrixXi f, Vector3f c, bool isTranslucent);
        
        bool intersects(Line ray, float &t);
        bool intersects(Line ray, Light light, float &t, Vector3f &normal, Vector3f &intersectionPoint) override;
        bool isTranslucent();

    private:
        void calculateTriangles(Eigen::MatrixXf vertices, Eigen::MatrixXi faces);
};
