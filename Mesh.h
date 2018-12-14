#include "Face.h"
#include "Object.h"
#include "Plane.h"
#include "Triangle.h"
#include "BoundingBox.h"
#include <Eigen/Core>
#include <vector>
#include <limits>
#include <math.h>
#include <stdlib.h> 

using namespace Eigen;
using namespace std;

class Mesh : public Object
{
    vector<Triangle> triangles;
    int noTriangles;
    BoundingBox boundingBox;

    bool translucent;

    public:
        Mesh(Eigen::MatrixXf v, Eigen::MatrixXi f, Eigen::MatrixXf vn, Eigen::MatrixXi fn, Vector3f c, bool isTranslucent);
        
        bool intersects(Line ray, float &t);
        bool intersects(Line ray, float &t, Vector3f &normal, Vector3f &intersectionPoint, int &face) override;
        float getBoundingBoxIntersect(Line ray) override;
        bool isTranslucent();
        Vector3f randomPoint(Vector3f &normal, int face = -1);

    private:
        void calculateTriangles(Eigen::MatrixXf vertices, Eigen::MatrixXi faces, Eigen::MatrixXf vn, Eigen::MatrixXi fn);
        float getTriangleArea(Vector3f point1, Vector3f point2, Vector3f point3);
};
