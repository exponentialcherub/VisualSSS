#include "Object.h"
#include "Plane.h"
#include "Triangle.h"
#include "BoundingBox.h"
#include "KdNode.h"
#include <Eigen/Core>
#include <vector>
#include <limits>
#include <math.h>
#include <stdlib.h> 

using namespace Eigen;
using namespace std;

/**
 * A 3D model defined by a list of vertices, faces and normals. Includes functions for finding intersections with the mesh and finding
 * random points.
 **/
class Mesh : public Object
{
    vector<Triangle> triangles;
    int noTriangles;
    BoundingBox boundingBox;
    KdNode* kdNode;

    bool translucent;

    public:
        Mesh(Eigen::MatrixXf v, Eigen::MatrixXi f, Eigen::MatrixXf vn, Eigen::MatrixXi fn, Vector3f c, bool isTranslucent);
        
        bool intersects(Line ray, float &t);
        bool intersects(Line ray, float &t, Vector3f &normal, Vector3f &intersectionPoint) override;
        float getBoundingBoxIntersect(Line ray) override;
        bool isTranslucent();
        Vector3f randomPoint(Vector3f &normal);

    private:
        void calculateTriangles(Eigen::MatrixXf vertices, Eigen::MatrixXi faces, Eigen::MatrixXf vn, Eigen::MatrixXi fn);
        float getTriangleArea(Vector3f point1, Vector3f point2, Vector3f point3);
};
