#include <vector>
#include <limits>
#include "BoundingBox.h"
#include "Triangle.h"
using namespace std;
class KdNode
{
    public:
     BoundingBox boundingBox;
     KdNode* left;
     KdNode* right;
     std::vector<Triangle*> triangles;


     KdNode* build(vector<Triangle*>& tris, int depth) const;
     bool hit(KdNode* node, const Line& ray, float& t, float& tmin, Vector3f &intersectionPoint, Vector3f &normal);
};