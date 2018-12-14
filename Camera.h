#include <Eigen/Dense>

using namespace Eigen;

/**
 * The view of the image, defined by directional vectors and position.
 **/
class Camera
{
    public:
        Vector3f focalPoint;
        float focalLength;
        float height;
        float width;
        Vector3f upVector;
        Vector3f rightVector;
        Vector3f forwardVector;
};