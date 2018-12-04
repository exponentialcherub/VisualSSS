#include "Object.h"
#include "Light.h"
#include <vector>

using namespace std;

class Scene
{
    vector<Object*> objects;
    Light light;

    public:
        Scene(Light theLight);
        void rayTrace(Line ray, float *pixel);
        void addObject(Object &object);
};