#pragma once
/**
 * Currently not in use, using a work around to get area in scene.cpp instead, just returning it directly..
 **/
#include <iostream>

class Face
{
    public:
        virtual float getArea() = 0;
};