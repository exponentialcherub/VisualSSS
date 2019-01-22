#include "KdNode.h"
#include <iostream>
#include <stdlib.h>

KdNode* KdNode::build(vector<Triangle*>& tris, int depth) const{
    KdNode* node = new KdNode();
    node->triangles = tris;
    node->left = NULL;
    node->right = NULL;
    node->boundingBox = BoundingBox();

    if(tris.size() == 0){
        return node;
    }

    if(tris.size() == 1){
        node->boundingBox = tris[0]->boundingBox;
        node->left = new KdNode();
        node->right = new KdNode();
        node->left->triangles = vector<Triangle*>();
        node->right->triangles = vector<Triangle*>();
        return node;
    }

    node->boundingBox = tris[0]->boundingBox;

    for(int i=1;i<tris.size();i++)
    {
        node->boundingBox.expand(tris[i]->boundingBox);
    }

    Vector3f midpt ={0,0,0};

    for(int i=0;i<tris.size();i++)
    {
        midpt = midpt + (tris[i]->midPoint() * (1.0/tris.size()));
    }
    vector<Triangle*> left_tris;
    vector<Triangle*> right_tris;
    int axis = node->boundingBox.longestAxis();
    for(int i=0;i<tris.size();i++)
    {
        switch(axis){
            case 0:
                midpt[0] >= tris[i]->midPoint()[0] ? right_tris.push_back(tris[i]) : left_tris.push_back(tris[i]);
                break;
            case 1:
                midpt[1] >= tris[i]->midPoint()[1] ? right_tris.push_back(tris[i]) : left_tris.push_back(tris[i]);
                break;
            case 2:
                midpt[2] >= tris[i]->midPoint()[2] ? right_tris.push_back(tris[i]) : left_tris.push_back(tris[i]);
                break;
        }
    }

    if(left_tris.size() == 0 && right_tris.size() > 0) left_tris = right_tris;
    if(right_tris.size() == 0 && left_tris.size() > 0) right_tris = left_tris;

    int matches = 0;
    for(int i=0;i<left_tris.size();i++)
    {
        for(int j=0;j<right_tris.size();j++)
        {
            if(left_tris[i] == right_tris[j])
            {
                matches++;
            }
        }
    }
    
    if((float)matches / left_tris.size() < 0.5 && (float)matches / right_tris.size() < 0.5)
    {
        node->left = build(left_tris, depth+1);
        node->right = build(right_tris, depth+1);
    }
    else
    {
        node->left = new KdNode();
        node->right = new KdNode();
        node->left->triangles = vector<Triangle*>();
        node->right->triangles = vector<Triangle*>();
    }

    return node;
}

bool KdNode::hit(KdNode* node, const Line& ray, float& t, float& tmin, Vector3f &intersectionPoint, Vector3f &normal)
{
    if(node->boundingBox.intersects(ray))
    {
        bool hit_tri = false;
        Vector3f hit_pt, local_hit_pt;

        if(node->left->triangles.size() > 0 || node->right->triangles.size() > 0)
        {
            float t1,t2;
            Vector3f a, b;
            Vector3f n1, n2;
            bool hitLeft = hit(node->left, ray, t1, tmin, a, n1);
            bool hitRight = hit(node->right, ray, t2, tmin, b, n2);
            if(hitLeft)
            {
                t = t1;
                normal = n1;
                intersectionPoint = a;
                if(hitRight && t2 < t1)
                {
                    t = t2;
                    normal = n2;
                    intersectionPoint = b;
                }
            }
            else if(hitRight)
            {
                t = t2;
                normal = n2;
                intersectionPoint = b;
            }
            //if(newT < t && newT > 0)
            //{
                
            //}
            return (hitLeft || hitRight);
        }
        else
        {
            float newT;
            Vector3f a;
            Vector3f newNormal;
            t = 100000;
            for(int i=0;i<node->triangles.size();i++)
            {
                if(node->triangles[i]->intersect(ray, a, newT, newNormal))
                {
                    if(t < newT)
                    {
                        continue;
                    }
                    t = newT;
                    normal = newNormal;
                    intersectionPoint = a;
                    hit_tri = true;
                    tmin = t;
                }
            }
            if(hit_tri)
            {
                return true;
            }
            return false;
        }
    }
}