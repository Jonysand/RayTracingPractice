//
// Created by goksu on 2/25/20.
//
#include "Scene.hpp"
#include <omp.h>

#pragma once
struct hit_payload
{
    float tNear;
    uint32_t index;
    Vector2f uv;
    Object *hit_obj;
};

class Renderer
{
public:
    void Render(const Scene &scene);
    void Render(const Scene &scene, bool MultiThread);
private:
    omp_lock_t lock;
    
    
};

inline void shadePixel(int spp, std::vector<Vector3f> framebuffer, int m, const Scene &scene, Vector3f eye_pos, Vector3f dir)
{
    for (int k = 0; k < spp; k++)
    {
        framebuffer[m] += scene.castRay(Ray(eye_pos, dir), 0) / spp;
    }
}