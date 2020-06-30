//
// Created by goksu on 2/25/20.
//

#include <fstream>
#include "Scene.hpp"
#include "Renderer.hpp"
#include <thread>

inline float deg2rad(const float &deg) { return deg * M_PI / 180.0; }

const float EPSILON = 0.01;

// The main render function. This where we iterate over all pixels in the image,
// generate primary rays and cast these rays into the scene. The content of the
// framebuffer is saved to a file.
void Renderer::Render(const Scene &scene)
{
    std::vector<Vector3f> framebuffer(scene.width * scene.height);

    float scale = tan(deg2rad(scene.fov * 0.5));
    float imageAspectRatio = scene.width / (float)scene.height;
    Vector3f eye_pos(278, 273, -800);
    int m = 0;

    // change the spp value to change sample ammount
    int spp = 1;
    std::cout << "Single Thread"
              << "\n";
    std::cout << "SPP: " << spp << "\n";
    for (uint32_t j = 0; j < scene.height; ++j)
    {
        for (uint32_t i = 0; i < scene.width; ++i)
        {
            // generate primary ray direction
            float x = (2 * (i + 0.5) / (float)scene.width - 1) *
                      imageAspectRatio * scale;
            float y = (1 - 2 * (j + 0.5) / (float)scene.height) * scale;

            Vector3f dir = normalize(Vector3f(-x, y, 1));
            m = j * scene.width + i;
            for (int k = 0; k < spp; k++)
            {
                framebuffer[m] += scene.castRaySmallTP(Ray(eye_pos, dir), 0) / spp;
            }
        }
        UpdateProgress(m, (float)scene.width, (float)scene.height);
    }
    int finished = 1;
    UpdateProgress(finished, (float)scene.width, (float)scene.height);

    // save framebuffer to file
    FILE *fp = fopen("binary.ppm", "wb");
    (void)fprintf(fp, "P6\n%d %d\n255\n", scene.width, scene.height);
    for (auto i = 0; i < scene.height * scene.width; ++i)
    {
        static unsigned char color[3];
        color[0] = (unsigned char)(255 * std::pow(clamp(0, 1, framebuffer[i].x), 0.6f));
        color[1] = (unsigned char)(255 * std::pow(clamp(0, 1, framebuffer[i].y), 0.6f));
        color[2] = (unsigned char)(255 * std::pow(clamp(0, 1, framebuffer[i].z), 0.6f));
        fwrite(color, 1, 3, fp);
    }
    fclose(fp);
}

/*----------
multi thread
----------*/
void Renderer::Render(const Scene &scene, bool MultiThread)
{
    if (!MultiThread)
    {
        Render(scene);
        return;
    }
    std::vector<Vector3f> framebuffer(scene.width * scene.height);

    float scale = tan(deg2rad(scene.fov * 0.5));
    float imageAspectRatio = scene.width / (float)scene.height;
    Vector3f eye_pos(278, 273, -800);
    int m = 0;

    // change the spp value to change sample ammount
    int spp = 256;
    std::cout << "Multiple Thread"
              << "\n";
    std::cout << "SPP: " << spp << "\n";
    int finished = 0;
#pragma omp parallel for num_threads(8)
    for (uint32_t m = 0; m < scene.height * scene.width; ++m)
    {
        int i = m % scene.width;
        int j = (m - i) / scene.height;
        // generate primary ray direction
        float x = (2 * (i + 0.5) / (float)scene.width - 1) * imageAspectRatio * scale;
        float y = (1 - 2 * (j + 0.5) / (float)scene.height) * scale;
        Vector3f dir = normalize(Vector3f(-x, y, 1));
        // #pragma omp parallel for reduction(+: finished)
        for (int k = 0; k < spp; k++)
        {
            framebuffer[m] += scene.castRay(Ray(eye_pos, dir), 0) / spp;
        }
        ++finished;
        UpdateProgress(finished, scene.height, scene.width);
    }

    // save framebuffer to file
    FILE *fp = fopen("binary.ppm", "wb");
    (void)fprintf(fp, "P6\n%d %d\n255\n", scene.width, scene.height);
    for (auto i = 0; i < scene.height * scene.width; ++i)
    {
        static unsigned char color[3];
        color[0] = (unsigned char)(255 * std::pow(clamp(0, 1, framebuffer[i].x), 0.45f));
        color[1] = (unsigned char)(255 * std::pow(clamp(0, 1, framebuffer[i].y), 0.45f));
        color[2] = (unsigned char)(255 * std::pow(clamp(0, 1, framebuffer[i].z), 0.45f));
        fwrite(color, 1, 3, fp);
    }
    fclose(fp);
    UpdateProgress(finished, scene.height, scene.width);
}