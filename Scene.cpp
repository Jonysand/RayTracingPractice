//
// Created by Göksu Güvendiren on 2019-05-14.
//

#include "Scene.hpp"

inline float GXX_Chi(float a)
{
    return a > 0 ? 1 : 0;
}

inline float GXX_Geometry(Vector3f v, Vector3f m, Vector3f n, float ag)
{
    float thetaV = acosf(dotProduct(v, n));
    return GXX_Chi(dotProduct(v, m) * dotProduct(v, n)) * 2.0 / (1.0 + sqrtf(1.0 + ag * ag * tanf(thetaV) * tanf(thetaV)));
}

void Scene::buildBVH()
{
    printf(" - Generating BVH...\n\n");
    this->bvh = new BVHAccel(objects, 1, BVHAccel::SplitMethod::NAIVE);
}

Intersection Scene::intersect(const Ray &ray) const
{
    return this->bvh->Intersect(ray);
}

void Scene::sampleLight(Intersection &pos, float &pdf) const
{
    float emit_area_sum = 0;
    for (uint32_t k = 0; k < objects.size(); ++k)
    {
        if (objects[k]->hasEmit())
        {
            emit_area_sum += objects[k]->getArea();
        }
    }
    float p = get_random_float() * emit_area_sum;
    emit_area_sum = 0;
    for (uint32_t k = 0; k < objects.size(); ++k)
    {
        if (objects[k]->hasEmit())
        {
            emit_area_sum += objects[k]->getArea();
            if (p <= emit_area_sum)
            {
                objects[k]->Sample(pos, pdf);
                break;
            }
        }
    }
}

bool Scene::trace(
    const Ray &ray,
    const std::vector<Object *> &objects,
    float &tNear, uint32_t &index, Object **hitObject)
{
    *hitObject = nullptr;
    for (uint32_t k = 0; k < objects.size(); ++k)
    {
        float tNearK = kInfinity;
        uint32_t indexK;
        Vector2f uvK;
        if (objects[k]->intersect(ray, tNearK, indexK) && tNearK < tNear)
        {
            *hitObject = objects[k];
            tNear = tNearK;
            index = indexK;
        }
    }

    return (*hitObject != nullptr);
}

// Implementation of Path Tracing
Vector3f Scene::castRay(const Ray &ray, int depth) const
{
    Vector3f L_dir;
    Vector3f L_indir_reflect;
    // TO DO Implement Path Tracing Algorithm here
    Intersection p = intersect(ray);
    if (!p.happened)
        return Vector3f();
    if (p.obj->hasEmit())
    {
        return p.m->getEmission();
    }

    Intersection x;
    float pdf_light = 0.0f;
    sampleLight(x, pdf_light);
    Vector3f wo = -ray.direction;
    Vector3f wi = normalize(x.coords - p.coords);
    Ray testRay(p.coords, wi);
    Intersection testInter = intersect(testRay);
    if (testInter.happened && (testInter.coords - x.coords).norm() < EPSILON)
    {
        L_dir = x.emit *
                p.m->eval(wi, wo, p.normal) *
                dotProduct(wi, p.normal) *
                dotProduct(-wi, x.normal) /
                pow((x.coords - p.coords).norm(), 2) /
                pdf_light;
    }
    float P_RR = get_random_float();
    if (P_RR < RussianRoulette || depth < 5)
    {
        wi = p.m->sample(wo, p.normal);
        Ray nextRay(p.coords, wi);
        Intersection nextInter = intersect(nextRay);
        if (nextInter.happened && !intersect(nextRay).obj->hasEmit())
        {
            float pdf = p.m->pdf(wi, wo, p.normal);
            pdf = pdf > EPSILON ? pdf : EPSILON;
            Vector3f f_r = p.m->eval(wi, wo, p.normal);
            L_indir_reflect = castRay(nextRay, ++depth) *
                              f_r *
                              dotProduct(wi, p.normal) /
                              pdf /
                              RussianRoulette;
        }
    }

    return L_dir + L_indir_reflect;
}

Vector3f Scene::castRaySmallTP(const Ray &ray, int depth) const
{
    Intersection p = intersect(ray); // intersection of this ray in the scene
    if (!p.happened)
        return Vector3f();
    Vector3f x = p.coords;                                       // intersection point coordinates
    Vector3f n = p.normal;                                       // intersection point surface normal
    Vector3f nl = dotProduct(n, ray.direction) < 0 ? n : n * -1; // oriented normal, same side with ray
    Vector3f f = p.m->Kd;                                        // surface color
    // float P = f.x > f.y && f.x > f.z ? f.x : f.y > f.z ? f.y : f.z; // Russian Roulete probability
    float P = 0.8; // Russian Roulete probability
    if (++depth > 5)
    {
        if (get_random_float() < P)
            f = f * (1 / P);
        else
            return Vector3f();
    }

    switch (p.m->m_type)
    {
    case DIFFUSE:
    {
        // sample next ray's direction
        float r1 = 2 * M_PI * get_random_float(), r2 = get_random_float(), r2s = sqrt(r2);
        Vector3f w = nl, u = normalize((fabs(w.x) > .1 ? Vector3f(0, 1, 0) : crossProduct(Vector3f(1, 0, 0), w))), v = crossProduct(w, u);
        Vector3f d = normalize(u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2));

        // sample light
        Vector3f e;
        Intersection light_sample;
        float pdf_light = 0.0f;
        sampleLight(light_sample, pdf_light);
        Ray LightRay(p.coords, light_sample.coords - x);
        Intersection LightInter = intersect(LightRay);
        Vector3f wi = normalize(light_sample.coords - p.coords);
        if (LightInter.happened && (LightInter.coords - light_sample.coords).norm() < EPSILON)
        {
            e = f * light_sample.emit *
                M_1_PI *
                dotProduct(wi, p.normal) *
                dotProduct(-wi, light_sample.normal) /
                pow((light_sample.coords - p.coords).norm(), 2) /
                pdf_light;
        }

        return e + f * castRaySmallTP(Ray(x, d), depth) + (depth < 2 ? p.m->getEmission() : Vector3f());
        break;
    }

    case METAL:
    {
        // sample next ray's direction
        float x_1 = get_random_float(), x_2 = get_random_float();
        float ag = (1.2 - 0.2 * sqrtf(fabs(dotProduct(ray.direction, n)))) * p.m->roughness;
        float theta = atanf(ag * sqrtf(x_1) / sqrtf(1 - x_1));
        float phi = 2 * M_PI * x_2;
        Vector3f sampled_n = Vector3f(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta));
        Ray reflRay(x, ray.direction - sampled_n * 2 * dotProduct(sampled_n, ray.direction)); // reflected ray
        bool into = dotProduct(n, nl) > 0;                                                    // Ray from outside going in?
        float nc = 1, nt = p.m->ior, nnt = into ? nc / nt : nt / nc, ddn = dotProduct(ray.direction, nl), cos2t;
        if ((cos2t = 1 - nnt * nnt * (1 - ddn * ddn)) < 0) // Total internal reflection
            return f * castRaySmallTP(reflRay, depth);
        Vector3f tdir = normalize(ray.direction * nnt - sampled_n * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t)))); // refracted ray
        float a = nt - nc, b = nt + nc, R0 = a * a / (b * b), c = 1 - (into ? -ddn : dotProduct(tdir, n));
        float Re = R0 + (1 - R0) * c * c * c * c * c, P = .25 + .5 * Re;
        float weight_term_1 = fabs(dotProduct(ray.direction, sampled_n)) * GXX_Geometry(ray.direction, sampled_n, n, ag) / fabs(dotProduct(ray.direction, n)) / fabs(dotProduct(sampled_n, n));
        float weight_term_r = weight_term_1 * GXX_Geometry(reflRay.direction, sampled_n, n, ag);
        float weight_term_t = weight_term_1 * GXX_Geometry(tdir, sampled_n, n, ag);
        if (get_random_float() < P)
            f = castRaySmallTP(reflRay, depth) / P * weight_term_r;
        else
            f = castRaySmallTP(Ray(x, tdir), depth) / P * weight_term_t;
        return f;
        break;
    }

    case MIRROR:
    {
        return castRaySmallTP(Ray(x, ray.direction - n * 2 * dotProduct(n, ray.direction)), depth) + p.m->getEmission();
        break;
    }

    case GLASS:
    {
        Ray reflRay(x, ray.direction - n * 2 * dotProduct(n, ray.direction)); // reflected ray
        bool into = dotProduct(n, nl) > 0;                                    // Ray from outside going in?
        float nc = 1, nt = p.m->ior, nnt = into ? nc / nt : nt / nc, ddn = dotProduct(ray.direction, nl), cos2t;
        if ((cos2t = 1 - nnt * nnt * (1 - ddn * ddn)) < 0) // Total internal reflection
            return castRaySmallTP(reflRay, depth);
        Vector3f tdir = normalize(ray.direction * nnt - n * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t)))); // refracted ray
        float a = nt - nc, b = nt + nc, R0 = a * a / (b * b), c = 1 - (into ? -ddn : dotProduct(tdir, n));
        float Re = R0 + (1 - R0) * c * c * c * c * c, Tr = 1 - Re, P = .25 + .5 * Re, RP = Re / P, TP = Tr / (1 - P);
        return (depth > 2 ? (get_random_float() < P ? // Russian roulette
                                 castRaySmallTP(reflRay, depth) * RP
                                                    : castRaySmallTP(Ray(x, tdir), depth) * TP)
                          : castRaySmallTP(reflRay, depth) * Re + castRaySmallTP(Ray(x, tdir), depth) * Tr) +
               p.m->getEmission();
        break;
    }

    default:
        return Vector3f();
        break;
    }
}