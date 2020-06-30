//
// Created by LEI XU on 5/16/19.
//

#ifndef RAYTRACING_MATERIAL_H
#define RAYTRACING_MATERIAL_H

#include "Vector.hpp"

enum MaterialType
{
    DIFFUSE,
    METAL,
    MIRROR,
    GLASS
};

class Material
{
private:
    // Compute reflection direction
    Vector3f reflect(const Vector3f &I, const Vector3f &N) const
    {
        return I + 2 * dotProduct(I, N) * N;
    }

    // Compute refraction direction using Snell's law
    //
    // We need to handle with care the two possible situations:
    //
    //    - When the ray is inside the object
    //
    //    - When the ray is outside.
    //
    // If the ray is outside, you need to make cosi positive cosi = -N.I
    //
    // If the ray is inside, you need to invert the refractive indices and negate the normal N
    Vector3f refract(const Vector3f &I, const Vector3f &N, const float &ior) const
    {
        float cosi = clamp(-1, 1, dotProduct(I, N));
        float etai = 1, etat = ior;
        Vector3f n = N;
        if (cosi < 0)
        {
            cosi = -cosi;
        }
        else
        {
            std::swap(etai, etat);
            n = -N;
        }
        float eta = etai / etat;
        float k = 1 - eta * eta * (1 - cosi * cosi);
        return k < 0 ? 0 : eta * I + (eta * cosi - sqrtf(k)) * n;
    }

    Vector3f refractEDIT(const Vector3f &I, const Vector3f &N, const float &ior)
    {
        float n1 = 1, n2 = ior;
        float cosTheta1 = dotProduct(I, N);
        bool into = cosTheta1 > 0;
        if (!into)
        {
            std::swap(n1, n2);
            cosTheta1 = -cosTheta1;
        }
        float N1_2 = n1 / n2;
        float det = 1 - N1_2 * N1_2 * (1 - cosTheta1 * cosTheta1);
        if (det < 0)
            return reflect(I, N);
        float cosTheta2 = sqrt(det);
        float sinTheta2 = sqrtf(1 - det);
        Vector3f B = -normalize(I - cosTheta1 * N * (into ? 1 : -1));
        return normalize(cosTheta2 * N * (into ? -1 : 1) + sinTheta2 * B);
    }
    Vector3f refractEDIT2(const Vector3f &I, const Vector3f &N, const float &ior)
    {
        Vector3f d = -I;
        Vector3f n = N;
        Vector3f nl = dotProduct(n, d) < 0 ? n : n * -1;
        bool into = dotProduct(n, nl) > 0;
        double nc = 1, nt = ior, nnt = into ? nc / nt : nt / nc, ddn = dotProduct(d, nl), cos2t;
        if ((cos2t = 1 - nnt * nnt * (1 - ddn * ddn)) < 0) // Total internal reflection
            return reflect(I, N);
        return normalize(d * nnt - n * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t))));
    }

    // Compute Fresnel equation
    //
    // \param I is the incident view direction
    //
    // \param N is the normal at the intersection point
    //
    // \param ior is the material refractive index
    //
    // \param[out] kr is the amount of light reflected
    void fresnel(const Vector3f &I, const Vector3f &N, const float &ior, float &kr, float &kt) const
    {
        float cosi = clamp(-1, 1, dotProduct(I, N));
        float etai = 1, etat = ior;
        if (cosi < EPSILON)
        {
            std::swap(etai, etat);
        }
        // Compute sini using Snell's law
        float sint = etai / etat * sqrtf(std::max(0.f, 1 - cosi * cosi));
        // Total internal reflection
        if (sint >= 1)
        {
            kr = 1;
        }
        else
        {
            float cost = sqrtf(std::max(0.f, 1 - sint * sint));
            cosi = fabsf(cosi);
            float Rs = ((etat * cosi) - (etai * cost)) / ((etat * cosi) + (etai * cost));
            float Rp = ((etai * cosi) - (etat * cost)) / ((etai * cosi) + (etat * cost));
            kr = (Rs * Rs + Rp * Rp) / 2;
        }
        // As a consequence of the conservation of energy, transmittance is given by:
        kt = 1 - kr;
    }
    void fresnelEDIT(const Vector3f &I, const Vector3f &N, const Vector3f &O, const float &ior, float &kr, float &kt) const
    {
        Vector3f n = N;
        float nc = 1, nt = ior;
        Vector3f nl = dotProduct(-I, n) < 0 ? n : n * -1;
        bool into = dotProduct(n, nl) > 0;
        float ddn = dotProduct(-I, nl);
        float a = nt - nc, b = nt + nc;
        float R0 = a * a / (b * b);                     // reflectance at normal incidence based on IOR
        float c = 1 - (into ? -ddn : dotProduct(O, n)); // 1-cos(theta)
        kr = R0 + (1 - R0) * c * c * c * c * c;         // fresnel reflectance
        kt = 1 - kr;
    }

    Vector3f toWorld(const Vector3f &a, const Vector3f &N)
    {
        Vector3f B, C;
        if (std::fabs(N.x) > std::fabs(N.y))
        {
            float invLen = 1.0f / std::sqrt(N.x * N.x + N.z * N.z);
            C = Vector3f(N.z * invLen, 0.0f, -N.x * invLen);
        }
        else
        {
            float invLen = 1.0f / std::sqrt(N.y * N.y + N.z * N.z);
            C = Vector3f(0.0f, N.z * invLen, -N.y * invLen);
        }
        B = crossProduct(C, N);
        return a.x * B + a.y * C + a.z * N;
    }

public:
    MaterialType m_type;
    //Vector3f m_color;
    Vector3f m_emission;
    float ior;
    Vector3f Kd = Vector3f(1.0), Ks = Vector3f(0.0);
    float roughness;
    //Texture tex;

    inline Material(MaterialType t = DIFFUSE, Vector3f e = Vector3f(0, 0, 0));
    inline MaterialType getType();
    //inline Vector3f getColor();
    inline Vector3f getColorAt(double u, double v);
    inline Vector3f getEmission();
    inline bool hasEmission();

    inline Vector3f sample(const Vector3f &wi, const Vector3f &N); // sample an incidence ray
    // given a ray, calculate the PdF of this ray
    inline float pdf(const Vector3f &wi, const Vector3f &wo, const Vector3f &N);
    // given a ray, calculate the contribution of this ray
    inline Vector3f eval(const Vector3f &wi, const Vector3f &wo, const Vector3f &N);

    // Functions for Microfacet
    // http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.443.6410&rep=rep1&type=pdf
    // http://www.codinglabs.net/article_physically_based_rendering_cook_torrance.aspx
    float chiGGX(float v)
    {
        return v > 0 ? 1 : 0;
    }
    float GGX_Distribution(Vector3f n, Vector3f h, float alpha)
    {
        float NoH = dotProduct(n, h);
        float alpha2 = alpha * alpha;
        float NoH2 = NoH * NoH;
        float den = NoH2 * alpha2 + (1 - NoH2);
        return (chiGGX(NoH) * alpha2) / (M_PI * den * den);
    }
    float GGX_PartialGeometryTerm(Vector3f v, Vector3f n, Vector3f h, float alpha)
    {
        float VoH2 = dotProduct(v, h);
        float chi = chiGGX(dotProduct(v, h) / dotProduct(v, n));
        VoH2 = VoH2 * VoH2;
        float tan2 = (1 - VoH2) / VoH2;
        return (chi * 2) / (1 + sqrt(1 + alpha * alpha * tan2));
    }
};

Material::Material(MaterialType t, Vector3f e)
{
    m_type = t;
    //m_color = c;
    m_emission = e;
    switch (t)
    {
    case METAL:
        ior = 1.5f;
        roughness = 0.1f;
        break;
    case GLASS:
        ior = 1.5f;
        roughness = 0.001f;
    default:
        break;
    }
}

MaterialType Material::getType() { return m_type; }
///Vector3f Material::getColor(){return m_color;}
Vector3f Material::getEmission() { return m_emission; }
bool Material::hasEmission()
{
    if (m_emission.norm() > EPSILON)
        return true;
    else
        return false;
}

Vector3f Material::getColorAt(double u, double v)
{
    return Vector3f();
}

Vector3f Material::sample(const Vector3f &wo, const Vector3f &N)
{
    switch (m_type)
    {
    case DIFFUSE:
    {
        // uniform sample on the hemisphere
        float x_1 = get_random_float(), x_2 = get_random_float();
        float z = std::fabs(1.0f - 2.0f * x_1);
        float r = std::sqrt(1.0f - z * z), phi = 2 * M_PI * x_2;
        Vector3f localRay(r * std::cos(phi), r * std::sin(phi), z);
        return toWorld(localRay, N);

        break;
    }
    case MIRROR:
    {
        return reflect(wo, N);
        break;
    }
    case GLASS:
    {
        // float FresnelKr = 1.0f;
        // float FresnelKt = 0.0f;
        // Vector3f reflDIR = reflect(wo, N);
        Vector3f refrDIR = refractEDIT2(wo, N, ior);
        // fresnelEDIT(refrDIR, N, wo, ior, FresnelKr, FresnelKt);
        // float P = .25 + .5 * FresnelKr;
        // if (get_random_float() < P)
        // {
        //     return reflDIR;
        // }
        // else
        // {
        //     return refrDIR;
        // }
        return refrDIR;
    }
    default:
    {
        float x_1 = get_random_float(), x_2 = get_random_float();
        float z = std::fabs(1.0f - 2.0f * x_1);
        float r = std::sqrt(1.0f - z * z), phi = 2 * M_PI * x_2;
        Vector3f localRay(r * std::cos(phi), r * std::sin(phi), z);
        return toWorld(localRay, N);
    }
    }
}

float Material::pdf(const Vector3f &wi, const Vector3f &wo, const Vector3f &N)
{
    switch (m_type)
    {
    case DIFFUSE:
    {
        // uniform sample probability 1 / (2 * PI)
        if (dotProduct(wi, N) > 0.0f)
            return 0.5f / M_PI;
        else
            return 0.0f;
        break;
    }
    case METAL:
    {
        if (dotProduct(wi, N) < 0.0f)
            return 0.0f;
        Vector3f h = normalize(wi + wo);
        return GGX_Distribution(h, N, roughness);
        break;
    }
    case MIRROR:
    {
        return 1.0f;
        break;
    }
    case GLASS:
    {
        float FresnelKr = 1.0f;
        float FresnelKt = 0.0f;
        fresnelEDIT(wi, N, wo, ior, FresnelKr, FresnelKt);
        float P = .25 + .5 * FresnelKr;
        bool isReflect = dotProduct(wi, N) > 0;
        if (isReflect)
            return P;
        else
            return (1 - P);
    }
    default:
    {
        if (dotProduct(wi, N) > 0.0f)
            return 0.5f / M_PI;
        else
            return 0.0f;
        break;
    }
    }
}

Vector3f Material::eval(const Vector3f &wi, const Vector3f &wo, const Vector3f &N)
{
    switch (m_type)
    {
    case DIFFUSE:
    {
        // calculate the contribution of diffuse model
        float cosalpha = dotProduct(N, wo);
        if (cosalpha > 0.0f)
        {
            Vector3f diffuse = Kd / M_PI;
            return diffuse;
        }
        else
            return Vector3f(0.0f);
        break;
    }
    case METAL:
    {
        float cosalpha = dotProduct(N, wo);
        if (cosalpha < 0.0f)
            return Vector3f(0.0f);
        Vector3f h = normalize(wi + wo);
        float FresnelKr = 1.0f;
        float FresnelKt = 0.0f;
        // fresnelEDIT(wi, N, wo, ior, FresnelKr, FresnelKt);
        fresnel(wi, N, ior, FresnelKr, FresnelKt);
        float ShadowingMaskingTerm = GGX_PartialGeometryTerm(N, N, h, roughness);
        float DSTerm = GGX_Distribution(h, N, roughness);
        return 0.2 * Kd / M_PI + 0.8 * FresnelKr * ShadowingMaskingTerm * DSTerm / (4.0f * fabs(dotProduct(wi, N) * dotProduct(wo, N)));
        break;
    }
    case GLASS:
    {
        // float FresnelKr = 1.0f;
        // float FresnelKt = 0.0f;
        // bool isReflect = dotProduct(wi, N) > 0;
        // fresnelEDIT(wi, N, wo, ior, FresnelKr, FresnelKt);
        // float P = .25 + .5 * FresnelKr;
        // float FresnelTerm = isReflect ? FresnelKr : FresnelKt;
        // return FresnelTerm;
        return 1;
        break;
    }
    case MIRROR:
    {
        float cosalpha = dotProduct(N, wo);
        bool isReflect = 1.0 - dotProduct(normalize(wi + wo), N) < EPSILON;
        if (cosalpha < 0.0f || !isReflect)
            return Vector3f(0.0f);
        else
            return Vector3f(1.0f);
        break;
    }
    default:
        // Diffuse
        float cosalpha = dotProduct(N, wo);
        if (cosalpha > 0.0f)
        {
            Vector3f diffuse = Kd / M_PI;
            return diffuse;
        }
        else
            return Vector3f(0.0f);
        break;
    }
}

#endif //RAYTRACING_MATERIAL_H
