// pbrt is Copyright(c) 1998-2020 Matt Pharr, Wenzel Jakob, and Greg Humphreys.
// The pbrt source code is licensed under the Apache License, Version 2.0.
// SPDX: Apache-2.0

#ifndef PBRT_GPU_WORKITEMS_H
#define PBRT_GPU_WORKITEMS_H

#include <pbrt/pbrt.h>

#include <pbrt/base/sampler.h>
#include <pbrt/film.h>
#include <pbrt/gpu/workqueue.h>
#include <pbrt/lightsamplers.h>
#include <pbrt/materials.h>
#include <pbrt/ray.h>
#include <pbrt/util/pstd.h>
#include <pbrt/util/soa.h>

namespace pbrt {

struct RaySamples {
    struct {
        Point2f u;
        Float uc;
    } direct;
    struct {
        Float uc, rr;
        Point2f u;
    } indirect;

    bool haveSubsurface;
    struct {
        Float uc;
        Point2f u;
    } subsurface;
};

template <>
struct SOA<RaySamples> {
  public:
    SOA() = default;

    SOA(int size, Allocator alloc) {
        direct = alloc.allocate_object<Float4>(size);
        indirect = alloc.allocate_object<Float4>(size);
        subsurface = alloc.allocate_object<Float4>(size);
    }

    PBRT_CPU_GPU
    RaySamples operator[](int i) const {
        RaySamples rs;
        Float4 dir = Load4(direct + i);
        rs.direct.u = Point2f(dir.v[0], dir.v[1]);
        rs.direct.uc = dir.v[2];

        Float4 ind = Load4(indirect + i);
        rs.indirect.uc = ind.v[0];
        rs.indirect.rr = ind.v[1];
        rs.indirect.u = Point2f(ind.v[2], ind.v[3]);

        rs.haveSubsurface = dir.v[3] != 0;
        if (rs.haveSubsurface) {
            Float4 ss = Load4(subsurface + i);
            rs.subsurface.uc = ss.v[0];
            rs.subsurface.u = Point2f(ss.v[1], ss.v[2]);
        }

        return rs;
    }

    struct GetSetIndirector {
        PBRT_CPU_GPU
        operator RaySamples() const { return (*(const SOA *)soa)[index]; }

        PBRT_CPU_GPU
        void operator=(RaySamples rs) {
            soa->direct[index] = Float4{rs.direct.u[0], rs.direct.u[1], rs.direct.uc,
                                        Float(rs.haveSubsurface)};
            soa->indirect[index] = Float4{rs.indirect.uc, rs.indirect.rr,
                                          rs.indirect.u[0], rs.indirect.u[1]};
            if (rs.haveSubsurface)
                soa->subsurface[index] =
                    Float4{rs.subsurface.uc, rs.subsurface.u.x, rs.subsurface.u.y, 0.f};
        }

        SOA *soa;
        int index;
    };

    PBRT_CPU_GPU
    GetSetIndirector operator[](int i) { return GetSetIndirector{this, i}; }

  private:
    Float4 *__restrict__ direct;
    Float4 *__restrict__ indirect;
    Float4 *__restrict__ subsurface;
};

struct PixelSampleState {
    Float filterWeight;
    Point2f pFilm;
    Point2i pPixel;
    SampledWavelengths lambda;
    SampledSpectrum L;
    SampledSpectrum cameraRayWeight;
    VisibleSurface visibleSurface;
    RaySamples samples;
};

struct RayWorkItem {
    Ray ray;
    int pixelIndex;
    SampledWavelengths lambda;
    SampledSpectrum beta, pdfUni, pdfNEE;
    Point3fi piPrev;
    Normal3f nPrev;
    Normal3f nsPrev;
    Float etaScale;
    int isSpecularBounce;
    int anyNonSpecularBounces;
};

struct EscapedRayWorkItem {
    SampledSpectrum beta, pdfUni, pdfNEE;
    SampledWavelengths lambda;
    Point3f rayo;
    Vector3f rayd;
    Point3fi piPrev;
    Normal3f nPrev, nsPrev;
    int specularBounce;
    int pixelIndex;
};

using EscapedRayQueue = WorkQueue<EscapedRayWorkItem>;

struct HitAreaLightWorkItem {
    LightHandle areaLight;
    SampledWavelengths lambda;
    SampledSpectrum beta, pdfUni, pdfNEE;
    Point3f p;
    Normal3f n;
    Point2f uv;
    Vector3f wo;
    Point3fi piPrev;
    Vector3f rayd;
    Float time;
    Normal3f nPrev, nsPrev;
    int isSpecularBounce;
    int pixelIndex;
};

using HitAreaLightQueue = WorkQueue<HitAreaLightWorkItem>;

struct ShadowRayWorkItem {
    Ray ray;
    Float tMax;
    SampledWavelengths lambda;
    SampledSpectrum Ld, pdfUni, pdfNEE;
    int pixelIndex;
};

template <typename Material>
struct MaterialEvalWorkItem {
    PBRT_CPU_GPU
    BumpEvalContext GetBumpEvalContext() const {
        BumpEvalContext ctx;
        ctx.p = Point3f(pi);
        ctx.uv = uv;
        ctx.shading.n = ns;
        ctx.shading.dpdu = dpdus;
        ctx.shading.dpdv = dpdvs;
        ctx.shading.dndu = dndus;
        ctx.shading.dndv = dndvs;
        return ctx;
    }

    PBRT_CPU_GPU
    MaterialEvalContext GetMaterialEvalContext(Normal3f ns, Vector3f dpdus) const {
        MaterialEvalContext ctx;
        ctx.wo = wo;
        ctx.n = n;
        ctx.ns = ns;
        ctx.dpdus = dpdus;
        ctx.p = Point3f(pi);
        ctx.uv = uv;
        return ctx;
    }

    const Material *material;
    SampledWavelengths lambda;
    SampledSpectrum beta, pdfUni;
    Point3fi pi;
    Normal3f n, ns;
    Vector3f dpdus, dpdvs;
    Normal3f dndus, dndvs;
    Vector3f wo;
    Point2f uv;
    Float time;
    int anyNonSpecularBounces;
    Float etaScale;
    MediumInterface mediumInterface;
    int pixelIndex;
};

struct GetBSSRDFAndProbeRayWorkItem {
    PBRT_CPU_GPU
    MaterialEvalContext GetMaterialEvalContext() const {
        MaterialEvalContext ctx;
        ctx.wo = wo;
        ctx.n = n;
        ctx.ns = ns;
        ctx.dpdus = dpdus;
        ctx.p = p;
        ctx.uv = uv;
        return ctx;
    }

    MaterialHandle material;
    SampledWavelengths lambda;
    SampledSpectrum beta, pdfUni;
    Point3f p;
    Vector3f wo;
    Normal3f n, ns;
    Vector3f dpdus;
    Point2f uv;
    MediumInterface mediumInterface;
    Float etaScale;
    int pixelIndex;
};

struct SubsurfaceScatterWorkItem {
    Point3f p0, p1;
    MaterialHandle material;
    TabulatedBSSRDF bssrdf;
    SampledWavelengths lambda;
    SampledSpectrum beta, pdfUni;
    Float weight;
    Float uLight;
    SubsurfaceInteraction ssi;
    MediumInterface mediumInterface;
    Float etaScale;
    int pixelIndex;
};

struct MediumTransitionWorkItem {
    Ray ray;
    SampledWavelengths lambda;
    SampledSpectrum beta, pdfUni, pdfNEE;
    Point3fi piPrev;
    Normal3f nPrev, nsPrev;
    int isSpecularBounce;
    int anyNonSpecularBounces;
    Float etaScale;
    int pixelIndex;
};

using MediumTransitionQueue = WorkQueue<MediumTransitionWorkItem>;

struct MediumSampleWorkItem {
    // Both enqueue types (have mtl and no hit)
    Ray ray;
    Float tMax;
    SampledWavelengths lambda;
    SampledSpectrum beta;
    SampledSpectrum pdfUni;
    SampledSpectrum pdfNEE;
    int pixelIndex;
    Point3fi piPrev;
    Normal3f nPrev;
    Normal3f nsPrev;
    int isSpecularBounce;
    int anyNonSpecularBounces;
    Float etaScale;

    // Have a hit material as well
    LightHandle areaLight;
    Point3fi pi;
    Normal3f n;
    Vector3f wo;
    Point2f uv;
    MaterialHandle material;
    Normal3f ns;
    Vector3f dpdus;
    Vector3f dpdvs;
    Normal3f dndus;
    Normal3f dndvs;
    MediumInterface mediumInterface;
};

struct MediumScatterWorkItem {
    Point3f p;
    SampledWavelengths lambda;
    SampledSpectrum beta, pdfUni;
    HGPhaseFunction phase;
    Vector3f wo;
    Float etaScale;
    MediumHandle medium;
    int pixelIndex;
};

#include "gpu_workitems_soa.h"

class RayQueue : public WorkQueue<RayWorkItem> {
  public:
    using WorkQueue::WorkQueue;

    PBRT_CPU_GPU
    int PushCameraRay(const Ray &ray, const SampledWavelengths &lambda, int pixelIndex) {
        int index = AllocateEntry();
        DCHECK(!ray.HasNaN());
        this->ray[index] = ray;
        this->pixelIndex[index] = pixelIndex;
        this->lambda[index] = lambda;
        this->beta[index] = SampledSpectrum(1.f);
        this->etaScale[index] = 1.f;
        this->anyNonSpecularBounces[index] = false;
        this->pdfUni[index] = SampledSpectrum(1.f);
        this->pdfNEE[index] = SampledSpectrum(1.f);
        this->isSpecularBounce[index] = false;
        return index;
    }

    PBRT_CPU_GPU
    int PushIndirect(const Ray &ray, const Point3fi &piPrev, const Normal3f &nPrev,
                     const Normal3f &nsPrev, const SampledSpectrum &beta,
                     const SampledSpectrum &pdfUni, const SampledSpectrum &pdfNEE,
                     const SampledWavelengths &lambda, Float etaScale,
                     bool isSpecularBounce, bool anyNonSpecularBounces, int pixelIndex) {
        int index = AllocateEntry();
        DCHECK(!ray.HasNaN());
        this->ray[index] = ray;
        this->pixelIndex[index] = pixelIndex;
        this->piPrev[index] = piPrev;
        this->nPrev[index] = nPrev;
        this->nsPrev[index] = nsPrev;
        this->beta[index] = beta;
        this->pdfUni[index] = pdfUni;
        this->pdfNEE[index] = pdfNEE;
        this->lambda[index] = lambda;
        this->anyNonSpecularBounces[index] = anyNonSpecularBounces;
        this->isSpecularBounce[index] = isSpecularBounce;
        this->etaScale[index] = etaScale;
        return index;
    }
};

class ShadowRayQueue : public WorkQueue<ShadowRayWorkItem> {
public:
    using WorkQueue::WorkQueue;

    PBRT_CPU_GPU
    void Push(const Ray &ray, Float tMax, SampledWavelengths lambda,
              SampledSpectrum Ld, SampledSpectrum pdfUni, SampledSpectrum pdfNEE,
              int pixelIndex) {
        WorkQueue<ShadowRayWorkItem>::Push(
            ShadowRayWorkItem{ray, tMax, lambda, Ld, pdfUni, pdfNEE, pixelIndex});
    }
};

class GetBSSRDFAndProbeRayQueue : public WorkQueue<GetBSSRDFAndProbeRayWorkItem> {
public:
    using WorkQueue::WorkQueue;

    PBRT_CPU_GPU
    int Push(MaterialHandle material, SampledWavelengths lambda,
             SampledSpectrum beta, SampledSpectrum pdfUni,
             Point3f p, Vector3f wo, Normal3f n, Normal3f ns,
             Vector3f dpdus, Point2f uv, MediumInterface mediumInterface,
             Float etaScale, int pixelIndex) {
        int index = AllocateEntry();
        this->material[index] = material;
        this->lambda[index] = lambda;
        this->beta[index] = beta;
        this->pdfUni[index] = pdfUni;
        this->p[index] = p;
        this->wo[index] = wo;
        this->n[index] = n;
        this->ns[index] = ns;
        this->dpdus[index] = dpdus;
        this->uv[index] = uv;
        this->mediumInterface[index] = mediumInterface;
        this->etaScale[index] = etaScale;
        this->pixelIndex[index] = pixelIndex;
        return index;
    }
};

class SubsurfaceScatterQueue : public WorkQueue<SubsurfaceScatterWorkItem> {
  public:
    using WorkQueue::WorkQueue;

    PBRT_CPU_GPU
    int Push(Point3f p0, Point3f p1, MaterialHandle material, TabulatedBSSRDF bssrdf,
             SampledWavelengths lambda, SampledSpectrum beta, SampledSpectrum pdfUni,
             MediumInterface mediumInterface, Float etaScale, int pixelIndex) {
        int index = AllocateEntry();
        this->p0[index] = p0;
        this->p1[index] = p1;
        this->material[index] = material;
        this->bssrdf[index] = bssrdf;
        this->lambda[index] = lambda;
        this->beta[index] = beta;
        this->pdfUni[index] = pdfUni;
        this->mediumInterface[index] = mediumInterface;
        this->etaScale[index] = etaScale;
        this->pixelIndex[index] = pixelIndex;
        return index;
    }
};

class MediumSampleQueue : public WorkQueue<MediumSampleWorkItem> {
  public:
    using WorkQueue::WorkQueue;

    using WorkQueue::Push;

    PBRT_CPU_GPU
    int Push(Ray ray, Float tMax, SampledWavelengths lambda, SampledSpectrum beta,
             SampledSpectrum pdfUni, SampledSpectrum pdfNEE, int pixelIndex,
             Point3fi piPrev, Normal3f nPrev, Normal3f nsPrev, int isSpecularBounce,
             int anyNonSpecularBounces, Float etaScale) {
        int index = AllocateEntry();
        this->ray[index] = ray;
        this->tMax[index] = tMax;
        this->lambda[index] = lambda;
        this->beta[index] = beta;
        this->pdfUni[index] = pdfUni;
        this->pdfNEE[index] = pdfNEE;
        this->pixelIndex[index] = pixelIndex;
        this->piPrev[index] = piPrev;
        this->nPrev[index] = nPrev;
        this->nsPrev[index] = nsPrev;
        this->isSpecularBounce[index] = isSpecularBounce;
        this->anyNonSpecularBounces[index] = anyNonSpecularBounces;
        this->etaScale[index] = etaScale;
        return index;
    }
};

using MediumScatterQueue = WorkQueue<MediumScatterWorkItem>;

using MaterialEvalQueue =
    MultiWorkQueue<MaterialEvalWorkItem, CoatedDiffuseMaterial, CoatedConductorMaterial,
                   ConductorMaterial, DielectricMaterial, DiffuseMaterial,
                   DiffuseTransmissionMaterial, HairMaterial, MeasuredMaterial,
                   SubsurfaceMaterial, ThinDielectricMaterial, MixMaterial>;

}  // namespace pbrt

#endif  // PBRT_GPU_WORKITEMS_H
