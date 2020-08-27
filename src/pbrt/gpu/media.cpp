// pbrt is Copyright(c) 1998-2020 Matt Pharr, Wenzel Jakob, and Greg Humphreys.
// The pbrt source code is licensed under the Apache License, Version 2.0.
// SPDX: Apache-2.0

#include <pbrt/gpu/pathintegrator.h>

#include <pbrt/gpu/accel.h>
#include <pbrt/gpu/launch.h>
#include <pbrt/media.h>

namespace pbrt {

// It's not unususal for these values to have very large or very small
// magnitudes after multiple (null) scattering events, even though in the
// end ratios like beta/pdfUni are generally around 1.  To avoid overflow,
// we rescale all three of them by the same factor when they become large.
PBRT_CPU_GPU
static inline void rescale(SampledSpectrum &beta, SampledSpectrum &pdfLight,
                           SampledSpectrum &pdfUni) {
    // Note that no precision is lost in the rescaling since we're always
    // multiplying by an exact power of 2.
    if (beta.MaxComponentValue() > 0x1p24f ||
        pdfLight.MaxComponentValue() > 0x1p24f ||
        pdfUni.MaxComponentValue() > 0x1p24f) {
        beta *= 1.f / 0x1p24f;
        pdfLight *= 1.f / 0x1p24f;
        pdfUni *= 1.f / 0x1p24f;
    } else if (beta.MaxComponentValue() < 0x1p-24f ||
               pdfLight.MaxComponentValue() < 0x1p-24f ||
               pdfUni.MaxComponentValue() < 0x1p-24f) {
        beta *= 0x1p24f;
        pdfLight *= 0x1p24f;
        pdfUni *= 0x1p24f;
    }
}

void GPUPathIntegrator::SampleMediumInteraction(int depth) {
    ForAllQueued(
        "Sample medium interaction", mediumSampleQueue, maxQueueSize,
        [=] PBRT_GPU(MediumSampleWorkItem ms, int index) {
            Ray ray = ms.ray;
            Float tMax = ms.tMax;

            PBRT_DBG("Sampling medium interaction ray index %d depth %d ray %f %f %f d %f %f "
                "%f tMax %f\n",
                ms.rayIndex, depth, ray.o.x, ray.o.y, ray.o.z, ray.d.x, ray.d.y, ray.d.z,
                tMax);

            SampledWavelengths lambda = ms.lambda;
            SampledSpectrum beta = ms.beta;
            SampledSpectrum pdfUni = ms.pdfUni;
            SampledSpectrum pdfNEE = ms.pdfNEE;
            SampledSpectrum L(0.f);
            RNG rng(Hash(tMax), Hash(ray.d));

            PBRT_DBG("Lambdas %f %f %f %f\n", lambda[0], lambda[1], lambda[2], lambda[3]);
            PBRT_DBG("Medium sample beta %f %f %f %f pdfUni %f %f %f %f pdfNEE %f %f %f %f\n",
                beta[0], beta[1], beta[2], beta[3], pdfUni[0], pdfUni[1], pdfUni[2],
                pdfUni[3], pdfNEE[0], pdfNEE[1], pdfNEE[2], pdfNEE[3]);

            // Sample the medium according to T_maj, the homogeneous
            // transmission function based on the majorant.
            bool scattered = false;
            ray.medium.SampleTmaj(
                ray, tMax, rng, lambda, [&](const MediumSample &mediumSample) {
                    if (!mediumSample.intr) {
                        // No interaction was sampled, but update the path
                        // throughput and unidirectional PDF to the end of
                        // the ray segment.
                        beta *= mediumSample.Tmaj;
                        pdfUni *= mediumSample.Tmaj;
                        PBRT_DBG("No intr: beta %f %f %f %f pdfUni %f %f %f %f\n", beta[0],
                            beta[1], beta[2], beta[3], pdfUni[0], pdfUni[1], pdfUni[2],
                            pdfUni[3]);
                        return false;
                    }

                    const MediumInteraction &intr = *mediumSample.intr;
                    const SampledSpectrum &sigma_a = intr.sigma_a;
                    const SampledSpectrum &sigma_s = intr.sigma_s;
                    const SampledSpectrum &Tmaj = mediumSample.Tmaj;

                    PBRT_DBG("Medium event Tmaj %f %f %f %f sigma_a %f %f %f %f sigma_s %f %f "
                        "%f %f\n",
                        Tmaj[0], Tmaj[1], Tmaj[2], Tmaj[3], sigma_a[0], sigma_a[1],
                        sigma_a[2], sigma_a[3], sigma_s[0], sigma_s[1], sigma_s[2],
                        sigma_s[3]);

                    // Add emission, if present.  Always do this and scale
                    // by sigma_a/sigma_maj rather than only doing it
                    // (without scaling) at absorption events.
                    if (depth < maxDepth && intr.Le)
                        L += beta * intr.Le * sigma_a /
                             (intr.sigma_maj[0] * pdfUni.Average());

                    // Compute probabilities for each type of scattering.
                    Float pAbsorb = sigma_a[0] / intr.sigma_maj[0];
                    Float pScatter = sigma_s[0] / intr.sigma_maj[0];
                    Float pNull = std::max<Float>(0, 1 - pAbsorb - pScatter);
                    PBRT_DBG("Medium scattering probabilities: %f %f %f\n", pAbsorb, pScatter,
                        pNull);

                    // And randomly choose one.
                    Float um = rng.Uniform<Float>();
                    int mode = SampleDiscrete({pAbsorb, pScatter, pNull}, um);

                    if (mode == 0) {
                        // Absorption--done.
                        PBRT_DBG("absorbed\n");
                        beta = SampledSpectrum(0.f);
                        // Tell the medium to stop traveral.
                        return false;
                    } else if (mode == 1) {
                        // Scattering.
                        PBRT_DBG("scattered\n");
                        beta *= Tmaj * sigma_s;
                        pdfUni *= Tmaj * sigma_s;

                        // TODO: don't hard code a phase function.
                        const HGPhaseFunction *phase =
                            intr.phase.CastOrNullptr<HGPhaseFunction>();
                        // Enqueue medium scattering work.
                        mediumScatterQueue->Push(MediumScatterWorkItem{
                            intr.p(), lambda, beta, pdfUni, ms.rayIndex, *phase, -ray.d,
                            ms.etaScale, ray.medium, ms.pixelIndex});
                        scattered = true;

                        return false;
                    } else {
                        // Null scattering.
                        PBRT_DBG("null-scattered\n");
                        SampledSpectrum sigma_n = intr.sigma_n();

                        beta *= Tmaj * sigma_n;
                        pdfUni *= Tmaj * sigma_n;
                        pdfNEE *= Tmaj * intr.sigma_maj;

                        rescale(beta, pdfUni, pdfNEE);

                        return true;
                    }
                });

            PBRT_DBG("Post ray medium sample L %f %f %f %f beta %f %f %f %f\n", L[0], L[1],
                L[2], L[3], beta[0], beta[1], beta[2], beta[3]);
            PBRT_DBG("Post ray medium sample pdfUni %f %f %f %f pdfNEE %f %f %f %f\n",
                pdfUni[0], pdfUni[1], pdfUni[2], pdfUni[3], pdfNEE[0], pdfNEE[1],
                pdfNEE[2], pdfNEE[3]);

            // Add any emission found to its pixel sample's L value.
            if (L) {
                SampledSpectrum Lp = pixelSampleState.L[ms.pixelIndex];
                pixelSampleState.L[ms.pixelIndex] = Lp + L;
                PBRT_DBG("Added emitted radiance %f %f %f %f at pixel index %d ray index %d\n",
                    L[0], L[1], L[2], L[3], ms.pixelIndex, ms.rayIndex);
            }

            // There's no more work to do if there was a scattering event in
            // the medium.
            if (scattered || !beta || depth == maxDepth)
                return;

            // Otherwise, enqueue bump and medium stuff...
            // FIXME: this is all basically duplicate code w/optix.cu
            if (ms.tMax == Infinity) {
                // no intersection
                if (escapedRayQueue) {
                    PBRT_DBG("Adding ray to escapedRayQueue pixel index %d depth %d\n",
                        ms.pixelIndex, depth);
                    escapedRayQueue->Push(EscapedRayWorkItem{
                        beta, pdfUni, pdfNEE, lambda, ray.o, ray.d, ms.piPrev, ms.nPrev,
                        ms.nsPrev, (int)ms.isSpecularBounce, ms.pixelIndex});
                }
            }

            MaterialHandle material = ms.material;
            if (!material) {
                Interaction intr(ms.pi, ms.n);
                intr.mediumInterface = &ms.mediumInterface;
                Ray newRay = intr.SpawnRay(ray.d);
                mediumTransitionQueue->Push(MediumTransitionWorkItem{
                    newRay, lambda, beta, pdfUni, pdfNEE, ms.piPrev, ms.nPrev, ms.nsPrev,
                    ms.isSpecularBounce, ms.anyNonSpecularBounces, ms.etaScale,
                    ms.pixelIndex});
#if 0
                // WHY NOT THIS?
                rayQueues[(depth + 1) & 1]->PushIndirect(newRay, ms.piPrev, ms.nPrev, ms.nsPrev,
                                                         beta, pdfUni, pdfNEE, lambda, ms.etaScale,
                                                         ms.isSpecularBounce, ms.anyNonSpecularBounces,
                                                         ms.pixelIndex);
#endif
                return;
            }

            if (ms.areaLight) {
                PBRT_DBG("Ray hit an area light: adding to hitAreaLightQueue pixel index %d "
                    "depth %d\n",
                    ms.pixelIndex, depth);
                // TODO: intr.wo == -ray.d?
                hitAreaLightQueue->Push(HitAreaLightWorkItem{
                    ms.areaLight, lambda, beta, pdfUni, pdfNEE, Point3f(ms.pi), ms.n,
                    ms.uv, -ray.d, ms.piPrev, ray.d, ray.time, ms.nPrev, ms.nsPrev,
                    ms.isSpecularBounce, ms.pixelIndex});
            }

            FloatTextureHandle displacement = material.GetDisplacement();

            MaterialEvalQueue *q =
                (material.CanEvaluateTextures(BasicTextureEvaluator()) &&
                 (!displacement ||
                  BasicTextureEvaluator().CanEvaluate({displacement}, {})))
                    ? basicEvalMaterialQueue
                    : universalEvalMaterialQueue;

            PBRT_DBG("Enqueuing for material eval, mtl tag %d", material.Tag());

            auto enqueue = [=](auto ptr) {
                using Material = typename std::remove_reference_t<decltype(*ptr)>;
                q->Push<Material>(MaterialEvalWorkItem<Material>{
                    ptr, lambda, beta, pdfUni, ms.pi, ms.n, ms.ns, ms.dpdus, ms.dpdvs,
                    ms.dndus, ms.dndvs, -ray.d, ms.uv, ray.time, ms.anyNonSpecularBounces,
                    ms.etaScale, ms.mediumInterface, ms.rayIndex, ms.pixelIndex});
            };
            material.Dispatch(enqueue);
        });

    if (depth == maxDepth)
        return;

    using PhaseFunction = HGPhaseFunction;
    std::string desc = std::string("Sample direct/indirect - Henyey Greenstein");
    ForAllQueued(
        desc.c_str(), mediumScatterQueue, maxQueueSize,
        [=] PBRT_GPU(MediumScatterWorkItem ms, int index) {
            RaySamples raySamples = rayQueues[depth & 1]->raySamples[ms.rayIndex];
            Float time = 0;  // TODO: FIXME
            Vector3f wo = ms.wo;

            // Sample direct lighting at medium scattering event.  First,
            // choose a light source.
            LightSampleContext ctx(Point3fi(ms.p), Normal3f(0, 0, 0), Normal3f(0, 0, 0));
            pstd::optional<SampledLight> sampledLight =
                lightSampler.Sample(ctx, raySamples.direct.uc);

            if (sampledLight) {
                LightHandle light = sampledLight->light;
                // And now sample a point on the light.
                LightLiSample ls = light.SampleLi(ctx, raySamples.direct.u, ms.lambda,
                                                  LightSamplingMode::WithMIS);
                if (ls && ls.L) {
                    Vector3f wi = ls.wi;
                    SampledSpectrum beta = ms.beta * ms.phase.p(wo, wi);

                    PBRT_DBG("Phase phase beta %f %f %f %f\n", beta[0], beta[1], beta[2],
                        beta[3]);

                    // Compute PDFs for direct lighting MIS calculation.
                    Float lightPDF = ls.pdf * sampledLight->pdf;
                    Float phasePDF =
                        IsDeltaLight(light.Type()) ? 0.f : ms.phase.PDF(wo, wi);
                    SampledSpectrum pdfUni = ms.pdfUni * phasePDF;
                    SampledSpectrum pdfNEE = ms.pdfUni * lightPDF;

                    SampledSpectrum Ld = beta * ls.L;
                    Ray ray(ms.p, ls.pLight.p() - ms.p, time, ms.medium);

                    // Enqueue shadow ray
                    shadowRayQueue->Push(ray, 1 - ShadowEpsilon, ms.lambda, Ld,
                                         pdfUni, pdfNEE, ms.pixelIndex);

                    PBRT_DBG("Enqueued medium shadow ray depth %d "
                        "Ld %f %f %f %f pdfUni %f %f %f %f "
                        "pdfNEE %f %f %f %f parent ray index %d parent pixel index %d\n",
                        depth, Ld[0], Ld[1], Ld[2], Ld[3], pdfUni[0], pdfUni[1],
                        pdfUni[2], pdfUni[3], pdfNEE[0], pdfNEE[1], pdfNEE[2], pdfNEE[3],
                        ms.rayIndex, ms.pixelIndex);
                }
            }

            // Sample indirect lighting.
            PhaseFunctionSample phaseSample =
                ms.phase.Sample_p(wo, raySamples.indirect.u);
            if (!phaseSample)
                return;

            SampledSpectrum beta = ms.beta * phaseSample.p;
            SampledSpectrum pdfUni = ms.pdfUni * phaseSample.pdf;
            SampledSpectrum pdfNEE = ms.pdfUni;

            // Russian roulette
            // TODO: should we even bother? Generally beta/pdfUni is one here,
            // due to the way scattering events are scattered and because we're
            // sampling exactly from the phase function's distribution...
            SampledSpectrum rrBeta = beta * ms.etaScale / pdfUni.Average();
            if (rrBeta.MaxComponentValue() < 1 && depth > 1) {
                Float q = std::max<Float>(0, 1 - rrBeta.MaxComponentValue());
                if (raySamples.indirect.rr < q) {
                    PBRT_DBG("RR terminated medium indirect with q %f ray index %d\n", q,
                        ms.rayIndex);
                    return;
                }
                pdfUni *= 1 - q;
                pdfNEE *= 1 - q;
            }

            Ray ray(ms.p, phaseSample.wi, time, ms.medium);
            bool isSpecularBounce = false;
            bool anyNonSpecularBounces = true;

            // Spawn indirect ray.
            rayQueues[(depth + 1) & 1]->PushIndirect(
                ray, Point3fi(ms.p), Normal3f(0, 0, 0), Normal3f(0, 0, 0), beta, pdfUni,
                pdfNEE, ms.lambda, ms.etaScale, isSpecularBounce, anyNonSpecularBounces,
                ms.pixelIndex);
            PBRT_DBG("Enqueuing indirect medium ray at depth %d ray index %d pixel index %d\n",
                depth + 1, ms.rayIndex, ms.pixelIndex);
        });
}

void GPUPathIntegrator::HandleMediumTransitions(int depth) {
    ForAllQueued(
        "Handle medium transitions", mediumTransitionQueue, maxQueueSize,
        [=] PBRT_GPU(MediumTransitionWorkItem mt, int index) {
            // Have to do this here, later, since we can't be writing into
            // the other ray queue in optix closest hit.  (Wait--really?
            // Why not? Basically boils down to current indirect enqueue (and other
            // places?))
            // TODO: figure this out...
            rayQueues[(depth + 1) & 1]->PushIndirect(
                mt.ray, mt.piPrev, mt.nPrev, mt.nsPrev, mt.beta, mt.pdfUni, mt.pdfNEE,
                mt.lambda, mt.etaScale, mt.isSpecularBounce, mt.anyNonSpecularBounces,
                mt.pixelIndex);
            PBRT_DBG("Enqueuied ray after medium transition at depth %d pixel index %d",
                depth + 1, mt.pixelIndex);
        });
}

}  // namespace pbrt
