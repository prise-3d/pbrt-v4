// pbrt is Copyright(c) 1998-2020 Matt Pharr, Wenzel Jakob, and Greg Humphreys.
// The pbrt source code is licensed under the Apache License, Version 2.0.
// SPDX: Apache-2.0

#ifndef PBRT_BASE_ESTIMATOR_H
#define PBRT_BASE_ESTIMATOR_H

#include <pbrt/pbrt.h>

#include <pbrt/base/filter.h>
#include <pbrt/util/pstd.h>
#include <pbrt/util/taggedptr.h>

#include <string>

namespace pbrt {

class MeanEstimator;
class MONEstimator;
class AlphaMONEstimator;
class AutoAlphaMONEstimator;
class PakMONEstimator;
class MeanOrMONEstimator;

class AtomicDouble;
struct PixelWindow;

// EstimatorHandle Definition
class EstimatorHandle : public TaggedPointer<MeanEstimator, MONEstimator, 
                                            AlphaMONEstimator, AutoAlphaMONEstimator, 
                                            PakMONEstimator, MeanOrMONEstimator> {
  public:
    // Estimator Interface
    
    PBRT_CPU_GPU inline void GetEstimation(const Point2i &pFilm, RGB &rgb, 
                                        Float &weightSum, 
                                        AtomicDouble* splatRGB) const;

    PBRT_CPU_GPU inline void Estimate(const PixelWindow &pixelWindow, RGB &rgb, 
                                        Float &weightSum, 
                                        AtomicDouble* splatRGB) const;

    PBRT_CPU_GPU inline void AddSample(const Point2i &pFilm, RGB &rgb, Float weight);

    inline void AddSplat(const Point2i &p, Float &wt, RGB &rgb);

    using TaggedPointer::TaggedPointer;

    static EstimatorHandle Create(const std::string &name, Bounds2i &pixelBounds, Allocator &alloc); 

    std::string ToString() const;
};

}  // namespace pbrt

#endif  // PBRT_BASE_ESTIMATOR_H
