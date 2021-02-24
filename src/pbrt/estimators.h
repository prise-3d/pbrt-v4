// pbrt is Copyright(c) 1998-2020 Matt Pharr, Wenzel Jakob, and Greg Humphreys.
// The pbrt source code is licensed under the Apache License, Version 2.0.
// SPDX: Apache-2.0

#ifndef PBRT_ESTIMATOR_H
#define PBRT_ESTIMATOR_H

// PhysLight code contributed by Anders Langlands and Luca Fascione
// Copyright (c) 2020, Weta Digital, Ltd.
// SPDX-License-Identifier: Apache-2.0

// P3D updates

#include <pbrt/pbrt.h>
#include <pbrt/util/pstd.h>
#include <pbrt/util/vecmath.h>
#include <pbrt/util/math.h>
#include <pbrt/util/color.h>

#include <pbrt/base/bxdf.h>
#include <pbrt/base/camera.h>
#include <pbrt/estimators.h>
#include <pbrt/bsdf.h>
#include <pbrt/util/colorspace.h>
#include <pbrt/util/parallel.h>
#include <pbrt/util/sampling.h>
#include <pbrt/util/spectrum.h>
#include <pbrt/util/transform.h>

#include <atomic>
#include <string>
#include <vector>
#include <math.h>


namespace pbrt {

// P3D update  mon parameter (if value of 1, no kmon use, hence classical mean)
const int nbuffers = 12;

// Need to externalise PixelWindow declaration for RGBFilm
// P3D Updates
struct PixelBuffer {
    PixelBuffer() = default;

    double rgbSum[3] = {0., 0., 0.};
    double squaredSum[3] = {0., 0., 0.};
    AtomicDouble splatRGB[3];
    double weightSum = 0.;
};

// P3D Updates
struct PixelWindow {
    PixelWindow() = default;

    PixelBuffer buffers[nbuffers];
    VarianceEstimator<Float> varianceEstimator;

    int windowSize = nbuffers; // number of buffers clusters
    int index = 0; // keep track of index used
    bool filled = false;
};

// Base Estimator class
// Methods:
// - Create: enable to create an instance of Estimator
// - AddSample: add sample inside image buffer (with use of window)
// - Estimate: estimate final pixel value using custom estimator
// - AddSplat: add splat inside each pixel buffer
class Estimator {

    public:

        ~Estimator() {};
        
        static std::unique_ptr<Estimator> Create(const std::string &name, Bounds2i &pixelBounds, Allocator &alloc); 

        PBRT_CPU_GPU
        virtual void Estimate(const Point2i &pFilm, RGB &rgb, Float &weightSum, AtomicDouble* splatRGB) const = 0;

        PBRT_CPU_GPU
        virtual void AddSample(const Point2i &pFilm, RGB &rgb, Float weight);

        void AddSplat(const Point2i &p, Float &wt, RGB &rgb) {
            PixelWindow &pixel = pixels[p];
            for (int i = 0; i < 3; ++i)
                pixel.buffers[pixel.index].splatRGB[i].Add(wt * rgb[i]); // add to current index
        }

        std::string ToString() const;

    protected:

        Estimator(const std::string &name, Bounds2i &pixelBounds, Allocator &alloc) 
            : name(name), pixelBounds(pixelBounds), pixels(pixelBounds, alloc) {};
       
        std::string name;
        Bounds2i pixelBounds;

        Array2D<PixelWindow> pixels;
};

// Mean Estimator class
class MeanEstimator : public Estimator {

    public:

        MeanEstimator(const std::string &name, Bounds2i &pixelBounds, Allocator &alloc) 
            : Estimator(name, pixelBounds, alloc) {}; 

        PBRT_CPU_GPU
        void Estimate(const Point2i &pFilm, RGB &rgb, Float &weightSum, AtomicDouble* splatRGB) const;
};

// MON Estimator class
// Median of meaNs: use of median value from available mean buffers
class MONEstimator : public Estimator {

    public:

        MONEstimator(const std::string &name, Bounds2i &pixelBounds, Allocator &alloc) 
            : Estimator(name, pixelBounds, alloc) {};

        PBRT_CPU_GPU
        void Estimate(const Point2i &pFilm, RGB &rgb, Float &weightSum, AtomicDouble* splatRGB) const;

};

// AlphaMON Estimator class
// Median of meaNs: use of median value from available mean buffers
// use of an alpha criterion for convergence
class AlphaMONEstimator : public Estimator {

    public:

        AlphaMONEstimator(const std::string &name, Bounds2i &pixelBounds, Allocator &alloc, Float alpha) 
            : Estimator(name, pixelBounds, alloc), alpha(alpha) {}; 

        PBRT_CPU_GPU
        void Estimate(const Point2i &pFilm, RGB &rgb, Float &weightSum, AtomicDouble* splatRGB) const;

        Float getAlpha() {
            return alpha;
        };

        void setAlpha(Float alpha) {
            this->alpha = alpha;
        };

    private:
        Float alpha; // confidence criterion in [0, 1]
};

// AlphaMON Estimator class
// Median of meaNs: use of median value from available mean buffers
// use of an alpha criterion for convergence
class AutoAlphaMONEstimator : public Estimator {

    public:

        AutoAlphaMONEstimator(const std::string &name, Bounds2i &pixelBounds, Allocator &alloc, unsigned n) 
            : Estimator(name, pixelBounds, alloc), n(n) {

            meanEstimator = std::make_unique<MeanEstimator>("mean", pixelBounds, alloc);
            monEstimator = std::make_unique<MONEstimator>("mon", pixelBounds, alloc);

            for (int i = 0; i < n + 1; i++) {

                // determine alpha parameter based
                Float alpha = (1.0 / Float(n)) * (i);

                std::unique_ptr<AlphaMONEstimator> amonSpecific = std::make_unique<AlphaMONEstimator>("amon", pixelBounds, alloc, alpha);

                alphaMonEstimators.push_back(std::move(amonSpecific)); 
            }
        }; 

        PBRT_CPU_GPU
        void Estimate(const Point2i &pFilm, RGB &rgb, Float &weightSum, AtomicDouble* splatRGB) const;

    private:
        unsigned n; // number of step when searching optimal alpha value
        std::vector<std::unique_ptr<AlphaMONEstimator>> alphaMonEstimators;
        std::unique_ptr<MeanEstimator> meanEstimator;
        std::unique_ptr<MONEstimator> monEstimator;
};

// PakMON Estimation class
// Based of MON Estimator but with confidence criteria of median's neighborhood buffers
class PakMONEstimator : public Estimator {

    public:

        PakMONEstimator(const std::string &name, Bounds2i &pixelBounds, Allocator &alloc) : Estimator(name, pixelBounds, alloc) {}; 

        PBRT_CPU_GPU
        void Estimate(const Point2i &pFilm, RGB &rgb, Float &weightSum, AtomicDouble* splatRGB) const;

    private:
        PBRT_CPU_GPU
        Float getEntropy(pstd::vector<Float> values) const;

};

// Mean or MON estimator
// Based on Confidence Interval criteria of pixel
// Final estimator is chosen (mean or MON)
class MeanOrMONEstimator : public Estimator {

    public:

        MeanOrMONEstimator(const std::string &name, Bounds2i &pixelBounds, Allocator &alloc) : Estimator(name, pixelBounds, alloc) { 
            meanEstimator = std::make_unique<MeanEstimator>("mean", pixelBounds, alloc);
            monEstimator = std::make_unique<MONEstimator>("mon", pixelBounds, alloc);
        }; 

        PBRT_CPU_GPU
        void Estimate(const Point2i &pFilm, RGB &rgb, Float &weightSum, AtomicDouble* splatRGB) const;

    private:
        std::unique_ptr<MeanEstimator> meanEstimator;
        std::unique_ptr<MONEstimator> monEstimator;
};

}  // namespace pbrt

#endif  // PBRT_ESTIMATOR_H