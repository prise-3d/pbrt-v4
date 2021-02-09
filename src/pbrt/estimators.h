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
const int nbuffers = 11;

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
class Estimator {

    public:

        ~Estimator() {};
        
        static std::unique_ptr<Estimator> Create(const std::string &name); 

        PBRT_CPU_GPU
        virtual void Estimate(const PixelWindow &pixelWindow, RGB &rgb, Float &weightSum, AtomicDouble* splatRGB) const = 0;

        std::string ToString() const;

    protected:

        Estimator(const std::string &name) : name(name) {};
       
        std::string name;
};

// Mean Estimator class
class MeanEstimator : public Estimator {

    public:

        MeanEstimator(const std::string &name) : Estimator(name) {}; 

        PBRT_CPU_GPU
        void Estimate(const PixelWindow &pixelWindow, RGB &rgb, Float &weightSum, AtomicDouble* splatRGB) const;

};

// MON Estimator class
// Median of meaNs: use of median value from available mean buffers
class MONEstimator : public Estimator {

    public:

        MONEstimator(const std::string &name) : Estimator(name) {}; 

        PBRT_CPU_GPU
        void Estimate(const PixelWindow &pixelWindow, RGB &rgb, Float &weightSum, AtomicDouble* splatRGB) const;

};

// PakMON Estimation class
// Based of MON Estimator but with confidence criteria of median's neighborhood buffers
class PakMONEstimator : public Estimator {

    public:

        PakMONEstimator(const std::string &name) : Estimator(name) {}; 

        PBRT_CPU_GPU
        void Estimate(const PixelWindow &pixelWindow, RGB &rgb, Float &weightSum, AtomicDouble* splatRGB) const;

    private:
        PBRT_CPU_GPU
        Float getEntropy(pstd::vector<Float> values) const;

};

// Mean or MON estimator
// Based on Confidence Interval criteria of pixel
// Final estimator is chosen (mean or MON)
class MeanOrMONEstimator : public Estimator {

    public:

        MeanOrMONEstimator(const std::string &name) : Estimator(name) {
            meanEstimator = Estimator::Create("mean");
            monEstimator = Estimator::Create("mon");
        }; 

        PBRT_CPU_GPU
        void Estimate(const PixelWindow &pixelWindow, RGB &rgb, Float &weightSum, AtomicDouble* splatRGB) const;

    private:
        std::unique_ptr<Estimator> meanEstimator;
        std::unique_ptr<Estimator> monEstimator;
};

}  // namespace pbrt

#endif  // PBRT_ESTIMATOR_H