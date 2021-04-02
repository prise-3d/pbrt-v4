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
    double cubicSum[3] = {0., 0., 0.};
    AtomicDouble splatRGB[3];
    double weightSum = 0.;

    void Clear() {

        for (int i = 0; i < 3; i++) {
            rgbSum[i] = 0.;
            squaredSum[i] = 0.;
            cubicSum[i] = 0.;
            splatRGB[i] = 0.;
        }

        weightSum = 0.;
    }
};

// P3D Updates
struct PixelWindow {
    PixelWindow() = default;

    PixelBuffer buffers[nbuffers];
    VarianceEstimator<Float> varianceEstimator;

    int windowSize = nbuffers; // number of buffers clusters
    int index = 0; // keep track of index used
    int nsamples = 0; // keep track of nsamples;
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

// approximated Bayesian Median of Means Estimator class
class aBMMEstimator : public Estimator {

    public:

        aBMMEstimator(const std::string &name) : Estimator(name) {}; 

        PBRT_CPU_GPU
        void Estimate(const PixelWindow &pixelWindow, RGB &rgb, Float &weightSum, AtomicDouble* splatRGB) const;

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

// AlphaMON Estimator class
// Median of meaNs: use of median value from available mean buffers
// use of an alpha criterion for convergence
class AlphaMONEstimator : public Estimator {

    public:

        AlphaMONEstimator(const std::string &name) : Estimator(name) {}; 

        PBRT_CPU_GPU
        void Estimate(const PixelWindow &pixelWindow, RGB &rgb, Float &weightSum, AtomicDouble* splatRGB) const;

        PBRT_CPU_GPU
        void Estimate(const PixelWindow &pixelWindow, RGB &rgb, Float &weightSum, AtomicDouble* splatRGB, Float alpha) const;
};

// AlphaDistMON Estimator class
// Median of meaNs: use of median value from available mean buffers
// use of an alpha criterion for convergence using whole package
class AlphaDistMONEstimator : public Estimator {

    public:

        AlphaDistMONEstimator(const std::string &name) : Estimator(name) {}; 

        PBRT_CPU_GPU
        void Estimate(const PixelWindow &pixelWindow, RGB &rgb, Float &weightSum, AtomicDouble* splatRGB) const;

        PBRT_CPU_GPU
        void Estimate(const PixelWindow &pixelWindow, RGB &rgb, Float &weightSum, AtomicDouble* splatRGB, Float alpha) const;
};

// GiniMON Estimator class
// Median of meaNs: use of median value from available mean buffers
// Use of Gini in order to well use \alpha criterion
class GiniDistMONEstimator : public Estimator {

    public:

        GiniDistMONEstimator(const std::string &name) : Estimator(name) {

            // default alpha value
            alphaDistMoNEstimator = std::make_unique<AlphaDistMONEstimator>("admon");
        }; 

        PBRT_CPU_GPU
        void Estimate(const PixelWindow &pixelWindow, RGB &rgb, Float &weightSum, AtomicDouble* splatRGB) const;

    protected:
        std::unique_ptr<AlphaDistMONEstimator> alphaDistMoNEstimator;

        PBRT_CPU_GPU
        Float getGini(pstd::vector<Float> values) const;
};

class GiniDistPartialMONEstimator : public GiniDistMONEstimator {

    public:

        GiniDistPartialMONEstimator(const std::string &name) : GiniDistMONEstimator(name) {};

        PBRT_CPU_GPU
        void Estimate(const PixelWindow &pixelWindow, RGB &rgb, Float &weightSum, AtomicDouble* splatRGB) const;
};

// GiniMON Estimator class
// Median of meaNs: use of median value from available mean buffers
// Use of Gini in order to well use \alpha criterion
class GiniMONEstimator : public Estimator {

    public:

        GiniMONEstimator(const std::string &name) : Estimator(name) {

            // default alpha value
            alphaMoNEstimator = std::make_unique<AlphaMONEstimator>("amon");
        }; 

        PBRT_CPU_GPU
        void Estimate(const PixelWindow &pixelWindow, RGB &rgb, Float &weightSum, AtomicDouble* splatRGB) const;

    protected:
        std::unique_ptr<AlphaMONEstimator> alphaMoNEstimator;

        PBRT_CPU_GPU
        Float getEntropy(pstd::vector<Float> values) const;

        PBRT_CPU_GPU
        Float getGini(pstd::vector<Float> values) const;
};

class GiniBinaryMONEstimator : public GiniMONEstimator {

    public:

        GiniBinaryMONEstimator(const std::string &name) : GiniMONEstimator(name) {};

        PBRT_CPU_GPU
        void Estimate(const PixelWindow &pixelWindow, RGB &rgb, Float &weightSum, AtomicDouble* splatRGB) const;
};

class GiniPartialMONEstimator : public GiniMONEstimator {

    public:

        GiniPartialMONEstimator(const std::string &name) : GiniMONEstimator(name) {};

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
            meanEstimator = std::make_unique<MeanEstimator>("mean");
            monEstimator = std::make_unique<MONEstimator>("mon");
        }; 

        PBRT_CPU_GPU
        void Estimate(const PixelWindow &pixelWindow, RGB &rgb, Float &weightSum, AtomicDouble* splatRGB) const;

    private:
        std::unique_ptr<MeanEstimator> meanEstimator;
        std::unique_ptr<MONEstimator> monEstimator;
};

}  // namespace pbrt

#endif  // PBRT_ESTIMATOR_H