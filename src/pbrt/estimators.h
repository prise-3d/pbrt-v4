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
#include <pbrt/options.h>

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

// P3D update : 6 as suggested by Zirr
const int maxnbuffers = 6;

// Need to externalise PixelWindow declaration for RGBFilm
// P3D Updates
struct PixelBuffer {
    PixelBuffer() = default;

    // here we store the rgbSum of pixel
    double rgbSum[3] = {0., 0., 0.};
    AtomicDouble splatRGB[3];
    double weightSum[3] = {0., 0., 0.};

    void Clear() {

        for (int i = 0; i < 3; i++) {
            rgbSum[i] = 0.;
            splatRGB[i] = 0.;
            weightSum[i] = 0.;
        }
    }
};

// P3D Updates
struct PixelWindow {
    PixelWindow() = default;

    PixelBuffer buffers[maxnbuffers];
    VarianceEstimator<Float> varianceEstimator;

    int windowSize = maxnbuffers; // number of buffers clusters
    int nsamples = 0; // keep track of real added nsamples;
    double rgbSum[3] = {0., 0., 0.}; 

    double weightSum = 0;
    AtomicDouble splatRGB[3];
    double cascadeStart = 1;
    double cascadeBase = std::pow(262144, 1. / (double)maxnbuffers); // change using b = j^sqrt(S_{max}) where S_{max} = 262144
    int kmin = 1;
    int k = 32;
    int N = *Options->pixelSamples * *Options->nimages;
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

// Jung implemented estimator Estimator class
class ZirrEstimator : public Estimator {

    public:

        ZirrEstimator(const std::string &name) : Estimator(name) {}; 

        PBRT_CPU_GPU
        void Estimate(const PixelWindow &pixelWindow, RGB &rgb, Float &weightSum, AtomicDouble* splatRGB) const;

};

}  // namespace pbrt

#endif  // PBRT_ESTIMATOR_H