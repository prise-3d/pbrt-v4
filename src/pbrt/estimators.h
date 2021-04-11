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

// P3D update : 17 as suggested by Jung
const int maxnbuffers = 17;

// Need to externalise PixelWindow declaration for RGBFilm
// P3D Updates
struct PixelBuffer {
    PixelBuffer() = default;

    // here we store the rgbSum of pixel
    double rgbSum[3] = {0., 0., 0.};
    AtomicDouble splatRGB[3];
    double weightSum = 0.;

    void Clear() {

        for (int i = 0; i < 3; i++) {
            rgbSum[i] = 0.;
            splatRGB[i] = 0.;
        }

        weightSum = 0.;
    }
};

// P3D Updates
struct PixelWindow {
    PixelWindow() = default;

    PixelBuffer buffers[maxnbuffers];
    VarianceEstimator<Float> varianceEstimator;

    int windowSize = maxnbuffers; // number of buffers clusters
    int index = 0; // keep track of index used
    int nsamples = 0; // keep track of real added nsamples;
    bool medianUse = false; // set if median is used or not
    double rgbSum[3] = {0., 0., 0.}; // store final expected sample using median
    double allrgbSum[3] = {0., 0., 0.}; // store final expected sample using median
    double squaredSum[3] = {0., 0., 0.}; // keep track of current squared sum values of all samples
    double mean[3] = {0., 0., 0.}; // keep track of current mean value with all samples (need to compute every `windowSize`)
    double currentStd = 0.;

    double weightSum = 0;
    double allWeightSum = 0;
    AtomicDouble splatRGB[3];

    /**
    *   Update current std value based on last encountered samples
    *   Update the current median expected value
    *   Clear each temp `windowSize` buffers for median computation
    */
    void update() {

        // compute here the median of the current buffers
        Float currentWeight = 0.;

        if (windowSize == 1) {
            // store channel information
            for (int i = 0; i < 3; i++) {
                currentWeight += buffers[0].weightSum;
                rgbSum[i] += buffers[0].rgbSum[i];
                splatRGB[i] = splatRGB[i] + buffers[0].splatRGB[i];
            }
        } 
        else
        {
            // based on channel numbers
            for (int i = 0; i < 3; i++) {

                // store channel information
                pstd::vector<Float> cvalues;
                pstd::vector<Float> weightsSum;
                pstd::vector<double> csplats;

                for (int j = 0; j < windowSize; j++) {
                    cvalues.push_back(buffers[j].rgbSum[i]);
                    // per channel management (but weight can be different depending of median buffer)
                    weightsSum.push_back(buffers[j].weightSum);
                    csplats.push_back(buffers[j].splatRGB[i]);             
                }

                // temp storage in order to sort values
                pstd::vector<Float> means(cvalues);
                pstd::vector<int> sortedIndices = means.sort();

                Float medianWeight, median = 0.;
                double medianSplat = 0;

                // compute median from means
                // find associated weightsum index and use it
                // Classical MON
                // if (windowSize % 2 == 1){
                unsigned unsortedIndex = sortedIndices[int(windowSize/2)];

                median = cvalues[unsortedIndex];
                medianWeight = weightsSum[unsortedIndex];
                medianSplat = csplats[unsortedIndex];
                // }
                // else{
                //     int k_mean = int(windowSize/2);
                //     unsigned firstIndex = sortedIndices[k_mean - 1];
                //     unsigned secondIndex = sortedIndices[k_mean];

                //     median = (cvalues[firstIndex] + cvalues[secondIndex]) / 2;
                //     medianWeight = (weightsSum[firstIndex] + weightsSum[secondIndex]) / 2;
                //     medianSplat = (csplats[firstIndex]  + csplats[secondIndex]) / 2;
                // }

                // store channel information
                currentWeight += medianWeight;
                rgbSum[i] += median;
                splatRGB[i] = splatRGB[i] + medianSplat;
            }
        }

        weightSum += (currentWeight / 3);

        // clear now the buffers data for the sliding window
        for (int i = 0; i < windowSize; i++) {
            buffers[i].Clear();
        }
    }
    // update the enable window size of the estimator dynamically
    // after computing stdScene, this method will be call to adapt windowSize in Film
    void updateSize(const Float &stdScene) {

        // std::cout << "Start updating window size" << std::endl;
        // Float stdSum = 0.;

        // // std::cout << "Update current pixelWindow (" << windowSize << ")" << std::endl;

        // // compute current mean and std
        // for (int i = 0; i < 3; i++) {
        //     mean[i] = allrgbSum[i]  / nsamples;
        //     stdSum += (squaredSum[i] / nsamples) - (mean[i] * mean[i]);
        // }

        // // divide per number of chanels and get current std
        // currentStd = std::sqrt(stdSum / 3);

        // Float stdRatio = currentStd / stdScene; 

        // if (stdRatio < 1.) {
        //     windowSize = 1;
        // } else {
        //     windowSize = 2 * (std::floor(std::log2(stdRatio)) + 1) + 1;
        // }

        // // if (windowSize < 0) {
        // //     std::cout << currentStd << " vs " << stdScene << std::endl;
        // //     std::cout << "Ratio is " << stdRatio << " => wsize " << windowSize << std::endl;
        // // }
        // // // TODO : check why windowSize is sometimes negative...
        // if (windowSize > maxnbuffers || windowSize < 1)
        //     windowSize = maxnbuffers;

        windowSize = 3;
    }
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
class JungEstimator : public Estimator {

    public:

        JungEstimator(const std::string &name) : Estimator(name) {}; 

        PBRT_CPU_GPU
        void Estimate(const PixelWindow &pixelWindow, RGB &rgb, Float &weightSum, AtomicDouble* splatRGB) const;

};

}  // namespace pbrt

#endif  // PBRT_ESTIMATOR_H