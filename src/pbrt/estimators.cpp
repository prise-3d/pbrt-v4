// pbrt is Copyright(c) 1998-2020 Matt Pharr, Wenzel Jakob, and Greg Humphreys.
// The pbrt source code is licensed under the Apache License, Version 2.0.
// SPDX: Apache-2.0

// PhysLight code contributed by Anders Langlands and Luca Fascione
// Copyright (c) 2020, Weta Digital, Ltd.
// SPDX-License-Identifier: Apache-2.0

// P3D updates

#include <pbrt/estimators.h>

#include <pbrt/pbrt.h>
#include <pbrt/film.h>
#include <pbrt/util/pstd.h>
#include <pbrt/util/vecmath.h>
#include <pbrt/util/math.h>
#include <pbrt/util/color.h>

namespace pbrt {

std::unique_ptr<Estimator> Estimator::Create(const std::string &name) {

    std::unique_ptr<Estimator> estimator;

    // TODO: Later use of paramset maybe..
    // if (name == "mean")
    //     estimator = std::make_unique<MeanEstimator>(name);
    // else if (name == "mon")
    //     estimator = std::make_unique<MONEstimator>(name);
    if (name == "djung")
        estimator = std::make_unique<JungEstimator>(name);
    // else if (name == "amon")
    //     estimator = std::make_unique<AlphaMONEstimator>(name);
    // else if (name == "admon")
    //     estimator = std::make_unique<AlphaDistMONEstimator>(name);
    // else if (name == "gini-mon")
    //     estimator = std::make_unique<GiniMONEstimator>(name);
    // else if (name == "gini-binary-mon")
    //     estimator = std::make_unique<GiniBinaryMONEstimator>(name);
    // else if (name == "gini-partial-mon")
    //     estimator = std::make_unique<GiniPartialMONEstimator>(name);
    // else if (name == "gini-dmon")
    //     estimator = std::make_unique<GiniDistMONEstimator>(name);
    // else if (name == "gini-partial-dmon")
    //     estimator = std::make_unique<GiniDistPartialMONEstimator>(name);
    // else if (name == "pakmon")
    //     estimator = std::make_unique<PakMONEstimator>(name);
    // else if (name == "mean-or-mon")
    //     estimator = std::make_unique<MeanOrMONEstimator>(name);
    // else if (name == "abmm")
    //     estimator = std::make_unique<ABMMEstimator>(name);
    // else if (name == "gabmm")
    //     estimator = std::make_unique<GABMMEstimator>(name);
    else {
        printf("%s: estimator type unknown. Use of default: jung", name.c_str());
        estimator = std::make_unique<JungEstimator>(name);
    }

    if (!estimator)
        printf("%s: unable to create estimator.", name.c_str());

    return estimator;
}

std::string Estimator::ToString() const {
    return name + "Estimator";
}

void JungEstimator::Estimate(const PixelWindow &pixelWindow, RGB &rgb, Float &weightSum, AtomicDouble* splatRGB) const
{

    Float currentWeight = 0.;
    
    if (pixelWindow.windowSize == 1) {
        // store channel information
        for (int i = 0; i < 3; i++) {
            currentWeight += pixelWindow.buffers[0].weightSum;
            rgb[i] = pixelWindow.buffers[0].rgbSum[i];
            splatRGB[i] = Float(pixelWindow.buffers[0].splatRGB[i]);
        }
    } 
    else
    {
        // based on channel numbers
        for (int i = 0; i < 3; i++) {

            pstd::vector<double> cvalues;
            pstd::vector<double> csplats;
            pstd::vector<double> sortedValues;
            pstd::vector<int> indices;
            pstd::vector<double> weightsSum;
            
            // std::cout << "Current values are:" << std::endl;
            // store channel information in available temp buffer
            for (int j = 0; j < pixelWindow.windowSize; j++) {

                // std::cout << "-- Channel (" << i << "), b[" << j << "] : " << pixelWindow.buffers[j].rgbSum[i] << std::endl; 
                cvalues.push_back(pixelWindow.buffers[j].rgbSum[i]);
                sortedValues.push_back(pixelWindow.buffers[j].rgbSum[i]);
                // per channel management (but weight can be different depending of median buffer)
                weightsSum.push_back(pixelWindow.buffers[j].weightSum);
                csplats.push_back(pixelWindow.buffers[j].splatRGB[i]);  
                indices.push_back(j);  
            }

            // Need now to sort data
            //std::sort(std::begin(pixelWindow.indices), std::begin(pixelWindow.indices) + pixelWindow.windowSize, 
            //    [&](int i, int j){return pixelWindow.cvalues[i] < pixelWindow.cvalues[j];
            //});

            int jj, kk, min, tempI;
            double temp;
            int n = pixelWindow.windowSize;

            for (jj = 0; jj < n - 1; jj++) {
                min = jj;
                for (kk = jj + 1; kk < n; kk++)
                    if (sortedValues[kk] < sortedValues[min])
                        min = kk;

                temp = sortedValues[jj];
                sortedValues[jj] = sortedValues[min];
                sortedValues[min] = temp;

                tempI = indices[jj];
                indices[jj] = indices[min];
                indices[min] = tempI;
            }

            //std::cout << "Sorted values are:" << std::endl;
            // for (int j = 0; j < pixelWindow.windowSize; j++) {
            //    std::cout << "-- Channel (" << i << "), b[" << pixelWindow.indices[j] << "] : " << pixelWindow.sortedValues[j] << std::endl;
            //}

            //std::cout << "Median index :" << pixelWindow.indices[int(pixelWindow.windowSize/2)] << std::endl;
            //std::cout << "Median value :" << pixelWindow.cvalues[int(pixelWindow.windowSize/2)] << std::endl;

            // need to find median value of pixelWindow.cavlues

            // Float medianWeight, median = 0.;
            // double medianSplat = 0;

            // // compute median from means
            // // find associated weightsum index and use it
            // // Classical MON
            // // if (windowSize % 2 == 1){
            unsigned unsortedIndex = indices[int(pixelWindow.windowSize/2)];

            // median = pixelWindow.cvalues[unsortedIndex];
            // medianWeight = pixelWindow.weightsSum[unsortedIndex];
            // medianSplat = pixelWindow.csplats[unsortedIndex];
            // // }
            // // else{
            // //     int k_mean = int(windowSize/2);
            // //     unsigned firstIndex = sortedIndices[k_mean - 1];
            // //     unsigned secondIndex = sortedIndices[k_mean];

            // //     median = (cvalues[firstIndex] + cvalues[secondIndex]) / 2;
            // //     medianWeight = (weightsSum[firstIndex] + weightsSum[secondIndex]) / 2;
            // //     medianSplat = (csplats[firstIndex]  + csplats[secondIndex]) / 2;
            // // }

            // // store channel information
            currentWeight += weightsSum[unsortedIndex];
            // pixelWindow.rgbSum[i] += pixelWindow.cvalues[unsortedIndex] * pixelWindow.windowSize;
            // pixelWindow.splatRGB[i] = pixelWindow.splatRGB[i] + pixelWindow.csplats[unsortedIndex];

            rgb[i] = cvalues[unsortedIndex];
            splatRGB[i] = Float(csplats[unsortedIndex]);
        }
    }

    weightSum = (currentWeight / 3);
};


}