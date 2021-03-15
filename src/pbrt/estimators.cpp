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

    if (name == "mean")
        estimator = std::make_unique<MeanEstimator>(name);
    else if (name == "mon")
        estimator = std::make_unique<MONEstimator>(name);
    else if (name == "amon")
        estimator = std::make_unique<AlphaMONEstimator>(name);
    else if (name == "gini_mon")
        estimator = std::make_unique<GiniMONEstimator>(name);
    else if (name == "gini_binary_mon")
        estimator = std::make_unique<GiniBinaryMONEstimator>(name);
    else if (name == "gini_partial_mon")
        estimator = std::make_unique<GiniPartialMONEstimator>(name);
    else if (name == "pakmon")
        estimator = std::make_unique<PakMONEstimator>(name);
    else if (name == "mean_or_mon")
        estimator = std::make_unique<MeanOrMONEstimator>(name);
    else {
        printf("%s: estimator type unknown. Use of default: mean", name.c_str());
        estimator = std::make_unique<MeanEstimator>(name);
    }

    if (!estimator)
        printf("%s: unable to create estimator.", name.c_str());

    return estimator;
}

std::string Estimator::ToString() const {
    return name + "Estimator";
}

void MeanEstimator::Estimate(const PixelWindow &pixelWindow, RGB &rgb, Float &weightSum, AtomicDouble* splatRGB) const
{
    
    weightSum = 0.;

    // get each weightSum of pixelMoN
    for (int j = 0; j < pixelWindow.windowSize; j++) {
        weightSum += pixelWindow.buffers[j].weightSum;
    }

    // based on channel numbers
    for (int i = 0; i < 3; i++) {

        // loop over pixels (used as means storage) for computing real channel value
        rgb[i] = 0.;
        splatRGB[i] = 0.;

        for (int j = 0; j < pixelWindow.windowSize; j++) {
            rgb[i] += pixelWindow.buffers[j].rgbSum[i];
            splatRGB[i] = splatRGB[i] + pixelWindow.buffers[j].splatRGB[i];
        }
    }
};

void MONEstimator::Estimate(const PixelWindow &pixelWindow, RGB &rgb, Float &weightSum, AtomicDouble* splatRGB) const
{
    
    weightSum = 0.;

    // based on channel numbers
    for (int i = 0; i < 3; i++) {

        // store channel information
        pstd::vector<Float> cvalues;
        pstd::vector<Float> weightsSum;
        pstd::vector<double> csplats;

        for (int j = 0; j < pixelWindow.windowSize; j++) {
            cvalues.push_back(pixelWindow.buffers[j].rgbSum[i]);
            // per channel management (but weight can be different depending of median buffer)
            weightsSum.push_back(pixelWindow.buffers[j].weightSum);
            csplats.push_back(pixelWindow.buffers[j].splatRGB[i]);
        }

        // temp storage in order to sort values
        pstd::vector<Float> means(cvalues);
        pstd::vector<int> sortedIndices = means.sort();

        Float monWeight, monMean = 0.;
        double monSplat = 0;

        // compute median from means
        // find associated weightsum index and use it
        // Classical MON
        if (nbuffers % 2 == 1){
            unsigned unsortedIndex = sortedIndices[int(nbuffers/2)];

            monMean = cvalues[unsortedIndex];
            monWeight = weightsSum[unsortedIndex];
            monSplat = csplats[unsortedIndex];
        }
        else{
            int k_mean = int(nbuffers/2);
            unsigned firstIndex = sortedIndices[k_mean - 1];
            unsigned secondIndex = sortedIndices[k_mean];

            monMean = (cvalues[firstIndex] + cvalues[secondIndex]) / 2;
            monWeight = (weightsSum[firstIndex] + weightsSum[secondIndex]) / 2;
            monSplat = (csplats[firstIndex]  + csplats[secondIndex]) / 2;
        }

        // store channel information
        weightSum += monWeight;
        rgb[i] = monMean;
        splatRGB[i] = monSplat;
    }

    // divide per number of channel the weightSum
    weightSum /= 3;
};

void PakMONEstimator::Estimate(const PixelWindow &pixelWindow, RGB &rgb, Float &weightSum, AtomicDouble* splatRGB) const
{
    weightSum = 0.;

    // based on channel numbers
    for (int i = 0; i < 3; i++) {

        // store channel information
        pstd::vector<Float> cvalues;
        pstd::vector<Float> weightsSum;
        pstd::vector<double> csplats;

        for (int j = 0; j < pixelWindow.windowSize; j++) {
            cvalues.push_back(pixelWindow.buffers[j].rgbSum[i]);
            // per channel management (but weight can be different depending of median buffer)
            weightsSum.push_back(pixelWindow.buffers[j].weightSum);
            csplats.push_back(pixelWindow.buffers[j].splatRGB[i]);
        }

        // temp storage in order to sort values
        pstd::vector<Float> means(cvalues);
        pstd::vector<int> sortedIndices = means.sort();

        // sum storage
        Float meansSum = 0;

        // PakMON expected output
        Float weight, mean = 0.;
        double csplat = 0;

        // by default classical MON values
        if (nbuffers % 2 == 1){
            unsigned unsortedIndex = sortedIndices[int(nbuffers/2)];

            mean = cvalues[unsortedIndex];
            weight = weightsSum[unsortedIndex];
            csplat = csplats[unsortedIndex];
        }
        else{
            int k_mean = int(nbuffers/2);
            unsigned firstIndex = sortedIndices[k_mean - 1];
            unsigned secondIndex = sortedIndices[k_mean];

            mean = (cvalues[firstIndex] + cvalues[secondIndex]) / 2;
            weight = (weightsSum[firstIndex] + weightsSum[secondIndex]) / 2;
            csplat = (csplats[firstIndex]  + csplats[secondIndex]) / 2;
        }

        for (int j = 0; j < cvalues.size(); j++)
            meansSum += cvalues[j];

        Float currentMean = meansSum / cvalues.size();
            
        // compute variance distance evolution
        pstd::vector<Float> distances;

        for (int j = 2; j < means.size(); j++) {
            
            // use of sorted means in order to compute variance evolution by step 1
            // compute variance of each elements
            Float var = 0;
            
            // use of previously sorted means
            for(int k = 0; k < j; k++)
            {
                var += (means[k] - currentMean) * (means[k] - currentMean);
            }
            var /= (j + 1);

            // add new 
            distances.push_back(var);
        }

        // use of variance evolution and compute entropy
        Float distancesEntropy = getEntropy(distances);

        // Computation of PakMON using \alpha and \rho value
        unsigned middleIndex = int(nbuffers / 2);

        // alpha and rho automatically set value
        Float alpha = distancesEntropy;

        if (alpha < 0.000000001) {
            alpha = 0.000000001;
        }
        int rho = (int)(middleIndex * distancesEntropy) - (int)(nbuffers * 0.15); // try using avoid 30% (total) of current nbuffers

        unsigned lowerIndex = 0;
        unsigned higherIndex = 0;
    
        // get current lower and higher index 
        if (nbuffers % 2 == 0) {
            
            lowerIndex = middleIndex - 1;
            higherIndex = middleIndex;
            
        } else {
            lowerIndex = middleIndex;
            higherIndex = middleIndex;
        }

        // use of sorted means and relative sorted indices
        for (int j = 1; j < rho + 1; j++) {

            // current neighbor multiple factor
            Float multFactor = pow(alpha, j);

            // add left and right neighbor contribution
            mean += means[lowerIndex - j] * multFactor;
            mean += means[higherIndex + j] * multFactor;
            
            // weighting contribution to take in account
            // use of this index to retrieve the associated weightsSum
            weight += weightsSum[sortedIndices[lowerIndex]] * multFactor;
            weight += weightsSum[sortedIndices[higherIndex]] * multFactor;

            csplat += csplats[sortedIndices[lowerIndex]] * multFactor;
            csplat += csplats[sortedIndices[higherIndex]] * multFactor;
        }

        // store channel information
        weightSum += weight;
        rgb[i] = mean;
        splatRGB[i] = csplat;
    }

    // divide per number of channel the weightSum
    weightSum /= 3;
};

Float PakMONEstimator::getEntropy(pstd::vector<Float> values) const {

    // computation of squared values
    Float sumEigenValues = 0;
    pstd::vector<Float> eigenValues(values.size());

    for (int i = 0; i < values.size(); i++) {
        Float sum = values[i] * values[i];
        eigenValues[i] = sum;
        sumEigenValues += sum;
    }

    // normalization the squared values
    pstd::vector<Float> v(values.size());

    for (int i = 0; i < values.size(); i++) {
        // add of epsilon value
        v[i] = eigenValues[i] / (sumEigenValues + 0.000000000001);
    }

    // computation of entropy
    Float entropy = 0;

    for (int i = 0; i < values.size(); i++) {
        if (v[i] > 0) {
            entropy += v[i] * log(v[i]);
        }
    }

    entropy *= -1;

    entropy /= log(values.size());

    return entropy;
};


void AlphaMONEstimator::Estimate(const PixelWindow &pixelWindow, RGB &rgb, Float &weightSum, AtomicDouble* splatRGB) const
{
   this->Estimate(pixelWindow, rgb, weightSum, splatRGB, 0.5); // default use of 0.5 confidence
};

void AlphaMONEstimator::Estimate(const PixelWindow &pixelWindow, RGB &rgb, Float &weightSum, AtomicDouble* splatRGB, Float alpha) const
{
    weightSum = 0.;

    // based on channel numbers
    for (int i = 0; i < 3; i++) {

        // store channel information
        pstd::vector<Float> cvalues;
        pstd::vector<Float> weightsSum;
        pstd::vector<double> csplats;

        for (int j = 0; j < pixelWindow.windowSize; j++) {
            cvalues.push_back(pixelWindow.buffers[j].rgbSum[i]);
            // per channel management (but weight can be different depending of median buffer)
            weightsSum.push_back(pixelWindow.buffers[j].weightSum);
            csplats.push_back(pixelWindow.buffers[j].splatRGB[i]);
        }

        // temp storage in order to sort values
        pstd::vector<Float> means(cvalues);
        pstd::vector<int> sortedIndices = means.sort();

        // sum storage
        Float meansSum = 0;

        // PakMON expected output
        Float weight, mean = 0.;
        double csplat = 0;

        // by default classical MON values
        if (nbuffers % 2 == 1){
            unsigned unsortedIndex = sortedIndices[int(nbuffers/2)];

            mean = cvalues[unsortedIndex];
            weight = weightsSum[unsortedIndex];
            csplat = csplats[unsortedIndex];
        }
        else{
            int k_mean = int(nbuffers/2);
            unsigned firstIndex = sortedIndices[k_mean - 1];
            unsigned secondIndex = sortedIndices[k_mean];

            mean = (cvalues[firstIndex] + cvalues[secondIndex]) / 2;
            weight = (weightsSum[firstIndex] + weightsSum[secondIndex]) / 2;
            csplat = (csplats[firstIndex]  + csplats[secondIndex]) / 2;
        }

        // Computation of PakMON using \alpha and \rho value
        unsigned middleIndex = int(nbuffers / 2);

        unsigned lowerIndex = 0;
        unsigned higherIndex = 0;
        unsigned until = int(alpha * middleIndex); // TODO check if correct way to get until
    
        // get current lower and higher index 
        if (nbuffers % 2 == 0) {
            
            lowerIndex = middleIndex - 1;
            higherIndex = middleIndex;
            until = middleIndex - 1;
            
        } else {
            lowerIndex = middleIndex;
            higherIndex = middleIndex;
        }

        // use of sorted means and relative sorted indices
        for (int j = 1; j < until + 1; j++) {

            // current neighbor multiple factor
            // Float multFactor = pow(alpha, j);

            // add left and right neighbor contribution
            mean += means[lowerIndex - j];
            mean += means[higherIndex + j];
            
            // weighting contribution to take in account
            // use of this index to retrieve the associated weightsSum
            weight += weightsSum[sortedIndices[lowerIndex]];
            weight += weightsSum[sortedIndices[higherIndex]];

            csplat += csplats[sortedIndices[lowerIndex]];
            csplat += csplats[sortedIndices[higherIndex]];
        }

        // store channel information
        weightSum += weight;
        rgb[i] = mean;
        splatRGB[i] = csplat;
    }

    // divide per number of channel the weightSum
    weightSum /= 3;
};

void MeanOrMONEstimator::Estimate(const PixelWindow &pixelWindow, RGB &rgb, Float &weightSum, AtomicDouble* splatRGB) const
{
    // Check use of mon or use of mean based on IC
    Float meanICSum = 0.;
    Float monICSum = 0.;

    // based on channel numbers
    for (int i = 0; i < 3; i++) {

        // loop over pixels (used as means storage) for computing real channel value
        pstd::vector<Float> cvalues;
        pstd::vector<Float> csquared;
        pstd::vector<Float> weightsSum;
        pstd::vector<double> csplats;

        for (int j = 0; j < pixelWindow.windowSize; j++) {
            cvalues.push_back(pixelWindow.buffers[j].rgbSum[i]);
            csquared.push_back(pixelWindow.buffers[j].squaredSum[i]);
            weightsSum.push_back(pixelWindow.buffers[j].weightSum);
            csplats.push_back(pixelWindow.buffers[j].splatRGB[i]);
        }

        pstd::vector<Float> means(cvalues); // copy of current channel values

        // depending of confidence interval of predictor, compute the required value
        // need to compute mean value, hence all values are needed
        Float meanSquaredValues = 0.;
        
        Float meanWeight, meanMean = 0.;
        double meanSplat = 0;

        for (int j = 0; j < nbuffers; j++) {
            
            meanMean += means[j];
            meanSquaredValues += csquared[j];
            meanWeight += weightsSum[j];
            meanSplat += csplats[j];
        }

        /////////////////////////////////
        // confidence interval of Mean //
        /////////////////////////////////
        Float currentMeanMean = meanMean / meanWeight;

        // use of weight as number of samples
        Float meanStdValue = (meanSquaredValues / meanWeight) - (currentMeanMean * currentMeanMean);
        Float meanIC = (1.96 * meanStdValue) / std::sqrt(meanWeight);

        /////////////////////////////////
        // confidence interval of MON  //
        /////////////////////////////////

        // sort means vector and get sorted indices as output
        Float monSquaredValue = 0;

        Float monWeight, monMean = 0.;
        double monSplat = 0;

        pstd::vector<int> sortedIndices = means.sort();

        // compute median from means
        // find associated weightsum index and use it
        // Classical MON
        if (nbuffers % 2 == 1){
            unsigned unsortedIndex = sortedIndices[int(nbuffers/2)];

            monWeight = weightsSum[unsortedIndex];
            monMean = cvalues[unsortedIndex];
            monSquaredValue = csquared[unsortedIndex];
            monSplat = csplats[unsortedIndex];
        }
        else{
            int k_mean = int(nbuffers/2);
            unsigned firstIndex = sortedIndices[k_mean - 1];
            unsigned secondIndex = sortedIndices[k_mean];

            monWeight = (weightsSum[firstIndex] + weightsSum[secondIndex]) / 2;
            monMean = (cvalues[firstIndex] + cvalues[secondIndex]) / 2;
            monSquaredValue = (csquared[firstIndex] + csquared[secondIndex]) / 2;
            monSplat = (csplats[firstIndex]  + csplats[secondIndex]) / 2;
        }

        Float currentMonMean = monMean / monWeight;
        Float monStdValue = (monSquaredValue / monWeight) - (currentMonMean * currentMonMean);
        Float monIC = (1.96 * monStdValue) / std::sqrt(monWeight);

        meanICSum += meanIC;
        monICSum += monIC;
    }

    // check use of mean or mon
    if ((meanICSum / 3)  <= (monICSum / 3)) {
        meanEstimator->Estimate(pixelWindow, rgb, weightSum, splatRGB);
    } else {
        monEstimator->Estimate(pixelWindow, rgb, weightSum, splatRGB);
    }
};


Float GiniMONEstimator::getEntropy(pstd::vector<Float> values) const {

    // computation of squared values
    Float sumEigenValues = 0;
    pstd::vector<Float> eigenValues(values.size());

    for (int i = 0; i < values.size(); i++) {
        Float sum = values[i] * values[i];
        eigenValues[i] = sum;
        sumEigenValues += sum;
    }

    // normalization the squared values
    pstd::vector<Float> v(values.size());

    for (int i = 0; i < values.size(); i++) {
        // add of epsilon value
        v[i] = eigenValues[i] / (sumEigenValues + 0.000000000001);
    }

    // computation of entropy
    Float entropy = 0;

    for (int i = 0; i < values.size(); i++) {
        if (v[i] > 0) {
            entropy += v[i] * log(v[i]);
        }
    }

    entropy *= -1;

    entropy /= log(values.size());

    return entropy;
};

Float GiniMONEstimator::getGini(pstd::vector<Float> values) const {

    // get indices of array
    int n = values.size();
    Float arraySum = 0;
    Float indexArraySum = 0;

    Float minValue = Infinity;

    // get min value
    for (int i = 0; i < n; i++)
        if (values[i] < minValue)
            minValue = values[i];

    // need to sort obtained values
    std::sort(values.begin(), values.end());

    // avoid 0 value and store index
    for (int i = 0; i < n; i++) {

        // avoid negative value
        if (minValue < 0)
            values[i] -= minValue; 

        values[i] += 0.00000000001; // epsilon value
        arraySum += values[i];
        indexArraySum += (2 * (i + 1) - n - 1) * values[i];
    }

    return indexArraySum / (n * arraySum);
}

void GiniMONEstimator::Estimate(const PixelWindow &pixelWindow, RGB &rgb, Float &weightSum, AtomicDouble* splatRGB) const
{
    // for each channel number
    Float giniSum = 0;

    for (int i = 0; i < 3; i++) {

        pstd::vector<Float> cvalues;

        for (int j = 0; j < pixelWindow.windowSize; j++) {
            cvalues.push_back(pixelWindow.buffers[j].rgbSum[i] / pixelWindow.buffers[j].weightSum);
        }

        giniSum += this->getGini(cvalues);
    }

    Float giniMean = giniSum / 3.;

    // Need opposite of GiniMean
    alphaMoNEstimator->Estimate(pixelWindow, rgb, weightSum, splatRGB, 1 - giniMean);
};


void GiniBinaryMONEstimator::Estimate(const PixelWindow &pixelWindow, RGB &rgb, Float &weightSum, AtomicDouble* splatRGB) const
{
    // for each channel number
    Float giniSum = 0;

    for (int i = 0; i < 3; i++) {

        pstd::vector<Float> cvalues;

        for (int j = 0; j < pixelWindow.windowSize; j++) {
            cvalues.push_back(pixelWindow.buffers[j].rgbSum[i] / pixelWindow.buffers[j].weightSum);
        }

        giniSum += this->getGini(cvalues);
    }

    Float giniMean = giniSum / 3.;

    if (giniMean < 0.25)
        alphaMoNEstimator->Estimate(pixelWindow, rgb, weightSum, splatRGB, 1.);
    else
        alphaMoNEstimator->Estimate(pixelWindow, rgb, weightSum, splatRGB, 0.);
};


void GiniPartialMONEstimator::Estimate(const PixelWindow &pixelWindow, RGB &rgb, Float &weightSum, AtomicDouble* splatRGB) const
{
    // for each channel number
    Float giniSum = 0;

    for (int i = 0; i < 3; i++) {

        pstd::vector<Float> cvalues;

        for (int j = 0; j < pixelWindow.windowSize; j++) {
            cvalues.push_back(pixelWindow.buffers[j].rgbSum[i] / pixelWindow.buffers[j].weightSum);
        }

        giniSum += this->getGini(cvalues);
    }

    Float giniMean = giniSum / 3.;

    // totally be confident with Mean, otherwise use of Gini indicator (check also 0.5)
    if (giniMean < 0.25)
        alphaMoNEstimator->Estimate(pixelWindow, rgb, weightSum, splatRGB, 1.);
    else
        alphaMoNEstimator->Estimate(pixelWindow, rgb, weightSum, splatRGB, 1 - giniMean);
};

}