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
    if (name == "zirr")
        estimator = std::make_unique<ZirrEstimator>(name);
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
        printf("%s: estimator type unknown. Use of default: zirr", name.c_str());
        estimator = std::make_unique<ZirrEstimator>(name);
    }

    if (!estimator)
        printf("%s: unable to create estimator.", name.c_str());

    return estimator;
}

std::string Estimator::ToString() const {
    return name + "Estimator";
}

void ZirrEstimator::Estimate(const PixelWindow &pixelWindow, RGB &rgb, Float &weightSum, AtomicDouble* splatRGB) const
{
    
    weightSum = pixelWindow.weightSum;

    // based on channel numbers
    for (int i = 0; i < 3; i++) {

        // loop over pixels (used as means storage) for computing real channel value
        rgb[i] = pixelWindow.rgbSum[i];
        splatRGB[i] = Float(pixelWindow.splatRGB[i]);
    }
};


}