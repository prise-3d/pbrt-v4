// pbrt is Copyright(c) 1998-2020 Matt Pharr, Wenzel Jakob, and Greg Humphreys.
// The pbrt source code is licensed under the Apache License, Version 2.0.
// SPDX: Apache-2.0

#ifndef PBRT_OPTIONS_H
#define PBRT_OPTIONS_H

#include <pbrt/pbrt.h>
#include <pbrt/util/log.h>
#include <pbrt/util/pstd.h>
#include <pbrt/util/vecmath.h>

#include <string>

namespace pbrt {

// RenderingCoordinateSystem Definition
enum class RenderingCoordinateSystem { Camera, CameraWorld, World };

// BasicOptions Definition
struct BasicOptions {
    int nThreads = 0;
    int seed = 0;
    bool quickRender = false;
    bool quiet = false;
    bool recordPixelStatistics = false;
    bool upgrade = false;
    bool disablePixelJitter = false, disableWavelengthJitter = false;
    bool forceDiffuse = false;
    bool useGPU = false;
    RenderingCoordinateSystem renderingSpace = RenderingCoordinateSystem::CameraWorld;
};

// PBRTOptions Definiton
struct PBRTOptions : BasicOptions {
    LogLevel logLevel = LogLevel::Error;
    pstd::optional<int> pixelSamples;
    pstd::optional<int> nimages = 1; // P3D updates: use of number of images to generate
    pstd::optional<int> startIndex = 0; // P3D updates: index of image to start with
    pstd::optional<int> ndigits = 6; // P3D updates: number of digits for file format
    pstd::optional<int> gpuDevice;
    pstd::optional<int> monk; // P3D update k mon parameter
    std::string imageFile;
    std::string folderName = "temp"; // P3D updates: use of temp or specific folder where to save the computed images
    std::string mseReferenceImage, mseReferenceOutput;
    std::string debugStart;
    std::string displayServer;
    pstd::optional<Bounds2f> cropWindow;
    pstd::optional<Bounds2i> pixelBounds;

    std::string ToString() const;
};

// Options Global Variable Declaration
extern PBRTOptions *Options;

#if defined(PBRT_BUILD_GPU_RENDERER) && defined(__CUDACC__)
extern __constant__ BasicOptions OptionsGPU;
#endif

PBRT_CPU_GPU inline const BasicOptions &GetOptions() {
#if defined(PBRT_IS_GPU_CODE)
    return OptionsGPU;
#else
    return *Options;
#endif
}

}  // namespace pbrt

#endif  // PBRT_OPTIONS_H
