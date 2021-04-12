// pbrt is Copyright(c) 1998-2020 Matt Pharr, Wenzel Jakob, and Greg Humphreys.
// The pbrt source code is licensed under the Apache License, Version 2.0.
// SPDX: Apache-2.0

#ifndef PBRT_FILM_H
#define PBRT_FILM_H

// PhysLight code contributed by Anders Langlands and Luca Fascione
// Copyright (c) 2020, Weta Digital, Ltd.
// SPDX-License-Identifier: Apache-2.0

#include <pbrt/pbrt.h>

#include <pbrt/base/bxdf.h>
#include <pbrt/base/camera.h>
#include <pbrt/base/film.h>
#include <pbrt/estimators.h>
#include <pbrt/bsdf.h>
#include <pbrt/util/color.h>
#include <pbrt/util/colorspace.h>
#include <pbrt/util/parallel.h>
#include <pbrt/util/pstd.h>
#include <pbrt/util/sampling.h>
#include <pbrt/util/spectrum.h>
#include <pbrt/util/transform.h>
#include <pbrt/util/vecmath.h>
#include <pbrt/util/math.h>

#include <atomic>
#include <map>
#include <string>
#include <thread>
#include <vector>
#include <math.h>

namespace pbrt {

// PixelSensor Definition
class PixelSensor {
  public:
    // PixelSensor Public Methods
    static PixelSensor *Create(const ParameterDictionary &parameters,
                               const RGBColorSpace *colorSpace, Float exposureTime,
                               const FileLoc *loc, Allocator alloc);

    static PixelSensor *CreateDefault(Allocator alloc = {});

    PixelSensor(Spectrum r, Spectrum g, Spectrum b, const RGBColorSpace *outputColorSpace,
                Float wbTemp, Float imagingRatio, Allocator alloc)
        : r_bar(r, alloc), g_bar(g, alloc), b_bar(b, alloc), imagingRatio(imagingRatio) {
        // Compute XYZ from camera RGB matrix
        // Compute _rgbCamera_ values for training swatches
        DenselySampledSpectrum sensorIllum = Spectra::D(wbTemp, alloc);
        Float rgbCamera[nSwatchReflectances][3];
        for (int i = 0; i < nSwatchReflectances; ++i) {
            RGB rgb = ProjectReflectance<RGB>(swatchReflectances[i], &sensorIllum, &r_bar,
                                              &g_bar, &b_bar);
            for (int c = 0; c < 3; ++c)
                rgbCamera[i][c] = rgb[c];
        }

        // Compute _xyzOutput_ values for training swatches
        Float xyzOutput[24][3];
        Float sensorWhiteG = InnerProduct(&sensorIllum, &g_bar);
        Float sensorWhiteY = InnerProduct(&sensorIllum, &Spectra::Y());
        for (size_t i = 0; i < nSwatchReflectances; ++i) {
            Spectrum s = swatchReflectances[i];
            XYZ xyz =
                ProjectReflectance<XYZ>(s, &outputColorSpace->illuminant, &Spectra::X(),
                                        &Spectra::Y(), &Spectra::Z()) *
                (sensorWhiteY / sensorWhiteG);
            for (int c = 0; c < 3; ++c)
                xyzOutput[i][c] = xyz[c];
        }

        // Initialize _XYZFromSensorRGB_ using linear least squares
        pstd::optional<SquareMatrix<3>> m =
            LinearLeastSquares<3>(rgbCamera, xyzOutput, nSwatchReflectances);
        if (!m)
            ErrorExit("Sensor XYZ from RGB matrix could not be solved.");
        XYZFromSensorRGB = *m;
    }

    PixelSensor(const RGBColorSpace *outputColorSpace, Float wbTemp, Float imagingRatio,
                Allocator alloc)
        : r_bar(&Spectra::X(), alloc),
          g_bar(&Spectra::Y(), alloc),
          b_bar(&Spectra::Z(), alloc),
          imagingRatio(imagingRatio) {
        // Compute white balancing matrix for XYZ _PixelSensor_
        if (wbTemp != 0) {
            auto whiteIlluminant = Spectra::D(wbTemp, alloc);
            Point2f sourceWhite = SpectrumToXYZ(&whiteIlluminant).xy();
            Point2f targetWhite = outputColorSpace->w;
            XYZFromSensorRGB = WhiteBalance(sourceWhite, targetWhite);
        }
    }

    PBRT_CPU_GPU
    RGB ToSensorRGB(SampledSpectrum L, const SampledWavelengths &lambda) const {
        L = SafeDiv(L, lambda.PDF());
        return imagingRatio * RGB((r_bar.Sample(lambda) * L).Average(),
                                  (g_bar.Sample(lambda) * L).Average(),
                                  (b_bar.Sample(lambda) * L).Average());
    }

    // PixelSensor Public Members
    SquareMatrix<3> XYZFromSensorRGB;

  private:
    // PixelSensor Private Methods
    template <typename Triplet>
    static Triplet ProjectReflectance(Spectrum r, Spectrum illum, Spectrum b1,
                                      Spectrum b2, Spectrum b3);

    // PixelSensor Private Members
    DenselySampledSpectrum r_bar, g_bar, b_bar;
    Float imagingRatio;
    static constexpr int nSwatchReflectances = 24;
    static Spectrum swatchReflectances[nSwatchReflectances];
};

// PixelSensor Inline Methods
template <typename Triplet>
inline Triplet PixelSensor::ProjectReflectance(Spectrum refl, Spectrum illum, Spectrum b1,
                                               Spectrum b2, Spectrum b3) {
    Triplet result;
    Float g_integral = 0;
    for (Float lambda = Lambda_min; lambda <= Lambda_max; ++lambda) {
        g_integral += b2(lambda) * illum(lambda);
        result[0] += b1(lambda) * refl(lambda) * illum(lambda);
        result[1] += b2(lambda) * refl(lambda) * illum(lambda);
        result[2] += b3(lambda) * refl(lambda) * illum(lambda);
    }
    return result / g_integral;
}

// VisibleSurface Definition
class VisibleSurface {
  public:
    // VisibleSurface Public Methods
    PBRT_CPU_GPU
    VisibleSurface(const SurfaceInteraction &si, const CameraTransform &cameraTransform,
                   const SampledSpectrum &albedo, const SampledWavelengths &lambda);

    PBRT_CPU_GPU
    operator bool() const { return set; }

    VisibleSurface() = default;

    std::string ToString() const;

    // VisibleSurface Public Members
    bool set = false;
    Point3f p;
    Normal3f n, ns;
    Float time = 0;
    Float dzdx = 0, dzdy = 0;
    SampledSpectrum albedo;
};

// FilmBaseParameters Definition
struct FilmBaseParameters {
    FilmBaseParameters(const ParameterDictionary &parameters, Filter filter,
                       const PixelSensor *sensor, const FileLoc *loc);
    FilmBaseParameters(Point2i fullResolution, Bounds2i pixelBounds, Filter filter,
                       Float diagonal, const PixelSensor *sensor, std::string filename)
        : fullResolution(fullResolution),
          pixelBounds(pixelBounds),
          filter(filter),
          diagonal(diagonal),
          sensor(sensor),
          filename(filename) {}

    Point2i fullResolution;
    Bounds2i pixelBounds;
    Filter filter;
    Float diagonal;
    const PixelSensor *sensor;
    std::string filename;
};

// FilmBase Definition
class FilmBase {
  public:
    // FilmBase Public Methods
    FilmBase(FilmBaseParameters p)
        : fullResolution(p.fullResolution),
          pixelBounds(p.pixelBounds),
          filter(p.filter),
          diagonal(p.diagonal * .001f),
          sensor(p.sensor),
          filename(p.filename) {
        CHECK(!pixelBounds.IsEmpty());
        CHECK_GE(pixelBounds.pMin.x, 0);
        CHECK_LE(pixelBounds.pMax.x, fullResolution.x);
        CHECK_GE(pixelBounds.pMin.y, 0);
        CHECK_LE(pixelBounds.pMax.y, fullResolution.y);
        LOG_VERBOSE("Created film with full resolution %s, pixelBounds %s",
                    fullResolution, pixelBounds);
    }

    PBRT_CPU_GPU
    Point2i FullResolution() const { return fullResolution; }
    PBRT_CPU_GPU
    Bounds2i PixelBounds() const { return pixelBounds; }
    PBRT_CPU_GPU
    Float Diagonal() const { return diagonal; }
    PBRT_CPU_GPU
    Filter GetFilter() const { return filter; }
    PBRT_CPU_GPU
    const PixelSensor *GetPixelSensor() const { return sensor; }
    std::string GetFilename() const { return filename; }

    // P3D Updates
    void SetFilename(std::string filename) { this->filename = filename; }

    // P3D Updates 
    // compute std of current film
    void computeStd() { }


    PBRT_CPU_GPU
    SampledWavelengths SampleWavelengths(Float u) const {
        return SampledWavelengths::SampleXYZ(u);
    }

    PBRT_CPU_GPU
    Bounds2f SampleBounds() const;

    std::string BaseToString() const;

  protected:
    // FilmBase Protected Members
    Point2i fullResolution;
    Bounds2i pixelBounds;
    Filter filter;
    Float diagonal;
    const PixelSensor *sensor;
    std::string filename;
};

// RGBFilm Definition
class RGBFilm : public FilmBase {
  public:

    // RGBFilm Public Methods
    PBRT_CPU_GPU
    void Clear() { 
        
        ParallelFor2D(pixelBounds, [&](Point2i p) {
            for (int i = 0; i < pixels[p].windowSize; i++) {
                pixels[p].buffers[i].Clear();
            }
            pixels[p].varianceEstimator = VarianceEstimator<Float>();
        }); 
    };

    PBRT_CPU_GPU
    RGB GetPixelRGB(const Point2i &p, Float splatScale = 1) const {

        // P3D Updates
        const PixelWindow &pixelWindow = pixels[p];
        // update estimated value for rgbSum and weightSum

        //RGB rgb(pixel.rgbSum[0], pixel.rgbSum[1], pixel.rgbSum[2]);
        // Normalize _rgb_ with weight sum

        AtomicDouble splatRGB[3];
        RGB rgb(0, 0, 0);
        Float weightSum = 0.;

        estimator->Estimate(pixelWindow, rgb, weightSum, splatRGB);
        
        if (weightSum != 0)
            rgb /= weightSum;

        // Add splat value at pixel
        for (int c = 0; c < 3; ++c)
            rgb[c] += splatScale * splatRGB[c] / filterIntegral;

        // Convert _rgb_ to output RGB color space
        rgb = outputRGBFromSensorRGB * rgb;

        return rgb;
    }

    PBRT_CPU_GPU
    void AddSample(const Point2i &pFilm, SampledSpectrum L,
                   const SampledWavelengths &lambda, const VisibleSurface *,
                   Float weight) {
        // Convert sample radiance to _PixelSensor_ RGB
        RGB rgb = sensor->ToSensorRGB(L, lambda);

        // Optionally clamp sensor RGB value
        Float m = std::max({rgb.r, rgb.g, rgb.b});
        if (m > maxComponentValue) {
            rgb *= maxComponentValue / m;
        }

        DCHECK(InsideExclusive(pFilm, pixelBounds));

        // Update pixel values with filtered sample contribution
        PixelWindow &pixelWindow = pixels[pFilm];

        // add to current Film
        rgbSum[0] += rgb[0];
        rgbSum[1] += rgb[1];
        rgbSum[2] += rgb[2];

        squaredSum[0] += (rgb[0] * rgb[0]);
        squaredSum[1] += (rgb[1] * rgb[1]);
        squaredSum[2] += (rgb[2] * rgb[2]);

        nsamples += 1;
        weightSum += weight;

        // add sample to current buffer with weight
        pixelWindow.buffers[pixelWindow.index].rgbSum[0] += rgb[0];
        pixelWindow.buffers[pixelWindow.index].rgbSum[1] += rgb[1];
        pixelWindow.buffers[pixelWindow.index].rgbSum[2] += rgb[2];

        pixelWindow.buffers[pixelWindow.index].weightSum += weight;
       
        // Keep all information about samples
        pixelWindow.allrgbSum[0] += rgb[0];
        pixelWindow.allrgbSum[1] += rgb[1];
        pixelWindow.allrgbSum[2] += rgb[2];

        pixelWindow.squaredSum[0] += (rgb[0] * rgb[0]);
        pixelWindow.squaredSum[1] += (rgb[1] * rgb[1]);
        pixelWindow.squaredSum[2] += (rgb[2] * rgb[2]);

        pixelWindow.allWeightSum += weight;

        pixelWindow.index += 1;
        pixelWindow.nsamples += 1;

        if (pixelWindow.index >= pixelWindow.windowSize) {
            Float currentWeight = 0.;

            if (pixelWindow.windowSize == 1) {
                // store channel information
                for (int i = 0; i < 3; i++) {
                    currentWeight += pixelWindow.buffers[0].weightSum;
                    pixelWindow.rgbSum[i] += pixelWindow.buffers[0].rgbSum[i];
                    pixelWindow.splatRGB[i] = pixelWindow.splatRGB[i] + pixelWindow.buffers[0].splatRGB[i];
                }
            } 
            else
            {
                // based on channel numbers
                for (int i = 0; i < 3; i++) {

                    // store channel information in available temp buffer
                    for (int j = 0; j < pixelWindow.windowSize; j++) {
                        pixelWindow.cvalues[i] = pixelWindow.buffers[j].rgbSum[i];
                        pixelWindow.sortedValues[i] = pixelWindow.buffers[j].rgbSum[i];
                        // per channel management (but weight can be different depending of median buffer)
                        pixelWindow.weightsSum[i] = pixelWindow.buffers[j].weightSum;
                        pixelWindow.csplats[i] = pixelWindow.buffers[j].splatRGB[i];  
                        pixelWindow.indices[j] = j; // restore indices
                    }

                    // Need now to sort data
                    //std::sort(std::begin(pixelWindow.indices), std::begin(pixelWindow.indices) + pixelWindow.windowSize, 
                    //    [&](int i, int j){return pixelWindow.cvalues[i] < pixelWindow.cvalues[j];
                    //});

                    int j, k, min, tempI;
                    Float temp;
                    int n = pixelWindow.windowSize;

                    for (j = 0; j < n - 1; j++) {
                        min = j;
                        for (k = j + 1; k < n; k++)
                            if (pixelWindow.sortedValues[k] < pixelWindow.sortedValues[min])
                                min = k;

                        temp = pixelWindow.sortedValues[j];
                        pixelWindow.sortedValues[j] = pixelWindow.sortedValues[min];
                        pixelWindow.sortedValues[min] = temp;

                        tempI = pixelWindow.indices[j];
                        pixelWindow.indices[j] = pixelWindow.indices[min];
                        pixelWindow.indices[min] = tempI;
                    }

                    // need to find median value of pixelWindow.cavlues

                    // Float medianWeight, median = 0.;
                    // double medianSplat = 0;

                    // // compute median from means
                    // // find associated weightsum index and use it
                    // // Classical MON
                    // // if (windowSize % 2 == 1){
                    unsigned unsortedIndex = pixelWindow.indices[int(pixelWindow.windowSize/2)];

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
                    currentWeight += pixelWindow.weightsSum[unsortedIndex];
                    pixelWindow.rgbSum[i] += pixelWindow.cvalues[unsortedIndex];
                    pixelWindow.splatRGB[i] = pixelWindow.splatRGB[i] + pixelWindow.csplats[unsortedIndex];
                }
            }

            pixelWindow.weightSum += (currentWeight / 3);

            // // clear now the buffers data for the sliding window
            for (int i = 0; i < pixelWindow.windowSize; i++) {

                for (int j = 0; j < 3; j++) {
                    pixelWindow.buffers[i].rgbSum[j] = 0.;
                    pixelWindow.buffers[i].splatRGB[j] = 0.;
                }

                pixelWindow.buffers[i].weightSum = 0.;
            }

            // TODO : add dynamic windowSize
            pixelWindow.windowSize = 1;
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


            // reset index to sliding over WindowSize
            pixelWindow.index = 0;
        }
    }

    PBRT_CPU_GPU
    bool UsesVisibleSurface() const { return false; }

    RGBFilm() = default;

    RGBFilm(FilmBaseParameters p, const RGBColorSpace *colorSpace,
            Float maxComponentValue = Infinity, bool writeFP16 = true,
            Allocator alloc = {});

    static RGBFilm *Create(const ParameterDictionary &parameters, Float exposureTime,
                           Filter filter, const RGBColorSpace *colorSpace,
                           const FileLoc *loc, Allocator alloc);

    PBRT_CPU_GPU
    void AddSplat(const Point2f &p, SampledSpectrum v, const SampledWavelengths &lambda);

    void WriteImage(ImageMetadata metadata, Float splatScale = 1, unsigned imageIndex = 1);
    Image GetImage(ImageMetadata *metadata, Float splatScale = 1);

    std::string ToString() const;

    PBRT_CPU_GPU
    RGB ToOutputRGB(SampledSpectrum L, const SampledWavelengths &lambda) const {
        RGB sensorRGB = sensor->ToSensorRGB(L, lambda);
        return outputRGBFromSensorRGB * sensorRGB;
    }

    PBRT_CPU_GPU
    void ComputeStd() {
        
        Float stdSum = 0;

        // ParallelFor2D(pixelBounds, [&](Point2i p) {
            
        //     stdSum += pixels[p].currentStd;
        // });

        // return stdSum / (fullResolution.x * fullResolution.y);

        // // compute current mean and std
        for (int i = 0; i < 3; i++) {
            Float mean = rgbSum[i] / nsamples;
            stdSum += (squaredSum[i] / nsamples) - (mean * mean);
        }

        // divide per number of chanels and get current film std
        currentStd = std::sqrt(stdSum / 3);
        // std::cout << "Scene std is " << currentStd << " with " << nsamples << std::endl;
    }

  private:
    // RGBFilm::Pixel Definition
    // struct Pixel {
    //     Pixel() = default;
    //     double rgbSum[3] = {0., 0., 0.};
    //     double weightSum = 0.;
    //     AtomicDouble splatRGB[3];
    //     VarianceEstimator<Float> varianceEstimator;
    // };


    // RGBFilm Private Members
    double rgbSum[3] = {0., 0., 0.};
    double squaredSum[3] = {0., 0., 0.};
    double weightSum = 0;
    Float currentStd = 0;
    int nsamples = 0;

    const RGBColorSpace *colorSpace;
    Float maxComponentValue;
    bool writeFP16;
    Float filterIntegral;
    SquareMatrix<3> outputRGBFromSensorRGB;
    Array2D<PixelWindow> pixels;
    std::unique_ptr<Estimator> estimator;
};

// GBufferFilm Definition
class GBufferFilm : public FilmBase {
  public:

    PBRT_CPU_GPU
    void Clear() {};
        
    // GBufferFilm Public Methods
    GBufferFilm(FilmBaseParameters p, const RGBColorSpace *colorSpace,
                Float maxComponentValue = Infinity, bool writeFP16 = true,
                Allocator alloc = {});

    static GBufferFilm *Create(const ParameterDictionary &parameters, Float exposureTime,
                               Filter filter, const RGBColorSpace *colorSpace,
                               const FileLoc *loc, Allocator alloc);

    PBRT_CPU_GPU
    void AddSample(const Point2i &pFilm, SampledSpectrum L,
                   const SampledWavelengths &lambda, const VisibleSurface *visibleSurface,
                   Float weight);

    PBRT_CPU_GPU
    void AddSplat(const Point2f &p, SampledSpectrum v, const SampledWavelengths &lambda);

    PBRT_CPU_GPU
    RGB ToOutputRGB(const SampledSpectrum &L, const SampledWavelengths &lambda) const {
        RGB cameraRGB = sensor->ToSensorRGB(L, lambda);
        return outputRGBFromSensorRGB * cameraRGB;
    }

    PBRT_CPU_GPU
    bool UsesVisibleSurface() const { return true; }

    PBRT_CPU_GPU
    RGB GetPixelRGB(const Point2i &p, Float splatScale = 1) const {
        const Pixel &pixel = pixels[p];
        RGB rgb(pixel.rgbSum[0], pixel.rgbSum[1], pixel.rgbSum[2]);

        // Normalize pixel with weight sum
        Float weightSum = pixel.weightSum;
        if (weightSum != 0)
            rgb /= weightSum;

        // Add splat value at pixel
        for (int c = 0; c < 3; ++c)
            rgb[c] += splatScale * pixel.splatRGB[c] / filterIntegral;

        rgb = outputRGBFromSensorRGB * rgb;

        return rgb;
    }

    PBRT_CPU_GPU
    void ComputeStd() {
        // return 0.;
    }

    void WriteImage(ImageMetadata metadata, Float splatScale = 1, unsigned imageIndex = 1);
    Image GetImage(ImageMetadata *metadata, Float splatScale = 1);

    std::string ToString() const;

  private:
    // GBufferFilm::Pixel Definition
    struct Pixel {
        Pixel() = default;
        double rgbSum[3] = {0., 0., 0.};
        double weightSum = 0.;
        AtomicDouble splatRGB[3];
        Point3f pSum;
        Float dzdxSum = 0, dzdySum = 0;
        Normal3f nSum, nsSum;
        double albedoSum[3] = {0., 0., 0.};
        VarianceEstimator<Float> varianceEstimator[3];
    };

    // GBufferFilm Private Members
    Array2D<Pixel> pixels;
    const RGBColorSpace *colorSpace;
    Float maxComponentValue;
    bool writeFP16;
    Float filterIntegral;
    SquareMatrix<3> outputRGBFromSensorRGB;
};

PBRT_CPU_GPU
inline SampledWavelengths Film::SampleWavelengths(Float u) const {
    auto sample = [&](auto ptr) { return ptr->SampleWavelengths(u); };
    return Dispatch(sample);
}

PBRT_CPU_GPU
inline Bounds2f Film::SampleBounds() const {
    auto sb = [&](auto ptr) { return ptr->SampleBounds(); };
    return Dispatch(sb);
}

PBRT_CPU_GPU
inline Bounds2i Film::PixelBounds() const {
    auto pb = [&](auto ptr) { return ptr->PixelBounds(); };
    return Dispatch(pb);
}

PBRT_CPU_GPU
inline Point2i Film::FullResolution() const {
    auto fr = [&](auto ptr) { return ptr->FullResolution(); };
    return Dispatch(fr);
}

PBRT_CPU_GPU
inline void Film::Clear() {
    auto cl = [&](auto ptr) { return ptr->Clear(); };
    return Dispatch(cl);
}

PBRT_CPU_GPU
inline Float Film::Diagonal() const {
    auto diag = [&](auto ptr) { return ptr->Diagonal(); };
    return Dispatch(diag);
}

PBRT_CPU_GPU
inline Filter Film::GetFilter() const {
    auto filter = [&](auto ptr) { return ptr->GetFilter(); };
    return Dispatch(filter);
}

PBRT_CPU_GPU
inline bool Film::UsesVisibleSurface() const {
    auto uses = [&](auto ptr) { return ptr->UsesVisibleSurface(); };
    return Dispatch(uses);
}

// P3D update
PBRT_CPU_GPU
inline void Film::ComputeStd() {
    auto std = [&](auto ptr) { return ptr->ComputeStd(); };
    return Dispatch(std);
}


PBRT_CPU_GPU
inline RGB Film::GetPixelRGB(const Point2i &p, Float splatScale) const {
    auto get = [&](auto ptr) { return ptr->GetPixelRGB(p, splatScale); };
    return Dispatch(get);
}

PBRT_CPU_GPU
inline RGB Film::ToOutputRGB(const SampledSpectrum &L,
                             const SampledWavelengths &lambda) const {
    auto out = [&](auto ptr) { return ptr->ToOutputRGB(L, lambda); };
    return Dispatch(out);
}

PBRT_CPU_GPU
inline void Film::AddSample(const Point2i &pFilm, SampledSpectrum L,
                            const SampledWavelengths &lambda,
                            const VisibleSurface *visibleSurface, Float weight) {
    auto add = [&](auto ptr) {
        return ptr->AddSample(pFilm, L, lambda, visibleSurface, weight);
    };
    return Dispatch(add);
}

PBRT_CPU_GPU
inline const PixelSensor *Film::GetPixelSensor() const {
    auto filter = [&](auto ptr) { return ptr->GetPixelSensor(); };
    return Dispatch(filter);
}

}  // namespace pbrt

#endif  // PBRT_FILM_H
