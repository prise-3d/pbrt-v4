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

// P3D update k mon parameter (default 1, no kmon use, hence classical mean)
const int kmon = 11;
const std::vector<std::string> estimators = {"mean", "mon", "pakmon", "mean_or_mon"};

// PixelSensor Definition
class PixelSensor {
  public:
    // PixelSensor Public Methods
    static PixelSensor *Create(const ParameterDictionary &parameters,
                               const RGBColorSpace *colorSpace, Float exposureTime,
                               const FileLoc *loc, Allocator alloc);

    static PixelSensor *CreateDefault(Allocator alloc = {});

    PixelSensor(SpectrumHandle r, SpectrumHandle g, SpectrumHandle b,
                const RGBColorSpace *outputColorSpace, Float wbTemp, Float imagingRatio,
                Allocator alloc)
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
            SpectrumHandle s = swatchReflectances[i];
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
    RGB ToSensorRGB(const SampledSpectrum &L, const SampledWavelengths &lambda) const {
        return imagingRatio * RGB((r_bar.Sample(lambda) * L).Average(),
                                  (g_bar.Sample(lambda) * L).Average(),
                                  (b_bar.Sample(lambda) * L).Average());
    }

    // PixelSensor Public Members
    SquareMatrix<3> XYZFromSensorRGB;

  private:
    // PixelSensor Private Methods
    template <typename Triplet>
    static Triplet ProjectReflectance(SpectrumHandle r, SpectrumHandle illum,
                                      SpectrumHandle b1, SpectrumHandle b2,
                                      SpectrumHandle b3);

    // PixelSensor Private Members
    DenselySampledSpectrum r_bar, g_bar, b_bar;
    Float imagingRatio;
    static constexpr int nSwatchReflectances = 24;
    static SpectrumHandle swatchReflectances[nSwatchReflectances];
};

// PixelSensor Inline Methods
template <typename Triplet>
inline Triplet PixelSensor::ProjectReflectance(SpectrumHandle refl, SpectrumHandle illum,
                                               SpectrumHandle b1, SpectrumHandle b2,
                                               SpectrumHandle b3) {
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
    FilmBaseParameters(const ParameterDictionary &parameters, FilterHandle filter,
                       const PixelSensor *sensor, const FileLoc *loc);
    FilmBaseParameters(Point2i fullResolution, Bounds2i pixelBounds, FilterHandle filter,
                       Float diagonal, const PixelSensor *sensor, std::string filename)
        : fullResolution(fullResolution),
          pixelBounds(pixelBounds),
          filter(filter),
          diagonal(diagonal),
          sensor(sensor),
          filename(filename) {}

    Point2i fullResolution;
    Bounds2i pixelBounds;
    FilterHandle filter;
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
    FilterHandle GetFilter() const { return filter; }
    PBRT_CPU_GPU
    const PixelSensor *GetPixelSensor() const { return sensor; }
    std::string GetFilename() const { return filename; }

    // P3D Updates
    void SetFilename(std::string filename) { this->filename = filename; }
    
    PBRT_CPU_GPU
    SampledWavelengths SampleWavelengths(Float u) const {
        return SampledWavelengths::SampleXYZ(u);
    }

    std::string BaseToString() const;

    PBRT_CPU_GPU
    Bounds2f SampleBounds() const;

  protected:
    // FilmBase Protected Members
    Point2i fullResolution;
    Bounds2i pixelBounds;
    FilterHandle filter;
    Float diagonal;
    const PixelSensor *sensor;
    std::string filename;
};

// RGBFilm Definition
class RGBFilm : public FilmBase {
  public:
    // RGBFilm Public Methods
    PBRT_CPU_GPU
    RGB GetPixelRGB(const Point2i &p, Float splatScale = 1) const {

        // P3D Updates
        const PixelMON &pixelMON = pixels[p];
        // update estimated value for rgbSum and weightSum

        //RGB rgb(pixel.rgbSum[0], pixel.rgbSum[1], pixel.rgbSum[2]);
        // Normalize _rgb_ with weight sum

        pstd::vector<Float> estimations;
        pstd::vector<Float> weights;
        pstd::vector<double> splats;

        // default estimator used
        std::string currentEstimator = estimator;

        // TODO : check if better to use IC comparisons over all channels
        // TODO : improve this part...
        if (estimator.compare("mean_or_mon") == 0) {
            
            // if true, use of mean
            if (MeanOrMon(p)) {
                currentEstimator = "mean";
            } else{
                currentEstimator = "mon";
            }
            currentEstimator = "mean";
        }

        // based on channel numbers
        for (int i = 0; i < 3; i++) {

            // loop over pixels (used as means storage) for computing real channel value
            pstd::vector<Float> cvalues;
            pstd::vector<Float> csquared;
            pstd::vector<Float> weightsSum;
            pstd::vector<double> csplats;

            for (int j = 0; j < pixelMON.k; j++) {
                cvalues.push_back(pixelMON.means[j].rgbSum[i]);
                csquared.push_back(pixelMON.means[j].squaredSum[i]);
                weightsSum.push_back(pixelMON.means[j].weightSum);
                csplats.push_back(pixelMON.means[j].splatRGB[i]);
            }

            auto cestimation = Estimate(currentEstimator, cvalues, csquared, weightsSum, csplats);

            // add associated estimations and weights for the current channel obtained
            // weight depend of the pixels used for computing MON or PakMON
            estimations.push_back(cestimation.x);
            weights.push_back(cestimation.y);
            splats.push_back(cestimation.z);
        }

        // computed filter weight sum based on each channel
        Float weightSum = (weights[0] + weights[1] + weights[2]) / 3.;
        
        RGB rgb(estimations[0], estimations[1], estimations[2]);

        AtomicDouble splatRGB[3];
        splatRGB[0] = splats[0];
        splatRGB[1] = splats[1];
        splatRGB[2] = splats[2];
        
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
        PixelMON &pixelMON = pixels[pFilm];
        
        // add sample value inside current package
        pixelMON.means[pixelMON.index].rgbSum[0] += rgb[0];
        pixelMON.means[pixelMON.index].rgbSum[1] += rgb[1];
        pixelMON.means[pixelMON.index].rgbSum[2] += rgb[2];

        // add of squared sum of new pixel value
        pixelMON.means[pixelMON.index].squaredSum[0] += rgb[0] * rgb[0];
        pixelMON.means[pixelMON.index].squaredSum[1] += rgb[1] * rgb[1];
        pixelMON.means[pixelMON.index].squaredSum[2] += rgb[2] * rgb[2];

        pixelMON.means[pixelMON.index].weightSum += weight;

        pixelMON.index += 1;

        if (pixelMON.index >= pixelMON.k)
            pixelMON.index = 0;
        // P3D Updates
    }

    PBRT_CPU_GPU
    bool UsesVisibleSurface() const { return false; }

    PBRT_CPU_GPU
    bool MeanOrMon(const Point2i &p) const {

        const PixelMON &pixelMON = pixels[p];

        // return true if use of mean
        // return false if use of MoN
        Float meanICSum = 0.;
        Float monICSum = 0.;

        // based on channel numbers
        for (int i = 0; i < 3; i++) {

            // loop over pixels (used as means storage) for computing real channel value
            pstd::vector<Float> cvalues;
            pstd::vector<Float> csquared;
            pstd::vector<Float> weightsSum;
            pstd::vector<double> csplats;

            for (int j = 0; j < pixelMON.k; j++) {
                cvalues.push_back(pixelMON.means[j].rgbSum[i]);
                csquared.push_back(pixelMON.means[j].squaredSum[i]);
                weightsSum.push_back(pixelMON.means[j].weightSum);
                csplats.push_back(pixelMON.means[j].splatRGB[i]);
            }

            pstd::vector<Float> means(cvalues); // copy of current channel values

            // depending of confidence interval of predictor, compute the required value
            // need to compute mean value, hence all values are needed
            Float meanSquaredValues = 0.;
            
            Float meanWeight, meanMean = 0.;
            double meanSplat = 0;

            for (int i = 0; i < kmon; i++) {
                
                meanMean += means[i];
                meanWeight += weightsSum[i];
                meanSplat += csplats[i];
                meanSquaredValues += csquared[i];
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
            if (kmon % 2 == 1){
                unsigned unsortedIndex = sortedIndices[int(kmon/2)];

                monWeight = weightsSum[unsortedIndex];
                monMean = cvalues[unsortedIndex];
                monSquaredValue = csquared[unsortedIndex];
                monSplat = csplats[unsortedIndex];
            }
            else{
                int k_mean = int(kmon/2);
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

        return (meanICSum / 3)  <= (monICSum / 3);
    }

    PBRT_CPU_GPU
    Point3f Estimate(std::string currentEstimator, pstd::vector<Float> cvalues, pstd::vector<Float> csquared, pstd::vector<Float> weightsSum, pstd::vector<double> csplats) const {
            
        // int nElements = cvalues.size();
        pstd::vector<Float> means(cvalues); // copy of current channel values

        // expected output
        Float weight, mean = 0.;
        double csplat = 0;

        // if mean as current estimator or unknown estimator passed as parameter
        // TODO: fix use of in_array function
        if (currentEstimator.compare("mean") == 0) { // || !pstd::in_array<std::string>(estimator, estimators)) {

            // classical mean, hence all values are needed
            for (int i = 0; i < kmon; i++) {
                
                mean += means[i];
                weight += weightsSum[i];
                csplat += csplats[i];
            }

            return Point3f((mean / kmon), (weight / kmon), (csplat / kmon));
        }

        if (currentEstimator.compare("mon") == 0 || currentEstimator.compare("pakmon") == 0) {

            // sort means vector and get sorted indices as output
            pstd::vector<int> sortedIndices = means.sort();

            // compute median from means
            // find associated weightsum index and use it
            // Classical MON
            if (kmon % 2 == 1){
                unsigned unsortedIndex = sortedIndices[int(kmon/2)];

                weight = weightsSum[unsortedIndex];
                mean = cvalues[unsortedIndex];
                csplat = csplats[unsortedIndex];
            }
            else{
                int k_mean = int(kmon/2);
                unsigned firstIndex = sortedIndices[k_mean - 1];
                unsigned secondIndex = sortedIndices[k_mean];

                weight = (weightsSum[firstIndex] + weightsSum[secondIndex]) / 2;
                mean = (cvalues[firstIndex] + cvalues[secondIndex]) / 2;
                csplat = (csplats[firstIndex]  + csplats[secondIndex]) / 2;
            }
            
            // Here use of PAKMON => compromise between Mean and MON
            if (currentEstimator.compare("pakmon") == 0) {
                
                Float meansSum = 0;

                for (int i = 0; i < cvalues.size(); i++)
                    meansSum += cvalues[i];

                Float currentMean = meansSum / cvalues.size();
                    
                // compute variance distance evolution
                pstd::vector<Float> distances;

                for (int i = 0; i < means.size(); i++) {
                    
                    // use of sorted means in order to compute variance evolution by step 1
                    if (i > 1) {

                        // compute variance of each elements
                        Float var = 0;
                        
                        // use of previously sorted means
                        for(int j = 0; j < i; j++)
                        {
                            var += (means[j] - currentMean) * (means[j] - currentMean);
                        }
                        var /= (i + 1);

                        // add new 
                        distances.push_back(var);
                    }
                }

                // use of variance evolution and compute entropy
                Float distancesEntropy = getEntropy(distances);

                // Computation of PakMON using \alpha and \rho value
                unsigned middleIndex = int(kmon / 2);

                // alpha and rho automatically set value
                Float alpha = distancesEntropy;

                if (alpha < 0.000000001) {
                    alpha = 0.000000001;
                }
                int rho = (int)(middleIndex * distancesEntropy) - (int)(kmon * 0.25); // try using avoid 30% (total) of current kmon

                unsigned lowerIndex = 0;
                unsigned higherIndex = 0;
            
                // get current lower and higher index 
                if (kmon % 2 == 0) {
                    
                    lowerIndex = middleIndex - 1;
                    higherIndex = middleIndex;
                    
                } else {
                    lowerIndex = middleIndex - 1;
                    higherIndex = middleIndex - 1;
                }

                // std::cout << "MON value is: " << (mean / weight) << std::endl;

                // use of `vp` which stores ordered mean and 
                for (int i = 1; i < rho + 1; i++) {

                    // current neighbor multiple factor
                    Float multFactor = pow(alpha, i);

                    // add left and right neighbor contribution
                    mean += means[lowerIndex - i] * multFactor;
                    mean += means[higherIndex + i] * multFactor;
                    
                    // weighting contribution to take in account
                    // use of this index to retrieve the associated weightsSum
                    weight += weightsSum[sortedIndices[lowerIndex]] * multFactor;
                    weight += weightsSum[sortedIndices[higherIndex]] * multFactor;

                    csplat += csplats[sortedIndices[lowerIndex]] * multFactor;
                    csplat += csplats[sortedIndices[higherIndex]] * multFactor;
                }
            }
        }

        return Point3f(mean, weight, csplat);
    }

    PBRT_CPU_GPU
    Float getEntropy(pstd::vector<Float> values) const {

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
    }

    RGBFilm() = default;
    RGBFilm(FilmBaseParameters p, const RGBColorSpace *colorSpace,
            Float maxComponentValue = Infinity, bool writeFP16 = true,
            Allocator alloc = {});

    static RGBFilm *Create(const ParameterDictionary &parameters, Float exposureTime,
                           FilterHandle filter, const RGBColorSpace *colorSpace,
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

  private:
    // RGBFilm::Pixel Definition
    // struct Pixel {
    //     Pixel() = default;
    //     double rgbSum[3] = {0., 0., 0.};
    //     double weightSum = 0.;
    //     AtomicDouble splatRGB[3];
    //     VarianceEstimator<Float> varianceEstimator;
    // };

    // P3D Updates
    struct Pixel {
        Pixel() = default;

        double rgbSum[3] = {0., 0., 0.};
        double squaredSum[3] = {0., 0., 0.};
        AtomicDouble splatRGB[3];
        double weightSum = 0.;
    };

    // P3D Updates
    struct PixelMON {
        PixelMON() = default;

        Pixel means[kmon];
        VarianceEstimator<Float> varianceEstimator;

        int k = kmon; // number of means clusters
        int index = 0; // keep track of index used
        bool filled = false;
    };

    // RGBFilm Private Members
    const RGBColorSpace *colorSpace;
    Float maxComponentValue;
    bool writeFP16;
    Float filterIntegral;
    SquareMatrix<3> outputRGBFromSensorRGB;
    Array2D<PixelMON> pixels;
    std::string estimator = Options->estimator;
};

// GBufferFilm Definition
class GBufferFilm : public FilmBase {
  public:
    // GBufferFilm Public Methods
    GBufferFilm(FilmBaseParameters p, const RGBColorSpace *colorSpace,
                Float maxComponentValue = Infinity, bool writeFP16 = true,
                Allocator alloc = {});

    static GBufferFilm *Create(const ParameterDictionary &parameters, Float exposureTime,
                               FilterHandle filter, const RGBColorSpace *colorSpace,
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
inline SampledWavelengths FilmHandle::SampleWavelengths(Float u) const {
    auto sample = [&](auto ptr) { return ptr->SampleWavelengths(u); };
    return Dispatch(sample);
}

PBRT_CPU_GPU
inline Bounds2f FilmHandle::SampleBounds() const {
    auto sb = [&](auto ptr) { return ptr->SampleBounds(); };
    return Dispatch(sb);
}

PBRT_CPU_GPU
inline Bounds2i FilmHandle::PixelBounds() const {
    auto pb = [&](auto ptr) { return ptr->PixelBounds(); };
    return Dispatch(pb);
}

PBRT_CPU_GPU
inline Point2i FilmHandle::FullResolution() const {
    auto fr = [&](auto ptr) { return ptr->FullResolution(); };
    return Dispatch(fr);
}

PBRT_CPU_GPU
inline Float FilmHandle::Diagonal() const {
    auto diag = [&](auto ptr) { return ptr->Diagonal(); };
    return Dispatch(diag);
}

PBRT_CPU_GPU
inline FilterHandle FilmHandle::GetFilter() const {
    auto filter = [&](auto ptr) { return ptr->GetFilter(); };
    return Dispatch(filter);
}

PBRT_CPU_GPU
inline bool FilmHandle::UsesVisibleSurface() const {
    auto uses = [&](auto ptr) { return ptr->UsesVisibleSurface(); };
    return Dispatch(uses);
}

PBRT_CPU_GPU
inline RGB FilmHandle::GetPixelRGB(const Point2i &p, Float splatScale) const {
    auto get = [&](auto ptr) { return ptr->GetPixelRGB(p, splatScale); };
    return Dispatch(get);
}

PBRT_CPU_GPU
inline RGB FilmHandle::ToOutputRGB(const SampledSpectrum &L,
                                   const SampledWavelengths &lambda) const {
    auto out = [&](auto ptr) { return ptr->ToOutputRGB(L, lambda); };
    return Dispatch(out);
}

PBRT_CPU_GPU
inline void FilmHandle::AddSample(const Point2i &pFilm, SampledSpectrum L,
                                  const SampledWavelengths &lambda,
                                  const VisibleSurface *visibleSurface, Float weight) {
    auto add = [&](auto ptr) {
        return ptr->AddSample(pFilm, L, lambda, visibleSurface, weight);
    };
    return Dispatch(add);
}

PBRT_CPU_GPU
inline const PixelSensor *FilmHandle::GetPixelSensor() const {
    auto filter = [&](auto ptr) { return ptr->GetPixelSensor(); };
    return Dispatch(filter);
}

}  // namespace pbrt

#endif  // PBRT_FILM_H
