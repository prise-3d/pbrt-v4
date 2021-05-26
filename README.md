[P3D] pbrt, Version 4 (Early Release)
=====================================

[<img src="https://github.com/mmp/pbrt-v4/workflows/cpu-build-and-test/badge.svg">](https://github.com/mmp/pbrt-v4/actions?query=workflow%3Acpu-build-and-test)
[<img src="https://github.com/mmp/pbrt-v4/workflows/gpu-build-only/badge.svg">](https://github.com/mmp/pbrt-v4/actions?query=workflow%3Agpu-build-only)

![Transparent Machines frame, via @beeple](images/teaser-transparent-machines.png)

This is a **custom version** of the early release of [pbrt-v4](https://github.com/mmp/pbrt-v4) for the **PrISE-3D** project, the rendering system that will be
described in the (eventually) forthcoming fourth edition of *Physically
Based Rendering: From Theory to Implementation*.  (We hope to have an
online version of the book posted a few months into 2021 and printed books available
in Summer 2021.)

We are making this code available for hardy adventurers; it's not yet
extensively documented, but if you're familiar with previous versions of
pbrt, you should be able to make your away around it.  Our hope is that the
system will be useful to some people in its current form and that any bugs
in the current implementation might be found now, allowing us to correct
them before the book is final.

A number of scenes for pbrt-v4 are [available in a git
repository](https://github.com/mmp/pbrt-v4-scenes).

Custom version [Parameters]
---------------------------

Current version is an extension of pbrt-v4 with use of some specific needs required during thesis:

Extended command line parameter:

- `--folder`: {string} -- output folder of current rendered scene ;
- `--nimages`: {unsigned} -- number of independent images of `spp` samples to generate ;
- `--independent`: {bool} -- save or not in an independant way (default 1, hence true) ;
- `--startindex`: {unsigned} -- start output index of first generated image for new run ;
- `--estimator`: {string} -- Expected estimator as output ["mean", "mon", "pakmon", "mean_or_mon"] using `src/pbrt/estimators.h` Estimator factory.

**Note:** current version enable the use `MoN` (Median of meaNs) estimator as output:
- `nbuffers:` set the number of buffers (M-estimators) to use. It is a constant value in order to work on GPU. Value can be update and available at the top of the `src/pbrt/estimators.h` file (default 11). You need to compile again the pbrt version. A value of `1`, is equivalent to classical mean estimator ;
  
Custom version [Features]
-------------------------

#### Stereoscopic, AutoStereoscopic, and OmnidirectionalStereoscopic (ODS) cameras

Stereocopic camera example into `.pbrt` file:
```
Camera "stereoscopic" "float fov" 50
        "string view" "right" 
#        "string view" "left" 
       "float eyeDistance" [0.065]  

```

**Note:** it is necessary to generate the right and left images separately.

AutoStereocopic camera example into `.pbrt` file:
```
Camera "autostereoscopic" "float fov" 50
       "integer view" [7]
       "integer nbView" [8]
       "float eyeDistance" [0.065]
```

ODS camera example into `.pbrt` file:
```
Camera "ODS" "float fov" 50
        "string view" "right" 
#        "string view" "left" 
       "float eyeDistance" [0.065]  

**Note:** it is necessary to generate the right and left images separately.


**Note:** only one view can be generated each time. The calculated view is based on the desired number of views but also on the distance fixed at eye level.

#### rgbcat2

`rgbcat2` in an executable (available in build folder) in order to merge `left` and `right` images obtained from Stereoscopic camera.

Usage:

```bash
./rgbcat2 <image-left.png> <image-right.png> <output-image.png>
```
