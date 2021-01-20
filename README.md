[P3D] pbrt, Version 4 (Early Release)
=====================================

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
- `--monk`: {unsigned} -- use of Median of meaNs estimator instead of classical mean ;
- `--pakmon`: {bool} -- specify if PakMoN extension is used or not [0, 1].


__TODO:__
- `--independent`: {bool} -- adapt to GPU ; 
- `--monk`: {unsigned} -- adapt to GPU ; 
- `--pakmon`: {unsigned} -- adapt to GPU.

Custom version [Features]
-------------------------

#### Stereoscopic and AutoStereoscopic cameras

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

**Note:** only one view can be generated each time. The calculated view is based on the desired number of views but also on the distance fixed at eye level.

#### rgbcat2

`rgbcat2` in an executable (available in build folder) in order to merge `left` and `right` images obtained from Stereoscopic camera.

Usage:
```bash
./rgbcat2 <image-left.png> <image-right.png> <output-image.png>
```