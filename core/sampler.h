
/*
    pbrt source code Copyright(c) 1998-2009 Matt Pharr and Greg Humphreys.

    This file is part of pbrt.

    pbrt is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.  Note that the text contents of
    the book "Physically Based Rendering" are *not* licensed under the
    GNU GPL.

    pbrt is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

 */

#ifndef PBRT_CORE_SAMPLER_H
#define PBRT_CORE_SAMPLER_H

// core/sampler.h*
#include "pbrt.h"
#include "geometry.h"
#include "rng.h"
#include "memory.h"

// Sampling Declarations
class Sampler {
public:
    // Sampler Interface
    virtual ~Sampler();
    Sampler(int xstart, int xend, int ystart, int yend,
            int spp, float sopen, float sclose);
    virtual int GetMoreSamples(Sample *sample) = 0;
    virtual int MaximumSampleCount() = 0;
    virtual bool ReportResults(Sample *samples, const RayDifferential *rays,
        const Spectrum *Ls, const Intersection *isects, int count);
    virtual Sampler *GetSubSampler(int num, int count) = 0;
    virtual int RoundSize(int size) const = 0;

    // Sampler Public Data
    int xPixelStart, xPixelEnd, yPixelStart, yPixelEnd;
    int samplesPerPixel;
    float shutterOpen, shutterClose;
protected:
    // Sampler Protected Methods
    void ComputeSubWindow(int num, int count, int *xstart, int *xend, int *ystart, int *yend) const;
};


struct CameraSample {
    float imageX, imageY;
    float lensU, lensV;
    float time;
};


struct Sample : public CameraSample {
    // Sample Public Methods
    Sample(Sampler *sampler, SurfaceIntegrator *surf, VolumeIntegrator *vol,
        const Scene *scene);
    u_int Add1D(u_int num) {
        n1D.push_back(num);
        return n1D.size()-1;
    }
    u_int Add2D(u_int num) {
        n2D.push_back(num);
        return n2D.size()-1;
    }
    ~Sample() {
        if (oneD != NULL) {
            FreeAligned(oneD[0]);
            FreeAligned(oneD);
        }
    }
    Sample *Duplicate(int count, RNG *rng) const;

    // Sample Public Data
    vector<u_int> n1D, n2D;
    float **oneD, **twoD;
    mutable RNG *rng;
private:
    // Sample Private Methods
    void AllocateSampleMemory();
    Sample() { oneD = twoD = NULL; }
};



#endif // PBRT_CORE_SAMPLER_H
