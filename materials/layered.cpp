
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


// materials/layered.cpp*
#include "materials/layered.h"
#include "materials/matte.h"
#include "spectrum.h"
#include "reflection.h"
#include "paramset.h"
#include "texture.h"

// LayeredMaterial Method Definitions
BSDF *LayeredMaterial::GetBSDF(const DifferentialGeometry &dgGeom,
      const DifferentialGeometry &dgShading,
      MemoryArena &arena) const {
    // Allocate _BSDF_, possibly doing bump mapping with _bumpMap_

    BSDF *b1 = m1->GetBSDF(dgGeom, dgShading, arena); // m1: coating layer
    BSDF *b2 = m2->GetBSDF(dgGeom, dgShading, arena); // m2: base layer

    float eta_i = 1.0; //ior1->Evaluate(dgs); // air
    float eta_t = ior->Evaluate(dgShading); // coating

    Fresnel *f12 = BSDF_ALLOC(arena, FresnelDielectric)(eta_i, eta_t);
    Fresnel *f21 = BSDF_ALLOC(arena, FresnelDielectric)(eta_t, eta_i);

    Spectrum a = absorption->Evaluate(dgShading).Clamp();
    float d = thickness->Evaluate(dgShading);

    // Create a layered material on top of base b2 considering the paramters of coating b1: Fresnel, absorption, depth
    int n2 = b2->NumComponents();
    for (int i = 0; i < n2; ++i)
        b1->Add(BSDF_ALLOC(arena, LayeredBxDF)(b2->bxdfs[i], f12, f21, a, d, eta_i, eta_t));

    return b1;
}


LayeredMaterial *CreateLayeredMaterial(const Transform &xform,
        const TextureParams &mp, const Reference<Material> &m1,
        const Reference<Material> &m2) {
    Reference<Texture<float> > ior = mp.GetFloatTexture("ior", float(1.5f)); // ior of m1: default 1.5 for glass
    Reference<Texture<float> > d = mp.GetFloatTexture("thickness", float(1.0f));
    Reference<Texture<Spectrum> > a = mp.GetSpectrumTexture("absorption", Spectrum(0.5f));
    return new LayeredMaterial(m1, m2, ior, d, a);
}

