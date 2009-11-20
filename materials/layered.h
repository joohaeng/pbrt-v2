
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

#ifndef PBRT_MATERIALS_LAYERED_H
#define PBRT_MATERIALS_LAYERED_H

// materials/layered.h*
#include "pbrt.h"
#include "material.h"

// LayeredMaterial Declarations
class LayeredMaterial : public Material {
public:
    // LayeredMaterial Public Methods
    LayeredMaterial(const Reference<Material> &mat1, const Reference<Material> &mat2,
                //const Reference<Texture<float> > &ior_,
                float ior_,
                const Reference<Texture<float> > &d,
                const Reference<Texture<Spectrum> > &a,
                const Reference<Texture<float> > &tir_,
                const Reference<Texture<float> > &mf_normal_,
                const Reference<Texture<float> > &base_only_
	) 
	{
        m1 = mat1;
        m2 = mat2;
        ior = ior_;
        thickness = d;
        absorption = a;
        tir = tir_;
        mf_normal = mf_normal_;
        base_only = base_only_;
    }
    BSDF *GetBSDF(const DifferentialGeometry &dgGeom,
                  const DifferentialGeometry &dgShading,
                  MemoryArena &arena) const;
private:
    // LayeredMaterial Private Data
    Reference<Material> m1, m2;
    //Reference<Texture<float> > 	ior;
    float ior;
    Reference<Texture<float> > 	tir; // 1 for TIR computation, otherwise for no consideration.
    Reference<Texture<float> > 	mf_normal; // 0 to select a surface normal rather than a MF normal
    Reference<Texture<float> > 	base_only; // 0 to select a surface normal rather than a MF normal
    Reference<Texture<float> > 	thickness;
    Reference<Texture<Spectrum> > absorption;
};

LayeredMaterial *CreateLayeredMaterial(const Transform &xform,
    const TextureParams &mp, const Reference<Material> &m1,
    const Reference<Material> &m2);

#endif // PBRT_MATERIALS_LAYERED_H
