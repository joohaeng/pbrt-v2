
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


// core/reflection.cpp*
#include "reflection.h"
#include "spectrum.h"
#include "sampler.h"
#include "montecarlo.h"
#include <stdarg.h>

// BxDF Local Definitions
struct IrregIsoProc {
    // IrregIsoProc Public Methods
    IrregIsoProc() { sumWeights = 0.f; nFound = 0; }
    void operator()(const Point &p, const IrregIsotropicBRDFSample &sample,
                    float d2, float &maxDist2) {
        float weight = expf(-100.f * d2);
        v += weight * sample.v;
        sumWeights += weight;
        ++nFound;
    }
    Spectrum v;
    float sumWeights;
    int nFound;
};



// BxDF Utility Functions
Spectrum FrDiel(float cosi, float cost, const Spectrum &etai,
                const Spectrum &etat) {
    Spectrum Rparl = ((etat * cosi) - (etai * cost)) /
                     ((etat * cosi) + (etai * cost));
    Spectrum Rperp = ((etai * cosi) - (etat * cost)) /
                     ((etai * cosi) + (etat * cost));
    return (Rparl*Rparl + Rperp*Rperp) / 2.f;
}


Spectrum FrCond(float cosi, const Spectrum &eta, const Spectrum &k) {
    Spectrum tmp = (eta*eta + k*k) * cosi*cosi;
    Spectrum Rparl2 = (tmp - (2.f * eta * cosi) + 1) /
                      (tmp + (2.f * eta * cosi) + 1);
    Spectrum tmp_f = eta*eta + k*k;
    Spectrum Rperp2 =
        (tmp_f - (2.f * eta * cosi) + cosi*cosi) /
        (tmp_f + (2.f * eta * cosi) + cosi*cosi);
    return (Rparl2 + Rperp2) / 2.f;
}



// BxDF Method Definitions
Spectrum BRDFToBTDF::f(const Vector &wo, const Vector &wi) const {
    return brdf->f(wo, otherHemisphere(wi));
}


Spectrum BRDFToBTDF::Sample_f(const Vector &wo, Vector *wi,
                              float u1, float u2, float *pdf) const {
    Spectrum f = brdf->Sample_f(wo, wi, u1, u2, pdf);
    *wi = otherHemisphere(*wi);
    return f;
}


Spectrum ScaledBxDF::f(const Vector &wo, const Vector &wi) const {
    return s * bxdf->f(wo, wi);
}


Spectrum ScaledBxDF::Sample_f(const Vector &wo, Vector *wi,
                              float u1, float u2, float *pdf) const {
    Spectrum f = bxdf->Sample_f(wo, wi, u1, u2, pdf);
    return s * f;
}

inline Vector SnellDir(const Vector &w, float etai, float etat) {

    // Compute transmitted ray direction
    float sini2 = SinTheta2(w);
    float eta = etai / etat;
    float sint2 = eta * eta * sini2;

    // Handle total internal reflection for transmission
    //if (sint2 >= 1.) Assert("TIR"); //return 0.;
    float cost = sqrtf(max(0.f, 1.f - sint2));
	
	return Vector(eta * w.x, eta * w.y, cost);
}

inline float G(const Vector &wo, const Vector &wi, const Vector &wh) {
    float NdotWh = AbsCosTheta(wh);
    float NdotWo = AbsCosTheta(wo);
    float NdotWi = AbsCosTheta(wi);
    float WOdotWh = AbsDot(wo, wh);
	//Info("%f, %f, %f, %f", NdotWh, NdotWo, NdotWi, WOdotWh);
    return min(1.f, min((2.f * NdotWh * NdotWo / WOdotWh), (2.f * NdotWh * NdotWi / WOdotWh)));
}
Spectrum LayeredBxDF::f(const Vector &wo, const Vector &wi) const {
	// wor, wir: refr dir of wi and wo inside a coating layer
	Vector wor = SnellDir(wo, etai, etat);
	Vector wir = SnellDir(wi, etai, etat);
    Vector whr = Normalize(wir + wor);

	float tmp =	depth * (1.0f/CosTheta(wir) + 1.0f/CosTheta(wor));
	Spectrum a = (tmp > 0 ? Exp(-alpha * tmp) : Spectrum(1.));

	float g = G(wor, wir, whr);

	Spectrum spectrum_1 = Spectrum(1.f);

	//Spectrum t = Spectrum(1.f) - f21->Evaluate(CosTheta(wor));
	//if ( tir ) t = Spectrum(1.f - g) + (t * g);
	//t = Spectrum(1.f - g) + (t * g);
	Spectrum t = f21->Evaluate(CosTheta(wor)); // wo will be computed again
	t = spectrum_1 - ( tir ? t * g : t);

    return (spectrum_1 - f12->Evaluate(CosTheta(wi))) \
		* bxdf->f(wor, wir) * a * t;
	//return t;
	//return Spectrum(g);
	//return spectrum_1;
}

Spectrum LayeredBxDF::Sample_f(const Vector &wo, Vector *wi,
                              float u1, float u2, float *pdf) const {
	Vector wor = SnellDir(wo, etai, etat);
	Vector wir;
    Spectrum f = bxdf->Sample_f(wor, &wir, u1, u2, pdf);
	*wi = SnellDir(wir, etat, etai);
    return this->f(wo, *wi);
}

Fresnel::~Fresnel() { }
Spectrum FresnelConductor::Evaluate(float cosi) const {
    return FrCond(fabsf(cosi), eta, k);
}


Spectrum FresnelDielectric::Evaluate(float cosi) const {
    // Compute Fresnel reflectance for dielectric
    cosi = Clamp(cosi, -1.f, 1.f);

    // Compute indices of refraction for dielectric
    bool entering = cosi > 0.;
    float ei = eta_i, et = eta_t;
    if (!entering)
        swap(ei, et);

    // Compute _sint_ using Snell's law
    float sint = ei/et * sqrtf(max(0.f, 1.f - cosi*cosi));
    if (sint >= 1.) {
        // Handle total internal reflection
        return 1.;
    }
    else {
        float cost = sqrtf(max(0.f, 1.f - sint*sint));
        return FrDiel(fabsf(cosi), cost, ei, et);
    }
}


Spectrum SpecularReflection::Sample_f(const Vector &wo,
        Vector *wi, float u1, float u2, float *pdf) const {
    // Compute perfect specular reflection direction
    *wi = Vector(-wo.x, -wo.y, wo.z);
    *pdf = 1.f;
    return fresnel->Evaluate(CosTheta(wo)) * R / AbsCosTheta(*wi);
}


Spectrum SpecularTransmission::Sample_f(const Vector &wo,
        Vector *wi, float u1, float u2, float *pdf) const {
    // Figure out which $\eta$ is incident and which is transmitted
    bool entering = CosTheta(wo) > 0.;
    float ei = etai, et = etat;
    if (!entering)
        swap(ei, et);

    // Compute transmitted ray direction
    float sini2 = SinTheta2(wo);
    float eta = ei / et;
    float sint2 = eta * eta * sini2;

    // Handle total internal reflection for transmission
    if (sint2 >= 1.) return 0.;
    float cost = sqrtf(max(0.f, 1.f - sint2));
    if (entering) cost = -cost;
    float sintOverSini = eta;
    *wi = Vector(sintOverSini * -wo.x, sintOverSini * -wo.y, cost);
    *pdf = 1.f;
    Spectrum F = fresnel.Evaluate(CosTheta(wo));
    return (et*et)/(ei*ei) * (Spectrum(1.)-F) * T /
        AbsCosTheta(*wi);
}


Spectrum Lambertian::f(const Vector &wo, const Vector &wi) const {
    return RoverPI;
}


Spectrum OrenNayar::f(const Vector &wo, const Vector &wi) const {
    float sinthetai = SinTheta(wi);
    float sinthetao = SinTheta(wo);
    // Compute cosine term of Oren-Nayar model
    float maxcos = 0.f;
    if (sinthetai > 1e-4 && sinthetao > 1e-4) {
        float sinphii = SinPhi(wi), cosphii = CosPhi(wi);
        float sinphio = SinPhi(wo), cosphio = CosPhi(wo);
        float dcos = cosphii * cosphio + sinphii * sinphio;
        maxcos = max(0.f, dcos);
    }

    // Compute sine and tangent terms of Oren-Nayar model
    float sinalpha, tanbeta;
    if (AbsCosTheta(wi) > AbsCosTheta(wo)) {
        sinalpha = sinthetao;
        tanbeta = sinthetai / AbsCosTheta(wi);
    }
    else {
        sinalpha = sinthetai;
        tanbeta = sinthetao / AbsCosTheta(wo);
    }
    return R * INV_PI * (A + B * maxcos * sinalpha * tanbeta);
}


Microfacet::Microfacet(const Spectrum &reflectance, Fresnel *f,
                       MicrofacetDistribution *d)
    : BxDF(BxDFType(BSDF_REFLECTION | BSDF_GLOSSY)),
     R(reflectance), distribution(d), fresnel(f) {
}


Spectrum Microfacet::f(const Vector &wo, const Vector &wi) const {
    float cosThetaO = AbsCosTheta(wo);
    float cosThetaI = AbsCosTheta(wi);
    if (cosThetaI == 0.f || cosThetaO == 0.f) return Spectrum(0.f);
    Vector wh = Normalize(wi + wo);
    float cosThetaH = Dot(wi, wh);
    Spectrum F = fresnel->Evaluate(cosThetaH);
    return R * distribution->D(wh) * G(wo, wi, wh) * F /
               (4.f * cosThetaI * cosThetaO);
}


FresnelBlend::FresnelBlend(const Spectrum &d, const Spectrum &s,
                           MicrofacetDistribution *dist)
    : BxDF(BxDFType(BSDF_REFLECTION | BSDF_GLOSSY)), Rd(d), Rs(s) {
    distribution = dist;
}


Spectrum FresnelBlend::f(const Vector &wo, const Vector &wi) const {
    Spectrum diffuse = (28.f/(23.f*M_PI)) * Rd *
        (Spectrum(1.f) - Rs) *
        (1.f - powf(1.f - .5f * AbsCosTheta(wi), 5)) *
        (1.f - powf(1.f - .5f * AbsCosTheta(wo), 5));
    Vector wh = Normalize(wi + wo);
    Spectrum specular = distribution->D(wh) /
        (4.f * AbsDot(wi, wh) * max(AbsCosTheta(wi), AbsCosTheta(wo))) *
        SchlickFresnel(Dot(wi, wh));
    return diffuse + specular;
}


Point BRDFRemap(const Vector &wo, const Vector &wi) {
    float cosi = CosTheta(wi), coso = CosTheta(wo);
    float sini = SinTheta(wi), sino = SinTheta(wo);
    float phii = SphericalPhi(wi), phio = SphericalPhi(wo);
    float dphi = phii - phio;
    if (dphi < 0.) dphi += 2.f * M_PI;
    if (dphi > 2.f * M_PI) dphi -= 2.f * M_PI;
    if (dphi > M_PI) dphi = 2.f * M_PI - dphi;
    return Point(sini * sino, dphi / M_PI, cosi * coso);
}


Spectrum IrregIsotropicBRDF::f(const Vector &wo, const Vector &wi) const {
    Point m = BRDFRemap(wo, wi);
    float lastMaxDist2 = .001f;
    while (true) {
        // Try to find enough BRDF samples around _m_ within search radius
        IrregIsoProc proc;
        float maxDist2 = lastMaxDist2;
        isoBRDFData->Lookup(m, proc, maxDist2);
        if (proc.nFound > 2 || lastMaxDist2 > 1.5f)
            return proc.v.Clamp() / proc.sumWeights;
        lastMaxDist2 *= 2.f;
    }
}


Spectrum RegularHalfangleBRDF::f(const Vector &wo, const Vector &wi) const {
    // Compute $\wh$ and transform $\wi$ to halfangle coordinate system
    Vector wh = Normalize(wi + wo);
    float whTheta = SphericalTheta(wh);
    float whCosPhi = CosPhi(wh), whSinPhi = SinPhi(wh);
    float whCosTheta = CosTheta(wh), whSinTheta = SinTheta(wh);
    Vector whx(whCosPhi * whCosTheta, whSinPhi * whCosTheta, -whSinTheta);
    Vector why(-whSinPhi, whCosPhi, 0);
    Vector wd(Dot(wi, whx), Dot(wi, why), Dot(wi, wh));

    // Compute _index_ into measured BRDF tables
    float wdTheta = SphericalTheta(wd), wdPhi = SphericalPhi(wd);
    if (wdPhi > M_PI) wdPhi -= M_PI;

    // Compute indices _whThetaIndex_, _wdThetaIndex_, _wdPhiIndex_
#define REMAP(V, MAX, COUNT) \
        Clamp(int((V) / (MAX) * (COUNT)), 0, (COUNT)-1)
    int whThetaIndex = REMAP(sqrtf(max(0.f, whTheta / (M_PI / 2.f))),
                             1.f, nThetaH);
    int wdThetaIndex = REMAP(wdTheta, M_PI / 2.f, nThetaD);
    int wdPhiIndex = REMAP(wdPhi, M_PI, nPhiD);
#undef REMAP
    int index = wdPhiIndex + nPhiD * (wdThetaIndex + whThetaIndex * nThetaD);
    return Spectrum::FromRGB(&brdf[3*index]);
}


Spectrum BxDF::Sample_f(const Vector &wo, Vector *wi,
                        float u1, float u2, float *pdf) const {
    // Cosine-sample the hemisphere, flipping the direction if necessary
    *wi = CosineSampleHemisphere(u1, u2);
    if (wo.z < 0.) wi->z *= -1.f;
    *pdf = Pdf(wo, *wi);
    return f(wo, *wi);
}


float BxDF::Pdf(const Vector &wo, const Vector &wi) const {
    return SameHemisphere(wo, wi) ? fabsf(wi.z) * INV_PI : 0.f;
}


float BRDFToBTDF::Pdf(const Vector &wo,
        const Vector &wi) const {
    return brdf->Pdf(wo, -wi);
}


Spectrum Microfacet::Sample_f(const Vector &wo, Vector *wi,
                              float u1, float u2, float *pdf) const {
    distribution->Sample_f(wo, wi, u1, u2, pdf);
    if (!SameHemisphere(wo, *wi)) return Spectrum(0.f);
    return f(wo, *wi);
}


float Microfacet::Pdf(const Vector &wo, const Vector &wi) const {
    if (!SameHemisphere(wo, wi)) return 0.f;
    return distribution->Pdf(wo, wi);
}


void Blinn::Sample_f(const Vector &wo, Vector *wi,
                     float u1, float u2, float *pdf) const {
    SampleBlinn(wo, wi, u1, u2, pdf, exponent);
}


float Blinn::Pdf(const Vector &wo, const Vector &wi) const {
    return BlinnPdf(wo, wi, exponent);
}


void Anisotropic::Sample_f(const Vector &wo, Vector *wi,
                           float u1, float u2, float *pdf) const {
    // Sample from first quadrant and remap to hemisphere to sample $\wh$
    float phi, costheta;
    if (u1 < .25f) {
        sampleFirstQuadrant(4.f * u1, u2, &phi, &costheta);
    } else if (u1 < .5f) {
        u1 = 4.f * (.5f - u1);
        sampleFirstQuadrant(u1, u2, &phi, &costheta);
        phi = M_PI - phi;
    } else if (u1 < .75f) {
        u1 = 4.f * (u1 - .5f);
        sampleFirstQuadrant(u1, u2, &phi, &costheta);
        phi += M_PI;
    } else {
        u1 = 4.f * (1.f - u1);
        sampleFirstQuadrant(u1, u2, &phi, &costheta);
        phi = 2.f * M_PI - phi;
    }
    float sintheta = sqrtf(max(0.f, 1.f - costheta*costheta));
    Vector wh = SphericalDirection(sintheta, costheta, phi);
    if (!SameHemisphere(wo, wh)) wh = -wh;

    // Compute incident direction by reflecting about $\wh$
    *wi = -wo + 2.f * Dot(wo, wh) * wh;

    // Compute PDF for $\wi$ from anisotropic distribution
    float anisotropic_pdf = D(wh) / (4.f * Dot(wo, wh));
    if (Dot(wo, wh) < 0.f) anisotropic_pdf = 0.f;
    *pdf = anisotropic_pdf;
}


void Anisotropic::sampleFirstQuadrant(float u1, float u2,
        float *phi, float *costheta) const {
    if (ex == ey)
        *phi = M_PI * u1 * 0.5f;
    else
        *phi = atanf(sqrtf((ex+1.f) / (ey+1.f)) *
                     tanf(M_PI * u1 * 0.5f));
    float cosphi = cosf(*phi), sinphi = sinf(*phi);
    *costheta = powf(u2, 1.f/(ex * cosphi * cosphi +
                              ey * sinphi * sinphi + 1));
}


float Anisotropic::Pdf(const Vector &wo, const Vector &wi) const {
    Vector wh = Normalize(wo + wi);
    if (!SameHemisphere(wo, wh)) wh = -wh;
    // Compute PDF for $\wi$ from anisotropic distribution
    float anisotropic_pdf = D(wh) / (4.f * Dot(wo, wh));
    if (Dot(wo, wh) < 0.f) anisotropic_pdf = 0.f;
    return anisotropic_pdf;
}


Spectrum FresnelBlend::Sample_f(const Vector &wo, Vector *wi,
                                float u1, float u2, float *pdf) const {
    if (u1 < .5) {
        u1 = 2.f * u1;
        // Cosine-sample the hemisphere, flipping the direction if necessary
        *wi = CosineSampleHemisphere(u1, u2);
        if (wo.z < 0.) wi->z *= -1.f;
    }
    else {
        u1 = 2.f * (u1 - .5f);
        distribution->Sample_f(wo, wi, u1, u2, pdf);
        if (!SameHemisphere(wo, *wi)) return Spectrum(0.f);
    }
    *pdf = Pdf(wo, *wi);
    return f(wo, *wi);
}


float FresnelBlend::Pdf(const Vector &wo, const Vector &wi) const {
    if (!SameHemisphere(wo, wi)) return 0.f;
    return .5f * (fabsf(wi.z) * INV_PI + distribution->Pdf(wo, wi));
}


Spectrum BxDF::rho(const Vector &w, int nSamples,
                   const float *samples) const {
    Spectrum r = 0.;
    for (int i = 0; i < nSamples; ++i) {
        // Estimate one term of $\rho_\roman{hd}$
        Vector wi;
        float pdf = 0.f;
        Spectrum f = Sample_f(w, &wi, samples[2*i], samples[2*i+1], &pdf);
        if (pdf > 0.) r += f * AbsCosTheta(wi) / pdf;
    }
    return r / float(nSamples);
}


Spectrum BxDF::rho(int nSamples, const float *samples1,
        const float *samples2) const {
    Spectrum r = 0.;
    for (int i = 0; i < nSamples; ++i) {
        // Estimate one term of $\rho_\roman{hh}$
        Vector wo, wi;
        wo = UniformSampleHemisphere(samples1[2*i], samples1[2*i+1]);
        float pdf_o = INV_TWOPI, pdf_i = 0.f;
        Spectrum f = Sample_f(wo, &wi, samples2[2*i], samples2[2*i+1], &pdf_i);
        if (pdf_i > 0.)
            r += f * AbsCosTheta(wi) * AbsCosTheta(wo) / (pdf_o * pdf_i);
    }
    return r / (M_PI*nSamples);
}



// BSDF Method Definitions
BSDFSampleOffsets::BSDFSampleOffsets(int count, Sample *sample) {
    nSamples = count;
    dirOffset = sample->Add2D(nSamples);
    componentOffset = sample->Add1D(nSamples);
}


BSDFSample::BSDFSample(const Sample *sample,
                       const BSDFSampleOffsets &offsets, uint32_t num) {
    Assert(num < sample->n2D[offsets.dirOffset]);
    Assert(num < sample->n1D[offsets.componentOffset]);
    uDir[0] = sample->twoD[offsets.dirOffset][2*num];
    uDir[1] = sample->twoD[offsets.dirOffset][2*num+1];
    uComponent = sample->oneD[offsets.componentOffset][num];
}


Spectrum BSDF::Sample_f(const Vector &wo, Vector *wi, RNG &rng, BxDFType flags,
    BxDFType *sampledType) const {
    float pdf;
    BSDFSample bsdfSample(rng);
    Spectrum f = Sample_f(wo, wi, bsdfSample, &pdf, flags, sampledType);
    if (!f.IsBlack() && pdf > 0.) f /= pdf;
    return f;
}


Spectrum BSDF::Sample_f(const Vector &woW, Vector *wiW,
                        const BSDFSample &bsdfSample, float *pdf,
                        BxDFType flags, BxDFType *sampledType) const {
    // Choose which _BxDF_ to sample
    int matchingComps = NumComponents(flags);
    if (matchingComps == 0) {
        *pdf = 0.f;
        return Spectrum(0.f);
    }
    int which = min(Floor2Int(bsdfSample.uComponent * matchingComps),
                    matchingComps-1);
    BxDF *bxdf = NULL;
    int count = which;
    for (int i = 0; i < nBxDFs; ++i)
        if (bxdfs[i]->MatchesFlags(flags) && count-- == 0) {
            bxdf = bxdfs[i];
            break;
        }
    Assert(bxdf);

    // Sample chosen _BxDF_
    Vector wo = WorldToLocal(woW);
    Vector wi;
    *pdf = 0.f;
    Spectrum f = bxdf->Sample_f(wo, &wi, bsdfSample.uDir[0],
                                bsdfSample.uDir[1], pdf);
    if (*pdf == 0.f) return 0.f;
    if (sampledType) *sampledType = bxdf->type;
    *wiW = LocalToWorld(wi);

    // Compute overall PDF with all matching _BxDF_s
    if (!(bxdf->type & BSDF_SPECULAR) && matchingComps > 1) {
        for (int i = 0; i < nBxDFs; ++i) {
            if (bxdfs[i] != bxdf &&
                bxdfs[i]->MatchesFlags(flags))
                *pdf += bxdfs[i]->Pdf(wo, wi);
        }
    }
    if (matchingComps > 1) *pdf /= matchingComps;

    // Compute value of BSDF for sampled direction
    if (!(bxdf->type & BSDF_SPECULAR)) {
        f = 0.;
        if (Dot(*wiW, ng) * Dot(woW, ng) > 0) // ignore BTDFs
            flags = BxDFType(flags & ~BSDF_TRANSMISSION);
        else // ignore BRDFs
            flags = BxDFType(flags & ~BSDF_REFLECTION);
        for (int i = 0; i < nBxDFs; ++i)
            if (bxdfs[i]->MatchesFlags(flags))
                f += bxdfs[i]->f(wo, wi);
    }
    return f;
}


float BSDF::Pdf(const Vector &woW, const Vector &wiW,
        BxDFType flags) const {
    if (nBxDFs == 0.) return 0.;
    Vector wo = WorldToLocal(woW), wi = WorldToLocal(wiW);
    float pdf = 0.f;
    int matchingComps = 0;
    for (int i = 0; i < nBxDFs; ++i)
        if (bxdfs[i]->MatchesFlags(flags)) {
            ++matchingComps;
            pdf += bxdfs[i]->Pdf(wo, wi);
        }
    return matchingComps > 0 ? pdf / matchingComps : 0.f;
}


BSDF::BSDF(const DifferentialGeometry &dg, const Normal &ngeom, float e)
    : dgShading(dg), eta(e) {
    ng = ngeom;
    nn = dgShading.nn;
    sn = Normalize(dgShading.dpdu);
    tn = Cross(nn, sn);
    nBxDFs = 0;
}


Spectrum BSDF::f(const Vector &woW, const Vector &wiW,
                 BxDFType flags) const {
    Vector wi = WorldToLocal(wiW), wo = WorldToLocal(woW);
    if (Dot(wiW, ng) * Dot(woW, ng) > 0) // ignore BTDFs
        flags = BxDFType(flags & ~BSDF_TRANSMISSION);
    else // ignore BRDFs
        flags = BxDFType(flags & ~BSDF_REFLECTION);
    Spectrum f = 0.;
    for (int i = 0; i < nBxDFs; ++i)
        if (bxdfs[i]->MatchesFlags(flags))
            f += bxdfs[i]->f(wo, wi);
    return f;
}


Spectrum BSDF::rho(int nSamples, const float *samples1,
    const float *samples2, BxDFType flags) const {
    Spectrum ret(0.);
    for (int i = 0; i < nBxDFs; ++i)
        if (bxdfs[i]->MatchesFlags(flags))
            ret += bxdfs[i]->rho(nSamples, samples1, samples2);
    return ret;
}


Spectrum BSDF::rho(const Vector &wo, int nSamples,
        const float *samples, BxDFType flags) const {
    Spectrum ret(0.);
    for (int i = 0; i < nBxDFs; ++i)
        if (bxdfs[i]->MatchesFlags(flags))
            ret += bxdfs[i]->rho(wo, nSamples, samples);
    return ret;
}


