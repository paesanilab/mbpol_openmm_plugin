
/* Portions copyright (c) 2006 Stanford University and Simbios.
 * Contributors: Pande Group
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */


#include "MBPolReferenceDispersionForce.h"
#include <algorithm>
#include <cctype>
#include <cmath>

using std::vector;
using OpenMM::RealVec;

MBPolReferenceDispersionForce::MBPolReferenceDispersionForce( ) : _nonbondedMethod(NoCutoff), _cutoff(1.0e+10) {

    _periodicBoxDimensions = RealVec( 0.0, 0.0, 0.0 );
}

MBPolReferenceDispersionForce::NonbondedMethod MBPolReferenceDispersionForce::getNonbondedMethod( void ) const {
    return _nonbondedMethod;
}

void MBPolReferenceDispersionForce::setNonbondedMethod( MBPolReferenceDispersionForce::NonbondedMethod nonbondedMethod ){
    _nonbondedMethod = nonbondedMethod;
}

void MBPolReferenceDispersionForce::setCutoff( double cutoff ){
    _cutoff  = cutoff;
}

double MBPolReferenceDispersionForce::getCutoff( void ) const {
    return _cutoff;
}

void MBPolReferenceDispersionForce::setPeriodicBox( const RealVec& box ){
    _periodicBoxDimensions = box;
}

RealVec MBPolReferenceDispersionForce::getPeriodicBox( void ) const {
    return _periodicBoxDimensions;
}

const double C6_HH = 2.009358600184719e+01; // kcal/mol * A^(-6)
const double C6_OH = 8.349556669872743e+01; // kcal/mol * A^(-6)
const double C6_OO = 2.373212214147944e+02; // kcal/mol * A^(-6)

const double d6_HH =  9.406475169954112e+00; // A^(-1)
const double d6_OH =  9.775202425217957e+00; // A^(-1)
const double d6_OO =  9.295485815062264e+00; // A^(-1)

template <int N>
struct Factorial
{
    enum { value = N * Factorial<N - 1>::value };
};

template <>
struct Factorial<0>
{
    enum { value = 1 };
};

double tang_toennies(int n, const double& x)
{
    assert(n >= 0);

    int nn = n;

    double sum = 1.0 + x/nn;
    while (--nn != 0)
        sum = 1.0 + sum*x/nn;

    double tt = 1.0 - sum*std::exp(-x);

    if (std::fabs(tt) < 1.0e-8) {

        double term(1);
        for (nn = n; nn != 0; --nn)
            term *= x/nn;

        sum = 0.0;
        for (nn = n + 1; nn < 1000; ++nn) {
            term *= x/nn;
            sum += term;

            if (std::fabs(term/sum) < 1.0e-8)
                break;
        }

        tt = sum*std::exp(-x);
    }

    return tt;
}

inline double x6(const double& C6, const double& d6,
                 const OpenMM::RealVec& p1, const OpenMM::RealVec& p2,
                 OpenMM::RealVec& g1,       OpenMM::RealVec& g2)
{
    const double dx = (p1[0] - p2[0])*nm_to_A;
    const double dy = (p1[1] - p2[1])*nm_to_A;
    const double dz = (p1[2] - p2[2])*nm_to_A;

    const double rsq = dx*dx + dy*dy + dz*dz;
    const double r = std::sqrt(rsq);

    const double d6r = d6*r;
    const double tt6 = tang_toennies(6, d6r);

    const double inv_rsq = 1.0/rsq;
    const double inv_r6 = inv_rsq*inv_rsq*inv_rsq;

    const double e6 = C6*tt6*inv_r6;

    const double if6 = 1.0/Factorial<6>::value;

    const double grd = 6*e6*inv_rsq - C6*std::pow(d6, 7)*if6*std::exp(-d6r)/r;

    g1[0] += dx*grd * cal2joule * -nm_to_A;
    g2[0] -= dx*grd * cal2joule * -nm_to_A;

    g1[1] += dy*grd * cal2joule * -nm_to_A;
    g2[1] -= dy*grd * cal2joule * -nm_to_A;

    g1[2] += dz*grd * cal2joule * -nm_to_A;
    g2[2] -= dz*grd * cal2joule * -nm_to_A;

    return - e6 * cal2joule;
}

RealOpenMM MBPolReferenceDispersionForce::calculatePairIxn( int siteI, int siteJ,
                                                      const std::vector<RealVec>& particlePositions,
                                                      const std::vector<std::vector<int> >& allParticleIndices,
                                                      vector<RealVec>& forces ) const {

        // siteI and siteJ are indices in a oxygen-only array, in order to get the position of an oxygen, we need:
        // allParticleIndices[siteI][0]
        // first hydrogen: allParticleIndices[siteI][1]
        // second hydrogen: allParticleIndices[siteI][2]
        // same for the second water molecule



        const OpenMM::RealVec& Oa  = particlePositions[allParticleIndices[siteI][0]];
        const OpenMM::RealVec& Ha1 = particlePositions[allParticleIndices[siteI][1]];
        const OpenMM::RealVec& Ha2 = particlePositions[allParticleIndices[siteI][2]];

        const OpenMM::RealVec& Ob  = particlePositions[allParticleIndices[siteJ][0]];
        const OpenMM::RealVec& Hb1 = particlePositions[allParticleIndices[siteJ][1]];
        const OpenMM::RealVec& Hb2 = particlePositions[allParticleIndices[siteJ][2]];

        OpenMM::RealVec& gOa  = forces[allParticleIndices[siteI][0]];
        OpenMM::RealVec& gHa1 = forces[allParticleIndices[siteI][1]];
        OpenMM::RealVec& gHa2 = forces[allParticleIndices[siteI][2]];

        OpenMM::RealVec& gOb  = forces[allParticleIndices[siteJ][0]];
        OpenMM::RealVec& gHb1 = forces[allParticleIndices[siteJ][1]];
        OpenMM::RealVec& gHb2 = forces[allParticleIndices[siteJ][2]];

        const double HH6 =
               x6(C6_HH, d6_HH, Ha1, Hb1, gHa1, gHb1)
             + x6(C6_HH, d6_HH, Ha1, Hb2, gHa1, gHb2)
             + x6(C6_HH, d6_HH, Ha2, Hb1, gHa2, gHb1)
             + x6(C6_HH, d6_HH, Ha2, Hb2, gHa2, gHb2);

           const double OH6 =
               x6(C6_OH, d6_OH, Oa, Hb1, gOa, gHb1)
             + x6(C6_OH, d6_OH, Oa, Hb2, gOa, gHb2)
             + x6(C6_OH, d6_OH, Ob, Ha1, gOb, gHa1)
             + x6(C6_OH, d6_OH, Ob, Ha2, gOb, gHa2);

           const double OO6 =
               x6(C6_OO, d6_OO, Oa, Ob, gOa, gOb);

    RealOpenMM energy= (HH6 + OH6 + OO6);

    return energy;

}

RealOpenMM MBPolReferenceDispersionForce::calculateForceAndEnergy( int numParticles,
                                                             const vector<RealVec>& particlePositions,
                                                             const std::vector<std::vector<int> >& allParticleIndices,
                                                             const NeighborList& neighborList,
                                                             vector<RealVec>& forces ) const {


    RealOpenMM energy = 0.;
    for( unsigned int ii = 0; ii < neighborList.size(); ii++ ){

        OpenMM::AtomPair pair       = neighborList[ii];
        int siteI                   = pair.first;
        int siteJ                   = pair.second;

        energy                     += calculatePairIxn( siteI, siteJ,
                particlePositions, allParticleIndices, forces );

    }

    return energy;
}
