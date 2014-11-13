
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

#include "openmm/OpenMMException.h"
#include "openmm/internal/MBPolDispersionForceImpl.h"

#include "MBPolReferenceForce.h"
#include "MBPolReferenceDispersionForce.h"
#include "MBPolReferenceTwoBodyForce.h"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <map>
#include <string>
#include <iostream>

using std::pair;
using std::string;
using std::map;
using std::make_pair;

using std::vector;
using OpenMM::RealVec;
using namespace MBPolPlugin;

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

void MBPolReferenceDispersionForce::setDispersionParameters( const c6d6Datatype& c6d6Data) {
    _c6d6Data = c6d6Data;
}

inline double x6(const double& C6, const double& d6,
                 const OpenMM::RealVec& p1, const OpenMM::RealVec& p2,
                 OpenMM::RealVec& g1,       OpenMM::RealVec& g2)
{
    const double dx = (p1[0] - p2[0]);
    const double dy = (p1[1] - p2[1]);
    const double dz = (p1[2] - p2[2]);

    const double rsq = dx*dx + dy*dy + dz*dz;
    const double r = std::sqrt(rsq);

    const double d6r = d6*r/nm_to_A;
    const double tt6 = tang_toennies(6, d6r);

    const double inv_rsq = 1.0/rsq;
    const double inv_r6 = inv_rsq*inv_rsq*inv_rsq;

    const double e6 = C6*tt6*inv_r6/kcal_permol_Aminus6_to_kJ_permol_nmminus6;

    const double if6 = 1.0/Factorial<6>::value;

    const double grd = 6*e6*inv_rsq - C6/kcal_permol_Aminus6_to_kJ_permol_nmminus6*std::pow(d6/nm_to_A, 7)*if6*std::exp(-d6r)/r;

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
                                                      const std::vector<string>& allParticleElements,
                                                      vector<RealVec>& forces ) const {

        // FIXME improve exclusion of virtual sites and same-molecule interactions

        if ((allParticleElements[siteI] == "M") or (allParticleElements[siteJ] == "M"))
            return 0;
        if ((std::abs(siteI-siteJ) == 1) and (allParticleElements[siteI] != "O"))
            return 0;
        if ((std::abs(siteI-siteJ) == 2) and (allParticleElements[siteJ] == "O"))
            return 0;

        std::vector<RealVec> allPositions;

        allPositions.push_back(particlePositions[siteI] * nm_to_A);
        allPositions.push_back(particlePositions[siteJ] * nm_to_A);

        if( _nonbondedMethod == CutoffPeriodic )
            imageParticles(_periodicBoxDimensions, allPositions[0], allPositions[1]);

        c6d6Datatype :: const_iterator entry =
                _c6d6Data.find(make_pair(allParticleElements[siteI], allParticleElements[siteJ]));
        if (entry == _c6d6Data.end())
            entry = _c6d6Data.find(make_pair(allParticleElements[siteJ], allParticleElements[siteI]));
        if (entry == _c6d6Data.end())
            throw OpenMMException("Dispersion force parameters need to be defined for all pairs of elements: " +
                    allParticleElements[siteI] + ", " + allParticleElements[siteJ] + " missing");

        pair<double, double> c6d6 = entry->second;

        RealOpenMM energy= x6(c6d6.first, c6d6.second, allPositions[0], allPositions[1], forces[siteI], forces[siteJ]);

        return energy;

}

RealOpenMM MBPolReferenceDispersionForce::calculateForceAndEnergy( int numParticles,
                                                             const vector<RealVec>& particlePositions,
                                                             const std::vector<string>& allParticleElements,
                                                             const NeighborList& neighborList,
                                                             vector<RealVec>& forces ) const {


    RealOpenMM energy = 0.;
    for( unsigned int ii = 0; ii < neighborList.size(); ii++ ){

        OpenMM::AtomPair pair       = neighborList[ii];
        int siteI                   = pair.first;
        int siteJ                   = pair.second;

        energy                     += calculatePairIxn( siteI, siteJ,
                particlePositions, allParticleElements, forces );

    }

    return energy;
}
