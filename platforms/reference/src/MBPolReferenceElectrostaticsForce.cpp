
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

#include "MBPolReferenceElectrostaticsForce.h"
#include <algorithm>

// In case we're using some primitive version of Visual Studio this will
// make sure that erf() and erfc() are defined.
#include "openmm/internal/MSVC_erfc.h"

// from mbpol
#include "gammq.h"

using std::vector;
using OpenMM::RealVec;

#undef MBPOL_DEBUG

MBPolReferenceElectrostaticsForce::MBPolReferenceElectrostaticsForce( ) :
                                                   _nonbondedMethod(NoCutoff),
                                                   _numParticles(0), 
                                                   _electric(138.9354558456),
                                                   _dielectric(1.0),
                                                   _mutualInducedDipoleConverged(0),
                                                   _mutualInducedDipoleIterations(0),
                                                   _maximumMutualInducedDipoleIterations(100),
                                                   _mutualInducedDipoleEpsilon(1.0e+50),
                                                   _mutualInducedDipoleTargetEpsilon(1.0e-04),
                                                   _polarSOR(0.55),
                                                   _debye(48.033324),
                                                   _includeChargeRedistribution(true)
{
    initialize();
}

MBPolReferenceElectrostaticsForce::MBPolReferenceElectrostaticsForce( NonbondedMethod nonbondedMethod ) :
                                                   _nonbondedMethod(NoCutoff),
                                                   _numParticles(0), 
                                                   _electric(138.9354558456),
                                                   _dielectric(1.0),
                                                   _mutualInducedDipoleConverged(0),
                                                   _mutualInducedDipoleIterations(0),
                                                   _maximumMutualInducedDipoleIterations(100),
                                                   _mutualInducedDipoleEpsilon(1.0e+50),
                                                   _mutualInducedDipoleTargetEpsilon(1.0e-04),
                                                   _polarSOR(0.55),
                                                   _debye(48.033324),
                                                   _includeChargeRedistribution(true)
{
    initialize();
}

void MBPolReferenceElectrostaticsForce::initialize( void )
{

    unsigned int index    = 0;
    _mScale[index++]      = 0.0;
    _mScale[index++]      = 0.0;
    _mScale[index++]      = 0.0;
    _mScale[index++]      = 0.4;
    _mScale[index++]      = 0.8;

    index                 = 0;
    _dScale[index++]      = 0.0;
    _dScale[index++]      = 1.0;
    _dScale[index++]      = 1.0;
    _dScale[index++]      = 1.0;
    _dScale[index++]      = 1.0;

    index                 = 0;
    _pScale[index++]      = 0.0;
    _pScale[index++]      = 0.0;
    _pScale[index++]      = 0.0;
    _pScale[index++]      = 1.0;
    _pScale[index++]      = 1.0;

    index                 = 0;
    _uScale[index++]      = 1.0;
    _uScale[index++]      = 1.0;
    _uScale[index++]      = 1.0;
    _uScale[index++]      = 1.0;
    _uScale[index++]      = 1.0;

    return;
}

MBPolReferenceElectrostaticsForce::NonbondedMethod MBPolReferenceElectrostaticsForce::getNonbondedMethod( void ) const 
{
    return _nonbondedMethod;
}

void MBPolReferenceElectrostaticsForce::setNonbondedMethod( MBPolReferenceElectrostaticsForce::NonbondedMethod nonbondedMethod )
{
    _nonbondedMethod = nonbondedMethod;
}

MBPolReferenceElectrostaticsForce::PolarizationType MBPolReferenceElectrostaticsForce::getPolarizationType( void ) const 
{
    return _polarizationType;
}

void MBPolReferenceElectrostaticsForce::setPolarizationType( MBPolReferenceElectrostaticsForce::PolarizationType polarizationType )
{
    _polarizationType = polarizationType;
}

void MBPolReferenceElectrostaticsForce::setIncludeChargeRedistribution( bool includeChargeRedistribution ) {
    _includeChargeRedistribution = includeChargeRedistribution;
}

bool MBPolReferenceElectrostaticsForce::getIncludeChargeRedistribution( void ) const
{
    return _includeChargeRedistribution;
}

int MBPolReferenceElectrostaticsForce::getMutualInducedDipoleConverged( void ) const 
{
    return _mutualInducedDipoleConverged;
}

void MBPolReferenceElectrostaticsForce::setMutualInducedDipoleConverged( int mutualInducedDipoleConverged )
{
    _mutualInducedDipoleConverged = mutualInducedDipoleConverged;
}

int MBPolReferenceElectrostaticsForce::getMutualInducedDipoleIterations( void ) const 
{
    return _mutualInducedDipoleIterations;
}

void MBPolReferenceElectrostaticsForce::setMutualInducedDipoleIterations( int mutualInducedDipoleIterations )
{
    _mutualInducedDipoleIterations = mutualInducedDipoleIterations;
}

RealOpenMM MBPolReferenceElectrostaticsForce::getMutualInducedDipoleEpsilon( void ) const 
{
    return _mutualInducedDipoleEpsilon;
}

void MBPolReferenceElectrostaticsForce::setMutualInducedDipoleEpsilon( RealOpenMM mutualInducedDipoleEpsilon )
{
    _mutualInducedDipoleEpsilon = mutualInducedDipoleEpsilon;
}

int MBPolReferenceElectrostaticsForce::getMaximumMutualInducedDipoleIterations( void ) const 
{
    return _maximumMutualInducedDipoleIterations;
}

void MBPolReferenceElectrostaticsForce::setMaximumMutualInducedDipoleIterations( int maximumMutualInducedDipoleIterations )
{
    _maximumMutualInducedDipoleIterations = maximumMutualInducedDipoleIterations;
}

RealOpenMM MBPolReferenceElectrostaticsForce::getMutualInducedDipoleTargetEpsilon( void ) const 
{
    return _mutualInducedDipoleTargetEpsilon;
}

void MBPolReferenceElectrostaticsForce::setMutualInducedDipoleTargetEpsilon( RealOpenMM mutualInducedDipoleTargetEpsilon )
{
    _mutualInducedDipoleTargetEpsilon = mutualInducedDipoleTargetEpsilon;
}

void MBPolReferenceElectrostaticsForce::setupScaleMaps( const std::vector< std::vector< std::vector<int> > >& multipoleParticleCovalentInfo )
{

    /* Setup for scaling maps:
     *
     *     _scaleMaps[particleIndex][ScaleType] = map, where map[covalentIndex] = scaleFactor 
     *     _maxScaleIndex[particleIndex]        = max covalent index for particleIndex
     *
     *     multipoleParticleCovalentInfo[ii][jj], jj =0,1,2,3 contains covalent indices (c12, c13, c14, c15)
     *     multipoleParticleCovalentInfo[ii][jj], jj =4,5,6,7 contains covalent indices (p11, p12, p13, p14)
     *
     *     only including covalent particles w/ index >= ii
     */

    _scaleMaps.resize( multipoleParticleCovalentInfo.size() );
    _maxScaleIndex.resize( multipoleParticleCovalentInfo.size() );

    for( unsigned int ii = 0; ii < multipoleParticleCovalentInfo.size(); ii++ ){

        _scaleMaps[ii].resize(LAST_SCALE_TYPE_INDEX);
        _maxScaleIndex[ii] = 0;
        const std::vector< std::vector<int> >& covalentInfo = multipoleParticleCovalentInfo[ii];
        const std::vector<int> covalentListP11              = covalentInfo[MBPolElectrostaticsForce::PolarizationCovalent11];

        // pScale & mScale

        for( unsigned jj = 0; jj < MBPolElectrostaticsForce::PolarizationCovalent11; jj++ ){
            const std::vector<int> covalentList    = covalentInfo[jj];
            for( unsigned int kk = 0; kk < covalentList.size(); kk++ ){
                unsigned int covalentIndex             = static_cast<unsigned int>(covalentList[kk]);
                if( covalentIndex < ii )continue;

                // handle 0.5 factor for p14

                int hit = 0;
                if( jj == MBPolElectrostaticsForce::Covalent14 ){
                    for( unsigned int mm = 0; mm < covalentListP11.size() && hit == 0; mm++ ){
                        if( covalentListP11[mm]  == covalentIndex ){
                            hit = 1;
                        }
                    }
                } 
               
                _scaleMaps[ii][P_SCALE][covalentIndex] = hit ? 0.5*_pScale[jj+1] : _pScale[jj+1];
                _scaleMaps[ii][M_SCALE][covalentIndex] = _mScale[jj+1];
                _maxScaleIndex[ii]                     = _maxScaleIndex[ii] < covalentIndex ? covalentIndex : _maxScaleIndex[ii];
            }
        }

        // dScale & uScale

        for( unsigned jj = MBPolElectrostaticsForce::PolarizationCovalent11; jj < covalentInfo.size(); jj++ ){
            const std::vector<int> covalentList = covalentInfo[jj];
            for( unsigned int kk = 0; kk < covalentList.size(); kk++ ){
                unsigned int covalentIndex             = static_cast<unsigned int>(covalentList[kk]);
                if( covalentIndex < ii )continue;
                _scaleMaps[ii][D_SCALE][covalentIndex] = _dScale[jj-4];
                _scaleMaps[ii][U_SCALE][covalentIndex] = _uScale[jj-4];
                _maxScaleIndex[ii]                     = _maxScaleIndex[ii] < covalentIndex ? covalentIndex : _maxScaleIndex[ii];
            }
        }
    }

    //showScaleMapForParticle( 2, stderr );
    //showScaleMapForParticle( 10, stderr );

    return;
}

void MBPolReferenceElectrostaticsForce::showScaleMapForParticle( unsigned int particleI, FILE* log ) const 
{

#ifdef MBPOL_DEBUG
    (void) fprintf( log, "Scale map particle %5u maxIndex=%u\n", particleI, _maxScaleIndex[particleI] );

    std::string scaleNames[LAST_SCALE_TYPE_INDEX] = { "D", "P", "M" }; 
    for( unsigned int ii = 0; ii < _scaleMaps[particleI].size(); ii++ ){
        MapIntRealOpenMM scaleMap = _scaleMaps[particleI][ii];
        (void) fprintf( log, "  %s scale ", scaleNames[ii].c_str() );
        for( MapIntRealOpenMMCI jj = scaleMap.begin(); jj != scaleMap.end(); jj++ ){
            //if( jj->first > particleI && jj->second < 1.0 )
            if( jj->second < 1.0 )
            (void) fprintf( log, "%4d=%5.2f ", jj->first, jj->second );
        }    
        (void) fprintf( log, "\n" );
    }    
    (void) fprintf( log, "\n" );
    (void) fflush( log );
#endif
}

RealOpenMM MBPolReferenceElectrostaticsForce::getElectrostaticsScaleFactor( unsigned int particleI, unsigned int particleJ, ScaleType scaleType ) const 
{

    MapIntRealOpenMM  scaleMap   = _scaleMaps[particleI][scaleType];
    MapIntRealOpenMMCI isPresent = scaleMap.find( particleJ );
    if( isPresent != scaleMap.end() ){
        return isPresent->second;
    } else {
        return 1.0;
    }
}

void MBPolReferenceElectrostaticsForce::getDScaleAndPScale( unsigned int particleI, unsigned int particleJ, RealOpenMM& dScale, RealOpenMM& pScale ) const 
{
    dScale = getElectrostaticsScaleFactor( particleI, particleJ, D_SCALE );
    pScale = getElectrostaticsScaleFactor( particleI, particleJ, P_SCALE );
}

void MBPolReferenceElectrostaticsForce::getElectrostaticsScaleFactors( unsigned int particleI, unsigned int particleJ, std::vector<RealOpenMM>& scaleFactors ) const 
{
    scaleFactors[D_SCALE] = getElectrostaticsScaleFactor( particleI, particleJ, D_SCALE );
    scaleFactors[P_SCALE] = getElectrostaticsScaleFactor( particleI, particleJ, P_SCALE );
    scaleFactors[M_SCALE] = getElectrostaticsScaleFactor( particleI, particleJ, M_SCALE );
    scaleFactors[U_SCALE] = getElectrostaticsScaleFactor( particleI, particleJ, U_SCALE );
}

RealOpenMM MBPolReferenceElectrostaticsForce::normalizeRealVec( RealVec& vectorToNormalize ) const 
{
    RealOpenMM norm = SQRT( vectorToNormalize.dot( vectorToNormalize ) );
    if( norm > 0.0 ){
        vectorToNormalize *= (1.0/norm);
    }
    return norm;
}

void MBPolReferenceElectrostaticsForce::initializeRealOpenMMVector( vector<RealOpenMM>& vectorToInitialize ) const 
{
    RealOpenMM zero = 0.0;
    vectorToInitialize.resize( _numParticles );
    std::fill( vectorToInitialize.begin(), vectorToInitialize.end(), zero );
}

void MBPolReferenceElectrostaticsForce::initializeRealVecVector( vector<RealVec>& vectorToInitialize ) const 
{
    vectorToInitialize.resize( _numParticles );
    RealVec zeroVec( 0.0, 0.0, 0.0 );
    std::fill( vectorToInitialize.begin(), vectorToInitialize.end(), zeroVec );
}

void MBPolReferenceElectrostaticsForce::copyRealVecVector( const std::vector<OpenMM::RealVec>& inputVector, std::vector<OpenMM::RealVec>& outputVector ) const 
{ 
    outputVector.resize( inputVector.size() );
    for( unsigned int ii = 0; ii < inputVector.size(); ii++ ){
        outputVector[ii] = inputVector[ii];
    }   
    return;
}

void MBPolReferenceElectrostaticsForce::loadParticleData( const std::vector<RealVec>& particlePositions,
                                                      const std::vector<RealOpenMM>& charges,
                                                      const std::vector<RealOpenMM>& dipoles,
                                                      const std::vector<RealOpenMM>& quadrupoles,
                                                      const std::vector<RealOpenMM>& tholes,
                                                      const std::vector<RealOpenMM>& dampingFactors,
                                                      const std::vector<RealOpenMM>& polarity,
                                                      const std::vector<int>& multipoleAtomZs,
                                                      const std::vector<int>& multipoleAtomXs,
                                                      const std::vector<int>& multipoleAtomYs,
                                                      std::vector<ElectrostaticsParticleData>& particleData ) const 
{
   
    particleData.resize( _numParticles );
    for( unsigned int ii = 0; ii < _numParticles; ii++ ){

        particleData[ii].particleIndex        = ii;

        particleData[ii].position             = particlePositions[ii];
        particleData[ii].charge               = charges[ii];

        particleData[ii].dipole[0]            = dipoles[3*ii+0];
        particleData[ii].dipole[1]            = dipoles[3*ii+1];
        particleData[ii].dipole[2]            = dipoles[3*ii+2];

        particleData[ii].quadrupole[QXX]      = quadrupoles[9*ii+0];
        particleData[ii].quadrupole[QXY]      = quadrupoles[9*ii+1];
        particleData[ii].quadrupole[QXZ]      = quadrupoles[9*ii+2];
        particleData[ii].quadrupole[QYY]      = quadrupoles[9*ii+4];
        particleData[ii].quadrupole[QYZ]      = quadrupoles[9*ii+5];
        particleData[ii].quadrupole[QZZ]      = quadrupoles[9*ii+8];

        particleData[ii].thole[TCC]           = tholes[5*ii+0];
        particleData[ii].thole[TCD]           = tholes[5*ii+1];
        particleData[ii].thole[TDD]          = tholes[5*ii+2];
        particleData[ii].thole[TDDOH]         = tholes[5*ii+3];
        particleData[ii].thole[TDDHH]         = tholes[5*ii+4];
        particleData[ii].dampingFactor        = dampingFactors[ii];
        particleData[ii].polarity             = polarity[ii];

        particleData[ii].multipoleAtomZs = multipoleAtomZs[ii];
        particleData[ii].multipoleAtomXs = multipoleAtomXs[ii];
        particleData[ii].multipoleAtomYs = multipoleAtomYs[ii];

    }
}

void MBPolReferenceElectrostaticsForce::zeroFixedElectrostaticsFields( void )
{
    initializeRealVecVector( _fixedElectrostaticsField );
    initializeRealVecVector( _fixedElectrostaticsFieldPolar );
}

void MBPolReferenceElectrostaticsForce::checkChiralCenterAtParticle( ElectrostaticsParticleData& particleI, int axisType,
                                                                 ElectrostaticsParticleData& particleZ, ElectrostaticsParticleData& particleX, 
                                                                 ElectrostaticsParticleData& particleY ) const 
{

    if( axisType == MBPolElectrostaticsForce::ZThenX ){
        return;
    }
        
    RealVec deltaAD   = particleI.position - particleY.position;
    RealVec deltaBD   = particleZ.position - particleY.position;
    RealVec deltaCD   = particleX.position - particleY.position;

    RealVec deltaC    = deltaBD.cross( deltaCD );
    RealOpenMM volume = deltaC.dot( deltaAD );

    if( volume < 0.0 ){
        particleI.dipole[1]         *= -1.0; // pole(3,i)
        particleI.quadrupole[QXY]   *= -1.0; // pole(6,i)  && pole(8,i)
        particleI.quadrupole[QYZ]   *= -1.0; // pole(10,i) && pole(12,i)
    }
    return;

}

void MBPolReferenceElectrostaticsForce::checkChiral( std::vector<ElectrostaticsParticleData>& particleData,
                                                 const std::vector<int>& multipoleAtomXs,
                                                 const std::vector<int>& multipoleAtomYs,
                                                 const std::vector<int>& multipoleAtomZs,
                                                 const std::vector<int>& axisTypes ) const 
{

    for( unsigned int ii = 0; ii < _numParticles; ii++ ){
        if( multipoleAtomYs[ii] > -1 ){
            checkChiralCenterAtParticle( particleData[ii], axisTypes[ii],
                                         particleData[multipoleAtomZs[ii]],
                                         particleData[multipoleAtomXs[ii]],
                                         particleData[multipoleAtomYs[ii]] );
        }
    }
    return;
}

void MBPolReferenceElectrostaticsForce::applyRotationMatrixToParticle(       ElectrostaticsParticleData& particleI,
                                                                   const ElectrostaticsParticleData& particleZ,
                                                                   const ElectrostaticsParticleData& particleX,
                                                                         ElectrostaticsParticleData* particleY,
                                                                         int axisType ) const 
{

    // handle case where rotation matrix is identity (e.g. single ion)

    // get coordinates of this atom and the z & x axis atoms
    // compute the vector between the atoms and 1/sqrt(d2), d2 is distance between
    // this atom and the axis atom

    RealVec vectorY;
    RealVec vectorZ = particleZ.position - particleI.position;
    RealVec vectorX = particleX.position - particleI.position;

    normalizeRealVec( vectorZ );
 
    // branch based on axis type
 
    if( axisType == MBPolElectrostaticsForce::Bisector ){
 
        // bisector
  
        // dx = dx1 + dx2 (in TINKER code)
       
        normalizeRealVec( vectorX );
        vectorZ      += vectorX;
        normalizeRealVec( vectorZ );
       
    } else if( axisType == MBPolElectrostaticsForce::ZBisect ){
 
        // z-bisect
  
        // dx = dx1 + dx2 (in TINKER code)
       
        normalizeRealVec( vectorX );

        vectorY  = particleY->position - particleI.position;
        normalizeRealVec( vectorY );

        vectorX += vectorY;
        normalizeRealVec( vectorX );
       
    } else if( axisType == MBPolElectrostaticsForce::ThreeFold ){
 
        // 3-fold
  
        // dx = dx1 + dx2 + dx3 (in TINKER code)
       
        normalizeRealVec( vectorX );

        vectorY   = particleY->position - particleI.position;
        normalizeRealVec( vectorY );

        vectorZ  += vectorX +  vectorY;
        normalizeRealVec( vectorZ );
       
    } else if( axisType == MBPolElectrostaticsForce::ZOnly ){
 
        // z-only
  
        vectorX = RealVec( 0.1, 0.1, 0.1 );

    }
 
    RealOpenMM dot      = vectorZ.dot( vectorX );
    vectorX            -= vectorZ*dot;

    normalizeRealVec( vectorX );
    vectorY = vectorZ.cross( vectorX );

    RealVec rotationMatrix[3];
    rotationMatrix[0] = vectorX;
    rotationMatrix[1] = vectorY;
    rotationMatrix[2] = vectorZ; 

    RealVec labDipole;
    for( int ii = 0; ii < 3; ii++ ){
        labDipole[ii] = particleI.dipole[0]*rotationMatrix[0][ii];
        for( int jj = 1; jj < 3; jj++ ){
            labDipole[ii] += particleI.dipole[jj]*rotationMatrix[jj][ii];
        }
    }
    particleI.dipole = labDipole;
 
    RealOpenMM mPole[3][3];
    RealOpenMM rPole[3][3] = { { 0.0, 0.0, 0.0 },
                               { 0.0, 0.0, 0.0 },
                               { 0.0, 0.0, 0.0 } };

    mPole[0][0] = particleI.quadrupole[QXX];
    mPole[0][1] = particleI.quadrupole[QXY];
    mPole[0][2] = particleI.quadrupole[QXZ];

    mPole[1][0] = particleI.quadrupole[QXY];
    mPole[1][1] = particleI.quadrupole[QYY];
    mPole[1][2] = particleI.quadrupole[QYZ];

    mPole[2][0] = particleI.quadrupole[QXZ];
    mPole[2][1] = particleI.quadrupole[QYZ];
    mPole[2][2] = particleI.quadrupole[QZZ];
 
    for( int ii = 0; ii < 3; ii++ ){
       for( int jj = ii; jj < 3; jj++ ){
          for( int kk = 0; kk < 3; kk++ ){
             for( int mm = 0; mm < 3; mm++ ){
                 rPole[ii][jj] += rotationMatrix[kk][ii]*rotationMatrix[mm][jj]*mPole[kk][mm];
             }
          }
       }
    }
 
    particleI.quadrupole[QXX] = rPole[0][0];
    particleI.quadrupole[QXY] = rPole[0][1];
    particleI.quadrupole[QXZ] = rPole[0][2];

    particleI.quadrupole[QYY] = rPole[1][1];
    particleI.quadrupole[QYZ] = rPole[1][2];
    particleI.quadrupole[QZZ] = rPole[2][2];

    return;

}

void MBPolReferenceElectrostaticsForce::applyRotationMatrix( std::vector<ElectrostaticsParticleData>& particleData,
                                                         const std::vector<int>& multipoleAtomXs,
                                                         const std::vector<int>& multipoleAtomYs,
                                                         const std::vector<int>& multipoleAtomZs,
                                                         const std::vector<int>& axisTypes ) const 
{

    for( unsigned int ii = 0; ii < _numParticles; ii++ ){
        if( multipoleAtomZs[ii] >= 0 && multipoleAtomXs[ii] >= 0 ){
            applyRotationMatrixToParticle( particleData[ii], particleData[multipoleAtomZs[ii]], particleData[multipoleAtomXs[ii]],
                                           multipoleAtomYs[ii] > -1 ? &particleData[multipoleAtomYs[ii]] : NULL, axisTypes[ii] );
        }
    }

    return;
}

RealOpenMM MBPolReferenceElectrostaticsForce::getAndScaleInverseRs(  const ElectrostaticsParticleData& particleI,
                                                                    const ElectrostaticsParticleData& particleK,
                                                          RealOpenMM r, bool justScale, int interactionOrder, int interactionType) const
{

    bool isSameWater = (particleI.multipoleAtomZs == particleK.particleIndex) or
                (particleI.multipoleAtomYs == particleK.particleIndex) or
                (particleI.multipoleAtomXs == particleK.particleIndex);

    // MB-Pol has additional charge-charge term:
    // rrI[1] = charge-charge (ts0 in mbpol)
    // rrI[3] = dipole-charge (ts1 in mbpol)
    // rrI[5] = dipole-dipole (ts2 in mbpol)
    // rrI[7] = quadrupole

    RealOpenMM rrI;
    if (justScale) {
        rrI         = 1.;
    } else {
        RealOpenMM rI             =  1.0/r;
        RealOpenMM r2I            =  rI*rI;
        rrI                    = rI;
        for( unsigned int ii  = 3; ii <= interactionOrder; ii=ii+2 ){
           rrI *= (ii-2)*r2I;
        }
    }

    RealOpenMM pgamma = std::min(particleI.thole[interactionType], particleK.thole[interactionType]);

    if (interactionType == TDD) {
        // dipole - dipole thole parameters is different for 3 cases:
        // 1) different water molecules : thole[TDD]
        // 2) same water hydrogen-hydrogen : thole[TDDHH]
        // 3) same water hydrogen-oxygen : thole[TDDOH]

        if (isSameWater) {
            // FIXME improve identification of oxygen, now relies only on particles order
            bool oneIsOxygen = ((particleI.multipoleAtomZs >= particleI.particleIndex) and
                                (particleI.multipoleAtomYs >= particleI.particleIndex) and
                                (particleI.multipoleAtomXs >= particleI.particleIndex)) or
                                ((particleK.multipoleAtomZs >= particleK.particleIndex) and
                                 (particleK.multipoleAtomYs >= particleK.particleIndex) and
                                 (particleK.multipoleAtomXs >= particleK.particleIndex));
            if (oneIsOxygen) {
                pgamma = std::min(particleI.thole[TDDOH],particleK.thole[TDDOH]);
            } else {
                pgamma = std::min(particleI.thole[TDDHH],particleK.thole[TDDHH]);
            }
        } else {
            pgamma = std::min(particleI.thole[TDD],particleK.thole[TDD]);
        }
    }

    RealOpenMM damp      = pow(particleI.dampingFactor*particleK.dampingFactor, 1/6.); // AA in MBPol

    if( (damp != 0.0) and ( damp > -50.0 ) ) { // damp or not

        RealOpenMM ratio       = pow(r/damp, 4); // rA4 in MBPol
        RealOpenMM dampForExp = -1 * pgamma * ratio;

        switch (interactionOrder) {

        case 1:
            return rrI * (1.0 - EXP(dampForExp) + pow(pgamma, 1.0/4.0)*(r/damp)*EXP(ttm::gammln(3.0/4.0))*ttm::gammq(3.0/4.0, -dampForExp));

        case 3:
            return rrI * ( 1.0 - EXP(dampForExp) );

        case 5:
            return rrI * ((1.0 - EXP(dampForExp)) - (4./3.) * pgamma * EXP(dampForExp) * ratio);

        case 7:
            return rrI * (((1.0 - EXP(dampForExp)) - (4./3.) * pgamma * EXP(dampForExp) * ratio) -
                                (4./15.) * pgamma * (4. * pgamma * ratio - 1.) * EXP(dampForExp) / pow(damp, 4) * pow(r, 4));
        }
    }

    return rrI;
}

void MBPolReferenceElectrostaticsForce::calculateFixedElectrostaticsFieldPairIxn( const ElectrostaticsParticleData& particleI,
                                                                         const ElectrostaticsParticleData& particleJ,
                                                                         RealOpenMM dScale, RealOpenMM pScale ) 
{

    if( particleI.particleIndex == particleJ.particleIndex )return;

    // in MBPol there is no contribution to the Fixed Electrostatics Field 
    // from atoms of the same water molecule. multipoleAtomZs is used for 
    // defining a reference frame for the water molecules and
    // contains the indices to the other 2 atoms in the same water molecule.

    bool isSameWater = (particleI.multipoleAtomZs == particleJ.particleIndex) or
            (particleI.multipoleAtomYs == particleJ.particleIndex) or
            (particleI.multipoleAtomXs == particleJ.particleIndex);
    if( isSameWater )return;

    RealVec deltaR    = particleJ.position - particleI.position;
    RealOpenMM r      = SQRT( deltaR.dot( deltaR ) );
 
    // get scaling factors, if needed
     
    // charge - charge
    RealOpenMM rr3 = getAndScaleInverseRs( particleI, particleJ,r,false,3,TCC);

    // field at particle I due multipoles at particle J

    RealOpenMM factor                           = rr3*particleJ.charge;

    RealVec field                               = deltaR*factor;

    unsigned int particleIndex                  = particleI.particleIndex;
    _fixedElectrostaticsField[particleIndex]        -= field*dScale;
    _fixedElectrostaticsFieldPolar[particleIndex]   -= field*pScale;
 
    // field at particle J due multipoles at particle I

    factor                                      = rr3*particleI.charge;
 
    field                                       = deltaR*factor;
    particleIndex                               = particleJ.particleIndex;
    _fixedElectrostaticsField[particleIndex]        += field*dScale;
    _fixedElectrostaticsFieldPolar[particleIndex]   += field*pScale;
 
    return;
}

void MBPolReferenceElectrostaticsForce::calculateFixedElectrostaticsField( const vector<ElectrostaticsParticleData>& particleData )
{

    // calculate fixed multipole fields

    // loop includes diagonal term ii == jj for GK ixn; other calculateFixedElectrostaticsFieldPairIxn() methods
    // skip calculations for this case

    for( unsigned int ii = 0; ii < _numParticles; ii++ ){
        for( unsigned int jj = ii; jj < _numParticles; jj++ ){

            // if site jj is less than max covalent scaling index then get/apply scaling constants
            // otherwise add unmodified field and fieldPolar to particle fields 

            RealOpenMM dScale, pScale;
//            if( jj <= _maxScaleIndex[ii] ){
//                getDScaleAndPScale( ii, jj, dScale, pScale );
//            } else {
                dScale = pScale = 1.0;
            //}
            calculateFixedElectrostaticsFieldPairIxn( particleData[ii], particleData[jj], dScale, pScale );
        }
    }
    return;
}

void MBPolReferenceElectrostaticsForce::initializeInducedDipoles( std::vector<UpdateInducedDipoleFieldStruct>& updateInducedDipoleFields )
{

    // initialize inducedDipoles

    _inducedDipole.resize( _numParticles );
    _inducedDipolePolar.resize( _numParticles );

    for( unsigned int ii = 0; ii < _numParticles; ii++ ){
        _inducedDipole[ii]       = _fixedElectrostaticsField[ii];
        _inducedDipolePolar[ii]  = _fixedElectrostaticsFieldPolar[ii];
    }

    return;
}

void MBPolReferenceElectrostaticsForce::calculateInducedDipolePairIxn( unsigned int particleI, 
                                                                   unsigned int particleJ,
                                                                   RealOpenMM rr3,
                                                                   RealOpenMM rr5,
                                                                   const RealVec& deltaR,
                                                                   const std::vector<RealVec>& inducedDipole,
                                                                   std::vector<RealVec>& field ) const 
{

    RealOpenMM dDotDelta            = rr5*(inducedDipole[particleJ].dot( deltaR ) );
    field[particleI]               += inducedDipole[particleJ]*rr3 + deltaR*dDotDelta;
    dDotDelta                       = rr5*(inducedDipole[particleI].dot( deltaR ) );
    field[particleJ]               += inducedDipole[particleI]*rr3 + deltaR*dDotDelta;

    return;
}

void MBPolReferenceElectrostaticsForce::calculateInducedDipolePairIxns( const ElectrostaticsParticleData& particleI, 
                                                                    const ElectrostaticsParticleData& particleJ,
                                                                    std::vector<UpdateInducedDipoleFieldStruct>& updateInducedDipoleFields )
{

   if( particleI.particleIndex == particleJ.particleIndex )return;

    RealVec deltaR       = particleJ.position - particleI.position;
    RealOpenMM r         =  SQRT( deltaR.dot( deltaR ) );
  
    RealOpenMM scale3 = getAndScaleInverseRs( particleI, particleJ, r, false, 3, TDD);
    RealOpenMM scale5 = getAndScaleInverseRs( particleI, particleJ, r, false, 5, TDD);
 
    for( unsigned int ii = 0; ii < updateInducedDipoleFields.size(); ii++ ){
        calculateInducedDipolePairIxn( particleI.particleIndex, particleJ.particleIndex, -scale3, scale5, deltaR,
                                       *(updateInducedDipoleFields[ii].inducedDipoles), updateInducedDipoleFields[ii].inducedDipoleField );
    }
    return;

}

void MBPolReferenceElectrostaticsForce::calculateInducedDipoleFields( const std::vector<ElectrostaticsParticleData>& particleData,
                                                                  std::vector<UpdateInducedDipoleFieldStruct>& updateInducedDipoleFields)
{

    for( unsigned int ii = 0; ii < particleData.size(); ii++ ){
        for( unsigned int jj = ii; jj < particleData.size(); jj++ ){
            calculateInducedDipolePairIxns( particleData[ii], particleData[jj], updateInducedDipoleFields );
        }
    }
    return;
}

RealOpenMM MBPolReferenceElectrostaticsForce::updateInducedDipoleFields( const std::vector<ElectrostaticsParticleData>& particleData,
                                                                     std::vector<UpdateInducedDipoleFieldStruct>& updateInducedDipoleFields)
{

    // (1) zero fields
    // (2) calculate induced dipole pair ixns 
    // (3) update induced dipoles based on pair ixns and calculate/return convergence factor, maxEpsilon

    RealVec zeroVec( 0.0, 0.0, 0.0 );
    for( unsigned int ii = 0; ii < updateInducedDipoleFields.size(); ii++ ){
        std::fill( updateInducedDipoleFields[ii].inducedDipoleField.begin(), updateInducedDipoleFields[ii].inducedDipoleField.end(), zeroVec );
    }

    calculateInducedDipoleFields( particleData, updateInducedDipoleFields);

    RealOpenMM maxEpsilon = 0.0;
    for( unsigned int kk = 0; kk < updateInducedDipoleFields.size(); kk++ ){
        RealOpenMM epsilon = updateInducedDipole( particleData,
                                                  *(updateInducedDipoleFields[kk].fixedElectrostaticsField),
                                                    updateInducedDipoleFields[kk].inducedDipoleField,
                                                  *(updateInducedDipoleFields[kk].inducedDipoles) );
   
        maxEpsilon = epsilon > maxEpsilon ? epsilon : maxEpsilon;
    }

    return maxEpsilon;
}

RealOpenMM MBPolReferenceElectrostaticsForce::updateInducedDipole( const std::vector<ElectrostaticsParticleData>& particleData,
                                                               const std::vector<RealVec>& fixedElectrostaticsField,
                                                               const std::vector<RealVec>& inducedDipoleField,
                                                               std::vector<RealVec>& inducedDipole )
{

    RealOpenMM epsilon                    = 0.0;
    for( unsigned int ii = 0; ii < particleData.size(); ii++ ){
        RealVec    oldValue               = inducedDipole[ii];
        RealVec    newValue               = fixedElectrostaticsField[ii] + inducedDipoleField[ii]*particleData[ii].polarity;
        RealVec    delta                  = newValue - oldValue;
        inducedDipole[ii]                 = oldValue + delta*_polarSOR;
        epsilon                          += delta.dot( delta );
    }
    return epsilon;
}

void MBPolReferenceElectrostaticsForce::convergeInduceDipoles( const std::vector<ElectrostaticsParticleData>& particleData,
                                                           std::vector<UpdateInducedDipoleFieldStruct>& updateInducedDipoleField)
{

    bool done                 = false;
    setMutualInducedDipoleConverged( false );
    int iteration             = 0;
    RealOpenMM currentEpsilon = 1.0e+50;

    // loop until (1) induced dipoles are converged or
    //            (2) iterations == max iterations or
    //            (3) convergence factor (spsilon) increases 

    while( !done ){

        RealOpenMM epsilon = updateInducedDipoleFields( particleData, updateInducedDipoleField);   
                   epsilon = _polarSOR*_debye*SQRT( epsilon/( static_cast<RealOpenMM>(_numParticles) ) );

        if( epsilon < getMutualInducedDipoleTargetEpsilon() ){
            setMutualInducedDipoleConverged( true );
            done = true;
        } else if( currentEpsilon < epsilon || iteration >= getMaximumMutualInducedDipoleIterations() ){
            done = true;
        }

        currentEpsilon = epsilon;
        iteration++;
    }
    setMutualInducedDipoleEpsilon( currentEpsilon );
    setMutualInducedDipoleIterations( iteration );

    return;
}

void MBPolReferenceElectrostaticsForce::calculateInducedDipoles( const std::vector<ElectrostaticsParticleData>& particleData )
{

    // calculate fixed electric fields

    zeroFixedElectrostaticsFields();
    calculateFixedElectrostaticsField( particleData );

    // initialize inducedDipoles
    // if polarization type is 'Direct', then return after initializing; otherwise 
    // converge induced dipoles.

    for( unsigned int ii = 0; ii < _numParticles; ii++ ){
        _fixedElectrostaticsField[ii]      *= particleData[ii].polarity;
        _fixedElectrostaticsFieldPolar[ii] *= particleData[ii].polarity; 
    }

    _inducedDipole.resize( _numParticles );
    _inducedDipolePolar.resize( _numParticles );
    std::vector<UpdateInducedDipoleFieldStruct> updateInducedDipoleField;
    updateInducedDipoleField.push_back( UpdateInducedDipoleFieldStruct( &_fixedElectrostaticsField,       &_inducedDipole ) );
    updateInducedDipoleField.push_back( UpdateInducedDipoleFieldStruct( &_fixedElectrostaticsFieldPolar,  &_inducedDipolePolar ) );

    initializeInducedDipoles( updateInducedDipoleField );

    if( getPolarizationType() == MBPolReferenceElectrostaticsForce::Direct ){
        setMutualInducedDipoleConverged( true );
        return;
    }

    // UpdateInducedDipoleFieldStruct contains induced dipole, fixed multipole fields and fields
    // due to other induced dipoles at each site

    convergeInduceDipoles( particleData, updateInducedDipoleField );

    return;
}

RealOpenMM MBPolReferenceElectrostaticsForce::calculateElectrostaticPairIxn( const std::vector<ElectrostaticsParticleData>& particleData,
                                                                         unsigned int iIndex,
                                                                         unsigned int kIndex,
                                                                         const std::vector<RealOpenMM>& scalingFactors,
                                                                         std::vector<RealVec>& forces,
                                                                         std::vector<RealVec>& torque ) const 
{
    RealOpenMM temp3,temp5,temp7;
    RealOpenMM gl[9],gli[7],glip[7];
    RealOpenMM sc[10],sci[8],scip[8];
    RealOpenMM gf[7],gfi[6],gti[6];
 
    ElectrostaticsParticleData particleI = particleData[iIndex];
    ElectrostaticsParticleData particleK = particleData[kIndex];

    RealVec delta       = particleK.position - particleI.position;
    RealOpenMM r2       = delta.dot( delta );

    // set conversion factor, cutoff and switching coefficients

    RealOpenMM f        = _electric/_dielectric;

    // set scale factors for permanent multipole and induced terms

    RealOpenMM pdi      = particleI.dampingFactor;

    // apply Thole polarization damping to scale factors

    RealOpenMM r        = SQRT(r2);
    RealOpenMM rr1      = 1.0/r;
    RealOpenMM rr3      = rr1/r2;
    RealOpenMM rr5      = 3.0*rr3/r2;
    RealOpenMM rr7      = 5.0*rr5/r2;
 
    // calculate scalar products for induced components

    sci[1] = _inducedDipole[iIndex][0]*_inducedDipole[kIndex][0] 
	   + _inducedDipole[iIndex][1]*_inducedDipole[kIndex][1] 
	   + _inducedDipole[iIndex][2]*_inducedDipole[kIndex][2];
    sci[2] = _inducedDipole[iIndex][0]*delta[0] 
	   + _inducedDipole[iIndex][1]*delta[1] 
	   + _inducedDipole[iIndex][2]*delta[2];
    sci[3] = _inducedDipole[kIndex][0]*delta[0] 
	   + _inducedDipole[kIndex][1]*delta[1] 
	   + _inducedDipole[kIndex][2]*delta[2];

    scip[1] = _inducedDipole[iIndex][0]*_inducedDipolePolar[kIndex][0]
	    + _inducedDipole[iIndex][1]*_inducedDipolePolar[kIndex][1]
	    + _inducedDipole[iIndex][2]*_inducedDipolePolar[kIndex][2]
	    + _inducedDipolePolar[iIndex][0]*_inducedDipole[kIndex][0] 
	    + _inducedDipolePolar[iIndex][1]*_inducedDipole[kIndex][1]
	    + _inducedDipolePolar[iIndex][2]*_inducedDipole[kIndex][2];
    scip[2] = _inducedDipolePolar[iIndex][0]*delta[0] 
	    + _inducedDipolePolar[iIndex][1]*delta[1] 
	    + _inducedDipolePolar[iIndex][2]*delta[2];
    scip[3] = _inducedDipolePolar[kIndex][0]*delta[0] 
	    + _inducedDipolePolar[kIndex][1]*delta[1] 
	    + _inducedDipolePolar[kIndex][2]*delta[2];

    // calculate the gl functions for permanent components

    gl[0] = particleI.charge*particleK.charge;

    // calculate the gl functions for induced components

    gli[0] = particleK.charge*sci[2] - particleI.charge*sci[3];

    glip[0] = particleK.charge*scip[2] - particleI.charge*scip[3];

    bool isSameWater = (particleI.multipoleAtomZs == particleK.particleIndex) or
                       (particleI.multipoleAtomYs == particleK.particleIndex) or
	   	       (particleI.multipoleAtomXs == particleK.particleIndex);
    // Same water atoms have no charge/charge interaction and 
    // no induced-dipole/charge interaction
    if( isSameWater ) {
        gl[0] = 0.;
        gli[0] = 0.;
        glip[0] = 0.;
    }

    // compute the energy contributions for this interaction

    RealOpenMM scale1CC = getAndScaleInverseRs( particleI, particleK, r, true, 1, TCC);
    RealOpenMM scale3CD = getAndScaleInverseRs( particleI, particleK, r, true, 3, TCD);

    RealOpenMM energy =       rr1* gl[0]*scale1CC  ; // charge-charge
    energy           += 0.5*( rr3*gli[0]*scale3CD ); // charge - induced dipole
    energy           *= f;

    RealOpenMM scale3CC = getAndScaleInverseRs( particleI, particleK, r, true, 3, TCC);
    RealOpenMM scale5CD = getAndScaleInverseRs( particleI, particleK, r, true, 5, TCD);
    RealOpenMM scale5DD = getAndScaleInverseRs( particleI, particleK, r, true, 5, TDD);
    RealOpenMM scale7DD = getAndScaleInverseRs( particleI, particleK, r, true, 7, TDD);

    // intermediate variables for the permanent components
    gf[0] = rr3*gl[0]*scale3CC ; // charge -charge

    // intermediate variables for the induced components

    gfi[0] = 0.5 * rr5 *  gli[0]*scale5CD + // charge - induced dipole
             0.5 * rr5 * glip[0]*scale5CD + // charge - induced dipole
             0.5 * rr5 * scip[1]*scale5DD + // induced dipole - induced dipole
           - 0.5 * rr7 * (sci[2]*scip[3] + scip[2]*sci[3])*scale7DD; // induced dipole - induced dipole

    // get the permanent force components

    RealVec ftm2 = delta*gf[0];

    // get the induced force components

    RealVec ftm2i  = delta*gfi[0];

    ftm2i += ( _inducedDipolePolar[iIndex] *  sci[3] + // iPdipole_i * idipole_k
                    _inducedDipole[iIndex] * scip[3] +
               _inducedDipolePolar[kIndex] *  sci[2] + // iPdipole_k * idipole_i
	            _inducedDipole[kIndex] * scip[2]  ) * 0.5 * rr5 * scale5DD;
    
    // Same water atoms have no induced-dipole/charge interaction
    if (not( isSameWater )) {
	ftm2i += ( 
                   ( _inducedDipole[iIndex] + 
		     _inducedDipolePolar[iIndex] )*-particleK.charge +
                   ( _inducedDipole[kIndex] + 
		     _inducedDipolePolar[kIndex] )* particleI.charge
                 ) * 0.5 * rr3 * scale3CD;
    }

    // account for partially excluded induced interactions
    // not needed for MB-pol/water, but will be necessary
    // for larger systems with 1-4 interactions
    // FIXME check how to disable this in the xml

    // MBPol charge derivative terms

//    gE_elec[ih1 + k] += GRDQ(0, 0, k)*phi[4*n + 1]  // phi(h1)
//                      + GRDQ(0, 1, k)*phi[4*n + 2]  // phi(h2)
//                      + GRDQ(0, 2, k)*phi[4*n + 3]; // phi(M)
//
//    gE_elec[ih2 + k] += GRDQ(1, 0, k)*phi[4*n + 1]  // phi(h1)
//                      + GRDQ(1, 1, k)*phi[4*n + 2]  // phi(h2)
//                      + GRDQ(1, 2, k)*phi[4*n + 3]; // phi(M)
//
//    gE_elec[io + k] += GRDQ(2, 0, k)*phi[4*n + 1]  // phi(h1)
//                     + GRDQ(2, 1, k)*phi[4*n + 2]  // phi(h2)
//                     + GRDQ(2, 2, k)*phi[4*n + 3]; // phi(M)

    if (getIncludeChargeRedistribution() and (not (isSameWater))){

        double distanceK, distanceI, 
	       scale1I, scale1K, scale3I, scale3K, 
	       inducedDipoleI, inducedDipoleK;
	RealVec deltaI, deltaK;

        for (size_t s = 0; s < 3; ++s) {

            // vsH1f, vsH2f, vsMf

            deltaI = particleData[particleI.otherSiteIndex[s]].position
		   - particleK.position;
            distanceI = SQRT(deltaI.dot(deltaI));
            deltaK = particleData[particleK.otherSiteIndex[s]].position
		   - particleI.position;
            distanceK = SQRT(deltaK.dot(deltaK));

            scale1I = getAndScaleInverseRs( particleData[particleI.otherSiteIndex[s]], particleK, distanceI, true, 1, TCC );
            scale3I = getAndScaleInverseRs( particleData[particleI.otherSiteIndex[s]], particleK, distanceI, true, 3, TCC );

            scale1K = getAndScaleInverseRs( particleData[particleK.otherSiteIndex[s]], particleI, distanceK, true, 1, TCC );
            scale3K = getAndScaleInverseRs( particleData[particleK.otherSiteIndex[s]], particleI, distanceK, true, 3, TCC );

            inducedDipoleI = _inducedDipole[kIndex].dot(deltaI);
            inducedDipoleK = _inducedDipole[iIndex].dot(deltaK);

            for (size_t i = 0; i < 3; ++i) {

                ftm2[i] +=  scale1I * (1.0/distanceI) * particleI.chargeDerivatives[s][i] * particleK.charge; // charge - charge
                ftm2[i] -=  scale1K * (1.0/distanceK) * particleK.chargeDerivatives[s][i] * particleI.charge; // charge - charge

                ftm2i[i] += scale3I * pow(1.0/distanceI,3) * particleI.chargeDerivatives[s][i] * inducedDipoleI;// charge - charge
                ftm2i[i] -= scale3K * pow(1.0/distanceK,3) * particleK.chargeDerivatives[s][i] * inducedDipoleK;// charge - charge

            }

        }
    }

    RealVec force   = ftm2*scalingFactors[M_SCALE] + ftm2i;
            force  *= f;

    forces[iIndex] -= force;
    forces[kIndex] += force;

    return energy;
}

void MBPolReferenceElectrostaticsForce::mapTorqueToForceForParticle( const ElectrostaticsParticleData& particleI,
                                                                 const ElectrostaticsParticleData& particleU,
                                                                 const ElectrostaticsParticleData& particleV,
                                                                       ElectrostaticsParticleData* particleW,
                                                                       int axisType, const Vec3& torque,
                                                                       std::vector<RealVec>& forces ) const 
{
 
    static const int U                  = 0;
    static const int V                  = 1;
    static const int W                  = 2;
    static const int R                  = 3;
    static const int S                  = 4;
    static const int UV                 = 5;
    static const int UW                 = 6;
    static const int VW                 = 7;
    static const int UR                 = 8;
    static const int US                 = 9;
    static const int VS                 = 10;
    static const int WS                 = 11;
    static const int LastVectorIndex    = 12;
    
    static const int X                  = 0;
    static const int Y                  = 1;
    static const int Z                  = 2;
    static const int I                  = 3;
    
    RealOpenMM norms[LastVectorIndex];
    RealOpenMM angles[LastVectorIndex][2];

    // ---------------------------------------------------------------------------------------
 
    // get coordinates of this atom and the z & x axis atoms
    // compute the vector between the atoms and 1/sqrt(d2), d2 is distance between
    // this atom and the axis atom

    if( axisType == MBPolElectrostaticsForce::NoAxisType ){
        return;
    }

    RealVec vectorU = particleU.position - particleI.position;
    norms[U]        = normalizeRealVec( vectorU );

    RealVec vectorV = particleV.position - particleI.position;
    norms[V]        = normalizeRealVec( vectorV );

    RealVec vectorW;
    if( particleW && (axisType == MBPolElectrostaticsForce::ZBisect || axisType == MBPolElectrostaticsForce::ThreeFold) ){
         vectorW = particleW->position - particleI.position;
    } else {
         vectorW = vectorU.cross( vectorV );
    }
    norms[W]  = normalizeRealVec( vectorW );
 
    RealVec vectorUV, vectorUW, vectorVW;
    vectorUV = vectorV.cross( vectorU );
    vectorUW = vectorW.cross( vectorU );
    vectorVW = vectorW.cross( vectorV );
    
    norms[UV]                     = normalizeRealVec( vectorUV );
    norms[UW]                     = normalizeRealVec( vectorUW );
    norms[VW]                     = normalizeRealVec( vectorVW );

    // angles[][0] is cosine of angle
    // angles[][1] is sine   of angle

    angles[UV][0]                 = vectorU.dot( vectorV );
    angles[UV][1]                 = SQRT( 1.0 - angles[UV][0]*angles[UV][0]);
    
    angles[UW][0]                 = vectorU.dot( vectorW );
    angles[UW][1]                 = SQRT( 1.0 - angles[UW][0]*angles[UW][0]);

    angles[VW][0]                 = vectorV.dot( vectorW );
    angles[VW][1]                 = SQRT( 1.0 - angles[VW][0]*angles[VW][0]);

    RealVec dphi;
    dphi[U]                       = vectorU.dot( torque );
    dphi[V]                       = vectorV.dot( torque );
    dphi[W]                       = vectorW.dot( torque );
    dphi                         *= -1.0;

    // branch based on axis type
 
    if( axisType == MBPolElectrostaticsForce::ZThenX || axisType == MBPolElectrostaticsForce::Bisector ){
 
        RealOpenMM factor1;
        RealOpenMM factor2;
        RealOpenMM factor3;
        RealOpenMM factor4;
        RealOpenMM half = 0.5;
    
        factor1                 =  dphi[V]/(norms[U]*angles[UV][1]);
        factor2                 =  dphi[W]/(norms[U]);
        factor3                 = -dphi[U]/(norms[V]*angles[UV][1]);
    
        if( axisType == MBPolElectrostaticsForce::Bisector ){ 
            factor2    *= half;
            factor4     = half*dphi[W]/(norms[V]);
        } else {
            factor4     = 0.0;
        }
 
        for( int ii = 0; ii < 3; ii++ ){
            double forceU                                        =  vectorUV[ii]*factor1 + factor2*vectorUW[ii];
            forces[particleU.particleIndex][ii]                 -=  forceU;

            double forceV                                        =  vectorUV[ii]*factor3 + factor4*vectorVW[ii];
            forces[particleV.particleIndex][ii]                 -=  forceV;

            forces[particleI.particleIndex][ii]                 +=  (forceU + forceV);
        }

    } else if( axisType == MBPolElectrostaticsForce::ZBisect ){

        RealVec vectorR           = vectorV + vectorW; 
        RealVec vectorS           = vectorU.cross( vectorR );

        norms[R]                  = normalizeRealVec( vectorR );
        norms[S]                  = normalizeRealVec( vectorS );

        RealVec vectorUR          =  vectorR.cross( vectorU );
        RealVec vectorUS          =  vectorS.cross( vectorU );
        RealVec vectorVS          =  vectorS.cross( vectorV );
        RealVec vectorWS          =  vectorS.cross( vectorW );

        norms[UR]                 = normalizeRealVec( vectorUR );
        norms[US]                 = normalizeRealVec( vectorUS );
        norms[VS]                 = normalizeRealVec( vectorVS );
        norms[WS]                 = normalizeRealVec( vectorWS );

        angles[UR][0]             = vectorU.dot( vectorR );
        angles[UR][1]             = SQRT( 1.0 - angles[UR][0]*angles[UR][0]);

        angles[US][0]             = vectorU.dot( vectorS );
        angles[US][1]             = SQRT( 1.0 - angles[US][0]*angles[US][0]);

        angles[VS][0]             = vectorV.dot( vectorS );
        angles[VS][1]             = SQRT( 1.0 - angles[VS][0]*angles[VS][0]);

        angles[WS][0]             = vectorW.dot( vectorS );
        angles[WS][1]             = SQRT( 1.0 - angles[WS][0]*angles[WS][0]);
 
        RealVec t1                = vectorV - vectorS*angles[VS][0];
        RealVec t2                = vectorW - vectorS*angles[WS][0];

        RealOpenMM notUsed        = normalizeRealVec( t1 );
              notUsed             = normalizeRealVec( t2 );

        RealOpenMM ut1cos         = vectorU.dot( t1 );
        RealOpenMM ut1sin         = SQRT( 1.0 - ut1cos*ut1cos);

        RealOpenMM ut2cos         = vectorU.dot( t2 );
        RealOpenMM ut2sin         = SQRT( 1.0 - ut2cos*ut2cos);

        RealOpenMM dphiR          = vectorR.dot( torque )*(-1.0);
        RealOpenMM dphiS          = vectorS.dot( torque )*(-1.0);

        RealOpenMM factor1        = dphiR/(norms[U]*angles[UR][1]);
        RealOpenMM factor2        = dphiS/(norms[U]);
        RealOpenMM factor3        = dphi[U]/(norms[V]*(ut1sin+ut2sin));
        RealOpenMM factor4        = dphi[U]/(norms[W]*(ut1sin+ut2sin));

        RealVec forceU            =  vectorUR*factor1 + vectorUS*factor2;
        forces[particleU.particleIndex]        -= forceU;

        RealVec forceV            = (vectorS*angles[VS][1] - t1*angles[VS][0])*factor3;
        forces[particleV.particleIndex]        -= forceV;

        RealVec forceW            = (vectorS*angles[WS][1] - t2*angles[WS][0])*factor4;
        forces[particleW->particleIndex]       -= forceW;

        forces[particleI.particleIndex]        += (forceU + forceV + forceW);

    } else if( axisType == MBPolElectrostaticsForce::ThreeFold ){

        // 3-fold

        for( int ii = 0; ii < 3; ii++ ){

            RealOpenMM du =  vectorUW[ii]*dphi[W]/(norms[U]*angles[UW][1]) +
                             vectorUV[ii]*dphi[V]/(norms[U]*angles[UV][1]) -
                             vectorUW[ii]*dphi[U]/(norms[U]*angles[UW][1]) -
                             vectorUV[ii]*dphi[U]/(norms[U]*angles[UV][1]);

            RealOpenMM dv =  vectorVW[ii]*dphi[W]/(norms[V]*angles[VW][1]) -
                             vectorUV[ii]*dphi[U]/(norms[V]*angles[UV][1]) -
                             vectorVW[ii]*dphi[V]/(norms[V]*angles[VW][1]) +
                             vectorUV[ii]*dphi[V]/(norms[V]*angles[UV][1]);

            RealOpenMM dw = -vectorUW[ii]*dphi[U]/(norms[W]*angles[UW][1]) -
                             vectorVW[ii]*dphi[V]/(norms[W]*angles[VW][1]) +
                             vectorUW[ii]*dphi[W]/(norms[W]*angles[UW][1]) +
                             vectorVW[ii]*dphi[W]/(norms[W]*angles[VW][1]);

            du /= 3.0;
            dv /= 3.0;
            dw /= 3.0;

            forces[particleU.particleIndex][ii] -= du;
            forces[particleV.particleIndex][ii] -= dv;
            if( particleW )
                forces[particleW->particleIndex][ii] -= dw;
            forces[particleI.particleIndex][ii] += (du + dv + dw);
        }

    } else if( axisType == MBPolElectrostaticsForce::ZOnly ){

        // z-only

        for( int ii = 0; ii < 3; ii++ ){
            RealOpenMM du                               = vectorUV[ii]*dphi[V]/(norms[U]*angles[UV][1]);
            forces[particleU.particleIndex][ii]        -= du;
            forces[particleI.particleIndex][ii]        += du;
        }
    }
 
    return;
 
}

void MBPolReferenceElectrostaticsForce::mapTorqueToForce( std::vector<ElectrostaticsParticleData>& particleData,
                                                      const std::vector<int>& multipoleAtomXs,
                                                      const std::vector<int>& multipoleAtomYs,
                                                      const std::vector<int>& multipoleAtomZs,
                                                      const std::vector<int>& axisTypes,
                                                      std::vector<RealVec>& torques,
                                                      std::vector<RealVec>& forces ) const 
{

    // map torques to forces

    for( unsigned int ii = 0; ii < particleData.size(); ii++ ){
        if( axisTypes[ii] != MBPolElectrostaticsForce::NoAxisType ){
             mapTorqueToForceForParticle( particleData[ii],
                                          particleData[multipoleAtomZs[ii]], particleData[multipoleAtomXs[ii]],
                                          multipoleAtomYs[ii] > -1 ? &particleData[multipoleAtomYs[ii]] : NULL,
                                          axisTypes[ii], torques[ii], forces ); 
        }
    }
    return;
}

RealOpenMM MBPolReferenceElectrostaticsForce::calculateElectrostatic( const std::vector<ElectrostaticsParticleData>& particleData,
                                                                  std::vector<RealVec>& torques,
                                                                  std::vector<RealVec>& forces )
{

    RealOpenMM energy = 0.0;
    std::vector<RealOpenMM> scaleFactors(LAST_SCALE_TYPE_INDEX);
    for( unsigned int kk = 0; kk < scaleFactors.size(); kk++ ){
        scaleFactors[kk] = 1.0;
    }   

    // main loop over particle pairs

    for( unsigned int ii = 0; ii < particleData.size(); ii++ ){
        for( unsigned int jj = ii+1; jj < particleData.size(); jj++ ){

            if( jj <= _maxScaleIndex[ii] ){
                getElectrostaticsScaleFactors( ii, jj, scaleFactors);
            }

            energy += calculateElectrostaticPairIxn( particleData, ii, jj, scaleFactors, forces, torques );

            if( jj <= _maxScaleIndex[ii] ){
                for( unsigned int kk = 0; kk < LAST_SCALE_TYPE_INDEX; kk++ ){
                    scaleFactors[kk] = 1.0;
                }
            }
        }
    }

    return energy;
}

void MBPolReferenceElectrostaticsForce::setup( const std::vector<RealVec>& particlePositions,
                                           const std::vector<RealOpenMM>& charges,
                                           const std::vector<RealOpenMM>& dipoles,
                                           const std::vector<RealOpenMM>& quadrupoles,
                                           const std::vector<RealOpenMM>& tholes,
                                           const std::vector<RealOpenMM>& dampingFactors,
                                           const std::vector<RealOpenMM>& polarity,
                                           const std::vector<int>& axisTypes,
                                           const std::vector<int>& multipoleAtomZs,
                                           const std::vector<int>& multipoleAtomXs,
                                           const std::vector<int>& multipoleAtomYs,
                                           const std::vector< std::vector< std::vector<int> > >& multipoleAtomCovalentInfo,
                                           std::vector<ElectrostaticsParticleData>& particleData )
{


    // load particle parameters into vector of ElectrostaticsParticleData
    // check for inverted chiral centers
    // apply rotation matrix to get lab frame dipole and quadrupoles
    // setup scaling factors
    // get induced dipoles
    // check if induced dipoles converged

    _numParticles = particlePositions.size();
    loadParticleData( particlePositions, charges, dipoles, quadrupoles,
                      tholes, dampingFactors, polarity, multipoleAtomZs, multipoleAtomXs, multipoleAtomYs, particleData );

    if (getIncludeChargeRedistribution())
    {
        for( unsigned int ii = 0; ii < _numParticles; ii=ii+4 ){ // FIXME this assumes only waters
            computeWaterCharge(particleData[ii], particleData[ii+1], particleData[ii+2], particleData[ii+3]);
        }
    }
    checkChiral( particleData, multipoleAtomXs, multipoleAtomYs, multipoleAtomZs, axisTypes );

    applyRotationMatrix( particleData, multipoleAtomXs, multipoleAtomYs, multipoleAtomZs, axisTypes );

    setupScaleMaps( multipoleAtomCovalentInfo );

    calculateInducedDipoles( particleData );

    if( !getMutualInducedDipoleConverged() ){
        std::stringstream message;
        message << "Induced dipoles did not converge: ";
        message << " iterations="      << getMutualInducedDipoleIterations();
        message << " eps="             << getMutualInducedDipoleEpsilon();
        throw OpenMMException(message.str());
    }

    return;
}

RealOpenMM MBPolReferenceElectrostaticsForce::calculateForceAndEnergy( const std::vector<RealVec>& particlePositions,
                                                                   const std::vector<RealOpenMM>& charges,
                                                                   const std::vector<RealOpenMM>& dipoles,
                                                                   const std::vector<RealOpenMM>& quadrupoles,
                                                                   const std::vector<RealOpenMM>& tholes,
                                                                   const std::vector<RealOpenMM>& dampingFactors,
                                                                   const std::vector<RealOpenMM>& polarity,
                                                                   const std::vector<int>& axisTypes,
                                                                   const std::vector<int>& multipoleAtomZs,
                                                                   const std::vector<int>& multipoleAtomXs,
                                                                   const std::vector<int>& multipoleAtomYs,
                                                                   const std::vector< std::vector< std::vector<int> > >& multipoleAtomCovalentInfo,
                                                                   std::vector<RealVec>& forces )
{

    // setup, including calculating induced dipoles
    // calculate electrostatic ixns including torques
    // map torques to forces
    
    std::vector<ElectrostaticsParticleData> particleData;
    setup( particlePositions, charges, dipoles, quadrupoles, tholes,
            dampingFactors, polarity, axisTypes, multipoleAtomZs, multipoleAtomXs, multipoleAtomYs,
            multipoleAtomCovalentInfo, particleData );

    std::vector<RealVec> torques;
    initializeRealVecVector( torques );
    RealOpenMM energy = calculateElectrostatic( particleData, torques, forces );

    mapTorqueToForce( particleData, multipoleAtomXs, multipoleAtomYs, multipoleAtomZs, axisTypes, torques, forces );

    return energy;
}

void MBPolReferenceElectrostaticsForce::calculateMBPolSystemElectrostaticsMoments( const std::vector<RealOpenMM>& masses,
                                                                           const std::vector<RealVec>& particlePositions,
                                                                           const std::vector<RealOpenMM>& charges,
                                                                           const std::vector<RealOpenMM>& dipoles,
                                                                           const std::vector<RealOpenMM>& quadrupoles,
                                                                           const std::vector<RealOpenMM>& tholes,
                                                                           const std::vector<RealOpenMM>& dampingFactors,
                                                                           const std::vector<RealOpenMM>& polarity,
                                                                           const std::vector<int>& axisTypes,
                                                                           const std::vector<int>& multipoleAtomZs,
                                                                           const std::vector<int>& multipoleAtomXs,
                                                                           const std::vector<int>& multipoleAtomYs,
                                                                           const std::vector< std::vector< std::vector<int> > >& multipoleAtomCovalentInfo,
                                                                           std::vector<RealOpenMM>& outputElectrostaticsMoments )
{

    // setup, including calculating induced dipoles
    // remove center of mass
    // calculate system moments

    std::vector<ElectrostaticsParticleData> particleData;
    setup( particlePositions, charges, dipoles, quadrupoles, tholes,
           dampingFactors, polarity, axisTypes, multipoleAtomZs, multipoleAtomXs, multipoleAtomYs,
           multipoleAtomCovalentInfo, particleData );

    RealOpenMM totalMass  = 0.0;
    RealVec centerOfMass  = RealVec( 0.0, 0.0, 0.0 );
    for( unsigned int ii  = 0; ii < _numParticles; ii++ ){
        RealOpenMM mass   = masses[ii];
        totalMass        += mass;
        centerOfMass     += particleData[ii].position*mass;
    }
    vector<RealVec> localPositions( _numParticles );
    if( totalMass > 0.0 ){
        centerOfMass  *= 1.0/totalMass;
    }
    for( unsigned int ii  = 0; ii < _numParticles; ii++ ){
        localPositions[ii] = particleData[ii].position - centerOfMass;
    }

    RealOpenMM netchg  = 0.0;

    RealVec dpl        = RealVec( 0.0, 0.0, 0.0 );

    RealOpenMM xxqdp   = 0.0;
    RealOpenMM xyqdp   = 0.0;
    RealOpenMM xzqdp   = 0.0;

    RealOpenMM yyqdp   = 0.0;
    RealOpenMM yzqdp   = 0.0;

    RealOpenMM zzqdp   = 0.0;

    for( unsigned int ii  = 0; ii < _numParticles; ii++ ){

        RealOpenMM charge         = particleData[ii].charge;
        RealVec position          = localPositions[ii];
        netchg                   += charge;

        RealVec netDipole         = (particleData[ii].dipole  + _inducedDipole[ii]);

        dpl                      += position*charge + netDipole;

        xxqdp                    += position[0]*position[0]*charge + 2.0*position[0]*netDipole[0];
        xyqdp                    += position[0]*position[1]*charge + position[0]*netDipole[1] + position[1]*netDipole[0];
        xzqdp                    += position[0]*position[2]*charge + position[0]*netDipole[2] + position[2]*netDipole[0];

        yyqdp                    += position[1]*position[1]*charge + 2.0*position[1]*netDipole[1];
        yzqdp                    += position[1]*position[2]*charge + position[1]*netDipole[2] + position[2]*netDipole[1];

        zzqdp                    += position[2]*position[2]*charge + 2.0*position[2]*netDipole[2];

    }

    // convert the quadrupole from traced to traceless form
 
    outputElectrostaticsMoments.resize( 13 );
    RealOpenMM qave                  = (xxqdp + yyqdp + zzqdp)/3.0;
    outputElectrostaticsMoments[4]        = 0.5*(xxqdp-qave);
    outputElectrostaticsMoments[5]        = 0.5*xyqdp;
    outputElectrostaticsMoments[6]        = 0.5*xzqdp;
    outputElectrostaticsMoments[8]        = 0.5*(yyqdp-qave);
    outputElectrostaticsMoments[9]        = 0.5*yzqdp;
    outputElectrostaticsMoments[12]       = 0.5*(zzqdp-qave);

    // add the traceless atomic quadrupoles to total quadrupole

    for( unsigned int ii = 0; ii < _numParticles; ii++ ){
        outputElectrostaticsMoments[4]  += particleData[ii].quadrupole[QXX];
        outputElectrostaticsMoments[5]  += particleData[ii].quadrupole[QXY];
        outputElectrostaticsMoments[6]  += particleData[ii].quadrupole[QXZ];
        outputElectrostaticsMoments[8]  += particleData[ii].quadrupole[QYY];
        outputElectrostaticsMoments[9]  += particleData[ii].quadrupole[QYZ];
        outputElectrostaticsMoments[12] += particleData[ii].quadrupole[QZZ];
    }
    outputElectrostaticsMoments[7]  = outputElectrostaticsMoments[5];
    outputElectrostaticsMoments[10] = outputElectrostaticsMoments[6];
    outputElectrostaticsMoments[11] = outputElectrostaticsMoments[9];
 
    RealOpenMM debye           = 4.80321;

    outputElectrostaticsMoments[0]  = netchg;

    dpl                       *= 10.0*debye;
    outputElectrostaticsMoments[1]  = dpl[0];
    outputElectrostaticsMoments[2]  = dpl[1];
    outputElectrostaticsMoments[3]  = dpl[2];
    
    debye *= 3.0;
    for( unsigned int ii = 4; ii < 13; ii++ ){
        outputElectrostaticsMoments[ii] *= 100.0*debye;
    }

    return;
}

RealOpenMM MBPolReferenceElectrostaticsForce::calculateElectrostaticPotentialForParticleGridPoint( const ElectrostaticsParticleData& particleI, const RealVec& gridPoint ) const 
{
  
    RealVec deltaR           = particleI.position - gridPoint;

    getPeriodicDelta( deltaR );

    RealOpenMM r2            = deltaR.dot( deltaR );
    RealOpenMM r             = SQRT( r2 );

    RealOpenMM rr1           = 1.0/r;
    RealOpenMM potential     = particleI.charge*rr1;

    RealOpenMM rr2           = rr1*rr1;
    RealOpenMM rr3           = rr1*rr2;

    RealOpenMM scd           = particleI.dipole.dot( deltaR );
    RealOpenMM scu           = _inducedDipole[particleI.particleIndex].dot( deltaR );
    potential               -= (scd + scu)*rr3;

    RealOpenMM rr5           = 3.0*rr3*rr2;
    RealOpenMM scq           = deltaR[0]*(particleI.quadrupole[QXX]*deltaR[0] + particleI.quadrupole[QXY]*deltaR[1] + particleI.quadrupole[QXZ]*deltaR[2]);
          scq               += deltaR[1]*(particleI.quadrupole[QXY]*deltaR[0] + particleI.quadrupole[QYY]*deltaR[1] + particleI.quadrupole[QYZ]*deltaR[2]);
          scq               += deltaR[2]*(particleI.quadrupole[QXZ]*deltaR[0] + particleI.quadrupole[QYZ]*deltaR[1] + particleI.quadrupole[QZZ]*deltaR[2]);
    potential               += scq*rr5;

    return potential;

}

void MBPolReferenceElectrostaticsForce::calculateElectrostaticPotential( const std::vector<RealVec>& particlePositions,
                                                                     const std::vector<RealOpenMM>& charges,
                                                                     const std::vector<RealOpenMM>& dipoles,
                                                                     const std::vector<RealOpenMM>& quadrupoles,
                                                                     const std::vector<RealOpenMM>& tholes,
                                                                     const std::vector<RealOpenMM>& dampingFactors,
                                                                     const std::vector<RealOpenMM>& polarity,
                                                                     const std::vector<int>& axisTypes,
                                                                     const std::vector<int>& multipoleAtomZs,
                                                                     const std::vector<int>& multipoleAtomXs,
                                                                     const std::vector<int>& multipoleAtomYs,
                                                                     const std::vector< std::vector< std::vector<int> > >& multipoleAtomCovalentInfo,
                                                                     const std::vector<RealVec>& grid,
                                                                     std::vector<RealOpenMM>& potential )
{

    // setup, including calculating induced dipoles
    // initialize potential 
    // calculate contribution of each particle to potential at grid point
    // apply prefactor

    std::vector<ElectrostaticsParticleData> particleData;
    setup( particlePositions, charges, dipoles, quadrupoles, tholes,
           dampingFactors, polarity, axisTypes, multipoleAtomZs, multipoleAtomXs, multipoleAtomYs,
           multipoleAtomCovalentInfo, particleData );

    potential.resize( grid.size() );
    for( unsigned int ii = 0; ii < grid.size(); ii++ ){
        potential[ii] = 0.0;
    }

    for( unsigned int ii = 0; ii < _numParticles; ii++ ){
        for( unsigned int jj = 0; jj < grid.size(); jj++ ){
            potential[jj] += calculateElectrostaticPotentialForParticleGridPoint( particleData[ii], grid[jj]  );
        }
    }

    RealOpenMM term = _electric/_dielectric;
    for( unsigned int ii = 0; ii < grid.size(); ii++ ){
        potential[ii] *= term;
    }

    return;
}

MBPolReferenceElectrostaticsForce::UpdateInducedDipoleFieldStruct::UpdateInducedDipoleFieldStruct( std::vector<OpenMM::RealVec>* inputFixed_E_Field, std::vector<OpenMM::RealVec>* inputInducedDipoles ) 
{ 
    fixedElectrostaticsField  = inputFixed_E_Field;
    inducedDipoles = inputInducedDipoles;
    inducedDipoleField.resize( fixedElectrostaticsField->size() );
}   

const int MBPolReferencePmeElectrostaticsForce::MBPOL_PME_ORDER = 5;

const RealOpenMM MBPolReferencePmeElectrostaticsForce::SQRT_PI = 1.77245385091;

MBPolReferencePmeElectrostaticsForce::MBPolReferencePmeElectrostaticsForce( void ) :
               MBPolReferenceElectrostaticsForce(PME),
               _cutoffDistance(0.9), _cutoffDistanceSquared(0.81),
               _pmeGridSize(0), _totalGridSize(0), _alphaEwald(0.0) 
{

    _fftplan = NULL;
    _pmeGrid = NULL;
    _pmeGridDimensions = IntVec( -1, -1, -1 );
} 

MBPolReferencePmeElectrostaticsForce::~MBPolReferencePmeElectrostaticsForce( )
{
    if( _fftplan ){
        fftpack_destroy( _fftplan );
    }
    if( _pmeGrid ){
        delete _pmeGrid;
    }
};
 
RealOpenMM MBPolReferencePmeElectrostaticsForce::getCutoffDistance( void ) const 
{
     return _cutoffDistance;
};
 
void MBPolReferencePmeElectrostaticsForce::setCutoffDistance( RealOpenMM cutoffDistance )
{
     _cutoffDistance        = cutoffDistance;
     _cutoffDistanceSquared = cutoffDistance*cutoffDistance;
};
 
RealOpenMM MBPolReferencePmeElectrostaticsForce::getAlphaEwald( void ) const 
{
     return _alphaEwald;
};

void MBPolReferencePmeElectrostaticsForce::setAlphaEwald( RealOpenMM alphaEwald )
{
     _alphaEwald = alphaEwald;
};

void MBPolReferencePmeElectrostaticsForce::getPmeGridDimensions( std::vector<int>& pmeGridDimensions ) const 
{

    pmeGridDimensions.resize( 3 );

    pmeGridDimensions[0] = _pmeGridDimensions[0];
    pmeGridDimensions[1] = _pmeGridDimensions[1];
    pmeGridDimensions[2] = _pmeGridDimensions[2];

    return;
};

void MBPolReferencePmeElectrostaticsForce::setPmeGridDimensions( std::vector<int>& pmeGridDimensions )
{

    if( (pmeGridDimensions[0] == _pmeGridDimensions[0]) && 
        (pmeGridDimensions[1] == _pmeGridDimensions[1]) &&
        (pmeGridDimensions[2] == _pmeGridDimensions[2]) )return;

    if( _fftplan ){
        fftpack_destroy(_fftplan);
    }
    fftpack_init_3d(&_fftplan,pmeGridDimensions[0], pmeGridDimensions[1], pmeGridDimensions[2]);

    _pmeGridDimensions[0] = pmeGridDimensions[0];
    _pmeGridDimensions[1] = pmeGridDimensions[1];
    _pmeGridDimensions[2] = pmeGridDimensions[2];

    initializeBSplineModuli( );
};

void MBPolReferencePmeElectrostaticsForce::setPeriodicBoxSize( RealVec& boxSize )
{

    if( boxSize[0] == 0.0 ||  boxSize[1] == 0.0 ||  boxSize[2] == 0.0 ){
        std::stringstream message;
        message << "Box size of zero is invalid.";
        throw OpenMMException(message.str());
    }

    _periodicBoxSize       = boxSize;

    _invPeriodicBoxSize[0] = 1.0/boxSize[0]; 
    _invPeriodicBoxSize[1] = 1.0/boxSize[1]; 
    _invPeriodicBoxSize[2] = 1.0/boxSize[2]; 

    return;
};

int compareInt2( const int2& v1, const int2& v2 )
{
    return v1[1] < v2[1];
}

void MBPolReferencePmeElectrostaticsForce::resizePmeArrays( void )
{

    _totalGridSize = _pmeGridDimensions[0]*_pmeGridDimensions[1]*_pmeGridDimensions[2];
    if( _pmeGridSize < _totalGridSize ){
        if( _pmeGrid ){
            delete _pmeGrid;
        }
        _pmeGrid      = new t_complex[_totalGridSize];
        _pmeGridSize  = _totalGridSize;
    }

    for( unsigned int ii = 0; ii < 3; ii++ ){
       _pmeBsplineModuli[ii].resize( _pmeGridDimensions[ii] );
       _thetai[ii].resize( MBPOL_PME_ORDER*_numParticles );
    }

    _iGrid.resize( _numParticles );
    _phi.resize( 20*_numParticles );
    _phid.resize( 10*_numParticles );
    _phip.resize( 10*_numParticles );
    _phidp.resize( 20*_numParticles );
    _pmeAtomRange.resize( _totalGridSize + 1);
    _pmeAtomGridIndex.resize( _numParticles );

    return;
}

void MBPolReferencePmeElectrostaticsForce::initializePmeGrid( void )
{
    if( _pmeGrid == NULL )return;
    //memset( _pmeGrid, 0, sizeof( t_complex )*_totalGridSize );

    for (int jj = 0; jj < _totalGridSize; jj++){
        _pmeGrid[jj].re = _pmeGrid[jj].im = 0.0;
    }
    return;
}    

void MBPolReferencePmeElectrostaticsForce::getPeriodicDelta( RealVec& deltaR ) const 
{
    deltaR[0]  -= FLOOR(deltaR[0]*_invPeriodicBoxSize[0]+0.5)*_periodicBoxSize[0];
    deltaR[1]  -= FLOOR(deltaR[1]*_invPeriodicBoxSize[1]+0.5)*_periodicBoxSize[1];
    deltaR[2]  -= FLOOR(deltaR[2]*_invPeriodicBoxSize[2]+0.5)*_periodicBoxSize[2];
}

void MBPolReferencePmeElectrostaticsForce::getPmeScale( RealVec& scale ) const
{
    scale[0] = static_cast<RealOpenMM>(_pmeGridDimensions[0])*_invPeriodicBoxSize[0];
    scale[1] = static_cast<RealOpenMM>(_pmeGridDimensions[1])*_invPeriodicBoxSize[1];
    scale[2] = static_cast<RealOpenMM>(_pmeGridDimensions[2])*_invPeriodicBoxSize[2];
}

void MBPolReferencePmeElectrostaticsForce::initializeBSplineModuli( void )
{

    // Initialize the b-spline moduli.

    int maxSize = -1;
    for( unsigned int ii = 0; ii < 3; ii++ ){
       _pmeBsplineModuli[ii].resize( _pmeGridDimensions[ii] );
        maxSize = maxSize  > _pmeGridDimensions[ii] ? maxSize : _pmeGridDimensions[ii];
    }

    RealOpenMM array[MBPOL_PME_ORDER];
    RealOpenMM x = 0.0;
    array[0]     = 1.0 - x;
    array[1]     = x;
    for( int k = 2; k < MBPOL_PME_ORDER; k++) {
        RealOpenMM denom = 1.0/k;
        array[k] = x*array[k-1]*denom;
        for (int i = 1; i < k; i++){
            array[k-i] = ((x+i)*array[k-i-1] + ((k-i+1)-x)*array[k-i])*denom;
        }
        array[0] = (1.0-x)*array[0]*denom;
    }

    vector<RealOpenMM> bsarray(maxSize+1, 0.0);
    for( int i = 2; i <= MBPOL_PME_ORDER+1; i++){
        bsarray[i] = array[i-2];
    }
    for( int dim = 0; dim < 3; dim++) {

        int size = _pmeGridDimensions[dim];

        // get the modulus of the discrete Fourier transform

        RealOpenMM factor = 2.0*M_PI/size;
        for (int i = 0; i < size; i++) {
            RealOpenMM sum1 = 0.0;
            RealOpenMM sum2 = 0.0;
            for (int j = 1; j <= size; j++) {
                RealOpenMM arg = factor*i*(j-1);
                sum1          += bsarray[j]*COS(arg);
                sum2          += bsarray[j]*SIN(arg);
            }
            _pmeBsplineModuli[dim][i] = (sum1*sum1 + sum2*sum2);
        }

        // fix for exponential Euler spline interpolation failure

        RealOpenMM eps = 1.0e-7;
        if (_pmeBsplineModuli[dim][0] < eps){
            _pmeBsplineModuli[dim][0] = 0.5*_pmeBsplineModuli[dim][1];
        }
        for (int i = 1; i < size-1; i++){
            if (_pmeBsplineModuli[dim][i] < eps){
                _pmeBsplineModuli[dim][i] = 0.5*(_pmeBsplineModuli[dim][i-1]+_pmeBsplineModuli[dim][i+1]);
            }
        }
        if (_pmeBsplineModuli[dim][size-1] < eps){
            _pmeBsplineModuli[dim][size-1] = 0.5*_pmeBsplineModuli[dim][size-2];
        }

        // compute and apply the optimal zeta coefficient

        int jcut = 50;
        for (int i = 1; i <= size; i++) {
            int k = i - 1;
            if (i > size/2)
                k = k - size;
            RealOpenMM zeta;
            if (k == 0){
                zeta = 1.0;
            } else {
                RealOpenMM sum1 = 1.0;
                RealOpenMM sum2 = 1.0;
                factor          = M_PI*k/size;
                for (int j = 1; j <= jcut; j++) {
                    RealOpenMM arg = factor/(factor+M_PI*j);
                    sum1           = sum1 + POW(arg,   MBPOL_PME_ORDER);
                    sum2           = sum2 + POW(arg, 2*MBPOL_PME_ORDER);
                }
                for (int j = 1; j <= jcut; j++) {
                    RealOpenMM arg  = factor/(factor-M_PI*j);
                    sum1           += POW(arg,   MBPOL_PME_ORDER);
                    sum2           += POW(arg, 2*MBPOL_PME_ORDER);
                }
                zeta = sum2/sum1;
            }
            _pmeBsplineModuli[dim][i-1] = _pmeBsplineModuli[dim][i-1]*(zeta*zeta);
        }
    }

    return;
}

void MBPolReferencePmeElectrostaticsForce::calculateFixedElectrostaticsFieldPairIxn( const ElectrostaticsParticleData& particleI,
                                                                            const ElectrostaticsParticleData& particleJ,
                                                                            RealOpenMM dscale, RealOpenMM pscale )
{

    // compute the real space portion of the Ewald summation

    if( particleI.particleIndex == particleJ.particleIndex )return;

    // in MBPol there is no contribution to the Fixed Multipole Field 
    // from atoms of the same water molecule. multipoleAtomZs is used for 
    // defining a reference frame for the water molecules and
    // contains the indices to the other 2 atoms in the same water molecule.

    bool isSameWater = (particleI.multipoleAtomZs == particleJ.particleIndex) or
                       (particleI.multipoleAtomYs == particleJ.particleIndex) or
                       (particleI.multipoleAtomXs == particleJ.particleIndex);

    RealVec deltaR    = particleJ.position - particleI.position;
    getPeriodicDelta( deltaR );
    RealOpenMM r2     = deltaR.dot( deltaR );

    if( r2 > _cutoffDistanceSquared )return;

    RealOpenMM r           = SQRT(r2);

    // calculate the error function damping terms

    RealOpenMM ralpha      = _alphaEwald*r;

    RealOpenMM bn0         = erfc(ralpha)/r;
    RealOpenMM alsq2       = 2.0*_alphaEwald*_alphaEwald;
    RealOpenMM alsq2n      = 1.0/(SQRT_PI*_alphaEwald);
    RealOpenMM exp2a       = EXP(-(ralpha*ralpha));
    alsq2n                *= alsq2;
    RealOpenMM bn1         = (bn0+alsq2n*exp2a)/r2;

    RealVec fim            = - deltaR * bn1 * particleJ.charge;
    RealVec fjm            = + deltaR * bn1 * particleI.charge;

// RealOpenMM rr3 = getAndScaleInverseRs( particleI, particleJ, r, false,3,TCC);
    // charge - charge
    RealOpenMM s3 = getAndScaleInverseRs( particleI, particleJ, r, true, 3,TCC);

    // FIXME verify this
    if( isSameWater ){
	s3 = 2;
    }
    RealOpenMM rr3 = (s3 - 1.)/(r2*r);

    RealVec fid            = - deltaR * rr3 * particleJ.charge;
    RealVec fjd            = + deltaR * rr3 * particleI.charge;

    RealVec fip            = - deltaR * rr3 * particleJ.charge;
    RealVec fjp            = + deltaR * rr3 * particleI.charge;

    // increment the field at each site due to this interaction

    unsigned int iIndex    = particleI.particleIndex;
    unsigned int jIndex    = particleJ.particleIndex;

    _fixedElectrostaticsField[iIndex]      += fim - fid;
    _fixedElectrostaticsField[jIndex]      += fjm - fjd;

    _fixedElectrostaticsFieldPolar[iIndex] += fim - fip;
    _fixedElectrostaticsFieldPolar[jIndex] += fjm - fjp;

    return;
}

void MBPolReferencePmeElectrostaticsForce::calculateFixedElectrostaticsField( const vector<ElectrostaticsParticleData>& particleData )
{

    // first calculate reciprocal space fixed multipole fields

    resizePmeArrays();
    computeMBPolBsplines( particleData );
    sort( _pmeAtomGridIndex.begin(), _pmeAtomGridIndex.end(), compareInt2 );
    findMBPolAtomRangeForGrid( particleData );
    initializePmeGrid();
    spreadFixedElectrostaticssOntoGrid( particleData );
    fftpack_exec_3d( _fftplan, FFTPACK_FORWARD, _pmeGrid, _pmeGrid);
    performMBPolReciprocalConvolution();
    fftpack_exec_3d( _fftplan, FFTPACK_BACKWARD, _pmeGrid, _pmeGrid);
    computeFixedPotentialFromGrid();
    recordFixedElectrostaticsField();

    // include self-energy portion of the multipole field
    // and initialize _fixedElectrostaticsFieldPolar to _fixedElectrostaticsField

    RealOpenMM term = (4.0/3.0)*(_alphaEwald*_alphaEwald*_alphaEwald)/SQRT_PI;
    for( unsigned int jj = 0; jj < _numParticles; jj++ ){
        RealVec selfEnergy             = particleData[jj].dipole*term;
        _fixedElectrostaticsField[jj]      += selfEnergy;
        _fixedElectrostaticsFieldPolar[jj]  = _fixedElectrostaticsField[jj];
    }

    // include direct space fixed multipole fields

    this->MBPolReferenceElectrostaticsForce::calculateFixedElectrostaticsField( particleData );

    return;
}

#define ARRAY(x,y) array[(x)-1+((y)-1)*MBPOL_PME_ORDER]

/**
 * This is called from computeBsplines().  It calculates the spline coefficients for a single atom along a single axis.
 */
void MBPolReferencePmeElectrostaticsForce::computeBSplinePoint( std::vector<RealOpenMM4>& thetai, RealOpenMM w  )
{

    RealOpenMM array[MBPOL_PME_ORDER*MBPOL_PME_ORDER];

    // initialization to get to 2nd order recursion

    ARRAY(2,2) = w;
    ARRAY(2,1) = 1.0 - w;

    // perform one pass to get to 3rd order recursion

    ARRAY(3,3) = 0.5 * w * ARRAY(2,2);
    ARRAY(3,2) = 0.5 * ((1.0+w)*ARRAY(2,1)+(2.0-w)*ARRAY(2,2));
    ARRAY(3,1) = 0.5 * (1.0-w) * ARRAY(2,1);

    // compute standard B-spline recursion to desired order

    for( int i = 4; i <= MBPOL_PME_ORDER; i++){
        int k = i - 1;
        RealOpenMM denom = 1.0 / k;
        ARRAY(i,i) = denom * w * ARRAY(k,k);
        for (int j = 1; j <= i-2; j++){
            ARRAY(i,i-j) = denom * ((w+j)*ARRAY(k,i-j-1)+(i-j-w)*ARRAY(k,i-j));
        }
        ARRAY(i,1) = denom * (1.0-w) * ARRAY(k,1);
    }

    // get coefficients for the B-spline first derivative

    int k = MBPOL_PME_ORDER - 1;
    ARRAY(k,MBPOL_PME_ORDER) = ARRAY(k,MBPOL_PME_ORDER-1);
    for (int i = MBPOL_PME_ORDER-1; i >= 2; i--)
        ARRAY(k,i) = ARRAY(k,i-1) - ARRAY(k,i);
    ARRAY(k,1) = -ARRAY(k,1);

    // get coefficients for the B-spline second derivative

    k = MBPOL_PME_ORDER - 2;
    ARRAY(k,MBPOL_PME_ORDER-1) = ARRAY(k,MBPOL_PME_ORDER-2);
    for (int i = MBPOL_PME_ORDER-2; i >= 2; i--)
        ARRAY(k,i) = ARRAY(k,i-1) - ARRAY(k,i);
    ARRAY(k,1) = -ARRAY(k,1);
    ARRAY(k,MBPOL_PME_ORDER) = ARRAY(k,MBPOL_PME_ORDER-1);
    for (int i = MBPOL_PME_ORDER-1; i >= 2; i--)
        ARRAY(k,i) = ARRAY(k,i-1) - ARRAY(k,i);
    ARRAY(k,1) = -ARRAY(k,1);

    // get coefficients for the B-spline third derivative

    k = MBPOL_PME_ORDER - 3;
    ARRAY(k,MBPOL_PME_ORDER-2) = ARRAY(k,MBPOL_PME_ORDER-3);
    for (int i = MBPOL_PME_ORDER-3; i >= 2; i--)
        ARRAY(k,i) = ARRAY(k,i-1) - ARRAY(k,i);
    ARRAY(k,1) = -ARRAY(k,1);
    ARRAY(k,MBPOL_PME_ORDER-1) = ARRAY(k,MBPOL_PME_ORDER-2);
    for (int i = MBPOL_PME_ORDER-2; i >= 2; i--)
        ARRAY(k,i) = ARRAY(k,i-1) - ARRAY(k,i);
    ARRAY(k,1) = -ARRAY(k,1);
    ARRAY(k,MBPOL_PME_ORDER) = ARRAY(k,MBPOL_PME_ORDER-1);
    for (int i = MBPOL_PME_ORDER-1; i >= 2; i--)
        ARRAY(k,i) = ARRAY(k,i-1) - ARRAY(k,i);
    ARRAY(k,1) = -ARRAY(k,1);

    // copy coefficients from temporary to permanent storage

    for (int i = 1; i <= MBPOL_PME_ORDER; i++){
        thetai[i-1] = RealOpenMM4(ARRAY(MBPOL_PME_ORDER,i), ARRAY(MBPOL_PME_ORDER-1,i), ARRAY(MBPOL_PME_ORDER-2,i), ARRAY(MBPOL_PME_ORDER-3,i));
    }

    return;
}

/**
 * Compute b-spline coefficients.
 */
void MBPolReferencePmeElectrostaticsForce::computeMBPolBsplines( const std::vector<ElectrostaticsParticleData>& particleData ) 
{

    //  get the B-spline coefficients for each multipole site

    for( unsigned int ii = 0; ii < _numParticles; ii++ ){
        RealVec position  = particleData[ii].position;
        getPeriodicDelta( position );
        IntVec igrid;
        for( unsigned int jj = 0; jj < 3; jj++ ){

            RealOpenMM w  = position[jj]*_invPeriodicBoxSize[jj];
            RealOpenMM fr = _pmeGridDimensions[jj]*(w-(int)(w+0.5)+0.5);
            int ifr       = static_cast<int>(fr);
            w             = fr - ifr;
            igrid[jj]     = ifr - MBPOL_PME_ORDER + 1;
            igrid[jj]    += igrid[jj] < 0 ? _pmeGridDimensions[jj] : 0;
            std::vector<RealOpenMM4> thetaiTemp(MBPOL_PME_ORDER);
            computeBSplinePoint( thetaiTemp, w);
            for( unsigned int kk = 0; kk < MBPOL_PME_ORDER; kk++ ){
                _thetai[jj][ii*MBPOL_PME_ORDER+kk] = thetaiTemp[kk];
            }
        }
    
        // Record the grid point.

        _iGrid[ii]               = igrid;
        _pmeAtomGridIndex[ii]    = int2( ii, igrid[0]*_pmeGridDimensions[1]*_pmeGridDimensions[2] + igrid[1]*_pmeGridDimensions[2] + igrid[2] );

    }

    return;
}

/**
 * For each grid point, find the range of sorted atoms associated with that point.
 */
void MBPolReferencePmeElectrostaticsForce::findMBPolAtomRangeForGrid( const vector<ElectrostaticsParticleData>& particleData ) 
{

    int last;
    int start = 0;
    last = (start == 0 ? -1 : _pmeAtomGridIndex[start-1][1]);
    for( unsigned int ii = start; ii < _numParticles; ++ii) {
        int2 atomData = _pmeAtomGridIndex[ii];
        int gridIndex = atomData[1];
        if (gridIndex != last)
        {
            for (int jj = last+1; jj <= gridIndex; ++jj){
                _pmeAtomRange[jj] = ii;
            }
            last = gridIndex;
        }
    }

    // Fill in values beyond the last atom.

    for (int j = last+1; j <= _totalGridSize; ++j){
        _pmeAtomRange[j] = _numParticles;
    }

    // The grid index won't be needed again.  Reuse that component to hold the z index, thus saving
    // some work in the multipole spreading.

    for( unsigned int ii = 0; ii < _numParticles; ii++) {
    
        RealOpenMM posz           = particleData[_pmeAtomGridIndex[ii][0]].position[2];
        posz                     -= FLOOR(posz*_invPeriodicBoxSize[2])*_periodicBoxSize[2];
        RealOpenMM w              = posz*_invPeriodicBoxSize[2];
        RealOpenMM fr             = _pmeGridDimensions[2]*(w-(int)(w+0.5)+0.5);
        int z                     = (static_cast<int>(fr)) - MBPOL_PME_ORDER + 1;
        _pmeAtomGridIndex[ii][1]  = z;
    }

    return;
}

void MBPolReferencePmeElectrostaticsForce::getGridPointGivenGridIndex( int gridIndex, IntVec& gridPoint ) const 
{

    gridPoint[0]  = gridIndex/(_pmeGridDimensions[1]*_pmeGridDimensions[2]);
    int remainder = gridIndex-gridPoint[0]*_pmeGridDimensions[1]*_pmeGridDimensions[2];
    gridPoint[1]  = remainder/_pmeGridDimensions[2];
    gridPoint[2]  = remainder-gridPoint[1]*_pmeGridDimensions[2];

    return;
}    

RealOpenMM MBPolReferencePmeElectrostaticsForce::computeFixedElectrostaticssGridValue( const vector<ElectrostaticsParticleData>& particleData,
                                                                              const int2& particleGridIndices, const RealVec& scale,
                                                                              int ix, int iy, const IntVec& gridPoint ) const 
{

    RealOpenMM gridValue = 0.0;
    for (int i = _pmeAtomRange[particleGridIndices[0]]; i < _pmeAtomRange[particleGridIndices[1]+1]; ++i) {
        int2 atomData = _pmeAtomGridIndex[i];
        int atomIndex = atomData[0];
        int z = atomData[1];
        int iz = gridPoint[2]-z+(gridPoint[2] >= z ? 0 : _pmeGridDimensions[2]);
        if( iz >= _pmeGridDimensions[2] ){
            iz -= _pmeGridDimensions[2];
        }

        RealOpenMM atomCharge       = particleData[atomIndex].charge;
        RealVec atomDipole          = RealVec( scale[0]*particleData[atomIndex].dipole[0],
                                               scale[1]*particleData[atomIndex].dipole[1],
                                               scale[2]*particleData[atomIndex].dipole[2] );

        RealOpenMM atomQuadrupoleXX =     scale[0]*scale[0]*particleData[atomIndex].quadrupole[QXX];
        RealOpenMM atomQuadrupoleXY = 2.0*scale[0]*scale[1]*particleData[atomIndex].quadrupole[QXY];
        RealOpenMM atomQuadrupoleXZ = 2.0*scale[0]*scale[2]*particleData[atomIndex].quadrupole[QXZ];
        RealOpenMM atomQuadrupoleYY =     scale[1]*scale[1]*particleData[atomIndex].quadrupole[QYY];
        RealOpenMM atomQuadrupoleYZ = 2.0*scale[1]*scale[2]*particleData[atomIndex].quadrupole[QYZ];
        RealOpenMM atomQuadrupoleZZ =     scale[2]*scale[2]*particleData[atomIndex].quadrupole[QZZ];

        RealOpenMM4 t = _thetai[0][atomIndex*MBPOL_PME_ORDER+ix];
        RealOpenMM4 u = _thetai[1][atomIndex*MBPOL_PME_ORDER+iy];
        RealOpenMM4 v = _thetai[2][atomIndex*MBPOL_PME_ORDER+iz];
        RealOpenMM term0 = atomCharge*u[0]*v[0] + atomDipole[1]*u[1]*v[0] + atomDipole[2]*u[0]*v[1] + atomQuadrupoleYY*u[2]*v[0] + atomQuadrupoleZZ*u[0]*v[2] + atomQuadrupoleYZ*u[1]*v[1];
        RealOpenMM term1 = atomDipole[0]*u[0]*v[0] + atomQuadrupoleXY*u[1]*v[0] + atomQuadrupoleXZ*u[0]*v[1];
        RealOpenMM term2 = atomQuadrupoleXX * u[0] * v[0];
        gridValue += term0*t[0] + term1*t[1] + term2*t[2];
    }
    return gridValue;
}

void MBPolReferencePmeElectrostaticsForce::spreadFixedElectrostaticssOntoGrid( const vector<ElectrostaticsParticleData>& particleData ) 
{

    RealVec scale;
    getPmeScale( scale );

    for (int gridIndex = 0; gridIndex < _totalGridSize; gridIndex++ ){

        IntVec gridPoint;
        getGridPointGivenGridIndex( gridIndex, gridPoint );

        RealOpenMM result = 0.0;
        for (int ix = 0; ix < MBPOL_PME_ORDER; ++ix)
        {
            int x = gridPoint[0]-ix+(gridPoint[0] >= ix ? 0 : _pmeGridDimensions[0]);
            for (int iy = 0; iy < MBPOL_PME_ORDER; ++iy)
            {
                int y  = gridPoint[1]-iy+(gridPoint[1] >= iy ? 0 : _pmeGridDimensions[1]);
                int z1 = gridPoint[2]-MBPOL_PME_ORDER+1;
                z1    += (z1 >= 0 ? 0 : _pmeGridDimensions[2]);
                int z2 = (z1 < gridPoint[2] ? gridPoint[2] : _pmeGridDimensions[2]-1);

                int2 particleGridIndices;
                particleGridIndices[0]  = x*_pmeGridDimensions[1]*_pmeGridDimensions[2]+y*_pmeGridDimensions[2]+z1;
                particleGridIndices[1]  = x*_pmeGridDimensions[1]*_pmeGridDimensions[2]+y*_pmeGridDimensions[2]+z2;
                result                 += computeFixedElectrostaticssGridValue( particleData, particleGridIndices, scale, ix, iy, gridPoint );

                if (z1 > gridPoint[2]){

                    particleGridIndices[0]  = x*_pmeGridDimensions[1]*_pmeGridDimensions[2]+y*_pmeGridDimensions[2];
                    particleGridIndices[1]  = x*_pmeGridDimensions[1]*_pmeGridDimensions[2]+y*_pmeGridDimensions[2]+gridPoint[2];
                    result                 += computeFixedElectrostaticssGridValue( particleData, particleGridIndices, scale, ix, iy, gridPoint );
                }
            }
        }
        _pmeGrid[gridIndex].re = result;
    }
    return;
}

void MBPolReferencePmeElectrostaticsForce::performMBPolReciprocalConvolution( void )
{

    RealOpenMM expFactor   = (M_PI*M_PI)/(_alphaEwald*_alphaEwald);
    RealOpenMM scaleFactor = 1.0/(M_PI*_periodicBoxSize[0]*_periodicBoxSize[1]*_periodicBoxSize[2]);

    for (int index = 0; index < _totalGridSize; index++)
    {
        int kx = index/(_pmeGridDimensions[1]*_pmeGridDimensions[2]);
        int remainder = index-kx*_pmeGridDimensions[1]*_pmeGridDimensions[2];
        int ky = remainder/_pmeGridDimensions[2];
        int kz = remainder-ky*_pmeGridDimensions[2];

        if (kx == 0 && ky == 0 && kz == 0){
            _pmeGrid[index].re = _pmeGrid[index].im = 0.0;
            continue;
        }

        int mx = (kx < (_pmeGridDimensions[0]+1)/2) ? kx : (kx-_pmeGridDimensions[0]);
        int my = (ky < (_pmeGridDimensions[1]+1)/2) ? ky : (ky-_pmeGridDimensions[1]);
        int mz = (kz < (_pmeGridDimensions[2]+1)/2) ? kz : (kz-_pmeGridDimensions[2]);

        RealOpenMM mhx = mx*_invPeriodicBoxSize[0];
        RealOpenMM mhy = my*_invPeriodicBoxSize[1];
        RealOpenMM mhz = mz*_invPeriodicBoxSize[2];

        RealOpenMM bx = _pmeBsplineModuli[0][kx];
        RealOpenMM by = _pmeBsplineModuli[1][ky];
        RealOpenMM bz = _pmeBsplineModuli[2][kz];

        RealOpenMM m2 = mhx*mhx+mhy*mhy+mhz*mhz;
        RealOpenMM denom = m2*bx*by*bz;
        RealOpenMM eterm = scaleFactor*EXP(-expFactor*m2)/denom;

        _pmeGrid[index].re *= eterm;
        _pmeGrid[index].im *= eterm;
    }
}

void MBPolReferencePmeElectrostaticsForce::computeFixedPotentialFromGrid( void )
{
    // extract the permanent multipole field at each site

    for (int m = 0; m < _numParticles; m++) {
        IntVec gridPoint = _iGrid[m];
        RealOpenMM tuv000 = 0.0;
        RealOpenMM tuv001 = 0.0;
        RealOpenMM tuv010 = 0.0;
        RealOpenMM tuv100 = 0.0;
        RealOpenMM tuv200 = 0.0;
        RealOpenMM tuv020 = 0.0;
        RealOpenMM tuv002 = 0.0;
        RealOpenMM tuv110 = 0.0;
        RealOpenMM tuv101 = 0.0;
        RealOpenMM tuv011 = 0.0;
        RealOpenMM tuv300 = 0.0;
        RealOpenMM tuv030 = 0.0;
        RealOpenMM tuv003 = 0.0;
        RealOpenMM tuv210 = 0.0;
        RealOpenMM tuv201 = 0.0;
        RealOpenMM tuv120 = 0.0;
        RealOpenMM tuv021 = 0.0;
        RealOpenMM tuv102 = 0.0;
        RealOpenMM tuv012 = 0.0;
        RealOpenMM tuv111 = 0.0;
        for (int iz = 0; iz < MBPOL_PME_ORDER; iz++) {
            int k = gridPoint[2]+iz-(gridPoint[2]+iz >= _pmeGridDimensions[2] ? _pmeGridDimensions[2] : 0);
            RealOpenMM4 v = _thetai[2][m*MBPOL_PME_ORDER+iz];
            RealOpenMM tu00 = 0.0;
            RealOpenMM tu10 = 0.0;
            RealOpenMM tu01 = 0.0;
            RealOpenMM tu20 = 0.0;
            RealOpenMM tu11 = 0.0;
            RealOpenMM tu02 = 0.0;
            RealOpenMM tu30 = 0.0;
            RealOpenMM tu21 = 0.0;
            RealOpenMM tu12 = 0.0;
            RealOpenMM tu03 = 0.0;
            for (int iy = 0; iy < MBPOL_PME_ORDER; iy++) {
                int j = gridPoint[1]+iy-(gridPoint[1]+iy >= _pmeGridDimensions[1] ? _pmeGridDimensions[1] : 0);
                RealOpenMM4 u = _thetai[1][m*MBPOL_PME_ORDER+iy];
                RealOpenMM4 t = RealOpenMM4(0.0, 0.0, 0.0, 0.0);
                for (int ix = 0; ix < MBPOL_PME_ORDER; ix++) {
                    int i = gridPoint[0]+ix-(gridPoint[0]+ix >= _pmeGridDimensions[0] ? _pmeGridDimensions[0] : 0);
                    int gridIndex = i*_pmeGridDimensions[1]*_pmeGridDimensions[2] + j*_pmeGridDimensions[2] + k;
                    RealOpenMM tq = _pmeGrid[gridIndex].re;
                    RealOpenMM4 tadd = _thetai[0][m*MBPOL_PME_ORDER+ix];
                    t[0] += tq*tadd[0];
                    t[1] += tq*tadd[1];
                    t[2] += tq*tadd[2];
                    t[3] += tq*tadd[3];
                }
                tu00 += t[0]*u[0];
                tu10 += t[1]*u[0];
                tu01 += t[0]*u[1];
                tu20 += t[2]*u[0];
                tu11 += t[1]*u[1];
                tu02 += t[0]*u[2];
                tu30 += t[3]*u[0];
                tu21 += t[2]*u[1];
                tu12 += t[1]*u[2];
                tu03 += t[0]*u[3];
            }
            tuv000 += tu00*v[0];
            tuv100 += tu10*v[0];
            tuv010 += tu01*v[0];
            tuv001 += tu00*v[1];
            tuv200 += tu20*v[0];
            tuv020 += tu02*v[0];
            tuv002 += tu00*v[2];
            tuv110 += tu11*v[0];
            tuv101 += tu10*v[1];
            tuv011 += tu01*v[1];
            tuv300 += tu30*v[0];
            tuv030 += tu03*v[0];
            tuv003 += tu00*v[3];
            tuv210 += tu21*v[0];
            tuv201 += tu20*v[1];
            tuv120 += tu12*v[0];
            tuv021 += tu02*v[1];
            tuv102 += tu10*v[2];
            tuv012 += tu01*v[2];
            tuv111 += tu11*v[1];
        }
        _phi[20*m] = tuv000;
        _phi[20*m+1] = tuv100;
        _phi[20*m+2] = tuv010;
        _phi[20*m+3] = tuv001;
        _phi[20*m+4] = tuv200;
        _phi[20*m+5] = tuv020;
        _phi[20*m+6] = tuv002;
        _phi[20*m+7] = tuv110;
        _phi[20*m+8] = tuv101;
        _phi[20*m+9] = tuv011;
        _phi[20*m+10] = tuv300;
        _phi[20*m+11] = tuv030;
        _phi[20*m+12] = tuv003;
        _phi[20*m+13] = tuv210;
        _phi[20*m+14] = tuv201;
        _phi[20*m+15] = tuv120;
        _phi[20*m+16] = tuv021;
        _phi[20*m+17] = tuv102;
        _phi[20*m+18] = tuv012;
        _phi[20*m+19] = tuv111;
    }
}

t_complex MBPolReferencePmeElectrostaticsForce::computeInducedDipoleGridValue( const int2& particleGridIndices, const RealVec& scale, int ix, int iy,
                                                                           const IntVec& gridPoint,
                                                                           const std::vector<RealVec>& inputInducedDipole,
                                                                           const std::vector<RealVec>& inputInducedDipolePolar ) const 
{

    
    // loop over particles contributing to this grid point

    t_complex gridValue;
    gridValue.re = gridValue.im = 0.0;
    
    for (int i = _pmeAtomRange[particleGridIndices[0]]; i < _pmeAtomRange[particleGridIndices[1]+1]; ++i){
        int2 atomData = _pmeAtomGridIndex[i];
        int atomIndex = atomData[0];
        int z = atomData[1];
        int iz = gridPoint[2]-z+(gridPoint[2] >= z ? 0 : _pmeGridDimensions[2]);
        if( iz >= _pmeGridDimensions[2] ){
            iz -= _pmeGridDimensions[2];
        }
        RealVec inducedDipole       = RealVec( scale[0]*inputInducedDipole[atomIndex][0],
                                               scale[1]*inputInducedDipole[atomIndex][1],
                                               scale[2]*inputInducedDipole[atomIndex][2] );

        RealVec inducedDipolePolar  = RealVec( scale[0]*inputInducedDipolePolar[atomIndex][0],
                                               scale[1]*inputInducedDipolePolar[atomIndex][1],
                                               scale[2]*inputInducedDipolePolar[atomIndex][2] );

        RealOpenMM4 t = _thetai[0][atomIndex*MBPOL_PME_ORDER+ix];
        RealOpenMM4 u = _thetai[1][atomIndex*MBPOL_PME_ORDER+iy];
        RealOpenMM4 v = _thetai[2][atomIndex*MBPOL_PME_ORDER+iz];

        RealOpenMM term01 = inducedDipole[1]*u[1]*v[0] + inducedDipole[2]*u[0]*v[1];
        RealOpenMM term11 = inducedDipole[0]*u[0]*v[0];
        RealOpenMM term02 = inducedDipolePolar[1]*u[1]*v[0] + inducedDipolePolar[2]*u[0]*v[1];
        RealOpenMM term12 = inducedDipolePolar[0]*u[0]*v[0];

        gridValue.re += term01*t[0] + term11*t[1];
        gridValue.im += term02*t[0] + term12*t[1];

    }
    return gridValue;
}

void MBPolReferencePmeElectrostaticsForce::spreadInducedDipolesOnGrid( const std::vector<RealVec>& inputInducedDipole,
                                                                   const std::vector<RealVec>& inputInducedDipolePolar )
{
    RealVec scale;
    getPmeScale( scale );

    for (int gridIndex = 0; gridIndex < _totalGridSize; gridIndex++ )
    {
        IntVec gridPoint;
        getGridPointGivenGridIndex( gridIndex, gridPoint );

        t_complex gridValue;
        gridValue.re  = gridValue.im = 0.0;

        for (int ix = 0; ix < MBPOL_PME_ORDER; ++ix)
        {
            int x = gridPoint[0]-ix+(gridPoint[0] >= ix ? 0 : _pmeGridDimensions[0]);
            for (int iy = 0; iy < MBPOL_PME_ORDER; ++iy)
            {
                int y   = gridPoint[1]-iy+(gridPoint[1] >= iy ? 0 : _pmeGridDimensions[1]);
                int z1  = gridPoint[2]-MBPOL_PME_ORDER+1;
                z1     += (z1 >= 0 ? 0 : _pmeGridDimensions[2]);
                int z2  = (z1 < gridPoint[2] ? gridPoint[2] : _pmeGridDimensions[2]-1);

                int2 particleGridIndices;
                particleGridIndices[0] = x*_pmeGridDimensions[1]*_pmeGridDimensions[2]+y*_pmeGridDimensions[2]+z1;
                particleGridIndices[1] = x*_pmeGridDimensions[1]*_pmeGridDimensions[2]+y*_pmeGridDimensions[2]+z2;
                gridValue             += computeInducedDipoleGridValue( particleGridIndices, scale, ix, iy, gridPoint, inputInducedDipole, inputInducedDipolePolar );

                if (z1 > gridPoint[2])
                {
                    particleGridIndices[0]  =  x*_pmeGridDimensions[1]*_pmeGridDimensions[2]+y*_pmeGridDimensions[2];
                    particleGridIndices[1]  =  x*_pmeGridDimensions[1]*_pmeGridDimensions[2]+y*_pmeGridDimensions[2]+gridPoint[2];
                    gridValue              +=  computeInducedDipoleGridValue( particleGridIndices, scale, ix, iy, gridPoint, inputInducedDipole, inputInducedDipolePolar );
                }
            }
        }
        _pmeGrid[gridIndex] = gridValue;
    }

    return;
}

void MBPolReferencePmeElectrostaticsForce::computeInducedPotentialFromGrid( void )
{
    // extract the induced dipole field at each site

    for (int m = 0; m < _numParticles; m++) {
        IntVec gridPoint = _iGrid[m];
        RealOpenMM tuv100_1 = 0.0;
        RealOpenMM tuv010_1 = 0.0;
        RealOpenMM tuv001_1 = 0.0;
        RealOpenMM tuv200_1 = 0.0;
        RealOpenMM tuv020_1 = 0.0;
        RealOpenMM tuv002_1 = 0.0;
        RealOpenMM tuv110_1 = 0.0;
        RealOpenMM tuv101_1 = 0.0;
        RealOpenMM tuv011_1 = 0.0;
        RealOpenMM tuv100_2 = 0.0;
        RealOpenMM tuv010_2 = 0.0;
        RealOpenMM tuv001_2 = 0.0;
        RealOpenMM tuv200_2 = 0.0;
        RealOpenMM tuv020_2 = 0.0;
        RealOpenMM tuv002_2 = 0.0;
        RealOpenMM tuv110_2 = 0.0;
        RealOpenMM tuv101_2 = 0.0;
        RealOpenMM tuv011_2 = 0.0;
        RealOpenMM tuv000 = 0.0;
        RealOpenMM tuv001 = 0.0;
        RealOpenMM tuv010 = 0.0;
        RealOpenMM tuv100 = 0.0;
        RealOpenMM tuv200 = 0.0;
        RealOpenMM tuv020 = 0.0;
        RealOpenMM tuv002 = 0.0;
        RealOpenMM tuv110 = 0.0;
        RealOpenMM tuv101 = 0.0;
        RealOpenMM tuv011 = 0.0;
        RealOpenMM tuv300 = 0.0;
        RealOpenMM tuv030 = 0.0;
        RealOpenMM tuv003 = 0.0;
        RealOpenMM tuv210 = 0.0;
        RealOpenMM tuv201 = 0.0;
        RealOpenMM tuv120 = 0.0;
        RealOpenMM tuv021 = 0.0;
        RealOpenMM tuv102 = 0.0;
        RealOpenMM tuv012 = 0.0;
        RealOpenMM tuv111 = 0.0;
        for (int iz = 0; iz < MBPOL_PME_ORDER; iz++) {
            int k = gridPoint[2]+iz-(gridPoint[2]+iz >= _pmeGridDimensions[2] ? _pmeGridDimensions[2] : 0);
            RealOpenMM4 v = _thetai[2][m*MBPOL_PME_ORDER+iz];
            RealOpenMM tu00_1 = 0.0;
            RealOpenMM tu01_1 = 0.0;
            RealOpenMM tu10_1 = 0.0;
            RealOpenMM tu20_1 = 0.0;
            RealOpenMM tu11_1 = 0.0;
            RealOpenMM tu02_1 = 0.0;
            RealOpenMM tu00_2 = 0.0;
            RealOpenMM tu01_2 = 0.0;
            RealOpenMM tu10_2 = 0.0;
            RealOpenMM tu20_2 = 0.0;
            RealOpenMM tu11_2 = 0.0;
            RealOpenMM tu02_2 = 0.0;
            RealOpenMM tu00 = 0.0;
            RealOpenMM tu10 = 0.0;
            RealOpenMM tu01 = 0.0;
            RealOpenMM tu20 = 0.0;
            RealOpenMM tu11 = 0.0;
            RealOpenMM tu02 = 0.0;
            RealOpenMM tu30 = 0.0;
            RealOpenMM tu21 = 0.0;
            RealOpenMM tu12 = 0.0;
            RealOpenMM tu03 = 0.0;
            for (int iy = 0; iy < MBPOL_PME_ORDER; iy++) {
                int j = gridPoint[1]+iy-(gridPoint[1]+iy >= _pmeGridDimensions[1] ? _pmeGridDimensions[1] : 0);
                RealOpenMM4 u = _thetai[1][m*MBPOL_PME_ORDER+iy];
                RealOpenMM t0_1 = 0.0;
                RealOpenMM t1_1 = 0.0;
                RealOpenMM t2_1 = 0.0;
                RealOpenMM t0_2 = 0.0;
                RealOpenMM t1_2 = 0.0;
                RealOpenMM t2_2 = 0.0;
                RealOpenMM t3 = 0.0;
                for (int ix = 0; ix < MBPOL_PME_ORDER; ix++) {
                    int i = gridPoint[0]+ix-(gridPoint[0]+ix >= _pmeGridDimensions[0] ? _pmeGridDimensions[0] : 0);
                    int gridIndex = i*_pmeGridDimensions[1]*_pmeGridDimensions[2] + j*_pmeGridDimensions[2] + k;
                    t_complex tq = _pmeGrid[gridIndex];
                    RealOpenMM4 tadd = _thetai[0][m*MBPOL_PME_ORDER+ix];
                    t0_1 += tq.re*tadd[0];
                    t1_1 += tq.re*tadd[1];
                    t2_1 += tq.re*tadd[2];
                    t0_2 += tq.im*tadd[0];
                    t1_2 += tq.im*tadd[1];
                    t2_2 += tq.im*tadd[2];
                    t3 += (tq.re+tq.im)*tadd[3];
                }
                tu00_1 += t0_1*u[0];
                tu10_1 += t1_1*u[0];
                tu01_1 += t0_1*u[1];
                tu20_1 += t2_1*u[0];
                tu11_1 += t1_1*u[1];
                tu02_1 += t0_1*u[2];
                tu00_2 += t0_2*u[0];
                tu10_2 += t1_2*u[0];
                tu01_2 += t0_2*u[1];
                tu20_2 += t2_2*u[0];
                tu11_2 += t1_2*u[1];
                tu02_2 += t0_2*u[2];
                RealOpenMM t0 = t0_1 + t0_2;
                RealOpenMM t1 = t1_1 + t1_2;
                RealOpenMM t2 = t2_1 + t2_2;
                tu00 += t0*u[0];
                tu10 += t1*u[0];
                tu01 += t0*u[1];
                tu20 += t2*u[0];
                tu11 += t1*u[1];
                tu02 += t0*u[2];
                tu30 += t3*u[0];
                tu21 += t2*u[1];
                tu12 += t1*u[2];
                tu03 += t0*u[3];
            }
            tuv100_1 += tu10_1*v[0];
            tuv010_1 += tu01_1*v[0];
            tuv001_1 += tu00_1*v[1];
            tuv200_1 += tu20_1*v[0];
            tuv020_1 += tu02_1*v[0];
            tuv002_1 += tu00_1*v[2];
            tuv110_1 += tu11_1*v[0];
            tuv101_1 += tu10_1*v[1];
            tuv011_1 += tu01_1*v[1];
            tuv100_2 += tu10_2*v[0];
            tuv010_2 += tu01_2*v[0];
            tuv001_2 += tu00_2*v[1];
            tuv200_2 += tu20_2*v[0];
            tuv020_2 += tu02_2*v[0];
            tuv002_2 += tu00_2*v[2];
            tuv110_2 += tu11_2*v[0];
            tuv101_2 += tu10_2*v[1];
            tuv011_2 += tu01_2*v[1];
            tuv000 += tu00*v[0];
            tuv100 += tu10*v[0];
            tuv010 += tu01*v[0];
            tuv001 += tu00*v[1];
            tuv200 += tu20*v[0];
            tuv020 += tu02*v[0];
            tuv002 += tu00*v[2];
            tuv110 += tu11*v[0];
            tuv101 += tu10*v[1];
            tuv011 += tu01*v[1];
            tuv300 += tu30*v[0];
            tuv030 += tu03*v[0];
            tuv003 += tu00*v[3];
            tuv210 += tu21*v[0];
            tuv201 += tu20*v[1];
            tuv120 += tu12*v[0];
            tuv021 += tu02*v[1];
            tuv102 += tu10*v[2];
            tuv012 += tu01*v[2];
            tuv111 += tu11*v[1];
        }
        _phid[10*m]   = 0.0;
        _phid[10*m+1] = tuv100_1;
        _phid[10*m+2] = tuv010_1;
        _phid[10*m+3] = tuv001_1;
        _phid[10*m+4] = tuv200_1;
        _phid[10*m+5] = tuv020_1;
        _phid[10*m+6] = tuv002_1;
        _phid[10*m+7] = tuv110_1;
        _phid[10*m+8] = tuv101_1;
        _phid[10*m+9] = tuv011_1;

        _phip[10*m]   = 0.0;
        _phip[10*m+1] = tuv100_2;
        _phip[10*m+2] = tuv010_2;
        _phip[10*m+3] = tuv001_2;
        _phip[10*m+4] = tuv200_2;
        _phip[10*m+5] = tuv020_2;
        _phip[10*m+6] = tuv002_2;
        _phip[10*m+7] = tuv110_2;
        _phip[10*m+8] = tuv101_2;
        _phip[10*m+9] = tuv011_2;

        _phidp[20*m] = tuv000;
        _phidp[20*m+1] = tuv100;
        _phidp[20*m+2] = tuv010;
        _phidp[20*m+3] = tuv001;
        _phidp[20*m+4] = tuv200;
        _phidp[20*m+5] = tuv020;
        _phidp[20*m+6] = tuv002;
        _phidp[20*m+7] = tuv110;
        _phidp[20*m+8] = tuv101;
        _phidp[20*m+9] = tuv011;
        _phidp[20*m+10] = tuv300;
        _phidp[20*m+11] = tuv030;
        _phidp[20*m+12] = tuv003;
        _phidp[20*m+13] = tuv210;
        _phidp[20*m+14] = tuv201;
        _phidp[20*m+15] = tuv120;
        _phidp[20*m+16] = tuv021;
        _phidp[20*m+17] = tuv102;
        _phidp[20*m+18] = tuv012;
        _phidp[20*m+19] = tuv111;
    }
    return;
}

RealOpenMM MBPolReferencePmeElectrostaticsForce::computeReciprocalSpaceFixedElectrostaticsForceAndEnergy( const std::vector<ElectrostaticsParticleData>& particleData,
                                                                                                 std::vector<RealVec>& forces, std::vector<RealVec>& torques ) const 
{
    RealOpenMM multipole[10];
    const int deriv1[] = {1, 4, 7, 8, 10, 15, 17, 13, 14, 19};
    const int deriv2[] = {2, 7, 5, 9, 13, 11, 18, 15, 19, 16};
    const int deriv3[] = {3, 8, 9, 6, 14, 16, 12, 19, 17, 18};
    RealVec scale;
    getPmeScale( scale );
    RealOpenMM energy = 0.0;
    for (int i = 0; i < _numParticles; i++ ) {

        // Compute the torque.

        multipole[0] = particleData[i].charge;

        multipole[1] =  particleData[i].dipole[0];
        multipole[2] =  particleData[i].dipole[1];
        multipole[3] =  particleData[i].dipole[2];

        multipole[4] = particleData[i].quadrupole[QXX];
        multipole[5] = particleData[i].quadrupole[QYY];
        multipole[6] = particleData[i].quadrupole[QZZ];

        multipole[7] = particleData[i].quadrupole[QXY]*2.0;
        multipole[8] = particleData[i].quadrupole[QXZ]*2.0;
        multipole[9] = particleData[i].quadrupole[QYZ]*2.0;

        const RealOpenMM* phi = &_phi[20*i];
        torques[i][0] += _electric*(multipole[3]*scale[1]*phi[2] - multipole[2]*scale[2]*phi[3]
                      + 2.0*(multipole[6]-multipole[5])*scale[1]*scale[2]*phi[9]
                      + multipole[8]*scale[0]*scale[1]*phi[7] + multipole[9]*scale[1]*scale[1]*phi[5]
                      - multipole[7]*scale[0]*scale[2]*phi[8] - multipole[9]*scale[2]*scale[2]*phi[6]);

        torques[i][1] += _electric*(multipole[1]*scale[2]*phi[3] - multipole[3]*scale[0]*phi[1]
                      + 2.0*(multipole[4]-multipole[6])*scale[0]*scale[2]*phi[8]
                      + multipole[7]*scale[1]*scale[2]*phi[9] + multipole[8]*scale[2]*scale[2]*phi[6]
                      - multipole[8]*scale[0]*scale[0]*phi[4] - multipole[9]*scale[0]*scale[1]*phi[7]);

        torques[i][2] += _electric*(multipole[2]*scale[0]*phi[1] - multipole[1]*scale[1]*phi[2]
                      + 2.0*(multipole[5]-multipole[4])*scale[0]*scale[1]*phi[7]
                      + multipole[7]*scale[0]*scale[0]*phi[4] + multipole[9]*scale[0]*scale[2]*phi[8]
                      - multipole[7]*scale[1]*scale[1]*phi[5] - multipole[8]*scale[1]*scale[2]*phi[9]);

        // Compute the force and energy.

        multipole[1] *= scale[0];
        multipole[2] *= scale[1];
        multipole[3] *= scale[2];
        multipole[4] *= scale[0]*scale[0];
        multipole[5] *= scale[1]*scale[1];
        multipole[6] *= scale[2]*scale[2];
        multipole[7] *= scale[0]*scale[1];
        multipole[8] *= scale[0]*scale[2];
        multipole[9] *= scale[1]*scale[2];

        RealVec f = RealVec( 0.0, 0.0, 0.0);
        for (int k = 0; k < 10; k++) {
            energy += multipole[k]*phi[k];
            f[0]   += multipole[k]*phi[deriv1[k]];
            f[1]   += multipole[k]*phi[deriv2[k]];
            f[2]   += multipole[k]*phi[deriv3[k]];
        }
        f[0]           *= scale[0];
        f[1]           *= scale[1];
        f[2]           *= scale[2];
        f              *= (_electric);
        forces[i]      -= f;

    }
    return (0.5*_electric*energy);
}

/**
 * Compute the forces due to the reciprocal space PME calculation for induced dipoles.
 */
RealOpenMM MBPolReferencePmeElectrostaticsForce::computeReciprocalSpaceInducedDipoleForceAndEnergy( MBPolReferenceElectrostaticsForce::PolarizationType polarizationType,
                                                                                                const std::vector<ElectrostaticsParticleData>& particleData,
                                                                                                std::vector<RealVec>& forces, std::vector<RealVec>& torques) const 
{

    RealOpenMM multipole[10];
    RealOpenMM inducedDipole[3];
    RealOpenMM inducedDipolePolar[3];
    RealOpenMM scales[3];
    const int deriv1[] = {1, 4, 7, 8, 10, 15, 17, 13, 14, 19};
    const int deriv2[] = {2, 7, 5, 9, 13, 11, 18, 15, 19, 16};
    const int deriv3[] = {3, 8, 9, 6, 14, 16, 12, 19, 17, 18};
    RealVec scale;
    getPmeScale( scale );
    RealOpenMM energy = 0.0;
    for (int i = 0; i < _numParticles; i++ ) {

        // Compute the torque.

        unsigned int iIndex = particleData[i].particleIndex;

        multipole[0] = particleData[i].charge;

        multipole[1] = particleData[i].dipole[0];
        multipole[2] = particleData[i].dipole[1];
        multipole[3] = particleData[i].dipole[2];

        multipole[4] = particleData[i].quadrupole[QXX];
        multipole[5] = particleData[i].quadrupole[QYY];
        multipole[6] = particleData[i].quadrupole[QZZ];
        multipole[7] = particleData[i].quadrupole[QXY]*2.0;
        multipole[8] = particleData[i].quadrupole[QXZ]*2.0;
        multipole[9] = particleData[i].quadrupole[QYZ]*2.0;

        torques[iIndex][0] += 0.5*_electric*(multipole[3]*scale[1]*_phidp[20*i+2] - multipole[2]*scale[2]*_phidp[20*i+3]
                      + 2.0*(multipole[6]-multipole[5])*scale[1]*scale[2]*_phidp[20*i+9]
                      + multipole[8]*scale[0]*scale[1]*_phidp[20*i+7] + multipole[9]*scale[1]*scale[1]*_phidp[20*i+5]
                      - multipole[7]*scale[0]*scale[2]*_phidp[20*i+8] - multipole[9]*scale[2]*scale[2]*_phidp[20*i+6]);

        torques[iIndex][1] += 0.5*_electric*(multipole[1]*scale[2]*_phidp[20*i+3] - multipole[3]*scale[0]*_phidp[20*i+1]
                      + 2.0*(multipole[4]-multipole[6])*scale[0]*scale[2]*_phidp[20*i+8]
                      + multipole[7]*scale[1]*scale[2]*_phidp[20*i+9] + multipole[8]*scale[2]*scale[2]*_phidp[20*i+6]
                      - multipole[8]*scale[0]*scale[0]*_phidp[20*i+4] - multipole[9]*scale[0]*scale[1]*_phidp[20*i+7]);

        torques[iIndex][2] += 0.5*_electric*(multipole[2]*scale[0]*_phidp[20*i+1] - multipole[1]*scale[1]*_phidp[20*i+2]
                      + 2.0*(multipole[5]-multipole[4])*scale[0]*scale[1]*_phidp[20*i+7]
                      + multipole[7]*scale[0]*scale[0]*_phidp[20*i+4] + multipole[9]*scale[0]*scale[2]*_phidp[20*i+8]
                      - multipole[7]*scale[1]*scale[1]*_phidp[20*i+5] - multipole[8]*scale[1]*scale[2]*_phidp[20*i+9]);

        // Compute the force and energy.

        multipole[1] *= scale[0];
        multipole[2] *= scale[1];
        multipole[3] *= scale[2];

        multipole[4] *= scale[0]*scale[0];
        multipole[5] *= scale[1]*scale[1];
        multipole[6] *= scale[2]*scale[2];
        multipole[7] *= scale[0]*scale[1];
        multipole[8] *= scale[0]*scale[2];
        multipole[9] *= scale[1]*scale[2];

        inducedDipole[0] = _inducedDipole[i][0];
        inducedDipole[1] = _inducedDipole[i][1];
        inducedDipole[2] = _inducedDipole[i][2];

        inducedDipolePolar[0] = _inducedDipolePolar[i][0];
        inducedDipolePolar[1] = _inducedDipolePolar[i][1];
        inducedDipolePolar[2] = _inducedDipolePolar[i][2];

        energy += scale[0]*inducedDipole[0]*_phi[20*i+1];
        energy += scale[1]*inducedDipole[1]*_phi[20*i+2];
        energy += scale[2]*inducedDipole[2]*_phi[20*i+3];

        RealVec f        = RealVec(0.0, 0.0, 0.0 );

        for (int k = 0; k < 3; k++) {

            int j1 = deriv1[k+1];
            int j2 = deriv2[k+1];
            int j3 = deriv3[k+1];

            f[0] += (inducedDipole[k]+inducedDipolePolar[k])*_phi[20*i+j1]*(scale[k]/scale[0]);
            f[1] += (inducedDipole[k]+inducedDipolePolar[k])*_phi[20*i+j2]*(scale[k]/scale[1]);
            f[2] += (inducedDipole[k]+inducedDipolePolar[k])*_phi[20*i+j3]*(scale[k]/scale[2]);
 
            if( polarizationType == MBPolReferenceElectrostaticsForce::Mutual )
            {
                f[0] += (inducedDipole[k]*_phip[10*i+j1] + inducedDipolePolar[k]*_phid[10*i+j1])*(scale[k]/scale[0]);
                f[1] += (inducedDipole[k]*_phip[10*i+j2] + inducedDipolePolar[k]*_phid[10*i+j2])*(scale[k]/scale[1]);
                f[2] += (inducedDipole[k]*_phip[10*i+j3] + inducedDipolePolar[k]*_phid[10*i+j3])*(scale[k]/scale[2]);
            }

        }

        f[0] *= scale[0];
        f[1] *= scale[1];
        f[2] *= scale[2];

        for (int k = 0; k < 10; k++) {
            f[0] += multipole[k]*_phidp[20*i+deriv1[k]];
            f[1] += multipole[k]*_phidp[20*i+deriv2[k]];
            f[2] += multipole[k]*_phidp[20*i+deriv3[k]];
        }

        f[0]           *= scale[0];
        f[1]           *= scale[1];
        f[2]           *= scale[2];
        f              *= (0.5*_electric);
        forces[iIndex] -= f;
    }
    return (0.5*_electric*energy);
}

void MBPolReferencePmeElectrostaticsForce::recordFixedElectrostaticsField( void )
{
    RealVec scale;
    getPmeScale( scale );
    scale *= -1.0;
    for (int i = 0; i < _numParticles; i++ ){
        _fixedElectrostaticsField[i][0] = scale[0]*_phi[20*i+1];
        _fixedElectrostaticsField[i][1] = scale[1]*_phi[20*i+2];
        _fixedElectrostaticsField[i][2] = scale[2]*_phi[20*i+3];
    }
    return;
}

void MBPolReferencePmeElectrostaticsForce::initializeInducedDipoles( std::vector<UpdateInducedDipoleFieldStruct>& updateInducedDipoleFields )
{

    this->MBPolReferenceElectrostaticsForce::initializeInducedDipoles( updateInducedDipoleFields );
    calculateReciprocalSpaceInducedDipoleField( updateInducedDipoleFields );
    return;
}

void MBPolReferencePmeElectrostaticsForce::recordInducedDipoleField( vector<RealVec>& field, vector<RealVec>& fieldPolar )
{
    RealVec scale;
    getPmeScale( scale );
    scale *= -1.0;
    for (int i = 0; i < _numParticles; i++ ) {

        field[i][0]      += scale[0]*_phid[10*i+1];
        field[i][1]      += scale[1]*_phid[10*i+2];
        field[i][2]      += scale[2]*_phid[10*i+3];

        fieldPolar[i][0] += scale[0]*_phip[10*i+1];
        fieldPolar[i][1] += scale[1]*_phip[10*i+2];
        fieldPolar[i][2] += scale[2]*_phip[10*i+3];
    }
    return;
}

void MBPolReferencePmeElectrostaticsForce::calculateReciprocalSpaceInducedDipoleField( std::vector<UpdateInducedDipoleFieldStruct>& updateInducedDipoleFields )
{
    // Perform PME for the induced dipoles.

    initializePmeGrid();
    spreadInducedDipolesOnGrid( *(updateInducedDipoleFields[0].inducedDipoles), *(updateInducedDipoleFields[1].inducedDipoles) );
    fftpack_exec_3d( _fftplan, FFTPACK_FORWARD, _pmeGrid, _pmeGrid);
    performMBPolReciprocalConvolution();
    fftpack_exec_3d( _fftplan, FFTPACK_BACKWARD, _pmeGrid, _pmeGrid);
    computeInducedPotentialFromGrid();
    recordInducedDipoleField( updateInducedDipoleFields[0].inducedDipoleField, updateInducedDipoleFields[1].inducedDipoleField );
}

void MBPolReferencePmeElectrostaticsForce::calculateInducedDipoleFields( const std::vector<ElectrostaticsParticleData>& particleData,
                                                                     std::vector<UpdateInducedDipoleFieldStruct>& updateInducedDipoleFields)
{

    // direct space ixns

    for( unsigned int ii = 0; ii < particleData.size(); ii++ ){
        for( unsigned int jj = ii + 1; jj < particleData.size(); jj++ ){
            calculateDirectInducedDipolePairIxns( particleData[ii], particleData[jj], updateInducedDipoleFields );
        }
    }

// FIXME segfault!   // reciprocal space ixns

    calculateReciprocalSpaceInducedDipoleField( updateInducedDipoleFields );

    // self ixn

    RealOpenMM term = (4.0/3.0)*(_alphaEwald*_alphaEwald*_alphaEwald)/SQRT_PI;
    for( unsigned int ii = 0; ii < updateInducedDipoleFields.size(); ii++ ){
        vector<RealVec>& inducedDipoles = *(updateInducedDipoleFields[ii].inducedDipoles);
        vector<RealVec>& field          = updateInducedDipoleFields[ii].inducedDipoleField;
        for( unsigned int jj = 0; jj < particleData.size(); jj++ ){
            field[jj] += inducedDipoles[jj]*term;
        }
    }

    return;
}

void MBPolReferencePmeElectrostaticsForce::calculateDirectInducedDipolePairIxn( unsigned int iIndex, unsigned int jIndex,
                                                                            RealOpenMM preFactor1, RealOpenMM preFactor2,
                                                                            const RealVec& delta,
                                                                            const std::vector<RealVec>& inducedDipole,
                                                                            std::vector<RealVec>& field ) const 
{

    // field at i due induced dipole at j

    RealOpenMM dur  = inducedDipole[jIndex].dot( delta );
    field[iIndex]  += delta*(dur*preFactor2) + inducedDipole[jIndex]*preFactor1;

    // field at j due induced dipole at i

               dur  = inducedDipole[iIndex].dot( delta );
    field[jIndex]  += delta*(dur*preFactor2) + inducedDipole[iIndex]*preFactor1;

    return;
}

void MBPolReferencePmeElectrostaticsForce::calculateDirectInducedDipolePairIxns( const ElectrostaticsParticleData& particleI,
                                                                             const ElectrostaticsParticleData& particleJ,
                                                                             std::vector<UpdateInducedDipoleFieldStruct>& updateInducedDipoleFields )
{

    // compute the real space portion of the Ewald summation
  
    RealOpenMM uscale = 1.0;
    RealVec deltaR    = particleJ.position - particleI.position;

    // periodic boundary conditions

    getPeriodicDelta( deltaR );
    RealOpenMM r2     = deltaR.dot( deltaR );

    if( r2 > _cutoffDistanceSquared )return;

    RealOpenMM r           = SQRT(r2);

    // calculate the error function damping terms

    RealOpenMM ralpha      = _alphaEwald*r;

    RealOpenMM bn0         = erfc(ralpha)/r;
    RealOpenMM alsq2       = 2.0*_alphaEwald*_alphaEwald;
    RealOpenMM alsq2n      = 1.0/(SQRT_PI*_alphaEwald);
    RealOpenMM exp2a       = EXP(-(ralpha*ralpha));
    alsq2n                *= alsq2;
    RealOpenMM bn1         = (bn0+alsq2n*exp2a)/r2;

    alsq2n                *= alsq2;
    RealOpenMM bn2         = (3.0*bn1+alsq2n*exp2a)/r2;

    // compute the error function scaled and unscaled terms

//    // TODO check if we can use get and scale
//    RealOpenMM scale3      = 1.0;
//    RealOpenMM scale5      = 1.0;
//    RealOpenMM damp        = particleI.dampingFactor*particleJ.dampingFactor;
//    if( damp != 0.0 ){
//
//        RealOpenMM ratio  = (r/damp);
//        ratio       = ratio*ratio*ratio;
//        // TODO implement variable thole in PME
//        RealOpenMM pgamma = particleI.thole[TCC] < particleJ.thole[TCC] ? particleI.thole[TCC] : particleJ.thole[TCC];
//        damp        = -pgamma*ratio;
//
//        if( damp > -50.0) {
//            RealOpenMM expdamp = expf(damp);
//            scale3        = 1.0 - expdamp;
//            scale5        = 1.0 - expdamp*(1.0-damp);
//        }
//    }

    RealOpenMM scale3 = getAndScaleInverseRs(particleI, particleJ, r, true, 3, TDD);
    RealOpenMM scale5 = getAndScaleInverseRs(particleI, particleJ, r, true, 5, TDD);

    RealOpenMM dsc3        = uscale*scale3;
    RealOpenMM dsc5        = uscale*scale5;

    RealOpenMM r3          = (r*r2);
    RealOpenMM r5          = (r3*r2);
    RealOpenMM rr3         = (1.0-dsc3)/r3;
    RealOpenMM rr5         = 3.0*(1.0-dsc5)/r5;

    RealOpenMM preFactor1  = rr3 - bn1;
    RealOpenMM preFactor2  = bn2 - rr5;

    for( unsigned int ii = 0; ii < updateInducedDipoleFields.size(); ii++ ){
        calculateDirectInducedDipolePairIxn( particleI.particleIndex, particleJ.particleIndex, preFactor1, preFactor2, deltaR,
                                            *(updateInducedDipoleFields[ii].inducedDipoles),
                                              updateInducedDipoleFields[ii].inducedDipoleField );
    }    

    return;
}

RealOpenMM MBPolReferencePmeElectrostaticsForce::calculatePmeSelfEnergy( const std::vector<ElectrostaticsParticleData>& particleData ) const 
{

    RealOpenMM cii = 0.0;
    RealOpenMM dii = 0.0;
    RealOpenMM qii = 0.0;
    for( unsigned int ii = 0; ii < _numParticles; ii++ ){

        const ElectrostaticsParticleData& particleI = particleData[ii];

        cii      +=  particleI.charge*particleI.charge;
        dii      +=  particleI.dipole.dot( particleI.dipole + _inducedDipole[ii] ) ;
    
        qii      +=  particleI.quadrupole[QXX]*particleI.quadrupole[QXX] +
                     particleI.quadrupole[QYY]*particleI.quadrupole[QYY] +
                     particleI.quadrupole[QZZ]*particleI.quadrupole[QZZ] +
                    (particleI.quadrupole[QXY]*particleI.quadrupole[QXY] +
                     particleI.quadrupole[QXZ]*particleI.quadrupole[QXZ] +
                     particleI.quadrupole[QYZ]*particleI.quadrupole[QYZ])*2.0;
    
    }

    RealOpenMM term      = 2.0*_alphaEwald*_alphaEwald;
    RealOpenMM energy    = (cii + term*(dii/3.0 + 2.0*term*qii/5.0));
               energy   *= -(_electric*_alphaEwald/(_dielectric*SQRT_PI));

    return energy;
}

void MBPolReferencePmeElectrostaticsForce::calculatePmeSelfTorque( const std::vector<ElectrostaticsParticleData>& particleData,
                                                               std::vector<RealVec>& torques ) const 
{

    RealOpenMM term = (2.0/3.0)*(_electric/_dielectric)*(_alphaEwald*_alphaEwald*_alphaEwald)/SQRT_PI;

    for( unsigned int ii = 0; ii < _numParticles; ii++ ){

        const ElectrostaticsParticleData& particleI = particleData[ii];
        RealVec ui                             = (_inducedDipole[ii] + _inducedDipolePolar[ii]);
        RealVec torque                         = particleI.dipole.cross( ui );
                torque                        *= term;
       torques[ii]                            += torque;
    }
    return;
}

RealOpenMM MBPolReferencePmeElectrostaticsForce::calculatePmeDirectElectrostaticPairIxn( const ElectrostaticsParticleData& particleI, 
                                                                                     const ElectrostaticsParticleData& particleJ,
                                                                                     std::vector<RealVec>& forces,
                                                                                     std::vector<RealVec>& torques ) const 
{

    unsigned int iIndex = particleI.particleIndex;
    unsigned int jIndex = particleJ.particleIndex;

    RealOpenMM energy;
    RealVec deltaR   = particleJ.position - particleI.position;
    getPeriodicDelta( deltaR );
    RealOpenMM r2    = deltaR.dot( deltaR );

    if( r2 > _cutoffDistanceSquared )return 0.0;

    RealOpenMM xr    = deltaR[0];
    RealOpenMM yr    = deltaR[1];
    RealOpenMM zr    = deltaR[2];

    RealOpenMM r      = SQRT(r2);
    RealOpenMM ck     = particleJ.charge;

    // set the permanent multipole and induced dipole values;

    RealOpenMM ci    = particleI.charge;

    // calculate the real space error function terms

    RealOpenMM ralpha = _alphaEwald*r;
    RealOpenMM bn0    = erfc(ralpha)/r;

    RealOpenMM alsq2  = 2.0*_alphaEwald*_alphaEwald;
    RealOpenMM alsq2n = 0.0;
    if( _alphaEwald > 0.0 ){
        alsq2n = 1.0/(SQRT_PI*_alphaEwald);
    }
    RealOpenMM exp2a  = EXP(-(ralpha*ralpha));

    alsq2n           *= alsq2;
    RealOpenMM bn1    = (bn0+alsq2n*exp2a)/r2;

    alsq2n           *= alsq2;
    RealOpenMM bn2    = (3.0*bn1+alsq2n*exp2a)/r2;

    alsq2n           *= alsq2;
    RealOpenMM bn3    = (5.0*bn2+alsq2n*exp2a)/r2;

    // apply Thole polarization damping to scale factors

    RealOpenMM rr1    = 1.0/r;
    RealOpenMM rr3    = rr1/r2;
    RealOpenMM rr5    = 3.0*rr3/r2;
    RealOpenMM rr7    = 5.0*rr5/r2;

    // calculate the scalar products for induced components

    RealOpenMM sci3  = _inducedDipole[iIndex][0]*xr 
	             + _inducedDipole[iIndex][1]*yr 
		     + _inducedDipole[iIndex][2]*zr;
    RealOpenMM sci4  = _inducedDipole[jIndex][0]*xr 
	             + _inducedDipole[jIndex][1]*yr 
		     + _inducedDipole[jIndex][2]*zr;

    RealOpenMM scip2 = _inducedDipole[iIndex][0]*_inducedDipolePolar[jIndex][0]
	             + _inducedDipole[iIndex][1]*_inducedDipolePolar[jIndex][1]
                     + _inducedDipole[iIndex][2]*_inducedDipolePolar[jIndex][2]
		     + _inducedDipolePolar[iIndex][0]*_inducedDipole[jIndex][0]
                     + _inducedDipolePolar[iIndex][1]*_inducedDipole[jIndex][1]
		     + _inducedDipolePolar[iIndex][2]*_inducedDipole[jIndex][2];

    RealOpenMM scip3 = _inducedDipolePolar[iIndex][0]*xr 
	             + _inducedDipolePolar[iIndex][1]*yr 
		     + _inducedDipolePolar[iIndex][2]*zr;
    RealOpenMM scip4 = _inducedDipolePolar[jIndex][0]*xr 
	             + _inducedDipolePolar[jIndex][1]*yr 
		     + _inducedDipolePolar[jIndex][2]*zr;

    // calculate the gl functions for permanent components

    RealOpenMM gl0           = ci*ck;

    // calculate the gl functions for induced components

    RealOpenMM gli1          = ck*sci3 - ci*sci4;

    RealOpenMM glip1         = ck*scip3 - ci*scip4;

    bool isSameWater = (particleI.multipoleAtomZs == particleJ.particleIndex) or
                       (particleI.multipoleAtomYs == particleJ.particleIndex) or
                       (particleI.multipoleAtomXs == particleJ.particleIndex);

    // in PME same water interactions are not excluded, 
    // but the scale factors are set to 0.
//    if( isSameWater ) {
////        gl0 = 0.;
////        gli1 = 0.;
////        glip1 = 0.;
//    }
    // compute the energy contributions for this interaction

    RealOpenMM e             = bn0*gl0;
    RealOpenMM ei            = 0.5 * (bn1*gli1);

    // get the real energy without any screening function

    RealOpenMM scale1CC =getAndScaleInverseRs(particleI,particleJ,r,true,1,TCC);
    RealOpenMM scale3CD =getAndScaleInverseRs(particleI,particleJ,r,true,3,TCD);

    if( isSameWater ) {
	scale1CC = scale3CD = 0.;
        //scale1CC = scale3CD = scale3DD = scale5DD = 0.;
    }
    RealOpenMM erl  =       rr1*gl0 *(1 - scale1CC) ; // charge-charge
    RealOpenMM erli = 0.5*( rr3*gli1*(1 - scale3CD)); // charge - induced dipole

    e                   = e - erl; // FIXME verify this
    ei                  = ei - erli;

    energy              = (e + ei);

#if 1

    RealOpenMM scale3CC =getAndScaleInverseRs(particleI,particleJ,r,true,3,TCC);
    RealOpenMM scale5CD =getAndScaleInverseRs(particleI,particleJ,r,true,5,TCD);
    RealOpenMM scale5DD =getAndScaleInverseRs(particleI,particleJ,r,true,5,TDD);
    RealOpenMM scale7DD =getAndScaleInverseRs(particleI,particleJ,r,true,7,TDD);

    if( isSameWater ) {
    	scale3CC = scale5CD = scale5DD = scale7DD = 0.;
    }

    // intermediate variables for permanent force terms

    RealOpenMM gf1 = bn1*gl0;

    RealOpenMM gfr1 = (1 - scale3CC) * rr3*gl0;

    // intermediate variables for induced force terms

    RealOpenMM gfi1  = 0.5*( bn2* ( gli1 
	                          + glip1
	                          + scip2 ) // inddip - inddip 
			   - bn3*(sci3*scip4+scip3*sci4));

    RealOpenMM gfi2 = -ck*bn1;
    RealOpenMM gfi3 =  ci*bn1;

    RealOpenMM gfri1 = 0.5*(rr5 * ( gli1  * (1 - scale5CD)   // charge - inddip 
			      	  + glip1 * (1 - scale5CD)   // charge - inddip 
				  + scip2 * (1 - scale5DD) ) // inddip - inddip
	    //FIXME Should there be an rr7 in front of sci3*scip4????!?
	                  - rr7 * (sci3*scip4+scip3*sci4)
			                  * (1 - scale7DD)   // inddip - inddip
			   );

    // get the permanent force with screening

    RealOpenMM ftm21 = gf1*xr;
    RealOpenMM ftm22 = gf1*yr;
    RealOpenMM ftm23 = gf1*zr;

    // get the permanent force without screening

    RealOpenMM ftm2r1 = gfr1*xr;
    RealOpenMM ftm2r2 = gfr1*yr;
    RealOpenMM ftm2r3 = gfr1*zr;

    // get the induced force with screening

    RealOpenMM ftm2i1 = gfi1*xr;
    // charge_i * inddip_j
    ftm2i1 += 0.5*( gfi2* (       _inducedDipole[iIndex][0]
			  +       _inducedDipolePolar[iIndex][0])
	          + gfi3* (       _inducedDipole[jIndex][0]
			  +       _inducedDipolePolar[jIndex][0])
    // inddipP_i* inddip_j
                  + bn2 * ( sci4 *_inducedDipolePolar[iIndex][0]
			  + scip4*_inducedDipole[iIndex][0]     ) 
	          + bn2 * ( sci3 *_inducedDipolePolar[jIndex][0]
			  + scip3*_inducedDipole[jIndex][0]     ));

    RealOpenMM ftm2i2 = gfi1*yr;
    // charge_i * inddip_j
    ftm2i2 += 0.5*( gfi2* (       _inducedDipole[iIndex][1]
			  +       _inducedDipolePolar[iIndex][1])
	          + gfi3* (       _inducedDipole[jIndex][1]
			  +       _inducedDipolePolar[jIndex][1])
    // inddipP_i* inddip_j
                  + bn2 * ( sci4 *_inducedDipolePolar[iIndex][1]
			  + scip4*_inducedDipole[iIndex][1]     ) 
	          + bn2 * ( sci3 *_inducedDipolePolar[jIndex][1]
			  + scip3*_inducedDipole[jIndex][1]     ));

    RealOpenMM ftm2i3 = gfi1*zr;
    // charge_i * inddip_j
    ftm2i3 += 0.5*( gfi2* (       _inducedDipole[iIndex][2]
			  +       _inducedDipolePolar[iIndex][2])
	          + gfi3* (       _inducedDipole[jIndex][2]
			  +       _inducedDipolePolar[jIndex][2])
    // inddipP_i* inddip_j
                  + bn2 * ( sci4 *_inducedDipolePolar[iIndex][2]
			  + scip4*_inducedDipole[iIndex][2]     ) 
	          + bn2 * ( sci3 *_inducedDipolePolar[jIndex][2]
			  + scip3*_inducedDipole[jIndex][2]     ));

    // get the induced force without screening

    RealOpenMM ftm2ri1 = gfri1*xr;
    ftm2ri1 += 0.5*rr5*(1 - scale5DD)*( sci4 * _inducedDipolePolar[iIndex][0]
				      + scip4* _inducedDipole[iIndex][0] 
			              + sci3 * _inducedDipolePolar[jIndex][0]
				      + scip3* _inducedDipole[jIndex][0]);

    // Same water atoms have no induced-dipole/charge interaction

    ftm2ri1 += 0.5*rr3*(1 - scale3CD)*( - ck *(_inducedDipole[iIndex][0]
		                              +_inducedDipolePolar[iIndex][0])
				        + ci *(_inducedDipole[jIndex][0]
				              +_inducedDipolePolar[jIndex][0]));

    RealOpenMM ftm2ri2 = gfri1*yr;
    ftm2ri2 += 0.5*rr5*(1 - scale5DD)*( sci4 * _inducedDipolePolar[iIndex][1]
				      + scip4* _inducedDipole[iIndex][1] 
			              + sci3 * _inducedDipolePolar[jIndex][1]
				      + scip3* _inducedDipole[jIndex][1]);

    ftm2ri2 += 0.5*rr3*(1 - scale3CD)*( - ck *(_inducedDipole[iIndex][1]
		                              +_inducedDipolePolar[iIndex][1])
				        + ci *(_inducedDipole[jIndex][1]
				              +_inducedDipolePolar[jIndex][1]));


    RealOpenMM ftm2ri3 = gfri1*zr;
    ftm2ri3 += 0.5*rr5*(1 - scale5DD)*( sci4 * _inducedDipolePolar[iIndex][2]
				      + scip4* _inducedDipole[iIndex][2] 
			              + sci3 * _inducedDipolePolar[jIndex][2]
				      + scip3* _inducedDipole[jIndex][2]);

    ftm2ri3 += 0.5*rr3*(1 - scale3CD)*( - ck *(_inducedDipole[iIndex][2]
		                              +_inducedDipolePolar[iIndex][2])
				        + ci *(_inducedDipole[jIndex][2]
				              +_inducedDipolePolar[jIndex][2]));

#endif


#if 0 // partially excluded + direct polarization correction
//    // account for partially excluded induced interactions
//
//    RealOpenMM temp3 = 0.5 * rr3 * ((gli1+gli6)*scalingFactors[P_SCALE]
//                               +(glip1+glip6)*scalingFactors[D_SCALE]);
//
//    RealOpenMM temp5 = 0.5 * rr5 * ((gli2+gli7)*scalingFactors[P_SCALE]
//                               +(glip2+glip7)*scalingFactors[D_SCALE]);
//
//    RealOpenMM temp7 = 0.5 * rr7 * (gli3*scalingFactors[P_SCALE]
//                               +glip3*scalingFactors[D_SCALE]);
//
//    RealOpenMM fridmp1 = temp3*ddsc31 + temp5*ddsc51 + temp7*ddsc71;
//    RealOpenMM fridmp2 = temp3*ddsc32 + temp5*ddsc52 + temp7*ddsc72;
//    RealOpenMM fridmp3 = temp3*ddsc33 + temp5*ddsc53 + temp7*ddsc73;
//
//    // find some scaling terms for induced-induced force
//
//    temp3         = 0.5 * rr3 * scalingFactors[U_SCALE] * scip2;
//    temp5         = -0.5 * rr5 * scalingFactors[U_SCALE] * (sci3*scip4+scip3*sci4);
//    RealOpenMM findmp1 = temp3*ddsc31 + temp5*ddsc51;
//    RealOpenMM findmp2 = temp3*ddsc32 + temp5*ddsc52;
//    RealOpenMM findmp3 = temp3*ddsc33 + temp5*ddsc53;

    // modify the forces for partially excluded interactions
    // FIXME check how to disable this in the xml

//    ftm2i1       -= (fridmp1 + findmp1);
//    ftm2i2       -= (fridmp2 + findmp2);
//    ftm2i3       -= (fridmp3 + findmp3);

//    // correction to convert mutual to direct polarization force
//
//    if( getPolarizationType() == MBPolReferenceElectrostaticsForce::Direct ){
//
//       RealOpenMM gfd     = 0.5 * (bn2*scip2 - bn3*(scip3*sci4+sci3*scip4));
//       ftm2i1       -= gfd*xr + 0.5*bn2*(sci4*_inducedDipolePolar[iIndex][0]+scip4*_inducedDipole[iIndex][0]+sci3*_inducedDipolePolar[jIndex][0]+scip3*_inducedDipole[jIndex][0]);
//       ftm2i2       -= gfd*yr + 0.5*bn2*(sci4*_inducedDipolePolar[iIndex][1]+scip4*_inducedDipole[iIndex][1]+sci3*_inducedDipolePolar[jIndex][1]+scip3*_inducedDipole[jIndex][1]);
//       ftm2i3       -= gfd*zr + 0.5*bn2*(sci4*_inducedDipolePolar[iIndex][2]+scip4*_inducedDipole[iIndex][2]+sci3*_inducedDipolePolar[jIndex][2]+scip3*_inducedDipole[jIndex][2]);
//
//       RealOpenMM gfdr    = 0.5 * (rr5*scip2*usc3 - rr7*(scip3*sci4 +sci3*scip4)*usc5);
//       RealOpenMM fdir1   = gfdr*xr + 0.5*usc5*rr5*(sci4*_inducedDipolePolar[iIndex][0]+scip4*_inducedDipole[iIndex][0] + sci3*_inducedDipolePolar[jIndex][0]+scip3*_inducedDipole[jIndex][0]);
//       RealOpenMM fdir2   = gfdr*yr + 0.5*usc5*rr5*(sci4*_inducedDipolePolar[iIndex][1]+scip4*_inducedDipole[iIndex][1] + sci3*_inducedDipolePolar[jIndex][1]+scip3*_inducedDipole[jIndex][1]);
//       RealOpenMM fdir3   = gfdr*zr + 0.5*usc5*rr5*(sci4*_inducedDipolePolar[iIndex][2]+scip4*_inducedDipole[iIndex][2] + sci3*_inducedDipolePolar[jIndex][2]+scip3*_inducedDipole[jIndex][2]);
//
//       ftm2i1       += fdir1 + findmp1;
//       ftm2i2       += fdir2 + findmp2;
//       ftm2i3       += fdir3 + findmp3;
//    }
#endif

#if 0 // torque
//    // intermediate variables for induced torque terms
//
//    RealOpenMM gti2  = 0.5 * bn2 * (sci4+scip4);
//    RealOpenMM gti3  = 0.5 * bn2 * (sci3+scip3);
//    RealOpenMM gti4  = gfi4;
//    RealOpenMM gti5  = gfi5;
//    RealOpenMM gti6  = gfi6;
//    RealOpenMM gtri2 = 0.5 * rr5 * (sci4*psc5+scip4*dsc5);
//    RealOpenMM gtri3 = 0.5 * rr5 * (sci3*psc5+scip3*dsc5);
//    RealOpenMM gtri4 = gfri4;
//    RealOpenMM gtri5 = gfri5;
//    RealOpenMM gtri6 = gfri6;
//
//    // get the permanent torque with screening
//
//    RealOpenMM ttm21 = -bn1*dixdk1 + gf2*dixr1
//        + gf4*(dixqkr1+dkxqir1+rxqidk1-2.0*qixqk1)
//        - gf5*rxqir1 - gf7*(rxqikr1+qkrxqir1);
//    RealOpenMM ttm22 = -bn1*dixdk2 + gf2*dixr2
//        + gf4*(dixqkr2+dkxqir2+rxqidk2-2.0*qixqk2)
//        - gf5*rxqir2 - gf7*(rxqikr2+qkrxqir2);
//    RealOpenMM ttm23 = -bn1*dixdk3 + gf2*dixr3
//        + gf4*(dixqkr3+dkxqir3+rxqidk3-2.0*qixqk3)
//        - gf5*rxqir3 - gf7*(rxqikr3+qkrxqir3);
//    RealOpenMM ttm31 = bn1*dixdk1 + gf3*dkxr1
//        - gf4*(dixqkr1+dkxqir1+rxqkdi1-2.0*qixqk1)
//        - gf6*rxqkr1 - gf7*(rxqkir1-qkrxqir1);
//    RealOpenMM ttm32 = bn1*dixdk2 + gf3*dkxr2
//        - gf4*(dixqkr2+dkxqir2+rxqkdi2-2.0*qixqk2)
//        - gf6*rxqkr2 - gf7*(rxqkir2-qkrxqir2);
//    RealOpenMM ttm33 = bn1*dixdk3 + gf3*dkxr3
//        - gf4*(dixqkr3+dkxqir3+rxqkdi3-2.0*qixqk3)
//        - gf6*rxqkr3 - gf7*(rxqkir3-qkrxqir3);
//
//    // get the permanent torque without screening
//
//    RealOpenMM ttm2r1 = -rr3*dixdk1 + gfr2*dixr1-gfr5*rxqir1
//        + gfr4*(dixqkr1+dkxqir1+rxqidk1-2.0*qixqk1)
//        - gfr7*(rxqikr1+qkrxqir1);
//    RealOpenMM ttm2r2 = -rr3*dixdk2 + gfr2*dixr2-gfr5*rxqir2
//        + gfr4*(dixqkr2+dkxqir2+rxqidk2-2.0*qixqk2)
//        - gfr7*(rxqikr2+qkrxqir2);
//    RealOpenMM ttm2r3 = -rr3*dixdk3 + gfr2*dixr3-gfr5*rxqir3
//        + gfr4*(dixqkr3+dkxqir3+rxqidk3-2.0*qixqk3)
//        - gfr7*(rxqikr3+qkrxqir3);
//    RealOpenMM ttm3r1 = rr3*dixdk1 + gfr3*dkxr1 -gfr6*rxqkr1
//        - gfr4*(dixqkr1+dkxqir1+rxqkdi1-2.0*qixqk1)
//        - gfr7*(rxqkir1-qkrxqir1);
//    RealOpenMM ttm3r2 = rr3*dixdk2 + gfr3*dkxr2 -gfr6*rxqkr2
//        - gfr4*(dixqkr2+dkxqir2+rxqkdi2-2.0*qixqk2)
//        - gfr7*(rxqkir2-qkrxqir2);
//    RealOpenMM ttm3r3 = rr3*dixdk3 + gfr3*dkxr3 -gfr6*rxqkr3
//        - gfr4*(dixqkr3+dkxqir3+rxqkdi3-2.0*qixqk3)
//        - gfr7*(rxqkir3-qkrxqir3);
//
//    // get the induced torque with screening
//
//    RealOpenMM ttm2i1 = -bn1*(dixuk1+dixukp1)*0.5
//        + gti2*dixr1 + gti4*(ukxqir1+rxqiuk1
//        + ukxqirp1+rxqiukp1)*0.5 - gti5*rxqir1;
//    RealOpenMM ttm2i2 = -bn1*(dixuk2+dixukp2)*0.5
//        + gti2*dixr2 + gti4*(ukxqir2+rxqiuk2
//        + ukxqirp2+rxqiukp2)*0.5 - gti5*rxqir2;
//    RealOpenMM ttm2i3 = -bn1*(dixuk3+dixukp3)*0.5
//        + gti2*dixr3 + gti4*(ukxqir3+rxqiuk3
//        + ukxqirp3+rxqiukp3)*0.5 - gti5*rxqir3;
//    RealOpenMM ttm3i1 = -bn1*(dkxui1+dkxuip1)*0.5
//        + gti3*dkxr1 - gti4*(uixqkr1+rxqkui1
//        + uixqkrp1+rxqkuip1)*0.5 - gti6*rxqkr1;
//    RealOpenMM ttm3i2 = -bn1*(dkxui2+dkxuip2)*0.5
//        + gti3*dkxr2 - gti4*(uixqkr2+rxqkui2
//        + uixqkrp2+rxqkuip2)*0.5 - gti6*rxqkr2;
//    RealOpenMM ttm3i3 = -bn1*(dkxui3+dkxuip3)*0.5
//        + gti3*dkxr3 - gti4*(uixqkr3+rxqkui3
//        + uixqkrp3+rxqkuip3)*0.5 - gti6*rxqkr3;
//
//    // get the induced torque without screening
//
//    RealOpenMM ttm2ri1 = -rr3*(dixuk1*psc3+dixukp1*dsc3)*0.5
//        + gtri2*dixr1 + gtri4*((ukxqir1+rxqiuk1)*psc5
//        +(ukxqirp1+rxqiukp1)*dsc5)*0.5 - gtri5*rxqir1;
//    RealOpenMM ttm2ri2 = -rr3*(dixuk2*psc3+dixukp2*dsc3)*0.5
//        + gtri2*dixr2 + gtri4*((ukxqir2+rxqiuk2)*psc5
//        +(ukxqirp2+rxqiukp2)*dsc5)*0.5 - gtri5*rxqir2;
//    RealOpenMM ttm2ri3 = -rr3*(dixuk3*psc3+dixukp3*dsc3)*0.5
//        + gtri2*dixr3 + gtri4*((ukxqir3+rxqiuk3)*psc5
//        +(ukxqirp3+rxqiukp3)*dsc5)*0.5 - gtri5*rxqir3;
//    RealOpenMM ttm3ri1 = -rr3*(dkxui1*psc3+dkxuip1*dsc3)*0.5
//        + gtri3*dkxr1 - gtri4*((uixqkr1+rxqkui1)*psc5
//        +(uixqkrp1+rxqkuip1)*dsc5)*0.5 - gtri6*rxqkr1;
//    RealOpenMM ttm3ri2 = -rr3*(dkxui2*psc3+dkxuip2*dsc3)*0.5
//        + gtri3*dkxr2 - gtri4*((uixqkr2+rxqkui2)*psc5
//        +(uixqkrp2+rxqkuip2)*dsc5)*0.5 - gtri6*rxqkr2;
//    RealOpenMM ttm3ri3 = -rr3*(dkxui3*psc3+dkxuip3*dsc3)*0.5
//        + gtri3*dkxr3 - gtri4*((uixqkr3+rxqkui3)*psc5
//        +(uixqkrp3+rxqkuip3)*dsc5)*0.5 - gtri6*rxqkr3;
//
#endif
    // handle the case where scaling is used

#if 1
    // it was (1.0 - -scalingFactors[M_SCALE]) in each term
    ftm21  = (ftm21-(1.0)*ftm2r1);
    ftm2i1 = (ftm2i1-ftm2ri1);
//    ttm21  = (ttm21-(1.0)*ttm2r1);
//    ttm2i1 = (ttm2i1-ttm2ri1);
//    ttm31  = (ttm31-(1.0)*ttm3r1);
//    ttm3i1 = (ttm3i1-ttm3ri1);

    ftm22  = (ftm22-(1.0)*ftm2r2);
    ftm2i2 = (ftm2i2-ftm2ri2);
//    ttm22  = (ttm22-(1.0)*ttm2r2);
//    ttm2i2 = (ttm2i2-ttm2ri2);
//    ttm32  = (ttm32-(1.0)*ttm3r2);
//    ttm3i2 = (ttm3i2-ttm3ri2);

    ftm23  = (ftm23-(1.0)*ftm2r3);
    ftm2i3 = (ftm2i3-ftm2ri3);
//    ttm23  = (ttm23-(1.0)*ttm2r3);
//    ttm2i3 = (ttm2i3-ttm2ri3);
//    ttm33  = (ttm33-(1.0)*ttm3r3);
//    ttm3i3 = (ttm3i3-ttm3ri3);

    // increment gradient due to force and torque on first site;
#endif
    RealOpenMM conversionFactor  = (_electric/_dielectric);

    energy                 *= conversionFactor;

#if 1
    forces[iIndex][0]      -= (ftm21 + ftm2i1)*conversionFactor;
    forces[iIndex][1]      -= (ftm22 + ftm2i2)*conversionFactor;
    forces[iIndex][2]      -= (ftm23 + ftm2i3)*conversionFactor;

    forces[jIndex][0]      += (ftm21 + ftm2i1)*conversionFactor;
    forces[jIndex][1]      += (ftm22 + ftm2i2)*conversionFactor;
    forces[jIndex][2]      += (ftm23 + ftm2i3)*conversionFactor;
#endif
//
//    torques[iIndex][0]     += (ttm21 + ttm2i1)*conversionFactor;
//    torques[iIndex][1]     += (ttm22 + ttm2i2)*conversionFactor;
//    torques[iIndex][2]     += (ttm23 + ttm2i3)*conversionFactor;
//
//    torques[jIndex][0]     += (ttm31 + ttm3i1)*conversionFactor;
//    torques[jIndex][1]     += (ttm32 + ttm3i2)*conversionFactor;
//    torques[jIndex][2]     += (ttm33 + ttm3i3)*conversionFactor;

    return energy;

}

RealOpenMM MBPolReferencePmeElectrostaticsForce::calculateElectrostatic( const std::vector<ElectrostaticsParticleData>& particleData,
                                                                     std::vector<RealVec>& torques, std::vector<RealVec>& forces )
{

    RealOpenMM energy = 0.0;
    std::vector<RealOpenMM> scaleFactors(LAST_SCALE_TYPE_INDEX);
    for( unsigned int kk = 0; kk < scaleFactors.size(); kk++ ){
        scaleFactors[kk] = 1.0;
    }   

    // loop over particle pairs for direct space interactions

    for( unsigned int ii = 0; ii < particleData.size(); ii++ ){
        for( unsigned int jj = ii+1; jj < particleData.size(); jj++ ){

            if( jj <= _maxScaleIndex[ii] ){
                getElectrostaticsScaleFactors( ii, jj, scaleFactors);
            }

            energy += calculatePmeDirectElectrostaticPairIxn( particleData[ii], particleData[jj], forces, torques );

            if( jj <= _maxScaleIndex[ii] ){
                for( unsigned int kk = 0; kk < LAST_SCALE_TYPE_INDEX; kk++ ){
                    scaleFactors[kk] = 1.0;
                }
            }
        }
    }

    calculatePmeSelfTorque( particleData, torques );
    energy += computeReciprocalSpaceInducedDipoleForceAndEnergy( getPolarizationType(), particleData, forces, torques );
    energy += computeReciprocalSpaceFixedElectrostaticsForceAndEnergy( particleData, forces, torques );
    energy += calculatePmeSelfEnergy( particleData );

    return energy;
}

void MBPolReferenceElectrostaticsForce::computeWaterCharge(
        ElectrostaticsParticleData& particleO, ElectrostaticsParticleData& particleH1,
        ElectrostaticsParticleData& particleH2,ElectrostaticsParticleData& particleM)
{
    const double Bohr_A = 0.52917721092; // CODATA 2010
    // M-site positioning (TTM2.1-F)
    const double gammaM = 0.426706882;

    const double gamma1 = 1.0 - gammaM;
    const double gamma2 = gammaM/2;
    const double ath0 = 1.82400520401572996557;
    const double costhe = -0.24780227221366464506;
    const double reoh = 0.958649;
    const double b1D = 1.0;
    const double a = 0.2999e0;
    const double b = -0.6932e0;
    const double c0 = 1.0099e0;
    const double c1 = -0.1801e0;
    const double c2 = 0.0892e0;


    const double e =  1.602176565e-19; // C CODATA 2010

    // interaction energy of 2 unit charges 1A apart
    const double E_cc = 1.0e-7*(c0*e*c0*e)/1.0e-10; // in J
    const double Na = 6.02214129e+23; // CODATA 2010
    const double kcal_J = 4184.0;

    const double CHARGECON = sqrt(E_cc*Na/kcal_J);

    const size_t idxD0[84] = {
           1, 1, 1, 2, 1, 1, 1, 2, 2, 3, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4,
           1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 4, 4, 5, 1, 1, 1, 1, 1,
           1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 5, 5, 6, 1, 1, 1, 1,
           1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5,
           5, 6, 6, 7
    };

    const size_t idxD1[84] = {
           1, 1, 2, 1, 1, 2, 3, 1, 2, 1, 1, 2, 3, 4, 1, 2, 3, 1, 2, 1,
           1, 2, 3, 4, 5, 1, 2, 3, 4, 1, 2, 3, 1, 2, 1, 1, 2, 3, 4, 5,
           6, 1, 2, 3, 4, 5, 1, 2, 3, 4, 1, 2, 3, 1, 2, 1, 1, 2, 3, 4,
           5, 6, 7, 1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 1, 2, 3, 4, 1, 2,
           3, 1, 2, 1
    };

    const size_t idxD2[84] = {
           1, 2, 1, 1, 3, 2, 1, 2, 1, 1, 4, 3, 2, 1, 3, 2, 1, 2, 1, 1,
           5, 4, 3, 2, 1, 4, 3, 2, 1, 3, 2, 1, 2, 1, 1, 6, 5, 4, 3, 2,
           1, 5, 4, 3, 2, 1, 4, 3, 2, 1, 3, 2, 1, 2, 1, 1, 7, 6, 5, 4,
           3, 2, 1, 6, 5, 4, 3, 2, 1, 5, 4, 3, 2, 1, 4, 3, 2, 1, 3, 2,
           1, 2, 1, 1
    };


    const double coefD[84] = {
          -2.1689686086730e-03, 1.4910379754728e-02, 5.3546078430060e-02,
          -7.4055995388666e-02,-3.7764333017616e-03, 1.4089887256484e-01,
          -6.2584207687264e-02,-1.1260393113022e-01,-5.7824159269319e-02,
           1.4360743650655e-02,-1.5469680141070e-02,-1.3036350092795e-02,
           2.7515837781556e-02, 1.4098478875076e-01,-2.7663168397781e-02,
          -5.2378176254797e-03,-1.0237198381792e-02, 8.9571999265473e-02,
           7.2920263098603e-03,-2.6873260551686e-01, 2.0220870325864e-02,
          -7.0764766270927e-02, 1.2140640273760e-01, 2.0978491966341e-02,
          -1.9443840512668e-01, 4.0826835370618e-02,-4.5365190474650e-02,
           6.2779900072132e-02,-1.3194351021000e-01,-1.4673032718563e-01,
           1.1894031277247e-01,-6.4952851564679e-03, 8.8503610374493e-02,
           1.4899437409291e-01, 1.3962841511565e-01,-2.6459446720450e-02,
          -5.0128914532773e-02, 1.8329676428116e-01,-1.5559089125095e-01,
          -4.0176879767592e-02, 3.6192059996636e-01, 1.0202887240343e-01,
           1.9318668580051e-01,-4.3435977107932e-01,-4.2080828803311e-02,
           1.9144626027273e-01,-1.7851138969948e-01, 1.0524533875070e-01,
          -1.7954071602185e-02, 5.2022455612120e-02,-2.8891891146828e-01,
          -4.7452036576319e-02,-1.0939400546289e-01, 3.5916564473568e-01,
          -2.0162789820172e-01,-3.5838629543696e-01, 5.6706523551202e-03,
           1.3849337488211e-01,-4.1733982195604e-01, 4.1641570764241e-01,
          -1.2243429796296e-01, 4.7141730971228e-02,-1.8224510249551e-01,
          -1.8880981556620e-01,-3.1992359561800e-01,-1.8567550546587e-01,
           6.1850530431280e-01,-6.1142756235141e-02,-1.6996135584933e-01,
           5.4252879499871e-01, 6.6128603899427e-01, 1.2107016404639e-02,
          -1.9633639729189e-01, 2.7652059420824e-03,-2.2684111109778e-01,
          -4.7924491598635e-01, 2.4287790137314e-01,-1.4296023329441e-01,
           8.9664665907006e-02,-1.4003228575602e-01,-1.3321543452254e-01,
          -1.8340983193745e-01, 2.3426707273520e-01, 1.5141050914514e-01
    };

    double ROH1[3], ROH2[3], RHH[3], dROH1(0), dROH2(0), dRHH(0);

    for (size_t i = 0; i < 3; ++i) {
        ROH1[i] = particleH1.position[i]*10. - particleO.position[i]*10.; // H1 - O
        ROH2[i] = particleH2.position[i]*10. - particleO.position[i]*10.; // H2 - O
        RHH[i] = particleH1.position[i]*10. - particleH2.position[i]*10.; // H1 - H2

        dROH1 += ROH1[i]*ROH1[i];
        dROH2 += ROH2[i]*ROH2[i];
        dRHH += RHH[i]*RHH[i];
    }

    dROH1 = std::sqrt(dROH1);
    dROH2 = std::sqrt(dROH2);
    dRHH = std::sqrt(dRHH);

    const double costh =
        (ROH1[0]*ROH2[0] + ROH1[1]*ROH2[1] + ROH1[2]*ROH2[2])/(dROH1*dROH2);

    const double efac = exp(-b1D*(std::pow((dROH1 - reoh), 2)
                                     + std::pow((dROH2 - reoh), 2)));

    const double x1 = (dROH1 - reoh)/reoh;
    const double x2 = (dROH2 - reoh)/reoh;
    const double x3 = costh - costhe;

    double fmat[3][16];

    for (size_t i = 0; i < 3; ++i) {
        fmat[i][0] = 0.0;
        fmat[i][1] = 1.0;
    }

    for (size_t j = 2; j < 16; ++j) {
        fmat[0][j] = fmat[0][j - 1]*x1;
        fmat[1][j] = fmat[1][j - 1]*x2;
        fmat[2][j] = fmat[2][j - 1]*x3;
    }

    // Calculate the dipole moment

    double p1(0), p2(0);
    double pl1 = costh;
    double pl2 = 0.5*(3*pl1*pl1 - 1.0);

    double dp1dr1(0);
    double dp1dr2(0);
    double dp1dcabc(0);
    double dp2dr1(0);
    double dp2dr2(0);
    double dp2dcabc(0);

    for (size_t j = 1; j < 84; ++j) {
        const size_t inI = idxD0[j];
        const size_t inJ = idxD1[j];
        const size_t inK = idxD2[j];

        p1 += coefD[j]*fmat[0][inI]*fmat[1][inJ]*fmat[2][inK];
        p2 += coefD[j]*fmat[0][inJ]*fmat[1][inI]*fmat[2][inK];

        dp1dr1 +=
            coefD[j]*(inI - 1)*fmat[0][inI - 1]*fmat[1][inJ]*fmat[2][inK];
        dp1dr2 +=
            coefD[j]*(inJ - 1)*fmat[0][inI]*fmat[1][inJ - 1]*fmat[2][inK];
        dp1dcabc +=
            coefD[j]*(inK - 1)*fmat[0][inI]*fmat[1][inJ]*fmat[2][inK - 1];
        dp2dr1 +=
            coefD[j]*(inJ - 1)*fmat[0][inJ - 1]*fmat[1][inI]*fmat[2][inK];
        dp2dr2 +=
            coefD[j]*(inI - 1)*fmat[0][inJ]*fmat[1][inI - 1]*fmat[2][inK];
        dp2dcabc +=
            coefD[j]*(inK - 1)*fmat[0][inJ]*fmat[1][inI]*fmat[2][inK - 1];
    }

    const double xx = Bohr_A;
    const double xx2 = xx*xx;

    dp1dr1 /= reoh/xx;
    dp1dr2 /= reoh/xx;
    dp2dr1 /= reoh/xx;
    dp2dr2 /= reoh/xx;

    const double pc0 =
        a*(std::pow(dROH1, b) + std::pow(dROH2, b))*(c0 + pl1*c1 + pl2*c2);

    const double dpc0dr1 =
        a*b*std::pow(dROH1, b - 1)*(c0 + pl1*c1 + pl2*c2)*xx2;
    const double dpc0dr2 =
        a*b*std::pow(dROH2, b - 1)*(c0 + pl1*c1 + pl2*c2)*xx2;
    double dpc0dcabc =
        a*(std::pow(dROH1, b) + std::pow(dROH2, b))*(c1 + 0.5*(6.0*pl1)*c2)*xx;

    const double defacdr1 = -2.0*b1D*(dROH1 - reoh)*efac*xx;
    const double defacdr2 = -2.0*b1D*(dROH2 - reoh)*efac*xx;

    dp1dr1 = dp1dr1*efac + p1*defacdr1 + dpc0dr1;
    dp1dr2 = dp1dr2*efac + p1*defacdr2 + dpc0dr2;
    dp1dcabc = dp1dcabc*efac + dpc0dcabc;
    dp2dr1 = dp2dr1*efac + p2*defacdr1 + dpc0dr1;
    dp2dr2 = dp2dr2*efac + p2*defacdr2 + dpc0dr2;
    dp2dcabc = dp2dcabc*efac + dpc0dcabc;

    p1 = coefD[0] + p1*efac + pc0*xx; // q^H1 in TTM2-F
    p2 = coefD[0] + p2*efac + pc0*xx; // q^H2 paper

    double chargeO= -(p1 + p2);  // Oxygen
    double chargeH1 = p1; // Hydrogen-1
    double chargeH2 = p2;  // Hydrogen-2
    double gamma2div1 = gamma2/gamma1;

    particleO.charge = 0.;
    particleH1.charge = chargeH1 + gamma2div1*(chargeH1 + chargeH2);
    particleH2.charge = chargeH2 + gamma2div1*(chargeH1 + chargeH2);
    particleM.charge = chargeO/gamma1;

    dp1dr1 /= xx;
    dp1dr2 /= xx;
    dp2dr1 /= xx;
    dp2dr2 /= xx;

    const double f1q1r13 = (dp1dr1 - (dp1dcabc*costh/dROH1))/dROH1;
    const double f1q1r23 = dp1dcabc/(dROH1*dROH2);
    const double f2q1r23 = (dp1dr2 - (dp1dcabc*costh/dROH2))/dROH2;
    const double f2q1r13 = dp1dcabc/(dROH2*dROH1);
    const double f1q2r13 = (dp2dr1 - (dp2dcabc*costh/dROH1))/dROH1;
    const double f1q2r23 = dp2dcabc/(dROH1*dROH2);
    const double f2q2r23 = (dp2dr2 - (dp2dcabc*costh/dROH2))/dROH2;
    const double f2q2r13 = dp2dcabc/(dROH2*dROH1);

    // first index is atom w.r.t. to which the derivative is
    // second index is the charge being differentiated

    enum ChargeDerivativesIndices { vsH1, vsH2, vsO };

    std:vector<RealVec> chargeDerivativesH1;
    chargeDerivativesH1.resize(3);

    //gradient of charge h1(second index) wrt displacement of h1(first index)
    for (size_t i = 0; i < 3; ++i) {
        chargeDerivativesH1[vsH1][i] = f1q1r13*ROH1[i] + f1q1r23*ROH2[i];
        chargeDerivativesH1[vsH2][i] = f2q1r13*ROH1[i] + f2q1r23*ROH2[i];
        chargeDerivativesH1[vsO][i] = -(chargeDerivativesH1[vsH1][i]+chargeDerivativesH1[vsH2][i]);
    }

    std::vector<RealVec> chargeDerivativesH2;
    chargeDerivativesH2.resize(3);

        //gradient of charge h1(second index) wrt displacement of h1(first index)
    for (size_t i = 0; i < 3; ++i) {
            chargeDerivativesH2[vsH1][i] = f1q2r13*ROH1[i] + f1q2r23*ROH2[i];
            chargeDerivativesH2[vsH2][i] = f2q2r13*ROH1[i] + f2q2r23*ROH2[i];
            chargeDerivativesH2[vsO][i] = -(chargeDerivativesH2[vsH1][i]+chargeDerivativesH2[vsH2][i]);
    }

    std::vector<RealVec> chargeDerivativesO;
    chargeDerivativesO.resize(3);

        //gradient of charge h1(second index) wrt displacement of h1(first index)
    for (size_t i = 0; i < 3; ++i) {
            chargeDerivativesO[vsH1][i] = -(chargeDerivativesH1[vsH1][i]+ chargeDerivativesH2[vsH1][i]);
            chargeDerivativesO[vsH2][i] =  -(chargeDerivativesH1[vsH2][i]+ chargeDerivativesH2[vsH2][i]);
            chargeDerivativesO[vsO][i] =  -(chargeDerivativesH1[vsO][i]+ chargeDerivativesH2[vsO][i]);
    }

    double sumH1, sumH2, sumO;

    for (size_t i = 0; i < 3; ++i) {
        particleM.chargeDerivatives[vsH1f][i] = 0.;
        particleM.chargeDerivatives[vsH2f][i] = 0.;
        particleM.chargeDerivatives[vsMf][i] = 0.;

        sumH1 = gamma2div1*(chargeDerivativesH1[vsH1][i]+chargeDerivativesH2[vsH1][i]);
        sumH2 = gamma2div1*(chargeDerivativesH1[vsH2][i]+chargeDerivativesH2[vsH2][i]);
        sumO = gamma2div1*(chargeDerivativesH1[vsO][i]+chargeDerivativesH2[vsO][i]);

        particleH1.chargeDerivatives[vsH1f][i] = chargeDerivativesH1[vsH1][i] + sumH1;
        particleH2.chargeDerivatives[vsH1f][i] = chargeDerivativesH1[vsH2][i] + sumH2;
        particleO.chargeDerivatives[vsH1f][i]  = chargeDerivativesH1[vsO][i] + sumO;

        particleH1.chargeDerivatives[vsH2f][i] = chargeDerivativesH2[vsH1][i] + sumH1;
        particleH2.chargeDerivatives[vsH2f][i] = chargeDerivativesH2[vsH2][i] + sumH2;
        particleO.chargeDerivatives[vsH2f][i]  = chargeDerivativesH2[vsO][i] +  sumO;

        particleH1.chargeDerivatives[vsMf][i] = chargeDerivativesO[vsH1][i] - 2*sumH1;
        particleH2.chargeDerivatives[vsMf][i] = chargeDerivativesO[vsH2][i] - 2*sumH2;
        particleO.chargeDerivatives[vsMf][i]  = chargeDerivativesO[vsO][i]  - 2*sumO;

    }

    // convert from q/A to q/nm
    for (unsigned int i = 0; i < 3; ++i) {
        for (unsigned int s = 0; s < 3; ++s) {
            particleH1.chargeDerivatives[s][i] *= 10;
            particleH2.chargeDerivatives[s][i] *= 10;
            particleO.chargeDerivatives[s][i] *= 10;
        }

    }

    // TODO implement as list

    particleH1.otherSiteIndex[vsH1f] = particleH1.particleIndex;
    particleH1.otherSiteIndex[vsH2f] = particleH2.particleIndex;
    particleH1.otherSiteIndex[vsMf]  = particleM.particleIndex;

    particleH2.otherSiteIndex[vsH1f] = particleH1.particleIndex;
    particleH2.otherSiteIndex[vsH2f] = particleH2.particleIndex;
    particleH2.otherSiteIndex[vsMf]  = particleM.particleIndex;

    particleM.otherSiteIndex[vsH1f] = particleH1.particleIndex;
    particleM.otherSiteIndex[vsH2f] = particleH2.particleIndex;
    particleM.otherSiteIndex[vsMf]  = particleM.particleIndex;

    particleO.otherSiteIndex[vsH1f] = particleH1.particleIndex;
    particleO.otherSiteIndex[vsH2f] = particleH2.particleIndex;
    particleO.otherSiteIndex[vsMf]  = particleM.particleIndex;

}
