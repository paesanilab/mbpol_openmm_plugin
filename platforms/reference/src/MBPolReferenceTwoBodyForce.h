
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

#ifndef __MBPolReferenceTwoBodyForce_H__
#define __MBPolReferenceTwoBodyForce_H__

#include "openmm/reference/RealVec.h"
#include "openmm/Vec3.h"
#include "openmm/reference/ReferenceNeighborList.h"
#include <string>
#include <vector>

using namespace  OpenMM;
enum TrimerAtomIndex { Oa, Ha1, Ha2, Ob, Hb1, Hb2, Oc, Hc1, Hc2 };
enum ClIndex { Oa_dummy, Ha1_dummy, Ha2_dummy, Cl, dummy1, dummy2, Cl2, dummy3, dummy4, Cl3 };
enum ExtraPointIndex { Xa1, Xa2, Xb1, Xb2, Xc1, Xc2 };

class MBPolReferenceTwoBodyForce;
typedef  RealOpenMM (MBPolReferenceTwoBodyForce::*CombiningFunction)( RealOpenMM x, RealOpenMM y) const;

void imageParticles(const RealVec& box, const RealVec & referenceParticle, RealVec& particleToImage);

void imageMolecules(const RealVec& box, std::vector<RealVec>& allPositions);

// ---------------------------------------------------------------------------------------

class MBPolReferenceTwoBodyForce {

public:

    /** 
     * This is an enumeration of the different methods that may be used for handling long range TwoBody forces.
     */
    enum NonbondedMethod {

        /**
         * No cutoff is applied to the interactions.  The full set of N^2 interactions is computed exactly.
         * This necessarily means that periodic boundary conditions cannot be used.  This is the default.
         */

        NoCutoff = 0,

        /**
         * Interactions beyond the cutoff distance are ignored.  
         */
        CutoffNonPeriodic = 1,
        /**
         * Periodic boundary conditions are used, so that each particle interacts only with the nearest periodic copy of
         * each other particle.  Interactions beyond the cutoff distance are ignored.  
         */
        CutoffPeriodic = 2,
    };
 
    /**---------------------------------------------------------------------------------------
       
       Constructor
       
       --------------------------------------------------------------------------------------- */
 
    MBPolReferenceTwoBodyForce( void );
 
    /**---------------------------------------------------------------------------------------
       
       Destructor
       
       --------------------------------------------------------------------------------------- */
 
    ~MBPolReferenceTwoBodyForce( ){};
 
    /**---------------------------------------------------------------------------------------
    
       Get nonbonded method
    
       @return nonbonded method
    
       --------------------------------------------------------------------------------------- */
    
    NonbondedMethod getNonbondedMethod( void ) const;

    /**---------------------------------------------------------------------------------------
    
       Set nonbonded method
    
       @param nonbonded method
    
       --------------------------------------------------------------------------------------- */
    
    void setNonbondedMethod( NonbondedMethod nonbondedMethod );

    /**---------------------------------------------------------------------------------------
    
       Get cutoff
    
       @return cutoff
    
       --------------------------------------------------------------------------------------- */
    
    double getCutoff( void ) const;

    /**---------------------------------------------------------------------------------------
    
       Set cutof
    
       @param cutoff
    
       --------------------------------------------------------------------------------------- */
    
    void setCutoff( double cutoff );


    void setPeriodicBox( const RealVec& box );

    /**---------------------------------------------------------------------------------------
    
       Get box dimensions
    
       @return box dimensions
    
       --------------------------------------------------------------------------------------- */
    
    RealVec getPeriodicBox( void ) const;

    /**---------------------------------------------------------------------------------------
    
       Calculate TwoBody ixn using neighbor list
    
       @param numParticles            number of particles
       @param particlePositions       Cartesian coordinates of particles
       @param indexIVs                position index for associated reducing particle
       @param sigmas                  particle sigmas 
       @param epsilons                particle epsilons
       @param reductions              particle reduction factors
       @param neighborList            neighbor list
       @param forces                  add forces to this vector
    
       @return energy
    
       --------------------------------------------------------------------------------------- */
    
    RealOpenMM calculateForceAndEnergy( int numParticles, const std::vector<OpenMM::RealVec>& particlePositions, 
                                        const std::vector<std::vector<int> >& allParticleIndices,
                                        const NeighborList& neighborList,
                                        std::vector<OpenMM::RealVec>& forces ) const;
         
private:

    NonbondedMethod _nonbondedMethod;
    double _cutoff;

    RealVec _periodicBoxDimensions;

    /**---------------------------------------------------------------------------------------

       Calculate pair ixn

       @param  combindedSigma       combined sigmas
       @param  combindedEpsilon     combined epsilons
       @param  particleIPosition    particle I position
       @param  particleJPosition    particle J position
       @param  force                output force

       @return energy for ixn

       --------------------------------------------------------------------------------------- */

    RealOpenMM calculatePairIxn( int siteI, int siteJ,
                                                          const std::vector<RealVec> & particlePositions,
                                                          const std::vector<std::vector<int> >& allParticleIndices,
                                                          std::vector<RealVec>& forces ) const;
};

// ---------------------------------------------------------------------------------------

#endif // _MBPolReferenceTwoBodyForce___
