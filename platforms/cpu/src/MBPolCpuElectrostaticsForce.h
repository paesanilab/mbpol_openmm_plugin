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

#ifndef __MBPolCpuElectrostaticsForce_H__
#define __MBPolCpuElectrostaticsForce_H__

#include <vector>
#include "openmm/internal/gmx_atomic.h"
//#include "openmm/cpu/RealVec.h"
//#include "openmm/cpu/SimTKOpenMMRealType.h"
#include "openmm/MBPolElectrostaticsForce.h"
#include "openmm/cpu/AlignedArray.h"
#include "openmm/internal/ThreadPool.h"
#include "openmm/internal/vectorize.h"
#include <map>
#include "openmm/reference/fftpack.h"
#include <complex>
#include <assert.h>

using std::vector;

typedef std::map< unsigned int, RealOpenMM> MapIntRealOpenMM;
typedef MapIntRealOpenMM::iterator MapIntRealOpenMMI;
typedef MapIntRealOpenMM::const_iterator MapIntRealOpenMMCI;


/**
 * 2-dimensional int vector
 */
class int2 {
public:
    /**
     * Create a int2 whose elements are all 0.
     */
    int2() {
        data[0] = data[1] = 0;
    }
    /**
     * Create a int2 with specified x, y components.
     */
    int2(int x, int y) {
        data[0] = x;
        data[1] = y;
    }
    int operator[](int index) const {
        assert(index >= 0 && index < 2);
        return data[index];
    }
    int& operator[](int index) {
        assert(index >= 0 && index < 2);
        return data[index];
    }

    // Arithmetic operators

    // unary plus
    int2 operator+() const {
        return int2(*this);
    }

    // plus
    int2 operator+(const int2& rhs) const {
        const int2& lhs = *this;
        return int2(lhs[0] + rhs[0], lhs[1] + rhs[1]);
    }

    int2& operator+=(const int2& rhs) {
        data[0] += rhs[0];
        data[1] += rhs[1];
        return *this;
    }

    int2& operator-=(const int2& rhs) {
        data[0] -= rhs[0];
        data[1] -= rhs[1];
        return *this;
    }

private:
    int data[2];
};

/**
 * 3-dimensional int vector
 */
class IntVec {
public:
    /**
     * Create a IntVec whose elements are all 0.
     */
    IntVec() {
        data[0] = data[1] = data[2] = 0;
    }
    /**
     * Create a IntVec with specified x, y, z, w components.
     */
    IntVec(int x, int y, int z) {
        data[0] = x;
        data[1] = y;
        data[2] = z;
    }
    int operator[](int index) const {
        assert(index >= 0 && index < 3);
        return data[index];
    }
    int& operator[](int index) {
        assert(index >= 0 && index < 3);
        return data[index];
    }

    // Arithmetic operators

    // unary plus
    IntVec operator+() const {
        return IntVec(*this);
    }

    // plus
    IntVec operator+(const IntVec& rhs) const {
        const IntVec& lhs = *this;
        return IntVec(lhs[0] + rhs[0], lhs[1] + rhs[1], lhs[2] + rhs[2] );
    }

    IntVec& operator+=(const IntVec& rhs) {
        data[0] += rhs[0];
        data[1] += rhs[1];
        data[2] += rhs[2];
        return *this;
    }

    IntVec& operator-=(const IntVec& rhs) {
        data[0] -= rhs[0];
        data[1] -= rhs[1];
        data[2] -= rhs[2];
        return *this;
    }

private:
    int data[3];
};

/**
 * 4-dimensional RealOpenMM vector
 */
class RealOpenMM4 {
public:
    /**
     * Create a RealOpenMM4 whose elements are all 0.
     */
    RealOpenMM4() {
        data[0] = data[1] = data[2] = data[3] = 0.0;
    }
    /**
     * Create a RealOpenMM4 with specified x, y, z, w components.
     */
    RealOpenMM4(RealOpenMM x, RealOpenMM y, RealOpenMM z, RealOpenMM w) {
        data[0] = x;
        data[1] = y;
        data[2] = z;
        data[3] = w;
    }
    RealOpenMM operator[](int index) const {
        assert(index >= 0 && index < 4);
        return data[index];
    }
    RealOpenMM& operator[](int index) {
        assert(index >= 0 && index < 4);
        return data[index];
    }

    // Arithmetic operators

    // unary plus
    RealOpenMM4 operator+() const {
        return RealOpenMM4(*this);
    }

    // plus
    RealOpenMM4 operator+(const RealOpenMM4& rhs) const {
        const RealOpenMM4& lhs = *this;
        return RealOpenMM4(lhs[0] + rhs[0], lhs[1] + rhs[1], lhs[2] + rhs[2],lhs[3] + rhs[3]);
    }

    RealOpenMM4& operator+=(const RealOpenMM4& rhs) {
        data[0] += rhs[0];
        data[1] += rhs[1];
        data[2] += rhs[2];
        data[3] += rhs[3];
        return *this;
    }

    RealOpenMM4& operator-=(const RealOpenMM4& rhs) {
        data[0] -= rhs[0];
        data[1] -= rhs[1];
        data[2] -= rhs[2];
        data[3] -= rhs[3];
        return *this;
    }

private:
    RealOpenMM data[4];
};

using namespace  OpenMM;
using namespace MBPolPlugin;

class MBPolCpuElectrostaticsForce {

   /**
    * MBPolCpuElectrostaticsForce is base class for ElectrostaticsForce calculations
    * MBPolCpuGeneralizedKirkwoodElectrostaticsForce  is derived class for Generalized Kirkwood calculations
    * MBPolCpuPmeElectrostaticsForce is derived class for PME calculations
    *
    * Below is a outline of the sequence of methods called to evaluate the force and energy
    * for each scenario: Generalized Kirkwood (GK) and PME.
    *
    * If 'virtual' appears before the method name, the method is overridden in one or more of the derived classes.
    *
    * calculateForceAndEnergy()                            calculate forces  and energy
    *
    *    setup()                                           rotate molecular multipole moments to lab frame
    *                                                      setup scaling maps and calculate induced dipoles (see calculateInducedDipoles below)
    *
    *    virtual calculateElectrostatic()                  calculate forces and torques
    *
    *                                                      GK case includes the following calls:
    *
    *                                                          MBPolCpuElectrostaticsForce::calculateElectrostatic()
    *                                                               loop over particle pairs: calculateElectrostaticPairIxn()
    *
    *                                                          TINKER's egk1a: calculateKirkwoodPairIxn()
    *
    *                                                          SASA force and energy: calculateCavityTermEnergyAndForces()
    *
    *                                                          TINKER's born1 (Born chain rule term): loop over particle pairs: calculateGrycukChainRulePairIxn()
    *
    *                                                          TINKER's ediff1: loop over particle pairs: calculateKirkwoodEDiffPairIxn()
    *
    *                                                      PME case includes the following calls:
    *
    *                                                          reciprocal [computeReciprocalSpaceInducedDipoleForceAndEnergy(),
    *                                                                      computeReciprocalSpaceFixedElectrostaticsForceAndEnergy]
    *
    *                                                          direct space calculations [calculatePmeDirectElectrostaticPairIxn()]
    *
    *                                                          self-energy [calculatePmeSelfEnergy()]
    *
    *                                                          torques [calculatePmeSelfTorque()]
    *
    *    mapTorqueToForce()                                map torques to forces
    *
    * setup()
    *    loadParticleData()                                load particle data (polarity, multipole moments, Thole factors, ...)
    *    checkChiral()                                     if needed, invert multipole moments at chiral centers
    *    applyRotationMatrix()                             rotate molecular multipole moments to lab frame
    *    calculateInducedDipoles()                         calculate induced dipoles
    *
    *
    * virtual calculateInducedDipoles()                    calculate induced dipoles:
    *                                                          field at each site due to fixed multipoles first calculated
    *                                                          if polarization type == Mutual, then loop until
    *                                                          induce dipoles converge.
    *                                                       For GK, include gkField in setup
    *                                                       For PME, base class method is used
    *
    *
    *     virtual zeroFixedElectrostaticsFields()                zero fixed multipole vectors; for GK includes zeroing of gkField vector
    *
    *     virtual calculateFixedElectrostaticsField()            calculate fixed multipole field -- particle pair loop
    *                                                       gkField also calculated for GK
    *                                                       for PME, reciprocal, direct space (particle pair loop) and self terms calculated
    *
    *
    *         virtual calculateFixedElectrostaticsFieldPairIxn() pair ixn for fixed multipole
    *                                                       gkField also calculated for GK
    *                                                       for PME, direct space ixn calculated here
    *
    *     virtual initializeInducedDipoles()                initialize induced dipoles; for PME, calculateReciprocalSpaceInducedDipoleField()
    *                                                       called in case polarization type == Direct
    *
    *     convergeInduceDipoles()                           loop until induced dipoles converge
    *
    *         updateInducedDipoleFields()                   update fields at each site due other induced dipoles
    *
    *           virtual calculateInducedDipoleFields()      calculate induced dipole field at each site by looping over particle pairs
    *                                                       for PME includes reciprocal space calculation calculateReciprocalSpaceInducedDipoleField(),
    *                                                       direct space calculateDirectInducedDipolePairIxns() and self terms
    *
    */

public:

    /**
     * This is an enumeration of the different methods that may be used for handling long range Electrostatics forces.
     */
    enum NonbondedMethod {

        /**
         * No cutoff is applied to the interactions.  The full set of N^2 interactions is computed exactly.
         * This necessarily means that periodic boundary conditions cannot be used.  This is the default.
         */
        NoCutoff = 0,

       /**
         * Periodic boundary conditions are used, and Particle-Mesh Ewald (PME) summation is used to compute the interaction of each particle
         * with all periodic copies of every other particle.
         */
        PME = 1
    };

    enum ChargeDerivativesIndicesFinal { vsH1f, vsH2f, vsMf };

    /**
     * Constructor
     *
     * @param nonbondedMethod nonbonded method
     */
    MBPolCpuElectrostaticsForce( NonbondedMethod nonbondedMethod, ThreadPool& threads);

    /**
     * Destructor
     *
     */
    virtual ~MBPolCpuElectrostaticsForce( ){};

    /**
     * Get nonbonded method.
     *
     * @return nonbonded method
     */
    NonbondedMethod getNonbondedMethod( void ) const;

    /**
     * Set nonbonded method.
     *
     * @param nonbondedMethod nonbonded method
     */
    void setNonbondedMethod( NonbondedMethod nonbondedMethod );

    /**
     * Get flag indicating if mutual induced dipoles are converged.
     *
     * @return nonzero if converged
     *
     */
    int getMutualInducedDipoleConverged( void ) const;

    /**
     * Get the number of iterations used in computing mutual induced dipoles.
     *
     * @return number of iterations
     *
     */
    int getMutualInducedDipoleIterations( void ) const;

    void setIncludeChargeRedistribution( bool includeChargeRedistribution );


    bool getIncludeChargeRedistribution( void ) const;

    void setTholeParameters( std::vector<RealOpenMM> tholeP) {
        _tholeParameters=tholeP;
    }

    std::vector<RealOpenMM> getTholeParameters( void ) const {
        assert(_tholeParameters.size() > 0);
        return _tholeParameters;
    }

    /**
     * Get the final epsilon for mutual induced dipoles.
     *
     *  @return epsilon
     *
     */
    RealOpenMM getMutualInducedDipoleEpsilon( void ) const;

    /**
     * Set the target epsilon for converging mutual induced dipoles.
     *
     * @param targetEpsilon target epsilon for converging mutual induced dipoles
     *
     */
    void setMutualInducedDipoleTargetEpsilon( RealOpenMM targetEpsilon );

    /**
     * Get the target epsilon for converging mutual induced dipoles.
     *
     * @return target epsilon for converging mutual induced dipoles
     *
     */
    RealOpenMM getMutualInducedDipoleTargetEpsilon( void ) const;

    /**
     * Set the maximum number of iterations to be executed in converging mutual induced dipoles.
     *
     * @param maximumMutualInducedDipoleIterations maximum number of iterations to be executed in converging mutual induced dipoles
     *
     */
    void setMaximumMutualInducedDipoleIterations( int maximumMutualInducedDipoleIterations );

    /**
     * Get the maximum number of iterations to be executed in converging mutual induced dipoles.
     *
     * @return maximum number of iterations to be executed in converging mutual induced dipoles
     *
     */
    int getMaximumMutualInducedDipoleIterations( void ) const;

    /**
     * Calculate force and energy.
     *
     * @param particlePositions         Cartesian coordinates of particles
     * @param charges                   scalar charges for each particle
     * @param dipoles                   molecular frame dipoles for each particle
     * @param quadrupoles               molecular frame quadrupoles for each particle
     * @param tholes                    Thole factors for each particle
     * @param dampingFactors            damping factors for each particle
     * @param polarity                  polarity for each particle
     * @param forces                    add forces to this vector
     *
     * @return energy
     */
    RealOpenMM calculateForceAndEnergy( const std::vector<OpenMM::RealVec>& particlePositions,
                                        const std::vector<RealOpenMM>& charges,
                                        const std::vector<int>& moleculeIndices,
                                        const std::vector<int>& atomTypes,
                                        const std::vector<RealOpenMM>& tholes,
                                        const std::vector<RealOpenMM>& dampingFactors,
                                        const std::vector<RealOpenMM>& polarity,
                                        std::vector<OpenMM::RealVec>& forces );

    /**
     * Calculate system multipole moments.
     *
     * @param masses                    particle masses
     * @param particlePositions         Cartesian coordinates of particles
     * @param charges                   scalar charges for each particle
     * @param tholes                    Thole factors for each particle
     * @param dampingFactors            dampling factors for each particle
     * @param polarity                  polarity for each particle
     * @param outputElectrostaticsMoments    output multipole moments
     */
    void calculateMBPolSystemElectrostaticsMoments( const std::vector<RealOpenMM>& masses,
                                                const std::vector<OpenMM::RealVec>& particlePositions,
                                                const std::vector<RealOpenMM>& charges,
                                                const std::vector<int>& moleculeIndices,
                                                const std::vector<int>& atomTypes,
                                                const std::vector<RealOpenMM>& tholes,
                                                const std::vector<RealOpenMM>& dampingFactors,
                                                const std::vector<RealOpenMM>& polarity,
                                                std::vector<RealOpenMM>& outputElectrostaticsMoments);

    /**
     * Calculate electrostatic potential at a set of grid points.
     *
     * @param particlePositions         Cartesian coordinates of particles
     * @param charges                   scalar charges for each particle
     * @param tholes                    Thole factors for each particle
     * @param dampingFactors            dampling factors for each particle
     * @param polarity                  polarity for each particle
     * @param input grid                input grid points to compute potential
     * @param outputPotential           output electrostatic potential
     */
    void calculateElectrostaticPotential( const std::vector<OpenMM::RealVec>& particlePositions,
                                          const std::vector<RealOpenMM>& charges,
                                          const std::vector<int>& moleculeIndices,
                                          const std::vector<int>& atomTypes,
                                          const std::vector<RealOpenMM>& tholes,
                                          const std::vector<RealOpenMM>& dampingFactors,
                                          const std::vector<RealOpenMM>& polarity,
                                          const std::vector<RealVec>& inputGrid,
                                          std::vector<RealOpenMM>& outputPotential );


	class MBPolCpuComputeForceTask;
    class ReduceThreadFieldTask;

protected:

    enum ElectrostaticsParticleDataEnum { PARTICLE_POSITION, PARTICLE_CHARGE,
                                     PARTICLE_THOLE, PARTICLE_DAMPING_FACTOR, PARTICLE_POLARITY, PARTICLE_FIELD,
                                     PARTICLE_FIELD_POLAR, GK_FIELD, PARTICLE_INDUCED_DIPOLE, PARTICLE_INDUCED_DIPOLE_POLAR };

    /*
     * Particle parameters and coordinates
     */
    class ElectrostaticsParticleData {
        public:
            unsigned int particleIndex;
            RealVec position;
            RealOpenMM charge;
            RealVec chargeDerivatives[3];
            unsigned int otherSiteIndex[3];
            RealOpenMM thole[5];
            RealOpenMM dampingFactor;
            RealOpenMM polarity;
            unsigned int moleculeIndex;
            unsigned int atomType;
    };

    /*
     * Helper class used in calculating induced dipoles
     */
    struct UpdateInducedDipoleFieldStruct {
            UpdateInducedDipoleFieldStruct( std::vector<OpenMM::RealVec>* inputFixed_E_Field, std::vector<OpenMM::RealVec>* inputInducedDipoles );
            std::vector<OpenMM::RealVec>* fixedElectrostaticsField;
            std::vector<OpenMM::RealVec>* inducedDipoles;
            std::vector<OpenMM::RealVec> inducedDipoleField;
    };

    unsigned int _numParticles;

    NonbondedMethod _nonbondedMethod;
    bool _includeChargeRedistribution;
    std::vector<RealOpenMM> _tholeParameters;
    RealOpenMM _electric;
    RealOpenMM _dielectric;

    std::vector<RealVec> _fixedElectrostaticsField;
    std::vector<RealVec> _fixedElectrostaticsFieldPolar;
    std::vector<RealVec> _inducedDipole;
    std::vector<RealVec> _inducedDipolePolar;

    int _mutualInducedDipoleConverged;
    int _mutualInducedDipoleIterations;
    int _maximumMutualInducedDipoleIterations;
    RealOpenMM  _mutualInducedDipoleEpsilon;
    RealOpenMM  _mutualInducedDipoleTargetEpsilon;
    RealOpenMM  _polarSOR;
    RealOpenMM  _debye;

	RealOpenMM * scale3;
	std::vector<int> scaleOffset;
	RealOpenMM * scale5;
    std::vector<UpdateInducedDipoleFieldStruct> * updateInducedDipoleFields;
    ThreadPool& threads;
	const std::vector<ElectrostaticsParticleData>* particleData;
    gmx_atomic_t* atomicCounter;
    // The following variables are used to make information accessible to the individual threads.
    //
    std::vector<AlignedArray<float> >* threadForce;
    std::vector<double> threadEnergy;

    std::vector<std::vector<std::vector<RealVec> > >* threadField;

    /**
     * Helper constructor method to centralize initialization of objects.
     *
     */
    void initialize( void );

    /**
     * Load particle data.
     *
     * @param particlePositions   particle coordinates
     * @param charges             charges
     * @param tholes              Thole parameters
     * @param dampingFactors      dampming factors
     * @param polarity            polarity
     * @param particleData        output data struct
     *
     */
    void loadParticleData( const std::vector<OpenMM::RealVec>& particlePositions,
                           const std::vector<RealOpenMM>& charges,
                           const std::vector<int>& moleculeIndices,
                           const std::vector<int>& atomTypes,
                           const std::vector<RealOpenMM>& tholes,
                           const std::vector<RealOpenMM>& dampingFactors,
                           const std::vector<RealOpenMM>& polarity,
                           std::vector<ElectrostaticsParticleData>& particleData ) const;

    void printPotential (std::vector<RealOpenMM> electrostaticPotential, RealOpenMM energy, std::string name, const std::vector<ElectrostaticsParticleData>& particleData );

    void computeWaterCharge(ElectrostaticsParticleData& particleO, ElectrostaticsParticleData& particleH1,
                   ElectrostaticsParticleData& particleH2,ElectrostaticsParticleData& particleM);

    /**
     * Calculate fixed multipole fields.
     *
     * @param particleData vector of particle data
     *
     */
    virtual void calculateFixedElectrostaticsField( const vector<ElectrostaticsParticleData>& particleData );

    /**
     * Set flag indicating if mutual induced dipoles are converged.
     *
     * @param converged nonzero if converged
     *
     */
    void setMutualInducedDipoleConverged( int converged );

    /**
     * Set number of iterations used in computing mutual induced dipoles.
     *
     * @param  number of iterations
     *
     */
    void setMutualInducedDipoleIterations( int iterations );

    /**
     * Set the final epsilon for mutual induced dipoles.
     *
     * @param epsilon
     *
     */
    void setMutualInducedDipoleEpsilon( RealOpenMM epsilon );

    /**
     * Calculate damped powers of 1/r.
     *
     * @param  particleI           index of particleI
     * @param  particleJ           index of particleJ
     * @param  dScale              output d-scale factor
     * @param  pScale              output p-scale factor
     */


    RealOpenMM getAndScaleInverseRs(  const ElectrostaticsParticleData& particleI,
                                                                        const ElectrostaticsParticleData& particleK,
                                                              RealOpenMM r, bool justScale, int interactionOrder, int interactionType) const;

    void getAndScaleInverseRs13justScaleTCC(  const ElectrostaticsParticleData& particleI,
                                                                    const ElectrostaticsParticleData& particleK,
                                                     const RealOpenMM pgamma, RealOpenMM r, RealOpenMM * scale3, RealOpenMM * scale5) const;

    /**
     * Zero fixed multipole fields.
     */
    virtual void zeroFixedElectrostaticsFields( void );

    /**
     * Calculate electric field at particle I due fixed multipoles at particle J and vice versa
     * (field at particle J due fixed multipoles at particle I).
     *
     * @param particleI               positions and parameters (charge, labFrame dipoles, quadrupoles, ...) for particle I
     * @param particleJ               positions and parameters (charge, labFrame dipoles, quadrupoles, ...) for particle J
     * @param dScale                  d-scale value for i-j interaction
     * @param pScale                  p-scale value for i-j interaction
     */
    virtual void calculateFixedElectrostaticsFieldPairIxn( const ElectrostaticsParticleData& particleI, const ElectrostaticsParticleData& particleJ);

    /**
     * Initialize induced dipoles
     *
     * @param updateInducedDipoleFields vector of UpdateInducedDipoleFieldStruct containing input induced dipoles and output fields
     */
    virtual void initializeInducedDipoles( std::vector<UpdateInducedDipoleFieldStruct>& updateInducedDipoleFields );

    /**
     * Calculate field at particle I due induced dipole at particle J and vice versa
     * (field at particle J due induced dipole at particle I).
     *
     * @param particleI               index of particle I
     * @param particleJ               index of particle J
     * @param rr3                     damped 1/r^3 factor
     * @param rr5                     damped 1/r^5 factor
     * @param delta                   delta of particle positions: particleJ.x - particleI.x, ...
     * @param inducedDipole           vector of induced dipoles
     * @param field                   vector of induced dipole fields
     */
    void calculateInducedDipolePairIxn( unsigned int particleI, unsigned int particleJ,
                                        const RealOpenMM & rr3, const RealOpenMM & rr5, const RealVec& delta,
                                        const std::vector<RealVec>& inducedDipole,
                                        std::vector<RealVec>& field ) const;

    virtual void precomputeScale35( const std::vector<ElectrostaticsParticleData>& particleData, RealOpenMM scale3[], RealOpenMM scale5[] );

    /**
     * Calculate induced dipole fields.
     *
     * @param particleData              vector of particle positions and parameters (charge, labFrame dipoles, quadrupoles, ...)
     * @param updateInducedDipoleFields vector of UpdateInducedDipoleFieldStruct containing input induced dipoles and output fields
     */
	virtual void calculateInducedDipoleFields( ThreadPool& threads, int threadIndex );
    void reduceThreadField( ThreadPool& threads, int threadIndex );
    /**
     * Converge induced dipoles.
     *
     * @param particleData              vector of particle positions and parameters (charge, labFrame dipoles, quadrupoles, ...)
     * @param updateInducedDipoleFields vector of UpdateInducedDipoleFieldStruct containing input induced dipoles and output fields
     */
    void convergeInduceDipoles( const std::vector<ElectrostaticsParticleData>& particleData,
                                std::vector<UpdateInducedDipoleFieldStruct>& calculateInducedDipoleField );

	void threadComputeForce(ThreadPool& threads, int threadIndex);
    /**
     * Update fields due to induced dipoles for each particle.
     *
     * @param particleData              vector of particle positions and parameters (charge, labFrame dipoles, quadrupoles, ...)
     * @param updateInducedDipoleFields vector of UpdateInducedDipoleFieldStruct containing input induced dipoles and output fields
     */
    RealOpenMM runUpdateInducedDipoleFields( const std::vector<ElectrostaticsParticleData>& particleData,
                                          std::vector<UpdateInducedDipoleFieldStruct>& calculateInducedDipoleField,
						 const RealOpenMM scale3[], const RealOpenMM scale5[]);

    /**
     * Update induced dipole for a particle given updated induced dipole field at the site.
     *
     * @param particleI                 positions and parameters (charge, labFrame dipoles, quadrupoles, ...) for particle I
     * @param fixedElectrostaticsField       fields due fixed multipoles at each site
     * @param inducedDipoleField        fields due induced dipoles at each site
     * @param inducedDipoles            output vector of updated induced dipoles
     */
    RealOpenMM updateInducedDipole( const std::vector<ElectrostaticsParticleData>& particleI,
                                    const std::vector<RealVec>& fixedElectrostaticsField,
                                    const std::vector<RealVec>& inducedDipoleField,
                                    std::vector<RealVec>& inducedDipoles);

    /**
     * Calculate induced dipoles.
     *
     * @param particleData      vector of particle positions and parameters (charge, labFrame dipoles, quadrupoles, ...)
     */
    virtual void calculateInducedDipoles( const std::vector<ElectrostaticsParticleData>& particleData );

    /**
     * Setup:
     *        if needed invert multipole moments at chiral centers
     *        rotate molecular multipole moments to lab frame
     *        setup scaling maps and
     *        calculate induced dipoles (see calculateInducedDipoles below)
     *
     * @param particlePositions         Cartesian coordinates of particles
     * @param charges                   scalar charges for each particle
     * @param tholes                    Thole factors for each particle
     * @param dampingFactors            dampling factors for each particle
     * @param polarity                  polarity for each particle
     * @param particleData              output vector of parameters (charge, labFrame dipoles, quadrupoles, ...) for particles
     *
     */
    void setup( const std::vector<OpenMM::RealVec>& particlePositions,
                const std::vector<RealOpenMM>& charges,
                const std::vector<int>& moleculeIndices,
                const std::vector<int>& atomTypes,
                const std::vector<RealOpenMM>& tholes,
                const std::vector<RealOpenMM>& dampingFactors,
                const std::vector<RealOpenMM>& polarity,
                std::vector<ElectrostaticsParticleData>& particleData );

    /**
     * Calculate electrostatic interaction between particles I and K.
     *
     * @param particleI         positions and parameters (charge, labFrame dipoles, quadrupoles, ...) for particle I
     * @param particleK         positions and parameters (charge, labFrame dipoles, quadrupoles, ...) for particle K
     * @param scalingFactors    scaling factors for interaction
     * @param forces            vector of particle forces to be updated
     */
    RealOpenMM calculateElectrostaticPairIxn(
            const std::vector<ElectrostaticsParticleData>& particleData,
                                                                                    unsigned int iIndex,
                                                                                    unsigned int kIndex,
                                                                                    std::vector<OpenMM::RealVec>& forces) const;

    /**
     * Calculate electrostatic forces
     *
     * @param particleData            vector of parameters (charge, labFrame dipoles, quadrupoles, ...) for particles
     * @param forces                  output forces
     *
     * @return energy
     */
    virtual RealOpenMM calculateElectrostatic( const std::vector<ElectrostaticsParticleData>& particleData,
                                               std::vector<OpenMM::RealVec>& forces );

    /**
     * Normalize a RealVec
     *
     * @param vectorToNormalize vector to normalize
     *
     * @return norm of vector on input
     *
     */
    RealOpenMM normalizeRealVec( RealVec& vectorToNormalize ) const;

    /**
     * Initialize vector of RealOpenMM (size=numParticles)
     *
     * @param vectorToInitialize vector to initialize
     *
     */
    void initializeRealOpenMMVector( vector<RealOpenMM>& vectorToInitialize ) const;

    /**
     * Initialize vector of RealVec (size=numParticles)
     *
     * @param vectorToInitialize vector to initialize
     *
     */
    void initializeRealVecVector( vector<RealVec>& vectorToInitialize ) const;

    /**
     * Copy vector of RealVec
     *
     * @param inputVector  vector to copy
     * @param outputVector output vector
     *
     */
    void copyRealVecVector( const std::vector<OpenMM::RealVec>& inputVector, std::vector<OpenMM::RealVec>& outputVector ) const;

    /**
     * Calculate potential at grid point due to a particle
     *
     * @param particleData            vector of parameters (charge, labFrame dipoles, quadrupoles, ...) for particles
     * @param gridPoint               grid point
     *
     * @return potential at grid point
     *
     */
    RealOpenMM calculateElectrostaticPotentialForParticleGridPoint( const ElectrostaticsParticleData& particleI, const RealVec& gridPoint ) const;

    /**
     * Apply periodic boundary conditions to difference in positions
     *
     * @param deltaR  difference in particle positions; modified on output after applying PBC
     *
     */
    virtual void getPeriodicDelta( RealVec& deltaR ) const {};
};

class MBPolCpuPmeElectrostaticsForce : public MBPolCpuElectrostaticsForce {

public:

    /**
     * Constructor
     *
     */
    MBPolCpuPmeElectrostaticsForce( ThreadPool& threads);

    /**
     * Destructor
     *
     */
    ~MBPolCpuPmeElectrostaticsForce( );

    /**
     * Get cutoff distance.
     *
     * @return cutoff distance
     *
     */
    RealOpenMM getCutoffDistance( void ) const;

    /**
     * Set cutoff distance.
     *
     * @return cutoff distance
     *
     */
    void setCutoffDistance( RealOpenMM cutoffDistance );

    /**
     * Get alpha used in Ewald summation.
     *
     * @return alpha
     *
     */
    RealOpenMM getAlphaEwald( void ) const;

    /**
     * Set alpha used in Ewald summation.
     *
     * @return alpha
     *
     */
    void setAlphaEwald( RealOpenMM alphaEwald );

    /**
     * Get PME grid dimensions.
     *
     * @param pmeGridDimensions contains PME grid dimensions upon return

     *
     */
    void getPmeGridDimensions( std::vector<int>& pmeGridDimensions ) const;

    /**
     * Set PME grid dimensions.
     *
     * @param pmeGridDimensions input PME grid dimensions
     *
     */
    void setPmeGridDimensions( std::vector<int>& pmeGridDimensions );

    /**
     * Set periodic box size.
     *
     * @param boxSize box dimensions
     */
     void setPeriodicBoxSize( RealVec& boxSize );

protected:

     /**
      * Resize PME arrays.
      *
      */
     void resizePmeArrays( void );

    RealOpenMM ewaldScalingReal( RealOpenMM r, int interactionOrder) const;


     /**
      * Calculate direct space electrostatic interaction between particles I and J.
      *
      * @param particleI         positions and parameters (charge, labFrame dipoles, quadrupoles, ...) for particle I
      * @param particleJ         positions and parameters (charge, labFrame dipoles, quadrupoles, ...) for particle J
      * @param scalingFactors    scaling factors for interaction
      * @param forces            vector of particle forces to be updated
      */
     RealOpenMM calculatePmeDirectElectrostaticPairIxn( const std::vector<ElectrostaticsParticleData>& particleData,
                            unsigned int iIndex, unsigned int jIndex,
                                                        std::vector<RealVec>& forces, std::vector<RealOpenMM>& electrostaticPotential ) const;


private:

    static const int MBPOL_PME_ORDER;
    static const RealOpenMM SQRT_PI;

    RealOpenMM _alphaEwald;
    RealOpenMM _cutoffDistance;
    RealOpenMM _cutoffDistanceSquared;

    RealVec _invPeriodicBoxSize;
    RealVec _periodicBoxSize;

    int _totalGridSize;
    IntVec _pmeGridDimensions;

    fftpack_t   _fftplan;

    unsigned int _pmeGridSize;
    t_complex* _pmeGrid;

    std::vector<RealOpenMM> _pmeBsplineModuli[3];
    std::vector<RealOpenMM4> _thetai[3];
    std::vector<IntVec> _iGrid;
    std::vector<RealOpenMM> _phi;
    std::vector<RealOpenMM> _phid;
    std::vector<RealOpenMM> _phip;
    std::vector<RealOpenMM> _phidp;
    std::vector<int> _pmeAtomRange;
    std::vector<int2> _pmeAtomGridIndex;
    std::vector<RealOpenMM4> _pmeBsplineTheta;
    std::vector<RealOpenMM4> _pmeBsplineDtheta;



    /**
     * Zero Pme grid.
     */
    void initializePmeGrid( void );

    /**
     * Modify input vector of differences in particle positions for periodic boundary conditions.
     *
     * @param delta                   input vector of difference in particle positios; on output adjusted for
     *                                periodic boundary conditions
     */
    void getPeriodicDelta( RealVec& deltaR ) const;

    /**
     * Get PME scale.
     *
     */
    void getPmeScale( RealVec& scale ) const;

    /**
     * Calculate damped inverse distances.
     *
     * @param particleI               positions and parameters (charge, labFrame dipoles, quadrupoles, ...) for particle I
     * @param particleJ               positions and parameters (charge, labFrame dipoles, quadrupoles, ...) for particle J
     * @param dScale                  d-scale value for i-j interaction
     * @param pScale                  p-scale value for i-j interaction
     * @param dampedDInverseDistances damped inverse distances (drr3,drr5,drr7 in udirect2a() in TINKER)
     * @param dampedPInverseDistances damped inverse distances (prr3,prr5,prr7 in udirect2a() in TINKER)
     */
    void getDampedInverseDistances( const ElectrostaticsParticleData& particleI,
                                                                  const ElectrostaticsParticleData& particleJ,
                                                                  RealOpenMM dscale, RealOpenMM pscale, RealOpenMM r,
                                                                  std::vector<RealOpenMM>& dampedDInverseDistances,
                                                                  std::vector<RealOpenMM>& dampedPInverseDistances ) const;

    /**
     * Initialize B-spline moduli.
     *
     */
    void initializeBSplineModuli( void );

    /**
     * Calculate direct-space field at site I due fixed multipoles at site J and vice versa.
     *
     * @param particleI               positions and parameters (charge, labFrame dipoles, quadrupoles, ...) for particle I
     * @param particleJ               positions and parameters (charge, labFrame dipoles, quadrupoles, ...) for particle J
     * @param dScale                  d-scale value for i-j interaction
     * @param pScale                  p-scale value for i-j interaction
     */
    void calculateFixedElectrostaticsFieldPairIxn( const ElectrostaticsParticleData& particleI, const ElectrostaticsParticleData& particleJ);

    /**
     * Calculate fixed multipole fields.
     *
     * @param particleData vector particle data
     *
     */
    void calculateFixedElectrostaticsField( const vector<ElectrostaticsParticleData>& particleData );

    /**
     * This is called from computeMBPolBsplines().  It calculates the spline coefficients for a single atom along a single axis.
     *
     * @param thetai output spline coefficients
     * @param w offset from grid point
     */
    void computeBSplinePoint(  std::vector<RealOpenMM4>& thetai, RealOpenMM w  );

    /**
     * Compute bspline coefficients.
     *
     * @param particleData   vector of particle positions and parameters (charge, labFrame dipoles, quadrupoles, ...)
     */
    void computeMBPolBsplines( const std::vector<ElectrostaticsParticleData>& particleData );

    /**
     * For each grid point, find the range of sorted atoms associated with that point.
     *
     * @param particleData              vector of particle positions and parameters (charge, labFrame dipoles, quadrupoles, ...)
     */
    void findMBPolAtomRangeForGrid( const vector<ElectrostaticsParticleData>& particleData );

    /**
     * Get grid point given grid index.
     *
     * @param gridIndex  input grid index
     * @param gridPoint  output grid point
     */
    void getGridPointGivenGridIndex( int gridIndex, IntVec& gridPoint ) const;

    /**
     * Compute induced dipole grid value.
     *
     * @param particleData            vector of particle positions and parameters (charge, labFrame dipoles, quadrupoles, ...)
     * @param particleGridIndices     particle grid indices
     * @param scale                   integer grid dimension/box size for each dimension
     * @param ix                      x-dimension offset value
     * @param iy                      y-dimension offset value
     * @param gridPoint               grid point for which value is to be computed
     * @param inputInducedDipole      induced dipole value
     * @param inputInducedDipolePolar induced dipole value
     */
     RealOpenMM computeFixedElectrostaticssGridValue( const vector<ElectrostaticsParticleData>& particleData,
                                                 const int2& particleGridIndices, const RealVec& scale, int ix, int iy, const IntVec& gridPoint ) const;

    /**
     * Spread fixed multipoles onto PME grid.
     *
     * @param particleData vector of particle positions and parameters (charge, labFrame dipoles, quadrupoles, ...)
     */
    void spreadFixedElectrostaticssOntoGrid( const vector<ElectrostaticsParticleData>& particleData );

    /**
     * Perform reciprocal convolution.
     *
     */
    void performMBPolReciprocalConvolution( void );

    /**
     * Compute reciprocal potential due fixed multipoles at each particle site.
     *
     */
    void computeFixedPotentialFromGrid(void );

    /**
     * Compute reciprocal potential due fixed multipoles at each particle site.
     *
     */
    void computeInducedPotentialFromGrid( void );

    /**
     * Calculate reciprocal space energy and force due to fixed multipoles.
     *
     * @param particleData    vector of particle positions and parameters (charge, labFrame dipoles, quadrupoles, ...)
     * @param forces          upon return updated vector of forces
     *
     * @return energy
     */
    RealOpenMM computeReciprocalSpaceFixedElectrostaticsForceAndEnergy( const std::vector<ElectrostaticsParticleData>& particleData,
                                                                   std::vector<RealVec>& forces, std::vector<RealOpenMM>& electrostaticPotential) const;

    /**
     * Set reciprocal space fixed multipole fields.
     *
     */
    void recordFixedElectrostaticsField( void );

    void precomputeScale35( const std::vector<ElectrostaticsParticleData>& particleData, RealOpenMM scale3[], RealOpenMM scale5[] );
    /**
     * Compute the potential due to the reciprocal space PME calculation for induced dipoles.
     *
     * @param updateInducedDipoleFields vector of UpdateInducedDipoleFieldStruct containing input induced dipoles and output fields
     */
    void calculateReciprocalSpaceInducedDipoleField( std::vector<UpdateInducedDipoleFieldStruct>& updateInducedDipoleFields );

    /**
     * Calculate field at particleI due to induced dipole at particle J and vice versa.
     *
     * @param iIndex        particle I index
     * @param jIndex        particle J index
     * @param preFactor1    first factor used in calculating field
     * @param preFactor2    second factor used in calculating field
     * @param deltaR        delta in particle positions after adjusting for periodic boundary conditions
     * @param inducedDipole vector of induced dipoles
     * @param field         vector of field at each particle due induced dipole of other particles
     */
    void calculateDirectInducedDipolePairIxn( unsigned int iIndex, unsigned int jIndex,
                                              RealOpenMM preFactor1, RealOpenMM preFactor2, const RealVec& delta,
                                              const std::vector<RealVec>& inducedDipole,
                                              std::vector<RealVec>& field ) const;

    /**
     * Calculate direct space field at particleI due to induced dipole at particle J and vice versa for
     * inducedDipole and inducedDipolePolar.
     *
     * @param particleI                 positions and parameters (charge, labFrame dipoles, quadrupoles, ...) for particle I
     * @param particleJ                 positions and parameters (charge, labFrame dipoles, quadrupoles, ...) for particle J
     * @param updateInducedDipoleFields vector of UpdateInducedDipoleFieldStruct containing input induced dipoles and output fields
     */
    void calculateDirectInducedDipolePairIxns( const ElectrostaticsParticleData& particleI,
                                               const ElectrostaticsParticleData& particleJ,
                                               std::vector<UpdateInducedDipoleFieldStruct>& updateInducedDipoleFields,
                                               RealOpenMM precomputedScale3,
                                               RealOpenMM precomputedScale5 );

    /**
     * Initialize induced dipoles
     *
     * @param updateInducedDipoleFields vector of UpdateInducedDipoleFieldStruct containing input induced dipoles and output fields
     */
    void initializeInducedDipoles( std::vector<UpdateInducedDipoleFieldStruct>& updateInducedDipoleFields );

    /**
     * Compute induced dipole grid value.
     *
     * @param atomIndices             indices of first and last atom contiputing to grid point value
     * @param scale                   integer grid dimension/box size for each dimension
     * @param ix                      x-dimension offset value
     * @param iy                      y-dimension offset value
     * @param gridPoint               grid point for which value is to be computed
     * @param inputInducedDipole      induced dipole value
     * @param inputInducedDipolePolar induced dipole polar value
     */
    t_complex computeInducedDipoleGridValue( const int2& atomIndices, const RealVec& scale, int ix, int iy, const IntVec& gridPoint,
                                             const std::vector<RealVec>& inputInducedDipole,
                                             const std::vector<RealVec>& inputInducedDipolePolar ) const;

    /**
     * Spread induced dipoles onto grid.
     *
     * @param inputInducedDipole      induced dipole value
     * @param inputInducedDipolePolar induced dipole polar value
     */
    void spreadInducedDipolesOnGrid( const std::vector<RealVec>& inputInducedDipole,
                                     const std::vector<RealVec>& inputInducedDipolePolar );

    /**
     * Calculate induced dipole fields.
     *
     * @param particleData              vector of particle positions and parameters (charge, labFrame dipoles, quadrupoles, ...)
     * @param updateInducedDipoleFields vector of UpdateInducedDipoleFieldStruct containing input induced dipoles and output fields
     */
	void calculateInducedDipoleFields( ThreadPool& threads, int threadIndex );

    /**
     * Set reciprocal space induced dipole fields.
     *
     * @param field       reciprocal space output induced dipole field value at each site
     * @param fieldPolar  reciprocal space output induced dipole polar field value at each site
     *
     */
    void recordInducedDipoleField( vector<RealVec>& field, vector<RealVec>& fieldPolar );

    /**
     * Compute Pme self energy.
     *
     * @param particleData            vector of parameters (charge, labFrame dipoles, quadrupoles, ...) for particles
     */
    RealOpenMM calculatePmeSelfEnergy( const std::vector<ElectrostaticsParticleData>& particleData, std::vector<RealVec>& forces, std::vector<RealOpenMM>& electrostaticPotential ) const;

    /**
     * Calculate reciprocal space energy/force for dipole interaction.
     *
     * @param particleData      vector of particle positions and parameters (charge, labFrame dipoles, quadrupoles, ...)
     * @param forces            vector of particle forces to be updated
     */
     RealOpenMM computeReciprocalSpaceInducedDipoleForceAndEnergy( const std::vector<ElectrostaticsParticleData>& particleData,
                                                                   std::vector<RealVec>& forces, std::vector<RealOpenMM>& electrostaticPotential) const;

    /**
     * Calculate electrostatic forces.
     *
     * @param particleData            vector of parameters (charge, labFrame dipoles, quadrupoles, ...) for particles
     * @param forces                  output forces
     *
     * @return energy
     */
    RealOpenMM calculateElectrostatic( const std::vector<ElectrostaticsParticleData>& particleData,
                                       std::vector<OpenMM::RealVec>& forces );

};

#endif // _MBPolCpuElectrostaticsForce___
