#ifndef OPENMM_MBPOL_MULTIPOLE_FORCE_H_
#define OPENMM_MBPOL_MULTIPOLE_FORCE_H_

/* -------------------------------------------------------------------------- *
 *                              OpenMMMBPol                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2012 Stanford University and the Authors.      *
 * Authors: Mark Friedrichs, Peter Eastman                                    *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "openmm/Force.h"
#include "openmm/OpenMMException.h"
#include "internal/windowsExportMBPol.h"
#include "openmm/Vec3.h"

#include <sstream>
#include <vector>

using namespace OpenMM;

namespace MBPolPlugin {

/**
 * This class implements the MBPol multipole interaction.
 *
 * To use it, create an MBPolElectrostaticsForce object then call addElectrostatics() once for each atom.  After
 * an entry has been added, you can modify its force field parameters by calling setElectrostaticsParameters().
 * This will have no effect on Contexts that already exist unless you call updateParametersInContext().
 */


class OPENMM_EXPORT_MBPOL MBPolElectrostaticsForce : public Force {

public:

    enum NonbondedMethod {

        /**
         * No cutoff is applied to nonbonded interactions.  The full set of N^2 interactions is computed exactly.
         * This necessarily means that periodic boundary conditions cannot be used.  This is the default.
         */
        NoCutoff = 0,

        /**
         * Periodic boundary conditions are used, and Particle-Mesh Ewald (PME) summation is used to compute the interaction of each particle
         * with all periodic copies of every other particle.
         */
        PME = 1
    };

    /**
     * Create an MBPolElectrostaticsForce.
     */
    MBPolElectrostaticsForce();

    /**
     * Get the number of particles in the potential function
     */
    int getNumElectrostatics() const {
        return multipoles.size();
    }

    /**
     * Get the method used for handling long-range nonbonded interactions.
     */
    NonbondedMethod getNonbondedMethod() const;

    /**
     * Set the method used for handling long-range nonbonded interactions.
     */
    void setNonbondedMethod(NonbondedMethod method);

    /**
     * Get the cutoff distance (in nm) being used for nonbonded interactions.  If the NonbondedMethod in use
     * is NoCutoff, this value will have no effect.
     *
     * @return the cutoff distance, measured in nm
     */
    double getCutoffDistance(void) const;

    /**
     * Set the cutoff distance (in nm) being used for nonbonded interactions.  If the NonbondedMethod in use
     * is NoCutoff, this value will have no effect.
     *
     * @param distance    the cutoff distance, measured in nm
     */
    void setCutoffDistance(double distance);

    void setIncludeChargeRedistribution( bool chargeRedistribution );

    bool getIncludeChargeRedistribution( void ) const;

    void setTholeParameters( std::vector<double> tholeP) {
        tholeParameters=tholeP;
    }

    std::vector<double> getTholeParameters( void ) const {
        return tholeParameters;
    }

    /**
     * Get the Ewald alpha parameter.  If this is 0 (the default), a value is chosen automatically
     * based on the Ewald error tolerance.
     *
     * @return the Ewald alpha parameter
     */
    double getAEwald() const;

    /**
     * Set the Ewald alpha parameter.  If this is 0 (the default), a value is chosen automatically
     * based on the Ewald error tolerance.
     *
     * @param Ewald alpha parameter
     */
    void setAEwald(double aewald);

    /**
     * Get the B-spline order to use for PME charge spreading
     *
     * @return the B-spline order
     */
    int getPmeBSplineOrder() const;

    /**
     * Get the PME grid dimensions.  If Ewald alpha is 0 (the default), this is ignored and grid dimensions
     * are chosen automatically based on the Ewald error tolerance.
     *
     * @return the PME grid dimensions
     */
    void getPmeGridDimensions(std::vector<int>& gridDimension) const;

   /**
     * Set the PME grid dimensions.  If Ewald alpha is 0 (the default), this is ignored and grid dimensions
     * are chosen automatically based on the Ewald error tolerance.
     *
     * @param the PME grid dimensions
     */
    void setPmeGridDimensions(const std::vector<int>& gridDimension);

    /**
     * Add multipole-related info for a particle
     *
     * @param charge               the particle's charge
     * @param multipoleAtomZ       index of first atom used in constructing lab<->molecular frames
     * @param multipoleAtomX       index of second atom used in constructing lab<->molecular frames
     * @param multipoleAtomY       index of second atom used in constructing lab<->molecular frames
     * @param dampingFactor        dampingFactor parameter
     * @param polarity             polarity parameter
     *
     * @return the index of the particle that was added
     */
    int addElectrostatics(double charge,
                     int moleculeIndex, int atomType, double dampingFactor, double polarity);

    /**
     * Get the multipole parameters for a particle.
     *
     * @param index                the index of the atom for which to get parameters
     * @param charge               the particle's charge
     * @param molecularDipole      the particle's molecular dipole (vector of size 3)
     * @param molecularQuadrupole  the particle's molecular quadrupole (vector of size 9)
     * @param axisType             the particle's axis type
     * @param multipoleAtomZ       index of first atom used in constructing lab<->molecular frames
     * @param multipoleAtomX       index of second atom used in constructing lab<->molecular frames
     * @param multipoleAtomY       index of second atom used in constructing lab<->molecular frames
     * @param dampingFactor        dampingFactor parameter
     * @param polarity             polarity parameter
     */
void getElectrostaticsParameters(int index, double& charge,
                     int& moleculeIndex, int& atomType, double& dampingFactor, double& polarity ) const;
    /**
     * Set the multipole parameters for a particle.
     *
     * @param index                the index of the atom for which to set parameters
     * @param charge               the particle's charge
     * @param molecularDipole      the particle's molecular dipole (vector of size 3)
     * @param molecularQuadrupole  the particle's molecular quadrupole (vector of size 9)
     * @param axisType             the particle's axis type
     * @param multipoleAtomZ       index of first atom used in constructing lab<->molecular frames
     * @param multipoleAtomX       index of second atom used in constructing lab<->molecular frames
     * @param multipoleAtomY       index of second atom used in constructing lab<->molecular frames
     * @param polarity             polarity parameter
     */

    void setElectrostaticsParameters(int index, double charge,
                     int moleculeIndex, int atomType, double dampingFactor, double polarity );

    /**
     * Get the max number of iterations to be used in calculating the mutual induced dipoles
     *
     * @return max number of iterations
     */
    int getMutualInducedMaxIterations(void) const;

    /**
     * Set the max number of iterations to be used in calculating the mutual induced dipoles
     *
     * @param max number of iterations
     */
    void setMutualInducedMaxIterations(int inputMutualInducedMaxIterations);

    /**
     * Get the target epsilon to be used to test for convergence of iterative method used in calculating the mutual induced dipoles
     *
     * @return target epsilon
     */
    double getMutualInducedTargetEpsilon(void) const;

    /**
     * Set the target epsilon to be used to test for convergence of iterative method used in calculating the mutual induced dipoles
     *
     * @param target epsilon
     */
    void setMutualInducedTargetEpsilon(double inputMutualInducedTargetEpsilon);

    /**
     * Get the error tolerance for Ewald summation.  This corresponds to the fractional error in the forces
     * which is acceptable.  This value is used to select the grid dimensions and separation (alpha)
     * parameter so that the average error level will be less than the tolerance.  There is not a
     * rigorous guarantee that all forces on all atoms will be less than the tolerance, however.
     * 
     * This can be overridden by explicitly setting an alpha parameter and grid dimensions to use.
     */
    double getEwaldErrorTolerance() const;
    /**
     * Get the error tolerance for Ewald summation.  This corresponds to the fractional error in the forces
     * which is acceptable.  This value is used to select the grid dimensions and separation (alpha)
     * parameter so that the average error level will be less than the tolerance.  There is not a
     * rigorous guarantee that all forces on all atoms will be less than the tolerance, however.
     * 
     * This can be overridden by explicitly setting an alpha parameter and grid dimensions to use.
     */
    void setEwaldErrorTolerance(double tol);

    /**
     * Get the electrostatic potential.
     *
     * @param inputGrid    input grid points over which the potential is to be evaluated
     * @param context      context
     * @param outputElectrostaticPotential output potential
     */

    void getElectrostaticPotential(const std::vector< Vec3 >& inputGrid,
                                    Context& context, std::vector< double >& outputElectrostaticPotential);

    /**
     * Update the multipole parameters in a Context to match those stored in this Force object.  This method
     * provides an efficient method to update certain parameters in an existing Context without needing to reinitialize it.
     * Simply call setElectrostaticsParameters() to modify this object's parameters, then call updateParametersInState() to
     * copy them over to the Context.
     * 
     * This method has several limitations.  The only information it updates is the parameters of multipoles.
     * All other aspects of the Force (the nonbonded method, the cutoff distance, etc.) are unaffected and can only be
     * changed by reinitializing the Context.  Furthermore, this method cannot be used to add new multipoles,
     * only to change the parameters of existing ones.
     */
    void updateParametersInContext(Context& context);

protected:
    ForceImpl* createImpl() const;
private:
    NonbondedMethod nonbondedMethod;
    double cutoffDistance;
    double aewald;
    int pmeBSplineOrder;
    std::vector<int> pmeGridDimension;
    int mutualInducedMaxIterations;
    double mutualInducedTargetEpsilon;
    double scalingDistanceCutoff;
    double electricConstant;
    double ewaldErrorTol;
    bool includeChargeRedistribution;
    std::vector<double> tholeParameters;
    class ElectrostaticsInfo;
    std::vector<ElectrostaticsInfo> multipoles;
};

class MBPolElectrostaticsForce::ElectrostaticsInfo {
public:

    int moleculeIndex, atomType;
    double charge, dampingFactor, polarity;

    ElectrostaticsInfo() {
        atomType = moleculeIndex = -1;
        charge = dampingFactor = polarity = 0.0;


    }

    ElectrostaticsInfo(double charge,
                   int moleculeIndex, int atomType, double dampingFactor, double polarity) :
        charge(charge), moleculeIndex(moleculeIndex), atomType(atomType),
        dampingFactor(dampingFactor), polarity(polarity) {

    }
};

enum TholeIndices { TCC, TCD, TDD, TDDOH, TDDHH };

} // namespace MBPolPlugin

#endif /*OPENMM_MBPOL_MULTIPOLE_FORCE_H_*/
