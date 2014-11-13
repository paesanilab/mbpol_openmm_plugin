#ifndef OPENMM_MBPOL_Dispersion_FORCE_H_
#define OPENMM_MBPOL_Dispersion_FORCE_H_

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
#include "internal/windowsExportMBPol.h"
#include <vector>
#include <string>
#include <map>


using std::string;
using std::vector;

using std::map;
using std::pair;
using namespace OpenMM;

typedef map< pair<const string, const string>, pair<double, double> > c6d6Datatype;


namespace MBPolPlugin {

/**
 * This class implements a buffered 14-7 potential used to model van der Waals forces.
 * 
 * To use it, create an MBPolDispersionForce object then call addParticle() once for each particle.  After
 * a particle has been added, you can modify its force field parameters by calling setParticleParameters().
 * This will have no effect on Contexts that already exist unless you call updateParametersInContext().
 * 
 * A unique feature of this class is that the interaction site for a particle does not need to be
 * exactly at the particle's location.  Instead, it can be placed a fraction of the distance from that
 * particle to another one.  This is typically done for hydrogens to place the interaction site slightly
 * closer to the parent atom.  The fraction is known as the "reduction factor", since it reduces the distance
 * from the parent atom to the interaction site.
 */

class OPENMM_EXPORT_MBPOL MBPolDispersionForce : public Force {
public:
    /**
     * This is an enumeration of the different methods that may be used for handling long range nonbonded forces.
     */
    enum NonbondedMethod {
        /**
         * No cutoff is applied to nonbonded interactions.  The full set of N^2 interactions is computed exactly.
         * This necessarily means that periodic boundary conditions cannot be used.  This is the default.
         */
        NoCutoff = 0,
        /**
         * Periodic boundary conditions are used, so that each particle interacts only with the nearest periodic copy of
         * each other particle.  Interactions beyond the cutoff distance are ignored.
         */
        CutoffPeriodic = 1,
        CutoffNonPeriodic = 2,
    };

    /**
     * Create an MBPol DispersionForce.
     */
    MBPolDispersionForce();

    /**
     * Get the number of particles
     */
    int getNumParticles() const {
        return parameters.size();
    }

    /**
     * Set the force field parameters for a vdw particle.
     * 
     * @param particleIndex   the particle index
     * @param parentIndex     the index of the parent particle
     * @param sigma           vdw sigma
     * @param epsilon         vdw epsilon
     * @param reductionFactor the fraction of the distance along the line from the parent particle to this particle
     *                        at which the interaction site should be placed
     */
    void setParticleParameters(int particleIndex, string atomElement);

    /**
     * Get the force field parameters for a vdw particle.
     * 
     * @param particleIndex   the particle index
     * @param parentIndex     the index of the parent particle
     * @param sigma           vdw sigma
     * @param epsilon         vdw epsilon
     * @param reductionFactor the fraction of the distance along the line from the parent particle to this particle
     *                        at which the interaction site should be placed
     */
    void getParticleParameters(int particleIndex, string & atomElement) const;


    /**
     * Add the force field parameters for a vdw particle.
     * 
     * @param parentIndex     the index of the parent particle
     * @param sigma           vdw sigma
     * @param epsilon         vdw epsilon
     * @param reductionFactor the fraction of the distance along the line from the parent particle to this particle
     *                        at which the interaction site should be placed
     * @return index of added particle
     */
    int addParticle(string atomElement);

    void addDispersionParameters(string firstElement, string secondElement, double c6, double d6);

    c6d6Datatype getDispersionParameters( void ) const;

    int getNumMolecules(void) const;
    /**
     * Set the cutoff distance.
     */
    void setCutoff(double cutoff);

    /**
     * Get the cutoff distance.
     */
    double getCutoff(void) const;

    /**
     * Get the method used for handling long range nonbonded interactions.
     */
    NonbondedMethod getNonbondedMethod() const;

    /**
     * Set the method used for handling long range nonbonded interactions.
     */
    void setNonbondedMethod(NonbondedMethod method);
    /**
     * Update the per-particle parameters in a Context to match those stored in this Force object.  This method provides
     * an efficient method to update certain parameters in an existing Context without needing to reinitialize it.
     * Simply call setParticleParameters() to modify this object's parameters, then call updateParametersInState()
     * to copy them over to the Context.
     * 
     * The only information this method updates is the values of per-particle parameters.  All other aspects of the Force
     * (the nonbonded method, the cutoff distance, etc.) are unaffected and can only be changed by reinitializing the Context.
     */
    void updateParametersInContext(Context& context);

    /**
     * Get whether to add a contribution to the energy that approximately represents the effect of VdW
     * interactions beyond the cutoff distance.  The energy depends on the volume of the periodic box, and is only
     * applicable when periodic boundary conditions are used.  When running simulations at constant pressure, adding
     * this contribution can improve the quality of results.
     */
    bool getUseDispersionCorrection() const {
        return useDispersionCorrection;
    }


    /**
     * Set whether to add a contribution to the energy that approximately represents the effect of VdW
     * interactions beyond the cutoff distance.  The energy depends on the volume of the periodic box, and is only
     * applicable when periodic boundary conditions are used.  When running simulations at constant pressure, adding
     * this contribution can improve the quality of results.
     */
    void setUseDispersionCorrection(bool useCorrection) {
        useDispersionCorrection = useCorrection;
    }


protected:
    ForceImpl* createImpl() const;
private:

    class DispersionInfo;
    NonbondedMethod nonbondedMethod;
    double cutoff;

    std::vector<DispersionInfo> parameters;
    c6d6Datatype c6d6Data;
    bool useDispersionCorrection;

};

class MBPolDispersionForce::DispersionInfo {
public:
    string atomElement;

    DispersionInfo() {

    }
    DispersionInfo( string atomElement) :
        atomElement(atomElement)  {
    }
};

} // namespace MBPolPlugin

#endif /*OPENMM_MBPOL_VDW_FORCE_H_*/

