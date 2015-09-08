#ifndef MBPOL_OPENMM_REFERENCE_KERNELS_H_
#define MBPOL_OPENMM_REFERENCE_KERNELS_H_

/* -------------------------------------------------------------------------- *
 *                              OpenMMMBPol                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
 * Authors:                                                                   *
 * Contributors:                                                              *
 *                                                                            *
 * This program is free software: you can redistribute it and/or modify       *
 * it under the terms of the GNU Lesser General Public License as published   *
 * by the Free Software Foundation, either version 3 of the License, or       *
 * (at your option) any later version.                                        *
 *                                                                            *
 * This program is distributed in the hope that it will be useful,            *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of             *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
 * GNU Lesser General Public License for more details.                        *
 *                                                                            *
 * You should have received a copy of the GNU Lesser General Public License   *
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.      *
 * -------------------------------------------------------------------------- */

#include "openmm/System.h"
#include "openmm/mbpolKernels.h"
#include "openmm/MBPolElectrostaticsForce.h"
#include "MBPolReferenceElectrostaticsForce.h"
#include "openmm/reference/ReferenceNeighborList.h"
#include "ReferenceThreeNeighborList.h"
#include "openmm/reference/SimTKOpenMMRealType.h"
#include <string>

using std::string;

namespace MBPolPlugin {
/**
 * This kernel is invoked by MBPolOneBodyForce to calculate the forces acting on the system and the energy of the system.
 */
class ReferenceCalcMBPolOneBodyForceKernel : public CalcMBPolOneBodyForceKernel {
public:
    ReferenceCalcMBPolOneBodyForceKernel(std::string name, const Platform& platform, const System& system);
    ~ReferenceCalcMBPolOneBodyForceKernel();
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the MBPolOneBodyForce this kernel will be used for
     */
    void initialize(const System& system, const MBPolOneBodyForce& force);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    double execute(ContextImpl& context, bool includeForces, bool includeEnergy);
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the MBPolOneBodyForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const MBPolOneBodyForce& force);
private:
    int numOneBodys;
    std::vector< std::vector<int> > allParticleIndices;
    const System& system;
    int usePBC;
};

/**
 * This kernel is invoked by MBPolElectrostaticsForce to calculate the forces acting on the system and the energy of the system.
 */
class ReferenceCalcMBPolElectrostaticsForceKernel : public CalcMBPolElectrostaticsForceKernel {
public:
    ReferenceCalcMBPolElectrostaticsForceKernel(std::string name, const Platform& platform, const System& system);
    ~ReferenceCalcMBPolElectrostaticsForceKernel();
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the MBPolElectrostaticsForce this kernel will be used for
     */
    void initialize(const System& system, const MBPolElectrostaticsForce& force);
    /**
     * Setup for MBPolReferenceElectrostaticsForce instance. 
     *
     * @param context        the current context
     *
     * @return pointer to initialized instance of MBPolReferenceElectrostaticsForce
     */
    MBPolReferenceElectrostaticsForce* setupMBPolReferenceElectrostaticsForce(ContextImpl& context );
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    double execute(ContextImpl& context, bool includeForces, bool includeEnergy);
    /** 
     * Calculate the electrostatic potential given vector of grid coordinates.
     *
     * @param context                      context
     * @param inputGrid                    input grid coordinates
     * @param outputElectrostaticPotential output potential 
     */
    void getElectrostaticPotential(ContextImpl& context, const std::vector< Vec3 >& inputGrid,
                                   std::vector< double >& outputElectrostaticPotential );

    /**
     * Get the system multipole moments.
     *
     * @param context                context 
     * @param outputElectrostaticsMonents vector of multipole moments:
                                     (charge,
                                      dipole_x, dipole_y, dipole_z,
                                      quadrupole_xx, quadrupole_xy, quadrupole_xz,
                                      quadrupole_yx, quadrupole_yy, quadrupole_yz,
                                      quadrupole_zx, quadrupole_zy, quadrupole_zz )
     */
    void getSystemElectrostaticsMoments(ContextImpl& context, std::vector< double >& outputElectrostaticsMoments);
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the MBPolElectrostaticsForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const MBPolElectrostaticsForce& force);

private:

    int numElectrostatics;
    MBPolElectrostaticsForce::NonbondedMethod nonbondedMethod;
    MBPolElectrostaticsForce::PolarizationType polarizationType;
    std::vector<RealOpenMM> charges;
    std::vector<RealOpenMM> dipoles;
    std::vector<RealOpenMM> quadrupoles;
    std::vector<RealOpenMM> tholes;
    std::vector<RealOpenMM> dampingFactors;
    std::vector<RealOpenMM> polarity;
    std::vector<int>   axisTypes;
    std::vector<int>   multipoleAtomZs;
    std::vector<int>   multipoleAtomXs;
    std::vector<int>   multipoleAtomYs;
    bool includeChargeRedistribution;
    std::vector< std::vector< std::vector<int> > > multipoleAtomCovalentInfo;

    int mutualInducedMaxIterations;
    RealOpenMM mutualInducedTargetEpsilon;

    bool usePme;
    RealOpenMM alphaEwald;
    RealOpenMM cutoffDistance;
    std::vector<int> pmeGridDimension;

    const System& system;
};

/**
 * This kernel is invoked to calculate the TwoBody forces acting on the system and the energy of the system.
 */
class ReferenceCalcMBPolTwoBodyForceKernel : public CalcMBPolTwoBodyForceKernel {
public:
    ReferenceCalcMBPolTwoBodyForceKernel(std::string name, const Platform& platform, const System& system);
    ~ReferenceCalcMBPolTwoBodyForceKernel();
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the MBPolTwoBodyForce this kernel will be used for
     */
    void initialize(const System& system, const MBPolTwoBodyForce& force);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    double execute(ContextImpl& context, bool includeForces, bool includeEnergy);
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the MBPolTwoBodyForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const MBPolTwoBodyForce& force);
private:
    int numParticles;
    int useCutoff;
    int usePBC;
    double cutoff;
    std::vector< std::vector<int> > allParticleIndices;
    const System& system;
    NeighborList* neighborList;
};

/**
 * This kernel is invoked to calculate the TwoBody forces acting on the system and the energy of the system.
 */
class ReferenceCalcMBPolThreeBodyForceKernel : public CalcMBPolThreeBodyForceKernel {
public:
    ReferenceCalcMBPolThreeBodyForceKernel(std::string name, const Platform& platform, const System& system);
    ~ReferenceCalcMBPolThreeBodyForceKernel();
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the MBPolThreeBodyForce this kernel will be used for
     */
    void initialize(const System& system, const MBPolThreeBodyForce& force);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    double execute(ContextImpl& context, bool includeForces, bool includeEnergy);
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the MBPolThreeBodyForce to copy the parameters from
     */
    void copyParametersToContext(ContextImpl& context, const MBPolThreeBodyForce& force);
private:
    int numParticles;
    int useCutoff;
    int usePBC;
    double cutoff;
    std::vector< std::vector<int> > allParticleIndices;
    const System& system;
    ThreeNeighborList* neighborList;
};

} // namespace MBPolPlugin

#endif /*MBPOL_OPENMM_REFERENCE_KERNELS_H*/
