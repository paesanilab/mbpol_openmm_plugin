#ifndef MBPOL_OPENMM_KERNELS_H_
#define MBPOL_OPENMM_KERNELS_H_

/* -------------------------------------------------------------------------- *
 *                             OpenMMMBPol                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2012 Stanford University and the Authors.      *
 * Authors:                                                                   *
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

#include "OpenMMMBPol.h"
#include "openmm/KernelImpl.h"
#include "openmm/System.h"
#include "openmm/Platform.h"

#include <set>
#include <string>
#include <vector>

namespace MBPolPlugin {

/**
 * This kernel is invoked by MBPolTorsionForce to calculate the forces acting on the system and the energy of the system.
 */
class CalcMBPolOneBodyForceKernel : public KernelImpl {

public:

    static std::string Name() {
        return "CalcMBPolOneBodyForce";
    }

    CalcMBPolOneBodyForceKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }

    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the OneBodyForce this kernel will be used for
     */
    virtual void initialize(const OpenMM::System& system, const MBPolOneBodyForce& force) = 0;

    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    virtual double execute(ContextImpl& context, bool includeForces, bool includeEnergy) = 0;
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the MBPolOneBodyForce to copy the parameters from
     */
    virtual void copyParametersToContext(ContextImpl& context, const MBPolOneBodyForce& force) = 0;
};

/**
 * This kernel is invoked by MBPolElectrostaticsForce to calculate the forces acting on the system and the energy of the system.
 */
class CalcMBPolElectrostaticsForceKernel : public KernelImpl {

public:

    static std::string Name() {
        return "CalcMBPolElectrostaticsForce";
    }

    CalcMBPolElectrostaticsForceKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }

    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the ElectrostaticsForce this kernel will be used for
     */
    virtual void initialize(const OpenMM::System& system, const MBPolElectrostaticsForce& force) = 0;

    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    virtual double execute(ContextImpl& context, bool includeForces, bool includeEnergy) = 0;

    virtual void getElectrostaticPotential( ContextImpl& context, const std::vector< Vec3 >& inputGrid,
                                            std::vector< double >& outputElectrostaticPotential ) = 0;

    virtual void getSystemElectrostaticsMoments( ContextImpl& context, std::vector< double >& outputElectrostaticsMonents ) = 0;
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the MBPolElectrostaticsForce to copy the parameters from
     */
    virtual void copyParametersToContext(ContextImpl& context, const MBPolElectrostaticsForce& force) = 0;
};

/**
 * This kernel is invoked by MBPolTwoBodyForce to calculate the TwoBody forces acting on the system and the TwoBody energy of the system.
 */
class CalcMBPolTwoBodyForceKernel : public KernelImpl {
public:

    static std::string Name() {
        return "CalcMBPolTwoBodyForce";
    }

    CalcMBPolTwoBodyForceKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }

    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the GBSAOBCForce this kernel will be used for
     */
    virtual void initialize(const OpenMM::System& system, const MBPolTwoBodyForce& force) = 0;

    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    virtual double execute(ContextImpl& context, bool includeForces, bool includeEnergy) = 0;
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the MBPolTwoBodyForce to copy the parameters from
     */
    virtual void copyParametersToContext(ContextImpl& context, const MBPolTwoBodyForce& force) = 0;
};

/**
 * This kernel is invoked by MBPolTwoBodyForce to calculate the TwoBody forces acting on the system and the TwoBody energy of the system.
 */
class CalcMBPolThreeBodyForceKernel : public KernelImpl {
public:

    static std::string Name() {
        return "CalcMBPolThreeBodyForce";
    }

    CalcMBPolThreeBodyForceKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }

    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the GBSAOBCForce this kernel will be used for
     */
    virtual void initialize(const OpenMM::System& system, const MBPolThreeBodyForce& force) = 0;

    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    virtual double execute(ContextImpl& context, bool includeForces, bool includeEnergy) = 0;
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the MBPolThreeBodyForce to copy the parameters from
     */
    virtual void copyParametersToContext(ContextImpl& context, const MBPolThreeBodyForce& force) = 0;
};

class CalcMBPolDispersionForceKernel : public KernelImpl {
public:

    static std::string Name() {
        return "CalcMBPolDispersionForce";
    }

    CalcMBPolDispersionForceKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
    }

    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param force      the GBSAOBCForce this kernel will be used for
     */

    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    virtual double execute(ContextImpl& context, bool includeForces, bool includeEnergy) = 0;
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the MBPolDispersionForce to copy the parameters from
     */
};

} // namespace MBPolPlugin

#endif /*MBPOL_OPENMM_KERNELS_H*/
