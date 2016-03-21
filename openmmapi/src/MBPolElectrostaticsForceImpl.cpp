/* -------------------------------------------------------------------------- *
 *                               OpenMMMBPol                                 *
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

#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/MBPolElectrostaticsForceImpl.h"
#include "openmm/mbpolKernels.h"
#include <stdio.h>
#include <cmath>

using namespace std;
using namespace  OpenMM;
using namespace MBPolPlugin;

using std::vector;

MBPolElectrostaticsForceImpl::MBPolElectrostaticsForceImpl(const MBPolElectrostaticsForce& owner) : owner(owner) {
}

MBPolElectrostaticsForceImpl::~MBPolElectrostaticsForceImpl() {
}

void MBPolElectrostaticsForceImpl::initialize(ContextImpl& context) {

    const System& system = context.getSystem();
    if (owner.getNumElectrostatics() != system.getNumParticles())
        throw OpenMMException("MBPolElectrostaticsForce must have exactly as many particles as the System it belongs to.");

    // check cutoff < 0.5*boxSize

    if (owner.getNonbondedMethod() == MBPolElectrostaticsForce::PME) {
        Vec3 boxVectors[3];
        system.getDefaultPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
        double cutoff = owner.getCutoffDistance();
        if (cutoff > 0.5*boxVectors[0][0] || cutoff > 0.5*boxVectors[1][1] || cutoff > 0.5*boxVectors[2][2])
            throw OpenMMException("MBPolElectrostaticsForce: The cutoff distance cannot be greater than half the periodic box size.");
    }   

    kernel = context.getPlatform().createKernel(CalcMBPolElectrostaticsForceKernel::Name(), context);
    kernel.getAs<CalcMBPolElectrostaticsForceKernel>().initialize(context.getSystem(), owner);
}

double MBPolElectrostaticsForceImpl::calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) {
    if ((groups&(1<<owner.getForceGroup())) != 0)
        return kernel.getAs<CalcMBPolElectrostaticsForceKernel>().execute(context, includeForces, includeEnergy);
    return 0.0;
}

std::vector<std::string> MBPolElectrostaticsForceImpl::getKernelNames() {
    std::vector<std::string> names;
    names.push_back(CalcMBPolElectrostaticsForceKernel::Name());
    return names;
}

void MBPolElectrostaticsForceImpl::getElectrostaticPotential( ContextImpl& context, const std::vector< Vec3 >& inputGrid,
                                                          std::vector< double >& outputElectrostaticPotential ){
    kernel.getAs<CalcMBPolElectrostaticsForceKernel>().getElectrostaticPotential(context, inputGrid, outputElectrostaticPotential);
}

void MBPolElectrostaticsForceImpl::updateParametersInContext(ContextImpl& context) {
    kernel.getAs<CalcMBPolElectrostaticsForceKernel>().copyParametersToContext(context, owner);
}
